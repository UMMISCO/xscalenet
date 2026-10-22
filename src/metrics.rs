//! Topological network metrics over a reconstructed [`Network`].
//!
//! Builds an undirected (optionally weighted) graph and computes graph-level and
//! node-level metrics: degree / strength, density, connected components, local &
//! global clustering, betweenness and closeness centrality (Brandes / BFS,
//! parallel over sources), and modularity-based communities (Louvain).

#[cfg(feature = "parallel")]
use rayon::prelude::*;
#[cfg(not(feature = "parallel"))]
use crate::par::prelude::*;
use std::collections::VecDeque;

use crate::network::Network;

/// Undirected graph with per-edge weights, built from a [`Network`].
pub struct Graph {
    pub n: usize,
    /// Adjacency list: `adj[i]` = `(neighbor, weight)` pairs.
    pub adj: Vec<Vec<(usize, f64)>>,
}

impl Graph {
    /// Build from a network over `n_nodes`. If `weighted`, edge weight is
    /// `|ecorr|` when present else `presence_score` else 1.0; otherwise all 1.0.
    pub fn from_network(net: &Network, n_nodes: usize, weighted: bool) -> Graph {
        let mut adj = vec![Vec::new(); n_nodes];
        for e in &net.edges {
            let w = if weighted {
                e.ecorr
                    .map(f64::abs)
                    .or(e.presence_score)
                    .filter(|v| *v > 0.0)
                    .unwrap_or(1.0)
            } else {
                1.0
            };
            adj[e.source].push((e.target, w));
            adj[e.target].push((e.source, w));
        }
        Graph { n: n_nodes, adj }
    }

    fn n_edges(&self) -> usize {
        self.adj.iter().map(|a| a.len()).sum::<usize>() / 2
    }
}

/// Per-node metrics. Centralities are normalized to `[0, 1]`.
#[derive(Clone, Debug, Default)]
pub struct NodeMetrics {
    pub degree: usize,
    pub strength: f64,
    pub betweenness: f64,
    pub closeness: f64,
    pub clustering: f64,
    pub community: usize,
    /// Degree within the node's own community (modular degree).
    pub modular_degree: usize,
    /// Betweenness within the node's community subgraph (normalized).
    pub modular_betweenness: f64,
}

/// Graph-level summary plus per-node metrics.
#[derive(Clone, Debug)]
pub struct GraphMetrics {
    pub n_nodes: usize,
    pub n_edges: usize,
    pub density: f64,
    pub avg_degree: f64,
    pub n_components: usize,
    pub largest_component: usize,
    /// Global clustering (transitivity): `3·triangles / connected triples`.
    pub global_clustering: f64,
    pub modularity: f64,
    pub n_communities: usize,
    pub nodes: Vec<NodeMetrics>,
}

/// Compute all metrics.
pub fn compute(g: &Graph) -> GraphMetrics {
    let n = g.n;
    let m = g.n_edges();

    let degree: Vec<usize> = g.adj.iter().map(|a| a.len()).collect();
    let strength: Vec<f64> = g.adj.iter().map(|a| a.iter().map(|&(_, w)| w).sum()).collect();

    let (n_components, largest_component) = components(g);
    let clustering = local_clustering(g);
    let global_clustering = transitivity(g);
    let betweenness = betweenness(g);
    let closeness = closeness(g);
    let (community, modularity) = louvain(g);
    let modular_degree: Vec<usize> = (0..n)
        .map(|i| g.adj[i].iter().filter(|&&(u, _)| community[u] == community[i]).count())
        .collect();
    let modular_betweenness = modular_betweenness(g, &community);
    let n_communities = {
        let mut c = community.clone();
        c.sort_unstable();
        c.dedup();
        c.len()
    };

    let nodes = (0..n)
        .map(|i| NodeMetrics {
            degree: degree[i],
            strength: strength[i],
            betweenness: betweenness[i],
            closeness: closeness[i],
            clustering: clustering[i],
            community: community[i],
            modular_degree: modular_degree[i],
            modular_betweenness: modular_betweenness[i],
        })
        .collect();

    GraphMetrics {
        n_nodes: n,
        n_edges: m,
        density: if n > 1 {
            2.0 * m as f64 / (n as f64 * (n as f64 - 1.0))
        } else {
            0.0
        },
        avg_degree: if n > 0 {
            2.0 * m as f64 / n as f64
        } else {
            0.0
        },
        n_components,
        largest_component,
        global_clustering,
        modularity,
        n_communities,
        nodes,
    }
}

/// Connected components: returns (count, size of the largest).
fn components(g: &Graph) -> (usize, usize) {
    let mut seen = vec![false; g.n];
    let mut count = 0;
    let mut largest = 0;
    for s in 0..g.n {
        if seen[s] {
            continue;
        }
        count += 1;
        let mut size = 0;
        let mut q = VecDeque::from([s]);
        seen[s] = true;
        while let Some(u) = q.pop_front() {
            size += 1;
            for &(v, _) in &g.adj[u] {
                if !seen[v] {
                    seen[v] = true;
                    q.push_back(v);
                }
            }
        }
        largest = largest.max(size);
    }
    (count, largest)
}

/// Local clustering coefficient per node (fraction of connected neighbor pairs).
fn local_clustering(g: &Graph) -> Vec<f64> {
    (0..g.n)
        .into_par_iter()
        .map(|i| {
            let neigh: Vec<usize> = g.adj[i].iter().map(|&(v, _)| v).collect();
            let k = neigh.len();
            if k < 2 {
                return 0.0;
            }
            let nset: std::collections::HashSet<usize> = neigh.iter().copied().collect();
            let mut links = 0usize;
            for (a, &u) in neigh.iter().enumerate() {
                for &w in &neigh[a + 1..] {
                    // count u-w edge once (u < w handled by iterating pairs)
                    if g.adj[u].iter().any(|&(x, _)| x == w) && nset.contains(&w) {
                        links += 1;
                    }
                }
            }
            2.0 * links as f64 / (k as f64 * (k as f64 - 1.0))
        })
        .collect()
}

/// Global transitivity = 3·triangles / connected triples.
fn transitivity(g: &Graph) -> f64 {
    let mut triangles = 0u64;
    let mut triples = 0u64;
    for i in 0..g.n {
        let neigh: Vec<usize> = g.adj[i].iter().map(|&(v, _)| v).collect();
        let k = neigh.len() as u64;
        triples += k.saturating_mul(k.saturating_sub(1)) / 2;
        for (a, &u) in neigh.iter().enumerate() {
            for &w in &neigh[a + 1..] {
                if g.adj[u].iter().any(|&(x, _)| x == w) {
                    triangles += 1;
                }
            }
        }
    }
    // Each triangle counted once per apex above; triples over-counts by design.
    if triples == 0 {
        0.0
    } else {
        3.0 * triangles as f64 / triples as f64
    }
}

/// Betweenness centrality (Brandes, unweighted, undirected), normalized to [0,1].
fn betweenness(g: &Graph) -> Vec<f64> {
    let n = g.n;
    // Brandes over every source; sources accumulate independently, so this is a
    // per-source fold + sum reduction. Rayon's two-closure fold/reduce has no
    // sequential analogue, so the single-thread build runs the plain fold.
    #[cfg(feature = "parallel")]
    let bc = (0..n)
        .into_par_iter()
        .fold(
            || vec![0.0f64; n],
            |mut acc, s| {
                brandes_source(g, s, &mut acc);
                acc
            },
        )
        .reduce(
            || vec![0.0f64; n],
            |mut a, b| {
                for i in 0..n {
                    a[i] += b[i];
                }
                a
            },
        );
    #[cfg(not(feature = "parallel"))]
    let bc = {
        let mut acc = vec![0.0f64; n];
        for s in 0..n {
            brandes_source(g, s, &mut acc);
        }
        acc
    };
    // Undirected: each pair counted twice; normalize by (n-1)(n-2).
    let norm = if n > 2 {
        (n as f64 - 1.0) * (n as f64 - 2.0)
    } else {
        1.0
    };
    bc.into_iter().map(|v| v / norm).collect()
}

fn brandes_source(g: &Graph, s: usize, delta_out: &mut [f64]) {
    let n = g.n;
    let mut sigma = vec![0.0f64; n];
    let mut dist = vec![-1i64; n];
    let mut pred: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut stack: Vec<usize> = Vec::new();
    sigma[s] = 1.0;
    dist[s] = 0;
    let mut q = VecDeque::from([s]);
    while let Some(v) = q.pop_front() {
        stack.push(v);
        for &(w, _) in &g.adj[v] {
            if dist[w] < 0 {
                dist[w] = dist[v] + 1;
                q.push_back(w);
            }
            if dist[w] == dist[v] + 1 {
                sigma[w] += sigma[v];
                pred[w].push(v);
            }
        }
    }
    let mut delta = vec![0.0f64; n];
    while let Some(w) = stack.pop() {
        for &v in &pred[w] {
            delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
        }
        if w != s {
            delta_out[w] += delta[w];
        }
    }
}

/// Betweenness computed within each community's induced subgraph, normalized
/// per community. Nodes bridging within their module score high; cross-module
/// structure is ignored (that is the point — modular vs total centrality).
fn modular_betweenness(g: &Graph, community: &[usize]) -> Vec<f64> {
    let n = g.n;
    let mut out = vec![0.0f64; n];
    let ncomm = community.iter().copied().max().map(|c| c + 1).unwrap_or(0);
    // Group nodes by community.
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); ncomm];
    for (i, &c) in community.iter().enumerate() {
        members[c].push(i);
    }
    for group in &members {
        let k = group.len();
        if k < 3 {
            continue; // no interior shortest paths to bridge
        }
        // Induced subgraph on `group`, relabeled 0..k.
        let local: std::collections::HashMap<usize, usize> =
            group.iter().enumerate().map(|(l, &gid)| (gid, l)).collect();
        let mut sub = Graph {
            n: k,
            adj: vec![Vec::new(); k],
        };
        for &gi in group {
            let li = local[&gi];
            for &(gj, w) in &g.adj[gi] {
                if let Some(&lj) = local.get(&gj) {
                    sub.adj[li].push((lj, w));
                }
            }
        }
        let bc = betweenness(&sub); // already normalized by (k-1)(k-2)
        for (l, &gid) in group.iter().enumerate() {
            out[gid] = bc[l];
        }
    }
    out
}

/// Closeness centrality (unweighted, per connected component), normalized.
fn closeness(g: &Graph) -> Vec<f64> {
    let n = g.n;
    (0..n)
        .into_par_iter()
        .map(|s| {
            let mut dist = vec![-1i64; n];
            dist[s] = 0;
            let mut q = VecDeque::from([s]);
            let mut sum = 0i64;
            let mut reached = 0i64;
            while let Some(v) = q.pop_front() {
                for &(w, _) in &g.adj[v] {
                    if dist[w] < 0 {
                        dist[w] = dist[v] + 1;
                        sum += dist[w];
                        reached += 1;
                        q.push_back(w);
                    }
                }
            }
            if sum == 0 {
                0.0
            } else {
                // Wasserman-Faust normalization for disconnected graphs.
                (reached as f64 / (n as f64 - 1.0)) * (reached as f64 / sum as f64)
            }
        })
        .collect()
}

/// Louvain community detection. Returns (community per node, modularity).
fn louvain(g: &Graph) -> (Vec<usize>, f64) {
    let n = g.n;
    if n == 0 {
        return (Vec::new(), 0.0);
    }
    // Work on a weighted supergraph that we coarsen each pass.
    let mut super_adj: Vec<Vec<(usize, f64)>> = g.adj.clone();
    // node -> original members map (which base nodes each super-node contains).
    let mut membership: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut comm = vec![0usize; n]; // final base-node -> community

    loop {
        let sn = super_adj.len();
        let two_m: f64 = super_adj
            .iter()
            .map(|a| a.iter().map(|&(_, w)| w).sum::<f64>())
            .sum();
        if two_m == 0.0 {
            break;
        }
        // node strengths, self-loops
        let k: Vec<f64> = super_adj
            .iter()
            .map(|a| a.iter().map(|&(_, w)| w).sum())
            .collect();
        let mut node_comm: Vec<usize> = (0..sn).collect();
        let mut comm_tot: Vec<f64> = k.clone(); // sum of strengths of nodes in community
        let mut improved_any = false;
        let mut moved = true;
        while moved {
            moved = false;
            for v in 0..sn {
                let cv = node_comm[v];
                // weight from v to each neighboring community
                let mut w_to: std::collections::HashMap<usize, f64> =
                    std::collections::HashMap::new();
                for &(u, w) in &super_adj[v] {
                    if u != v {
                        *w_to.entry(node_comm[u]).or_insert(0.0) += w;
                    }
                }
                // remove v from its community
                comm_tot[cv] -= k[v];
                let w_to_cv = w_to.get(&cv).copied().unwrap_or(0.0);
                // best community by modularity gain
                let mut best_c = cv;
                let mut best_gain = w_to_cv - k[v] * comm_tot[cv] / two_m;
                for (&c, &wc) in &w_to {
                    let gain = wc - k[v] * comm_tot[c] / two_m;
                    if gain > best_gain + 1e-12 {
                        best_gain = gain;
                        best_c = c;
                    }
                }
                comm_tot[best_c] += k[v];
                if best_c != cv {
                    node_comm[v] = best_c;
                    moved = true;
                    improved_any = true;
                }
            }
        }
        // Relabel communities 0..c and update base-node membership.
        let mut labels: Vec<usize> = node_comm.clone();
        let mut uniq = labels.clone();
        uniq.sort_unstable();
        uniq.dedup();
        let remap: std::collections::HashMap<usize, usize> =
            uniq.iter().enumerate().map(|(i, &c)| (c, i)).collect();
        for l in &mut labels {
            *l = remap[l];
        }
        for (sv, base_nodes) in membership.iter().enumerate() {
            for &b in base_nodes {
                comm[b] = labels[sv];
            }
        }
        if !improved_any || uniq.len() == sn {
            break;
        }
        // Coarsen: aggregate nodes in the same community into a super-node.
        let nc = uniq.len();
        let mut new_membership: Vec<Vec<usize>> = vec![Vec::new(); nc];
        for (sv, base_nodes) in membership.into_iter().enumerate() {
            new_membership[labels[sv]].extend(base_nodes);
        }
        let mut new_adj: Vec<std::collections::HashMap<usize, f64>> =
            vec![std::collections::HashMap::new(); nc];
        for v in 0..sn {
            let cv = labels[v];
            for &(u, w) in &super_adj[v] {
                *new_adj[cv].entry(labels[u]).or_insert(0.0) += w;
            }
        }
        super_adj = new_adj
            .into_iter()
            .map(|hm| hm.into_iter().collect())
            .collect();
        membership = new_membership;
    }

    let q = modularity(g, &comm);
    (comm, q)
}

/// Newman modularity of a partition on the base graph.
fn modularity(g: &Graph, comm: &[usize]) -> f64 {
    let two_m: f64 = g.adj.iter().map(|a| a.iter().map(|&(_, w)| w).sum::<f64>()).sum();
    if two_m == 0.0 {
        return 0.0;
    }
    let k: Vec<f64> = g.adj.iter().map(|a| a.iter().map(|&(_, w)| w).sum()).collect();
    let ncomm = comm.iter().copied().max().map(|c| c + 1).unwrap_or(0);
    let mut inside = vec![0.0f64; ncomm]; // sum of edge weights inside community (counts both directions)
    let mut tot = vec![0.0f64; ncomm]; // sum of strengths in community
    for i in 0..g.n {
        tot[comm[i]] += k[i];
        for &(j, w) in &g.adj[i] {
            if comm[j] == comm[i] {
                inside[comm[i]] += w; // each undirected edge counted twice
            }
        }
    }
    let mut q = 0.0;
    for c in 0..ncomm {
        q += inside[c] / two_m - (tot[c] / two_m).powi(2);
    }
    q
}

// (placeholder helper removed)

#[cfg(test)]
mod tests {
    use super::*;

    fn from_edges(n: usize, edges: &[(usize, usize)]) -> Graph {
        let mut adj = vec![Vec::new(); n];
        for &(a, b) in edges {
            adj[a].push((b, 1.0));
            adj[b].push((a, 1.0));
        }
        Graph { n, adj }
    }

    #[test]
    fn two_triangles_bridge() {
        // {0,1,2} and {3,4,5} triangles joined by the bridge 2-3.
        let g = from_edges(
            6,
            &[(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (2, 3)],
        );
        let m = compute(&g);
        assert_eq!(m.n_edges, 7);
        assert_eq!(m.n_components, 1);
        assert_eq!(m.largest_component, 6);
        // Two clear communities.
        assert_eq!(m.n_communities, 2);
        assert!(m.modularity > 0.3, "modularity {}", m.modularity);
        assert_eq!(m.nodes[0].community, m.nodes[1].community);
        assert_ne!(m.nodes[0].community, m.nodes[3].community);
        // A pure triangle node is fully clustered; the bridge node is not.
        assert!((m.nodes[0].clustering - 1.0).abs() < 1e-9);
        assert!(m.nodes[2].clustering < 1.0);
        // The bridge nodes carry the most betweenness.
        assert!(m.nodes[2].betweenness > m.nodes[0].betweenness);
        assert!(m.nodes[3].betweenness > m.nodes[5].betweenness);
    }

    #[test]
    fn path_graph_center_is_most_between() {
        // 0-1-2-3-4 path: node 2 is the most central.
        let g = from_edges(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
        let m = compute(&g);
        assert_eq!(m.n_components, 1);
        let bc: Vec<f64> = m.nodes.iter().map(|x| x.betweenness).collect();
        assert!(bc[2] > bc[1] && bc[1] > bc[0]);
        assert!((m.nodes[0].betweenness - 0.0).abs() < 1e-12); // endpoints
    }
}
