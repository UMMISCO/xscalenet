//! The end-to-end reconstruction pipeline as a library function — no I/O.
//!
//! [`run`] takes an already-loaded [`Dataset`] and a [`RunConfig`] and returns a
//! typed [`PipelineResult`] (per-method and consensus [`Network`]s). The CLI,
//! tests, and any embedder format that result themselves.

use crate::data::Dataset;
use crate::network::{Consensus, Network, NetworkKind, PipelineResult, RankedEdge, SpectralStats};
use crate::scalenet::{self, EdgeList, Method, Params};
use crate::scs::{self, MethodInfo};

/// What to reconstruct: which methods, with which parameters.
#[derive(Clone)]
pub struct RunConfig {
    pub methods: Vec<Method>,
    pub params: Params,
}

/// Run the spectral reconstruction + consensus pipeline.
///
/// `progress` receives coarse stage messages (e.g. `"# Reconstructing with
/// aracne..."`); pass a no-op closure to run silently.
pub fn run(
    data: &Dataset,
    config: &RunConfig,
    mut progress: impl FnMut(&str),
) -> Result<PipelineResult, String> {
    let prep = scalenet::prepare(data, &config.params)?;
    progress(&format!(
        "# Kept k={} eigenvectors, {} subgraphs, m={} vars/subgraph",
        prep.spectral.k,
        prep.subgraphs.len(),
        prep.subset_m
    ));

    // Precompute per-column ranks once (reused for every edge's ecorr).
    let ranks = data.rank_matrix();

    // Reconstruct each method (one edge list per threshold), then fill ecorr.
    let mut method_lists: Vec<Vec<EdgeList>> = Vec::with_capacity(config.methods.len());
    for &method in &config.methods {
        progress(&format!("# Reconstructing with {}...", method.name()));
        let mut lists = scalenet::run_method(data, &prep, method, &config.params);
        for el in &mut lists {
            scalenet::fill_ecorr(data, &ranks, el);
        }
        method_lists.push(lists);
    }

    let per_method: Vec<Vec<Network>> = method_lists
        .iter()
        .map(|lists| {
            lists
                .iter()
                .map(|el| Network::from_edge_list(el, NetworkKind::PerMethod))
                .collect()
        })
        .collect();

    // Consensus per threshold.
    let p = data.ncol();
    let mut consensus = Vec::with_capacity(config.params.pres_freq_thresh.len());
    for (ti, thr) in config.params.pres_freq_thresh.iter().enumerate() {
        let pairs: Vec<(MethodInfo, &EdgeList)> = config
            .methods
            .iter()
            .enumerate()
            .map(|(mi, &m)| (MethodInfo::for_method(m), &method_lists[mi][ti]))
            .collect();
        let (cons_el, raw) = scs::consensus(data, &pairs);
        let network = Network::from_edge_list(&cons_el, NetworkKind::Consensus);
        progress(&format!(
            "# Consensus @ presFreq {thr}: {} edges",
            network.edges.len()
        ));
        let ranked: Vec<RankedEdge> = raw.iter().map(|r| RankedEdge::from_raw(r, p)).collect();
        consensus.push(Consensus { network, ranked });
    }

    Ok(PipelineResult {
        nodes: data.names.clone(),
        methods: config.methods.clone(),
        thresholds: config.params.pres_freq_thresh.clone(),
        per_method,
        consensus,
        spectral: SpectralStats {
            k: prep.spectral.k,
            subset_m: prep.subset_m,
            n_subgraphs: prep.subgraphs.len(),
        },
        subgraphs: prep.subgraphs,
    })
}
