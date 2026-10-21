//! The reconstruction must be a function of the variable *set*, not of the
//! order the variables happen to appear in the input file.
//!
//! Before the `Lsym` fix (see `src/spectral.rs`) it was not: `Lrw = I - D^-1 W`
//! is asymmetric, and decomposing its lower triangle — which is what R's
//! `eigen(symmetric = TRUE)` does — divides each off-diagonal entry by the
//! degree of whichever endpoint has the larger index. Permuting the columns
//! changed the matrix, hence the eigenvectors, the spectral subsets, and every
//! edge. These tests pin the fixed behaviour and document the legacy one.

use std::collections::BTreeSet;

use xscalenet::data::Dataset;
use xscalenet::scalenet::{self, Method, Params};
use xscalenet::{discretize, spectral};

fn params(legacy: bool) -> Params {
    Params {
        var_perc: 0.2,
        eigen_perc: None,
        pres_freq_thresh: vec![0.3],
        epsilon: 0.001,
        restart: 21,
        seed: 6196,
        max_iter: 10_000,
        iss: 1.0,
        maxp: usize::MAX,
        iterative_eigen: false,
        legacy_lrw: legacy,
    }
}

fn pop2mat() -> Dataset {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/pop2mat.txt");
    let text = std::fs::read_to_string(path).unwrap();
    let (names, raw, nrow) = Dataset::parse_tsv_raw(&text).unwrap();
    discretize::discretize_to_dataset(names, &raw, nrow, 5)
}

/// A deterministic presence/absence matrix with real block structure but no
/// duplicate (or mutual-information-equivalent) variables, so its spectrum is
/// non-degenerate and the `G+`/`G-` cuts do not land inside tied blocks.
///
/// The bundled `pop2mat` fixture cannot be used for invariance: it has 19
/// duplicate columns, which make several eigenvalues exactly degenerate and put
/// three of four cuts inside tied blocks. There the subsets are genuinely
/// ill-defined — see `pop2mat_ambiguity_is_reported`.
fn synthetic() -> Dataset {
    let (blocks, per_block, nrow) = (3usize, 12usize, 120usize);
    let p = blocks * per_block;
    let mut state = 0x9E3779B97F4A7C15u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    // Per-block latent presence pattern, then per-variable noisy copies with a
    // variable-specific flip rate so no two columns coincide.
    let latents: Vec<Vec<u32>> = (0..blocks)
        .map(|_| (0..nrow).map(|_| (next() % 2) as u32).collect())
        .collect();
    let mut cols: Vec<Vec<u32>> = Vec::with_capacity(p);
    let mut names: Vec<String> = Vec::with_capacity(p);
    for (b, latent) in latents.iter().enumerate() {
        for v in 0..per_block {
            let flip = 5 + 2 * v as u64; // 5%..27% disagreement with the latent
            let col: Vec<u32> = (0..nrow)
                .map(|i| {
                    if next() % 100 < flip {
                        1 - latent[i]
                    } else {
                        latent[i]
                    }
                })
                .collect();
            names.push(format!("v{:02}_{}", b * per_block + v, b));
            cols.push(col);
        }
    }
    let levels = cols
        .iter()
        .map(|c| c.iter().max().copied().unwrap_or(0) + 1)
        .collect();
    let d = Dataset { names, cols, levels, nrow };
    for i in 0..p {
        for j in (i + 1)..p {
            assert_ne!(d.cols[i], d.cols[j], "synthetic fixture has duplicate columns");
        }
    }
    d
}

/// Guards the fixture: if this ever fails the invariance tests below are testing
/// an input whose subsets are not well defined, and their failure would be
/// meaningless.
#[test]
fn synthetic_fixture_is_well_conditioned() {
    let data = synthetic();
    let prep = scalenet::prepare(&data, &params(false)).unwrap();
    assert_eq!(
        prep.spectral.degenerate_pairs, 0,
        "fixture has degenerate retained eigenvalues"
    );
    assert_eq!(
        prep.boundary_ties, 0,
        "fixture has a G+/G- cut inside a tied block"
    );
    assert!(prep.spectral.k >= 3, "fixture is too trivial: k = {}", prep.spectral.k);
}

/// The ambiguity that duplicate variables create must be *reported*, not hidden:
/// no tie-break rule can make an ill-posed selection reproducible, so the run
/// tells the user instead.
#[test]
fn pop2mat_ambiguity_is_reported() {
    let data = pop2mat();
    let prep = scalenet::prepare(&data, &params(false)).unwrap();
    assert!(
        prep.boundary_ties > 0,
        "pop2mat has 19 duplicate columns; its G+/G- cuts should be flagged as tied"
    );
}

/// A deterministic shuffle of the column order (no `rand` dependency).
fn permutation(p: usize, seed: u64) -> Vec<usize> {
    let mut state = seed | 1;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut perm: Vec<usize> = (0..p).collect();
    for i in (1..p).rev() {
        perm.swap(i, (next() % (i as u64 + 1)) as usize);
    }
    perm
}

fn permute(d: &Dataset, perm: &[usize]) -> Dataset {
    Dataset {
        names: perm.iter().map(|&i| d.names[i].clone()).collect(),
        cols: perm.iter().map(|&i| d.cols[i].clone()).collect(),
        levels: perm.iter().map(|&i| d.levels[i]).collect(),
        nrow: d.nrow,
    }
}

/// Subgraphs as a set of name-sets: the identity of each subgraph, with both the
/// order of the subgraphs and the order within them discarded.
fn subset_names(d: &Dataset, subs: &[Vec<usize>]) -> BTreeSet<BTreeSet<String>> {
    subs.iter()
        .map(|s| s.iter().map(|&g| d.names[g].clone()).collect())
        .collect()
}

/// Present edges as unordered name pairs.
fn edge_names(d: &Dataset, method: Method, legacy: bool) -> BTreeSet<(String, String)> {
    let prep = scalenet::prepare(d, &params(legacy)).unwrap();
    let lists = scalenet::run_method(d, &prep, method, &params(legacy));
    lists[0]
        .edges
        .keys()
        .map(|&idx| {
            let (a, b) = xscalenet::util::pair_at(idx, d.ncol());
            let (x, y) = (d.names[a].clone(), d.names[b].clone());
            if x < y { (x, y) } else { (y, x) }
        })
        .collect()
}

#[test]
fn spectral_subsets_are_column_order_invariant() {
    let data = synthetic();
    let base = {
        let prep = scalenet::prepare(&data, &params(false)).unwrap();
        subset_names(&data, &prep.subgraphs)
    };
    assert!(!base.is_empty(), "no subgraphs to compare");

    for seed in [1u64, 2, 3, 12345, 99991] {
        let perm = permutation(data.ncol(), seed);
        let pd = permute(&data, &perm);
        let prep = scalenet::prepare(&pd, &params(false)).unwrap();
        assert_eq!(
            subset_names(&pd, &prep.subgraphs),
            base,
            "spectral subsets changed under column permutation (seed {seed})"
        );
    }
}

#[test]
fn aracne_edges_are_column_order_invariant() {
    let data = synthetic();
    let base = edge_names(&data, Method::Aracne, false);
    assert!(!base.is_empty(), "no ARACNE edges to compare");

    for seed in [1u64, 2, 3, 12345, 99991] {
        let perm = permutation(data.ncol(), seed);
        let pd = permute(&data, &perm);
        assert_eq!(
            edge_names(&pd, Method::Aracne, false),
            base,
            "ARACNE edge set changed under column permutation (seed {seed})"
        );
    }
}

#[test]
fn bayes_hc_edges_are_column_order_invariant() {
    let data = synthetic();
    let base = edge_names(&data, Method::BayesHc, false);
    assert!(!base.is_empty(), "no hill-climbing edges to compare");

    for seed in [1u64, 2, 3] {
        let perm = permutation(data.ncol(), seed);
        let pd = permute(&data, &perm);
        assert_eq!(
            edge_names(&pd, Method::BayesHc, false),
            base,
            "hill-climbing edge set changed under column permutation (seed {seed})"
        );
    }
}

/// The eigenvectors we hand to `prepare` must be eigenvectors of the *random-walk*
/// Laplacian `Lrw = I - D^-1 W` — the asymmetric matrix ScaleNet nominally
/// decomposes — not of the symmetric `Lsym` we actually decompose. This checks
/// the `v_Lrw = D^-1/2 u_Lsym` back-mapping against a direct multiplication by
/// the asymmetric `Lrw`.
#[test]
fn patched_eigenvectors_are_eigenvectors_of_lrw() {
    let data = pop2mat();
    let cols: Vec<&[u32]> = data.cols.iter().map(|c| c.as_slice()).collect();
    let sp = spectral::compute(
        &cols,
        &data.levels,
        &data.names,
        data.ncol(),
        None,
        (0.2 * data.ncol() as f64).ceil() as usize,
        false,
        6196,
        false,
    )
    .unwrap();

    // Rebuild Lrw over the surviving nodes, from the MI matrix the run used.
    let name_to_global: std::collections::HashMap<&str, usize> = data
        .names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let keep: Vec<usize> = sp
        .node_names
        .iter()
        .map(|n| name_to_global[n.as_str()])
        .collect();
    let q = keep.len();
    let degree: Vec<f64> = keep
        .iter()
        .map(|&i| sp.mim[i].iter().sum::<f64>())
        .collect();
    let lrw = |a: usize, b: usize| {
        if a == b {
            1.0
        } else {
            -sp.mim[keep[a]][keep[b]] / degree[a]
        }
    };

    for c in 0..sp.k {
        let v: Vec<f64> = (0..q).map(|r| sp.eigenvectors[r][c]).collect();
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(norm > 1e-12, "eigenvector {c} is numerically zero");
        let lambda = sp.eigenvalues[c];
        // ||Lrw v - lambda v||_inf, relative to the vector norm.
        let mut worst = 0.0f64;
        for a in 0..q {
            let row: f64 = (0..q).map(|b| lrw(a, b) * v[b]).sum();
            worst = worst.max((row - lambda * v[a]).abs());
        }
        assert!(
            worst / norm < 1e-6,
            "eigenvector {c} (lambda = {lambda}) is not an Lrw eigenvector: \
             relative residual {:.3e}",
            worst / norm
        );
    }
}

/// The Lanczos fast path (`--iterative-eigen`) and the dense path must select the
/// same spectral subsets for a fixed `k`.
#[test]
fn lanczos_and_dense_paths_agree() {
    let data = synthetic();
    let eigen_perc = Some(0.25);

    let mut dense = params(false);
    dense.eigen_perc = eigen_perc;
    let mut iter = dense.clone();
    iter.iterative_eigen = true;

    let pd = scalenet::prepare(&data, &dense).unwrap();
    let pi = scalenet::prepare(&data, &iter).unwrap();

    assert_eq!(pd.spectral.k, pi.spectral.k, "k differs between solvers");
    for c in 0..pd.spectral.k {
        assert!(
            (pd.spectral.eigenvalues[c] - pi.spectral.eigenvalues[c]).abs() < 1e-8,
            "eigenvalue {c} differs between solvers: {} vs {}",
            pd.spectral.eigenvalues[c],
            pi.spectral.eigenvalues[c]
        );
    }
    assert_eq!(
        subset_names(&data, &pi.subgraphs),
        subset_names(&data, &pd.subgraphs),
        "Lanczos and dense paths selected different spectral subsets"
    );
}

/// Documents the bug `--legacy-lrw` deliberately preserves. If this ever starts
/// failing, the legacy path has stopped reproducing pre-fix results and the
/// published analyses can no longer be regenerated.
#[test]
fn legacy_path_stays_column_order_dependent() {
    let data = synthetic();
    let base = edge_names(&data, Method::Aracne, true);
    let mut any_differs = false;
    for seed in [1u64, 2, 3, 12345, 99991] {
        let perm = permutation(data.ncol(), seed);
        let pd = permute(&data, &perm);
        if edge_names(&pd, Method::Aracne, true) != base {
            any_differs = true;
        }
    }
    assert!(
        any_differs,
        "the legacy path is no longer order-dependent — it is supposed to \
         reproduce the pre-fix behaviour exactly"
    );
}
