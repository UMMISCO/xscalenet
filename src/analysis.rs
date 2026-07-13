//! Functional-layer analysis of a reconstructed network as a single typed call —
//! the contract a web frontend, bindings, or the CLI consume.
//!
//! [`analyze`] runs the topological metrics (incl. modular centrality), and —
//! when a node annotation is supplied — ATC propagation, per-module Fisher
//! enrichment, and the multi-layer theme-proximity network.

use std::collections::BTreeMap;

use crate::atc::{self, AtcConfig};
use crate::enrichment::{self, EnrichmentTerm};
use crate::metrics::{self, Graph, GraphMetrics};
use crate::multilayer::{self, ThemeConfig, ThemeProximity};
use crate::network::Network;

/// What to compute.
#[derive(Clone, Copy, Debug, Default)]
pub struct AnalysisConfig {
    /// Weight the graph by `|ecorr|` (metrics/ATC/theme propagation).
    pub weighted: bool,
    /// Build the multi-layer theme-proximity network (needs an annotation).
    pub theme_network: bool,
}

/// The full functional-layer result for one network.
#[derive(Clone, Debug)]
pub struct NetworkAnalysis {
    pub metrics: GraphMetrics,
    /// ATC scores per node (`None` when no annotation was given).
    pub atc: Option<Vec<f64>>,
    /// Alternative ATC (row sums of the converged config).
    pub atc_last: Option<Vec<f64>>,
    /// Per-community Fisher enrichment (empty when no annotation).
    pub enrichment: Vec<(usize, Vec<EnrichmentTerm>)>,
    /// The theme-proximity (second-layer) network (`None` unless requested).
    pub theme_network: Option<ThemeProximity>,
    /// Annotation theme names, and each node's theme index (for export).
    pub themes: Vec<String>,
    pub labels: Vec<Option<usize>>,
}

/// Run the functional-layer analysis. `annotation` is `(labels, theme_names)`
/// where `labels[i]` is node *i*'s theme index (or `None`).
pub fn analyze(
    network: &Network,
    n_nodes: usize,
    annotation: Option<(&[Option<usize>], &[String])>,
    cfg: &AnalysisConfig,
) -> NetworkAnalysis {
    let g = Graph::from_network(network, n_nodes, cfg.weighted);
    let metrics = metrics::compute(&g);

    let (atc, atc_last) = match annotation {
        Some((labels, themes)) => {
            let t0 = atc::membership_from_labels(labels, themes.len());
            let acfg = AtcConfig {
                weighted: cfg.weighted,
                ..Default::default()
            };
            let r = atc::compute(&g, &t0, &acfg);
            (Some(r.atc), Some(r.atc_last_config))
        }
        None => (None, None),
    };

    let enrichment = match annotation {
        Some((labels, themes)) => {
            let reference: Vec<String> =
                labels.iter().filter_map(|l| l.map(|t| themes[t].clone())).collect();
            let mut by_comm: BTreeMap<usize, Vec<String>> = BTreeMap::new();
            for (i, nm) in metrics.nodes.iter().enumerate() {
                if let Some(t) = labels[i] {
                    by_comm.entry(nm.community).or_default().push(themes[t].clone());
                }
            }
            by_comm
                .into_iter()
                .map(|(c, interest)| (c, enrichment::enrichment_analysis(&interest, &reference)))
                .collect()
        }
        None => Vec::new(),
    };

    let theme_network = if cfg.theme_network {
        annotation.map(|(labels, themes)| {
            let t0 = atc::membership_from_labels(labels, themes.len());
            let tcfg = ThemeConfig {
                weighted: cfg.weighted,
                ..Default::default()
            };
            multilayer::theme_proximity(&g, &t0, themes, &tcfg)
        })
    } else {
        None
    };

    NetworkAnalysis {
        metrics,
        atc,
        atc_last,
        enrichment,
        theme_network,
        themes: annotation.map(|(_, t)| t.to_vec()).unwrap_or_default(),
        labels: annotation.map(|(l, _)| l.to_vec()).unwrap_or_default(),
    }
}
