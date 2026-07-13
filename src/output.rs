//! Writers for the scaleNet / consensus edge-list files (same byte layout as the
//! R package). They format the typed [`Network`] / [`RankedEdge`] domain model.

use std::collections::HashMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::enrichment::EnrichmentTerm;
use crate::metrics::GraphMetrics;
use crate::multilayer::ThemeProximity;
use crate::network::{Network, RankedEdge};
use crate::util::{fmt, pair_index};

fn fmt_opt(x: Option<f64>) -> String {
    match x {
        Some(v) => fmt(v),
        None => "NA".to_string(),
    }
}

/// Write a scaleNet edge list. By default every node pair is emitted (present
/// edges from the network, absent pairs as NA / 0 on the fly — the full-pair R
/// layout). With `edges_only`, only present edges are written.
pub fn write_edge_list(
    path: &Path,
    nodes: &[String],
    net: &Network,
    edges_only: bool,
) -> std::io::Result<()> {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir)?;
    }
    let p = nodes.len();
    let by_idx: HashMap<usize, &_> = net
        .edges
        .iter()
        .map(|e| (pair_index(e.source, e.target, p), e))
        .collect();

    let mut f = BufWriter::new(fs::File::create(path)?);
    writeln!(f, "x\ty\tepresenceScore\tepresence\teorientScore\teorient\tecorr")?;
    let mut idx = 0usize;
    for a in 0..p {
        let xa = &nodes[a];
        for b in (a + 1)..p {
            let yb = &nodes[b];
            match by_idx.get(&idx) {
                Some(e) => writeln!(
                    f,
                    "{xa}_<<_{yb}\t{xa}\t{yb}\t{}\t1\t{}\t{}\t{}",
                    fmt_opt(e.presence_score),
                    fmt_opt(e.orient_score),
                    fmt(e.orientation.code(net.kind)),
                    fmt_opt(e.ecorr),
                )?,
                None if !edges_only => {
                    writeln!(f, "{xa}_<<_{yb}\t{xa}\t{yb}\tNA\t0\tNA\tNA\tNA")?
                }
                None => {}
            }
            idx += 1;
        }
    }
    f.flush()?;
    Ok(())
}

/// Write a per-node analysis table: topological metrics and, if provided, ATC.
pub fn write_node_analysis(
    path: &Path,
    nodes: &[String],
    gm: &GraphMetrics,
    atc: Option<(&[f64], &[f64])>,
) -> std::io::Result<()> {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir)?;
    }
    let mut f = BufWriter::new(fs::File::create(path)?);
    writeln!(
        f,
        "node\tdegree\tstrength\tbetweenness\tcloseness\tclustering\tcommunity\tatc\tatc_last"
    )?;
    for (i, node) in nodes.iter().enumerate() {
        let n = &gm.nodes[i];
        let (atc_s, atc_last_s) = match atc {
            Some((a, al)) => (fmt(a[i]), fmt(al[i])),
            None => ("NA".to_string(), "NA".to_string()),
        };
        writeln!(
            f,
            "{node}\t{}\t{}\t{}\t{}\t{}\t{}\t{atc_s}\t{atc_last_s}",
            n.degree,
            fmt(n.strength),
            fmt(n.betweenness),
            fmt(n.closeness),
            fmt(n.clustering),
            n.community,
        )?;
    }
    f.flush()?;
    Ok(())
}

/// Write per-module annotation enrichment (one block per community).
pub fn write_enrichment(
    path: &Path,
    per_module: &[(usize, Vec<EnrichmentTerm>)],
) -> std::io::Result<()> {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir)?;
    }
    let mut f = BufWriter::new(fs::File::create(path)?);
    writeln!(f, "community\tterm\tinterest\treference\tpvalue\todds\tqvalue")?;
    for (c, terms) in per_module {
        for t in terms {
            writeln!(
                f,
                "{c}\t{}\t{}\t{}\t{}\t{}\t{}",
                t.term,
                t.interest_count,
                t.reference_count,
                fmt(t.pvalue),
                fmt(t.odds),
                fmt(t.qvalue),
            )?;
        }
    }
    f.flush()?;
    Ok(())
}

/// Write the theme-proximity (second-layer) network as an edge list over theme
/// pairs (upper triangle), sorted by descending proximity.
pub fn write_theme_network(path: &Path, tp: &ThemeProximity) -> std::io::Result<()> {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir)?;
    }
    let m = tp.themes.len();
    let mut rows: Vec<(usize, usize)> = Vec::new();
    for a in 0..m {
        for b in (a + 1)..m {
            rows.push((a, b));
        }
    }
    rows.sort_by(|&(a1, b1), &(a2, b2)| {
        tp.proximity[a2][b2].total_cmp(&tp.proximity[a1][b1])
    });
    let mut f = BufWriter::new(fs::File::create(path)?);
    writeln!(f, "theme1\ttheme2\tproximity\tdistance")?;
    for (a, b) in rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}",
            tp.themes[a],
            tp.themes[b],
            fmt(tp.proximity[a][b]),
            fmt(tp.distance[a][b]),
        )?;
    }
    f.flush()?;
    Ok(())
}

/// Write the ranked consensus "rawAvg" output (best first).
pub fn write_raw_consensus(
    path: &Path,
    nodes: &[String],
    ranked: &[RankedEdge],
) -> std::io::Result<()> {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir)?;
    }
    let mut f = BufWriter::new(fs::File::create(path)?);
    writeln!(f, "x\ty\tavg.rank\tavg.ort\tecorr")?;
    for e in ranked {
        let xa = &nodes[e.source];
        let yb = &nodes[e.target];
        writeln!(
            f,
            "{xa}_<<_{yb}\t{xa}\t{yb}\t{}\t{}\t{}",
            fmt(e.avg_rank),
            fmt(e.avg_ort),
            fmt_opt(e.ecorr),
        )?;
    }
    f.flush()?;
    Ok(())
}
