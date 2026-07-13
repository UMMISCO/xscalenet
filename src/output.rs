//! Writers for the scaleNet / consensus edge-list files (same byte layout as the
//! R package). They format the typed [`Network`] / [`RankedEdge`] domain model.

use std::collections::HashMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::analysis::NetworkAnalysis;
use crate::enrichment::EnrichmentTerm;
use crate::metrics::GraphMetrics;
use crate::multilayer::ThemeProximity;
use crate::network::{Network, Orientation, RankedEdge};
use crate::util::{fmt, pair_index};

fn json_str(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn json_num(x: Option<f64>) -> String {
    match x {
        Some(v) if v.is_finite() => format!("{v}"),
        _ => "null".to_string(),
    }
}

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
        "node\tdegree\tstrength\tbetweenness\tcloseness\tclustering\tcommunity\tmodular_degree\tmodular_betweenness\tatc\tatc_last"
    )?;
    for (i, node) in nodes.iter().enumerate() {
        let n = &gm.nodes[i];
        let (atc_s, atc_last_s) = match atc {
            Some((a, al)) => (fmt(a[i]), fmt(al[i])),
            None => ("NA".to_string(), "NA".to_string()),
        };
        writeln!(
            f,
            "{node}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{atc_s}\t{atc_last_s}",
            n.degree,
            fmt(n.strength),
            fmt(n.betweenness),
            fmt(n.closeness),
            fmt(n.clustering),
            n.community,
            n.modular_degree,
            fmt(n.modular_betweenness),
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

/// Write a self-contained JSON graph (nodes + metrics + annotation, edges, and
/// the optional theme-proximity layer) for a web frontend / IntegraNet.
pub fn write_json(
    path: &Path,
    nodes: &[String],
    net: &Network,
    an: &NetworkAnalysis,
) -> std::io::Result<()> {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir)?;
    }
    let mut f = BufWriter::new(fs::File::create(path)?);
    writeln!(f, "{{")?;

    // Nodes.
    writeln!(f, "  \"nodes\": [")?;
    for (i, name) in nodes.iter().enumerate() {
        let m = &an.metrics.nodes[i];
        let annot = an
            .labels
            .get(i)
            .and_then(|l| *l)
            .and_then(|t| an.themes.get(t))
            .map(|s| json_str(s))
            .unwrap_or_else(|| "null".to_string());
        let atc = an.atc.as_ref().map(|a| a[i]);
        let atc_last = an.atc_last.as_ref().map(|a| a[i]);
        let comma = if i + 1 < nodes.len() { "," } else { "" };
        writeln!(
            f,
            "    {{\"id\": {}, \"community\": {}, \"annotation\": {}, \"degree\": {}, \"strength\": {}, \"betweenness\": {}, \"closeness\": {}, \"clustering\": {}, \"modular_degree\": {}, \"modular_betweenness\": {}, \"atc\": {}, \"atc_last\": {}}}{}",
            json_str(name),
            m.community,
            annot,
            m.degree,
            json_num(Some(m.strength)),
            json_num(Some(m.betweenness)),
            json_num(Some(m.closeness)),
            json_num(Some(m.clustering)),
            m.modular_degree,
            json_num(Some(m.modular_betweenness)),
            json_num(atc),
            json_num(atc_last),
            comma,
        )?;
    }
    writeln!(f, "  ],")?;

    // Edges.
    writeln!(f, "  \"edges\": [")?;
    for (k, e) in net.edges.iter().enumerate() {
        let dir = match e.orientation {
            Orientation::Forward => "forward",
            Orientation::Backward => "backward",
            Orientation::Undirected => "undirected",
        };
        let comma = if k + 1 < net.edges.len() { "," } else { "" };
        writeln!(
            f,
            "    {{\"source\": {}, \"target\": {}, \"orientation\": {}, \"score\": {}, \"ecorr\": {}}}{}",
            json_str(&nodes[e.source]),
            json_str(&nodes[e.target]),
            json_str(dir),
            json_num(e.presence_score),
            json_num(e.ecorr),
            comma,
        )?;
    }
    write!(f, "  ]")?;

    // Optional theme-proximity layer.
    if let Some(tp) = &an.theme_network {
        writeln!(f, ",")?;
        writeln!(f, "  \"theme_network\": {{")?;
        let themes = tp
            .themes
            .iter()
            .map(|t| json_str(t))
            .collect::<Vec<_>>()
            .join(", ");
        writeln!(f, "    \"themes\": [{themes}],")?;
        writeln!(f, "    \"edges\": [")?;
        let m = tp.themes.len();
        let mut first = true;
        for a in 0..m {
            for b in (a + 1)..m {
                if !first {
                    writeln!(f, ",")?;
                }
                first = false;
                write!(
                    f,
                    "      {{\"source\": {}, \"target\": {}, \"proximity\": {}}}",
                    json_str(&tp.themes[a]),
                    json_str(&tp.themes[b]),
                    json_num(Some(tp.proximity[a][b])),
                )?;
            }
        }
        writeln!(f, "\n    ]")?;
        write!(f, "  }}")?;
    }
    writeln!(f, "\n}}")?;
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
