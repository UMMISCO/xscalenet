//! Input helpers shared by the CLI and language bindings: transparent gzip
//! reading, features-in-rows transposition, and node-annotation parsing.

use std::io::Read;

/// Read a file, transparently decompressing gzip (`.gz` extension or gzip magic).
pub fn read_input(path: &str) -> Result<String, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    if path.ends_with(".gz") || bytes.starts_with(&[0x1f, 0x8b]) {
        let mut s = String::new();
        flate2::read::GzDecoder::new(&bytes[..])
            .read_to_string(&mut s)
            .map_err(|e| format!("cannot decompress {path}: {e}"))?;
        Ok(s)
    } else {
        String::from_utf8(bytes).map_err(|e| format!("{path} is not valid UTF-8: {e}"))
    }
}

/// Transpose a features-in-rows matrix (first column = feature id, header =
/// observation ids) into the features-in-columns layout xscalenet expects.
pub fn transpose_tsv(text: &str) -> Result<String, String> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header = lines.next().ok_or("empty input")?;
    let nobs = header.split('\t').count().saturating_sub(1);
    let mut feat_names: Vec<&str> = Vec::new();
    let mut rows: Vec<Vec<&str>> = Vec::new();
    for (r, line) in lines.enumerate() {
        let mut it = line.split('\t');
        let fname = it.next().unwrap_or("");
        let vals: Vec<&str> = it.collect();
        if vals.len() != nobs {
            return Err(format!(
                "feature row {} has {} values, expected {nobs}",
                r + 1,
                vals.len()
            ));
        }
        feat_names.push(fname);
        rows.push(vals);
    }
    let mut out = String::with_capacity(text.len() + feat_names.len());
    out.push_str(&feat_names.join("\t"));
    out.push('\n');
    for j in 0..nobs {
        for (i, row) in rows.iter().enumerate() {
            if i > 0 {
                out.push('\t');
            }
            out.push_str(row[j]);
        }
        out.push('\n');
    }
    Ok(out)
}

/// Parse a node-annotation TSV (header, node id in column 1) and map each network
/// node to a categorical theme index for a chosen column (default: column 2).
/// Returns per-node `Some(theme)`/`None` plus the distinct theme names.
pub fn read_annotation_labels(
    text: &str,
    node_names: &[String],
    col: Option<&str>,
) -> Result<(Vec<Option<usize>>, Vec<String>), String> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header: Vec<&str> = lines
        .next()
        .ok_or("empty annotation file")?
        .split('\t')
        .map(str::trim)
        .collect();
    let col_idx = match col {
        Some(name) => header
            .iter()
            .position(|h| *h == name)
            .ok_or_else(|| format!("annotation column '{name}' not found"))?,
        None => 1.min(header.len().saturating_sub(1)),
    };
    let mut value_of: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').map(str::trim).collect();
        if f.len() > col_idx {
            value_of.insert(f[0].to_string(), f[col_idx].to_string());
        }
    }
    let mut themes: Vec<String> = value_of.values().cloned().collect();
    themes.sort();
    themes.dedup();
    let theme_idx: std::collections::HashMap<&String, usize> =
        themes.iter().enumerate().map(|(i, t)| (t, i)).collect();
    let labels = node_names
        .iter()
        .map(|n| value_of.get(n).and_then(|v| theme_idx.get(v).copied()))
        .collect();
    Ok((labels, themes))
}
