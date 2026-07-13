//! Writers for the scaleNet / consensus edge-list files (same layout as the R
//! package so outputs can be compared directly).

use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::data::Dataset;
use crate::scalenet::EdgeList;
use crate::scs::RawEdge;
use crate::util::fmt;

/// Write a full scaleNet edge list (all node pairs, combn order). Present edges
/// come from the sparse map; absent pairs are emitted as NA / 0 on the fly, so
/// the output is byte-identical to the dense version without materializing it.
pub fn write_edge_list(path: &Path, data: &Dataset, el: &EdgeList) -> std::io::Result<()> {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir)?;
    }
    let mut f = BufWriter::new(fs::File::create(path)?);
    writeln!(f, "x\ty\tepresenceScore\tepresence\teorientScore\teorient\tecorr")?;
    let p = data.ncol();
    let mut idx = 0usize;
    for a in 0..p {
        let xa = &data.names[a];
        for b in (a + 1)..p {
            let yb = &data.names[b];
            match el.get(idx) {
                Some(rec) => writeln!(
                    f,
                    "{xa}_<<_{yb}\t{xa}\t{yb}\t{}\t1\t{}\t{}\t{}",
                    fmt(rec.epresence_score),
                    fmt(rec.eorient_score),
                    fmt(rec.eorient),
                    fmt(rec.ecorr),
                )?,
                None => writeln!(f, "{xa}_<<_{yb}\t{xa}\t{yb}\tNA\t0\tNA\tNA\tNA")?,
            }
            idx += 1;
        }
    }
    f.flush()?;
    Ok(())
}

/// Write the ranked consensus "rawAvg" output.
pub fn write_raw_consensus(
    path: &Path,
    data: &Dataset,
    raw: &[RawEdge],
) -> std::io::Result<()> {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir)?;
    }
    let p = data.ncol();
    let mut f = BufWriter::new(fs::File::create(path)?);
    writeln!(f, "x\ty\tavg.rank\tavg.ort\tecorr")?;
    for e in raw {
        let (a, b) = crate::util::pair_at(e.idx, p);
        let xa = &data.names[a];
        let yb = &data.names[b];
        writeln!(
            f,
            "{xa}_<<_{yb}\t{xa}\t{yb}\t{}\t{}\t{}",
            fmt(e.avg_rank),
            fmt(e.avg_ort),
            fmt(e.ecorr),
        )?;
    }
    f.flush()?;
    Ok(())
}
