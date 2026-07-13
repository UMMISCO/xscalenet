# xscalenet documentation

`xscalenet` is an efficient Rust reimplementation of the
[ScaleNet](https://github.com/UMMISCO/scalenet) R package: reconstruction of
large biological networks from high-dimensional data via a **spectral
framework**, plus the **Spectral Consensus Strategy (SCS)** that fuses several
embedded reconstruction methods into one robust network.

> Affeldt S, Sokolovska N, Prifti E, Zucker J-D. *Spectral consensus strategy
> for accurate reconstruction of large biological networks.* BMC Bioinformatics
> 2016;17(Suppl 16):493. [PMC5249011](https://pmc.ncbi.nlm.nih.gov/articles/PMC5249011/)

## Read in this order

1. **[algorithm.md](algorithm.md)** — the method: what problem it solves, the
   spectral framework, and the consensus strategy, with the formulas as
   implemented.
2. **[architecture.md](architecture.md)** — how the code is organized: the
   end-to-end pipeline, module-by-module responsibilities, data model, and where
   parallelism lives.
3. **[r-mapping.md](r-mapping.md)** — how each Rust module maps to the original R
   package, what is faithfully reproduced, and what is intentionally out of scope.
4. **[performance.md](performance.md)** — the cost model per stage, the known
   bottlenecks, and the optimization direction.

## Quick start

```sh
cargo build --release
./target/release/xscalenet <input.tsv> <output_dir> [options]
```

`input.tsv` is tab-separated with **features in columns, observations in rows**
and a header of feature names. See [architecture.md](architecture.md#cli) for the
full option list and output layout.

## Project status

The core spectral + SCS path is implemented and validated against the R package's
fixtures (ARACNE reconstruction is numerically exact). Active work is tracked in
the modernization roadmap; the current priorities are a scalable discretizer and
a sparse edge model to remove the O(p²) cost inherited from the R design. See
[performance.md](performance.md) for details.
