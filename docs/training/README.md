# xscalenet training

A hands-on course for the xscalenet / IntegraNet stack, in **two tracks**:

| Track | For | Time | Prerequisite |
|---|---|---|---|
| **[Track A — Use it](track-user.md)** | analysts, biologists, anyone who reconstructs and explores networks | ~3 h | a terminal, no Rust |
| **[Track B — Extend it](track-dev.md)** | contributors to the engine | ~4 h | Rust basics, Track A modules A0–A2 |

Every exercise runs on data **bundled in the repo** (`tests/data/pop2mat.txt`, 45
variables × 25 observations), so nobody waits on a download or on client data.
[solutions.md](solutions.md) has the expected answers.

## Before you start — the shared setup (module A0)

Both tracks begin with [A0](track-user.md#a0--setup-15-min). Do it once:

```sh
cargo build --release              # CLI + library
cargo test  --release              # 16 tests, ~1 s — must be green
cargo build --release --features web --bin xscalenet-web
./target/release/xscalenet-web web # → http://localhost:8788/network.html
```

If `cargo test` is green and `curl -s localhost:8788/health` prints `ok`, you are
ready for either track.

## The one-paragraph summary

xscalenet reconstructs a network over `p` variables from `n` observations
**without ever learning one big graph**. It scores every variable pair by mutual
information, embeds that affinity matrix with a random-walk Laplacian, cuts the
spectrum into `2·(k−1)` overlapping variable subsets, reconstructs each subset
independently (ARACNE and/or Bayesian hill-climbing), and keeps the edges that
recur often enough across subsets. Several methods are then fused by rank
consensus (SCS). Downstream, the same crate computes topology metrics, Louvain
modules, Fisher enrichment, the FunNet theme layer, and ATC annotation
propagation — and ships them, plus the network, as one `graph.json` that the web
browser consumes.

## How the course maps to the reference docs

The training is **task-first**; the reference docs are **subject-first**. Each
module links to the reference section it drills:

| Module | Reference |
|---|---|
| A1, A4 | [algorithm.md](../algorithm.md) |
| A5 | [algorithm.md §7](../algorithm.md#7-downstream-analysis-metrics--atc) |
| A6 | [web/README.md](../../web/README.md) |
| A7, B6 | [suite.md](../suite.md) — the `graph.json` spec |
| B0–B3 | [architecture.md](../architecture.md) |
| B4, B7 | [r-mapping.md](../r-mapping.md) |
| B8 | [performance.md](../performance.md) |

## Teaching it

- **Half-day, mixed audience:** A0 → A2 → A5 → A6, then the A6 capstone. The
  browser work is what people remember; leave 45 min for it.
- **Full day, analysts:** all of Track A. Run A4's parameter sweep as a group and
  compare tables — the disagreement is the lesson.
- **Onboarding a contributor:** A0–A2 for context, then all of Track B. B4 ("add a
  method") is the exercise that actually transfers ownership of the code.
- **Ask for a real dataset in advance.** Every module works on `pop2mat`, but
  learners engage far more when module A6 runs on *their* matrix. It needs to be
  a plain TSV/`.gz`, variables in columns (or `--transpose`).

## Known caveats the course states out loud

Trainers should not be caught out by these; each is covered where it bites.

- **Louvain module detection is not reproducible run-to-run.** Reconstruction is
  deterministic, but community count and modularity Q vary between identical runs
  (observed on `pop2mat` @0.3: 3–4 communities, Q 0.329–0.343). See
  [A5](track-user.md#a5--the-functional-layer-30-min) and
  [B7](track-dev.md#b7--tests-and-validation-against-r-30-min).
- **`--json` from the CLI is *not* the rich `graph.json`** the web browser wants.
  See [B6](track-dev.md#b6--the-graphjson-contract-20-min).
- **`bayes_hc` differs from R's `bnlearn`** on ~20 % of edges by design; ARACNE is
  numerically exact. See [r-mapping.md](../r-mapping.md#fidelity-validated-against-the-r-fixtures-on-pop2mat).
