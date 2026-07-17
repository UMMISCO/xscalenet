#!/usr/bin/env bash
# Build the in-browser reconstruction module (wasm) + wasm-bindgen JS glue.
#
# Output → web/wasm/{xscalenet.js, xscalenet_bg.wasm}. The web frontend imports
# these to reconstruct a network 100% client-side, no server (single-thread: the
# `parallel`/Rayon feature is OFF, the sequential shim in src/par.rs is used).
#
# Requirements:
#   * a rustup nightly toolchain with the wasm32 target:
#       rustup toolchain install nightly
#       rustup target add wasm32-unknown-unknown --toolchain nightly
#   * wasm-bindgen-cli matching the pinned wasm-bindgen crate (=0.2.126):
#       cargo install wasm-bindgen-cli --version 0.2.126 --locked
#
# NOTE: on a Homebrew rust + rustup setup, `cargo`/`rustc` on PATH may be the
# Homebrew ones (older). This script pins RUSTC to the rustup nightly explicitly.
# NOTE: some nightlies ship a rust-lld whose @rpath can't find libLLVM.dylib
# ("Library not loaded: @rpath/libLLVM.dylib"), which SIGABRTs at link. The lib
# IS in the toolchain, just off rust-lld's search path — symlink it in:
#     NT=$(rustup which --toolchain nightly rustc | xargs dirname)/..
#     ln -sf "$NT/lib/libLLVM.dylib" \
#            "$NT/lib/rustlib/aarch64-apple-darwin/lib/libLLVM.dylib"
# Then --release links fine (~0.5 MB wasm). Shrink further with `wasm-opt -Oz`.
set -euo pipefail
cd "$(dirname "$0")/.."

NT="$(rustup which --toolchain nightly rustc 2>/dev/null | xargs dirname || true)"
[ -n "$NT" ] || { echo "nightly toolchain not found; see requirements in this script"; exit 1; }

PROFILE="${1:-debug}"           # debug (reliable) | release (smaller, may hit an lld bug)
FLAGS="--lib --no-default-features --features wasm --target wasm32-unknown-unknown"
[ "$PROFILE" = release ] && FLAGS="$FLAGS --release"

echo "▸ cargo build ($PROFILE) → wasm32"
RUSTC="$NT/rustc" "$NT/cargo" build $FLAGS

echo "▸ wasm-bindgen → web/wasm/"
mkdir -p web/wasm
wasm-bindgen --target web --out-dir web/wasm --out-name xscalenet \
  "target/wasm32-unknown-unknown/$PROFILE/xscalenet.wasm"

# optional size shrink if binaryen is installed
if command -v wasm-opt >/dev/null 2>&1; then
  echo "▸ wasm-opt -Oz"
  wasm-opt -Oz web/wasm/xscalenet_bg.wasm -o web/wasm/xscalenet_bg.wasm
fi
ls -la web/wasm/
echo "✓ done — open web/wasm/test.html over http to verify"
