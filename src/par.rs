//! Sequential fallback for the small subset of the Rayon API this crate uses.
//!
//! Active only when the `parallel` feature is **off** — e.g. the `wasm32`
//! single-thread build, where spawning threads is unavailable. Native and
//! default builds use real Rayon (via `rayon::prelude`) and are unaffected, so
//! the numerics stay byte-identical; only the *scheduling* changes here.
#![cfg(not(feature = "parallel"))]
// This module mirrors a slice of the Rayon API so the rest of the crate
// compiles unchanged with `parallel` off. Which of those methods a given build
// actually calls depends on the feature set, so an unused one here is the shim
// doing its job, not dead code to delete.
#![allow(dead_code)]

/// `into_par_iter()` → sequential `into_iter()` (ranges, `Vec`, …).
pub trait IntoParIter: IntoIterator + Sized {
    fn into_par_iter(self) -> Self::IntoIter {
        self.into_iter()
    }
}
impl<T: IntoIterator> IntoParIter for T {}

/// `par_iter()` / `par_iter_mut()` on slices & `Vec` → `iter()` / `iter_mut()`.
pub trait ParSlice {
    type Item;
    fn par_iter(&self) -> core::slice::Iter<'_, Self::Item>;
    fn par_iter_mut(&mut self) -> core::slice::IterMut<'_, Self::Item>;
}
impl<T> ParSlice for [T] {
    type Item = T;
    fn par_iter(&self) -> core::slice::Iter<'_, T> {
        self.iter()
    }
    fn par_iter_mut(&mut self) -> core::slice::IterMut<'_, T> {
        self.iter_mut()
    }
}
impl<T> ParSlice for Vec<T> {
    type Item = T;
    fn par_iter(&self) -> core::slice::Iter<'_, T> {
        self.iter()
    }
    fn par_iter_mut(&mut self) -> core::slice::IterMut<'_, T> {
        self.iter_mut()
    }
}

/// `par_iter()` / `par_iter_mut()` on `HashMap` → `iter()` / `iter_mut()`.
pub trait ParMap {
    type K;
    type V;
    fn par_iter(&self) -> std::collections::hash_map::Iter<'_, Self::K, Self::V>;
    fn par_iter_mut(&mut self) -> std::collections::hash_map::IterMut<'_, Self::K, Self::V>;
}
impl<K, V> ParMap for std::collections::HashMap<K, V> {
    type K = K;
    type V = V;
    fn par_iter(&self) -> std::collections::hash_map::Iter<'_, K, V> {
        self.iter()
    }
    fn par_iter_mut(&mut self) -> std::collections::hash_map::IterMut<'_, K, V> {
        self.iter_mut()
    }
}

/// `map_init(init, f)` → build the per-worker state once, thread it through a
/// plain sequential map (there is a single "worker").
pub struct MapInit<I, S, F> {
    it: I,
    state: S,
    f: F,
}
impl<I: Iterator, S, B, F: FnMut(&mut S, I::Item) -> B> Iterator for MapInit<I, S, F> {
    type Item = B;
    fn next(&mut self) -> Option<B> {
        let f = &mut self.f;
        let st = &mut self.state;
        self.it.next().map(|x| f(st, x))
    }
}
pub trait ParSeqExt: Iterator + Sized {
    fn map_init<S, B, INIT, F>(self, mut init: INIT, f: F) -> MapInit<Self, S, F>
    where
        INIT: FnMut() -> S,
        F: FnMut(&mut S, Self::Item) -> B,
    {
        MapInit {
            state: init(),
            it: self,
            f,
        }
    }
}
impl<I: Iterator> ParSeqExt for I {}

/// Mirrors `rayon::prelude::*` for the subset this crate uses.
pub mod prelude {
    pub use super::{IntoParIter, ParMap, ParSeqExt, ParSlice};
}
