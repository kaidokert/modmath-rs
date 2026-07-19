#![cfg_attr(not(test), no_std)]
#![cfg_attr(feature = "nightly", feature(const_trait_impl, const_ops, const_cmp))]
#![cfg_attr(feature = "nightly", allow(incomplete_features))]

//! Modular math implemented with traits.
//!
//! Provides modular arithmetic against any type implementing a minimal
//! set of `core::ops::` and `const_num_traits::` traits — primitive integers
//! or any bigint backend.
//!
//! Schoolbook surface lives in three bound-flavor modules:
//!
//! - [`basic`] — requires `Copy`. Operands by value.
//! - [`constrained`] — requires `Clone` instead. Mixed value / reference operands.
//! - [`strict`] — reference-based throughout. [`Overflowing`] instead of
//!   `Wrapping` arithmetic.
//!
//! Montgomery surface: [`Field`] / [`FieldCt`] / [`FieldNct`] are the
//! high-level type with precomputed parameters and lifetime-branded
//! [`Residue`] values. Each flavor's `montgomery` submodule has the
//! free-function building blocks (`compute_params`, `to_mont`, `mul`,
//! `mod_mul`, …); `<flavor>::montgomery::wide::*` exposes the wide-REDC
//! primitives (`redc`, `mul`, plus `ct::*` siblings).
//!
//! Constant-time variants — branchless finalize via [`subtle`] — sit
//! alongside their non-CT siblings ([`FieldCt`], `Field<T, Ct>`, `_ct`
//! function suffixes, `<flavor>::montgomery::wide::ct::*`).
//!
//! Backend-agnostic: any integer type implementing the required traits
//! works — built-in integers or third-party bigints. The `basic`
//! flavor's `Copy` bound rules out heap-allocated backends; those use
//! `constrained` or `strict`.
//!
//! ## Which surface to use (width safety)
//!
//! On a **runtime-width carrier** (a fixed-capacity, runtime-length bignum such
//! as `HeaplessBigInt`), a value's stored width reflects its magnitude, not the
//! ring it lives in — so a modular value must be established at the modulus
//! ("ring") width or width-sensitive ops fire at the wrong bit and truncate.
//! The [`FieldOps`] surface — Montgomery [`Field`]/[`FieldView`] and the
//! schoolbook [`SchoolbookField`] strategy — carries the ring width as a **type
//! invariant** of its [`Residue`], so a computation written against it cannot
//! seed a narrow value. **Prefer it for anything on a runtime-width carrier, and
//! for building your own reducer.**
//!
//! The raw-`T` schoolbook free functions ([`basic`]/[`constrained`]/[`strict`])
//! are the low-level convenience layer: they reduce *their own* operands to the
//! ring width internally, so a single `basic::mul(a, b, m)` is correct — but
//! they operate on ring-membership-untyped `T`, so a **hand-rolled reducer** that
//! seeds accumulators/coefficients from `T::zero()`/`one()` and grows them
//! *outside* these calls re-opens the width-seed hazard. Such reducers belong on
//! [`SchoolbookField`] (identities come from `field.zero()`/`one()`, which are
//! ring-width by construction), not on raw `T`.
//!
//! [`Overflowing`]: https://docs.rs/num-traits/latest/num_traits/ops/overflowing
//! [`subtle`]: https://crates.io/crates/subtle

mod add;
mod exp;
mod mul;
mod nonct;
mod parity;
mod sub;
#[cfg(test)]
mod test_carrier;
mod wide_mul;

mod field;
mod inv;
pub mod montgomery;

// Bound-flavor module hubs. Re-export the per-flavor schoolbook functions
// under short names (`modmath::basic::add` for `basic_mod_add`, etc.) —
// the `*_mod_` prefix becomes redundant once you're inside `modmath::basic`
// or sibling.
pub mod basic;
pub mod constrained;
pub mod strict;

pub use field::{
    Field, FieldCt, FieldFor, FieldNct, FieldOps, FieldView, MontStorage, Residue, ResidueCt,
    ResidueNct, SchoolbookField, SchoolbookResidue,
};
pub use nonct::NonCt;
pub use parity::Parity;
pub use wide_mul::WideMul;

#[cfg(feature = "nightly")]
pub use add::const_mod_add;
#[cfg(feature = "nightly")]
pub use exp::const_mod_exp;
// Flavor-neutral Montgomery primitives at the crate root: `NPrimeMethod`
// (algorithm selector for R>N `*_with_method` computations, used through
// `modmath::{basic,constrained,strict}::montgomery::compute_params_with_method`)
// plus the precompute helpers (`compute_n_prime_newton`, `compute_r_mod_n`,
// `compute_r2_mod_n`). The flavor-keyed wide-REDC
// wrappers live under `modmath::{basic,constrained,strict}::montgomery::wide::*`.
#[rustfmt::skip]
pub use montgomery::{
    NPrimeMethod,
    compute_n_prime_newton,
    compute_r_mod_n,
    compute_r_mod_n_ct,
    compute_r2_mod_n,
    compute_r2_mod_n_ct,
};
pub use montgomery::{CiosMontMul, CiosMontMulCt};

// Re-exported so the `Field::new_odd` / `basic_montgomery_mod_*_odd`
// surface is reachable as `modmath::Odd` without forcing downstreams to
// directly depend on `const-num-traits` for the typestate wrapper alone.
pub use const_num_traits::Odd;

// Same rationale for the divide-by-zero deletion path: the `*_nz`
// surface in `add`/`sub`/`mul`/`exp` is bounded on these capability
// traits; re-exporting them lets downstreams construct the proof
// (`m.into_nonzero()?`) and drive the infallible reductions without an
// extra direct dep on const-num-traits.
pub use const_num_traits::{DivNonZero, HasNonZero};
#[cfg(feature = "nightly")]
pub use mul::const_mod_mul;
#[cfg(feature = "nightly")]
pub use sub::const_mod_sub;

#[cfg(test)]
#[macro_export]
macro_rules! maybe_test {
    (on, $e:expr) => {
        $e
    };
    (off, $e:expr) => {
        ()
    };
}
