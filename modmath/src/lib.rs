#![cfg_attr(not(test), no_std)]
#![cfg_attr(feature = "nightly", feature(const_trait_impl, const_ops, const_cmp))]
#![cfg_attr(feature = "nightly", allow(incomplete_features))]

//! Modular math implemented with traits.
//!
//! Provides modular arithmetic against any type implementing a minimal
//! set of `core::ops::` and `num_traits::` traits — primitive integers
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
//! Tested against built-in integers, [`num-bigint`], [`crypto-bigint`],
//! [`bnum`], [`ibig`], and [`fixed-bigint`]. The `basic` flavor's `Copy`
//! bound rules out heap-allocated backends (`num-bigint`, `ibig`); those
//! use `constrained` or `strict`.
//!
//! [`Overflowing`]: https://docs.rs/num-traits/latest/num_traits/ops/overflowing
//! [`subtle`]: https://crates.io/crates/subtle
//! [`num-bigint`]: https://crates.io/crates/num-bigint
//! [`crypto-bigint`]: https://crates.io/crates/crypto-bigint
//! [`bnum`]: https://crates.io/crates/bnum
//! [`ibig`]: https://crates.io/crates/ibig
//! [`fixed-bigint`]: https://crates.io/crates/fixed-bigint

mod add;
mod exp;
mod mul;
mod parity;
mod sub;
mod wide_mul;

#[cfg(feature = "wide-mul")]
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

#[cfg(feature = "wide-mul")]
pub use field::{Field, FieldCt, FieldNct, MontStorage, Residue, ResidueCt, ResidueNct};
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
// `compute_r2_mod_n`) and `type_bit_width`. The flavor-keyed wide-REDC
// wrappers live under `modmath::{basic,constrained,strict}::montgomery::wide::*`.
#[rustfmt::skip]
pub use montgomery::{
    NPrimeMethod,
    type_bit_width,
    compute_n_prime_newton,
    compute_r_mod_n,
    compute_r2_mod_n,
};
#[cfg(feature = "wide-mul")]
pub use montgomery::{CiosMontMul, CiosMontMulCt};
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
