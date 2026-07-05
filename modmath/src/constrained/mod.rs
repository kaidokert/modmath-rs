//! Schoolbook modular arithmetic for `Clone`-bound types with mixed
//! by-value and by-reference parameters.
//!
//! Use this module when the operand type satisfies `T: Clone + ...`
//! but not necessarily `Copy` — heap-allocated bigints like
//! `num_bigint::BigUint`, `ibig::UBig`, and `bnum` types fit here.
//!
//! See [`modmath::basic`](crate::basic) for the simpler `Copy`-bound
//! surface and [`modmath::strict`](crate::strict) for the
//! reference-only surface (works with backends that don't expose
//! `WrappingMul` / `WrappingSub` on owned values).
//!
//! ## Pre-reduced variants
//!
//! [`pre_reduced`] re-exports `add` / `sub` / `mul` / `exp` siblings that
//! assume inputs are already in `[0, m)` and skip the `% m` reduction
//! step. Works for backends without `core::ops::Rem` and is the right
//! entry point inside loops where reduction is handled separately.
//! Note: `inv` has no pre-reduced sibling — modular inverse via the
//! extended Euclidean algorithm inside `inv` performs its own
//! magnitude-dependent loop, with no useful "pre-reduced" shortcut.

#[doc(inline)]
pub use crate::add::constrained_mod_add as add;
#[doc(inline)]
pub use crate::exp::constrained_mod_exp as exp;
#[doc(inline)]
pub use crate::inv::constrained_mod_inv as inv;
#[doc(inline)]
pub use crate::mul::constrained_mod_mul as mul;
#[doc(inline)]
pub use crate::sub::constrained_mod_sub as sub;

pub mod montgomery;

/// Pre-reduced variants. Caller guarantees inputs are in `[0, m)`;
/// no `Rem` bound.
pub mod pre_reduced {
    #[doc(inline)]
    pub use crate::add::constrained_mod_add_pr as add;
    #[doc(inline)]
    pub use crate::exp::constrained_mod_exp_pr as exp;
    #[doc(inline)]
    pub use crate::mul::constrained_mod_mul_pr as mul;
    #[doc(inline)]
    pub use crate::sub::constrained_mod_sub_pr as sub;
}

/// Non-zero-modulus variants. See [`basic::nonzero`](crate::basic::nonzero).
pub mod nonzero {
    #[doc(inline)]
    pub use crate::add::constrained_mod_add_nz as add;
    #[doc(inline)]
    pub use crate::exp::constrained_mod_exp_nz as exp;
    #[doc(inline)]
    pub use crate::mul::constrained_mod_mul_nz as mul;
    #[doc(inline)]
    pub use crate::sub::constrained_mod_sub_nz as sub;
}
