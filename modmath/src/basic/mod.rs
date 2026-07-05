//! Schoolbook modular arithmetic for `Copy`-bound types.
//!
//! Use this module when the operand type satisfies `T: Copy + ...` —
//! primitive integers (`u8`, `u16`, `u32`, `u64`, `u128`) and
//! stack-allocated bigints that impl `Copy`.
//!
//! For backends that don't impl `Copy`, see
//! [`modmath::constrained`](crate::constrained) (`Clone`-bound,
//! mixed value/reference parameters) or [`modmath::strict`](crate::strict)
//! (reference-based throughout).
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
pub use crate::add::basic_mod_add as add;
#[doc(inline)]
pub use crate::exp::basic_mod_exp as exp;
#[doc(inline)]
pub use crate::inv::basic_mod_inv as inv;
#[doc(inline)]
pub use crate::mul::basic_mod_mul as mul;
#[doc(inline)]
pub use crate::sub::basic_mod_sub as sub;

pub mod montgomery;

/// Pre-reduced variants. Caller guarantees inputs are in `[0, m)`;
/// no `Rem` bound.
pub mod pre_reduced {
    #[doc(inline)]
    pub use crate::add::basic_mod_add_pr as add;
    #[doc(inline)]
    pub use crate::exp::basic_mod_exp_pr as exp;
    #[doc(inline)]
    pub use crate::mul::basic_mod_mul_pr as mul;
    #[doc(inline)]
    pub use crate::sub::basic_mod_sub_pr as sub;
}

/// Non-zero-modulus variants. Caller provides `m: T::NonZero` (one
/// boundary-time `m.into_nonzero()?` proof); every `% m` reduction
/// inside collapses to `rem_nonzero` with no divide-by-zero panic
/// path when the carrier's `DivNonZero` impl elides the underlying
/// zero-check.
pub mod nonzero {
    #[doc(inline)]
    pub use crate::add::basic_mod_add_nz as add;
    #[doc(inline)]
    pub use crate::exp::basic_mod_exp_nz as exp;
    #[doc(inline)]
    pub use crate::mul::basic_mod_mul_nz as mul;
    #[doc(inline)]
    pub use crate::sub::basic_mod_sub_nz as sub;
}
