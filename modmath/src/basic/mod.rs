//! Schoolbook modular arithmetic for `Copy`-bound types.
//!
//! Use this module when the operand type satisfies `T: Copy + ...` —
//! primitive integers (`u8`, `u16`, `u32`, `u64`, `u128`) and
//! `Copy`-impl'ing bigints (`fixed_bigint::FixedUInt`,
//! `crypto_bigint::Uint`).
//!
//! For backends that don't impl `Copy`, see
//! [`modmath::constrained`](crate::constrained) (`Clone`-bound,
//! mixed value/reference parameters) or [`modmath::strict`](crate::strict)
//! (reference-based throughout).
//!
//! ## Pre-reduced variants
//!
//! [`pre_reduced`] mirrors this module's surface but assumes inputs
//! are already in `[0, m)`, skipping the `% m` reduction step. Works
//! for backends without `core::ops::Rem` and is the right entry point
//! inside loops where reduction is handled separately.

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
