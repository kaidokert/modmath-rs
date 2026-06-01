//! Schoolbook modular arithmetic with reference-based parameters
//! throughout — neither `Copy` nor `Clone` on owned `T` is required.
//!
//! Uses `OverflowingAdd` / `OverflowingSub` rather than
//! `WrappingAdd` / `WrappingSub`, the strictest bound profile in the
//! toolbox. Suitable for backends that don't expose the wrapping
//! operations on owned values, or when avoiding `Clone` is necessary
//! for the operand type.
//!
//! See [`modmath::basic`](crate::basic) (`Copy`-bound, simplest) and
//! [`modmath::constrained`](crate::constrained) (`Clone`-bound, mixed
//! value/reference) for less strict surfaces.
//!
//! ## Pre-reduced variants
//!
//! [`pre_reduced`] mirrors this module's surface but assumes inputs
//! are already in `[0, m)`, skipping the `% m` reduction step. Works
//! for backends without `core::ops::Rem` and is the right entry point
//! inside loops where reduction is handled separately.

#[doc(inline)]
pub use crate::add::strict_mod_add as add;
#[doc(inline)]
pub use crate::exp::strict_mod_exp as exp;
#[doc(inline)]
pub use crate::inv::strict_mod_inv as inv;
#[doc(inline)]
pub use crate::mul::strict_mod_mul as mul;
#[doc(inline)]
pub use crate::sub::strict_mod_sub as sub;

pub mod montgomery;

/// Pre-reduced variants. Caller guarantees inputs are in `[0, m)`;
/// no `Rem`-family bound.
pub mod pre_reduced {
    #[doc(inline)]
    pub use crate::add::strict_mod_add_pr as add;
    #[doc(inline)]
    pub use crate::exp::strict_mod_exp_pr as exp;
    #[doc(inline)]
    pub use crate::mul::strict_mod_mul_pr as mul;
    #[doc(inline)]
    pub use crate::sub::strict_mod_sub_pr as sub;
}
