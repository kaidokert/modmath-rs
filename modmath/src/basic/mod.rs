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
//! A pre-reduced value (already in `[0, m)`) is a residue in a fixed ring;
//! establish it once through [`SchoolbookField`](crate::SchoolbookField) /
//! [`FieldOps`](crate::FieldOps) rather than threading raw `T` past the `% m`
//! boundary — that surface carries the ring width as a type invariant and
//! cannot seed a narrow value.

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
