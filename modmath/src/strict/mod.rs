//! Schoolbook modular arithmetic with reference-based parameters
//! throughout — `Clone`, never `Copy`, so heap-allocated bigints
//! qualify.
//!
//! Uses `OverflowingAdd` / `OverflowingSub` rather than
//! `WrappingAdd` / `WrappingSub`, the strictest bound profile in the
//! toolbox. Suitable for backends that don't wrap on overflow —
//! arbitrary-precision types have no modulus-2^k to wrap at, but can
//! report overflow against a fixed working width.
//!
//! See [`modmath::basic`](crate::basic) (`Copy`-bound, simplest) and
//! [`modmath::constrained`](crate::constrained) (`Clone`-bound, mixed
//! value/reference) for less strict surfaces.
//!
//! A pre-reduced value (already in `[0, m)`) is a residue in a fixed ring;
//! establish it once through [`SchoolbookField`](crate::SchoolbookField) /
//! [`FieldOps`](crate::FieldOps) rather than threading raw `T` past the `% m`
//! boundary — that surface carries the ring width as a type invariant and
//! cannot seed a narrow value.

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

/// Non-zero-modulus variants. See [`basic::nonzero`](crate::basic::nonzero).
pub mod nonzero {
    #[doc(inline)]
    pub use crate::add::strict_mod_add_nz as add;
    #[doc(inline)]
    pub use crate::exp::strict_mod_exp_nz as exp;
    #[doc(inline)]
    pub use crate::mul::strict_mod_mul_nz as mul;
    #[doc(inline)]
    pub use crate::sub::strict_mod_sub_nz as sub;
}
