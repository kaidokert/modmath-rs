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
//! A pre-reduced value (already in `[0, m)`) is a residue in a fixed ring;
//! establish it once through [`SchoolbookField`](crate::SchoolbookField) /
//! [`FieldOps`](crate::FieldOps) rather than threading raw `T` past the `% m`
//! boundary — that surface carries the ring width as a type invariant and
//! cannot seed a narrow value.

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
