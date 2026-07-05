//! Sealed marker for "this type is not the `Ct` personality" — safe
//! to flow through variable-time algorithms.
//!
//! The non-`_pr` schoolbook surface excludes `Ct` carriers
//! structurally: it bounds on `Rem`/`Div`, which the `Ct` personality
//! doesn't implement. The `_pr` entries drop those bounds (the caller
//! pre-reduces), so nothing would stop a `Ct` operand from entering
//! their variable-time bodies and leaking through loop counts and
//! per-bit branches. [`NonCt`] closes that gap.
//!
//! Distinct from `const_num_traits::Personality` on purpose:
//! `Personality` describes a type's CT shape, while `NonCt` gates
//! algorithm choice on it — and the gate works by *absence* of an
//! impl, since Rust has no negative bounds.

use const_num_traits::{HasPersonality, Nct};

pub(crate) mod sealed {
    pub trait Sealed {}
}

/// Sealed marker: the type is not the `Ct` personality and may flow
/// through modmath's variable-time algorithms.
///
/// Blanket-implemented for every `T: HasPersonality<P = Nct>` — the
/// primitive integers (whose projection impls live in
/// `const-num-traits`) and any backend carrier declaring the `Nct`
/// personality. A `Ct` carrier projects `P = Ct`, misses the blanket
/// impl, and is rejected at compile time by every variable-time entry
/// bounded on `NonCt`. Sealed, so the only way into the gate is the
/// carrier's own personality declaration.
pub trait NonCt: sealed::Sealed {}

impl<T: HasPersonality<P = Nct>> sealed::Sealed for T {}
impl<T: HasPersonality<P = Nct>> NonCt for T {}

#[cfg(test)]
mod tests {
    use super::NonCt;

    /// Compile-check: primitives and Nct-personality carriers satisfy
    /// the gate. The Ct-personality non-impl belongs in a compile-fail
    /// fixture.
    #[test]
    fn nonct_impl_check() {
        use fixed_bigint::FixedUInt;
        fn assert_nonct<T: NonCt>() {}
        assert_nonct::<u8>();
        assert_nonct::<u16>();
        assert_nonct::<u32>();
        assert_nonct::<u64>();
        assert_nonct::<u128>();
        assert_nonct::<usize>();
        assert_nonct::<FixedUInt<u32, 4>>();
        assert_nonct::<FixedUInt<u32, 8>>();
    }
}
