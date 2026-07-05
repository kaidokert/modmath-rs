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

pub(crate) mod sealed {
    pub trait Sealed {}
}

/// Sealed marker: the type is not the `Ct` personality and may flow
/// through modmath's variable-time algorithms.
///
/// Implemented for the primitive integers and, behind the
/// `fixed-bigint` feature (on by default), for
/// `fixed_bigint::FixedUInt<W, N, Nct>`. The `Ct` personality
/// deliberately has no impl, so variable-time entries bounded on
/// `NonCt` reject `Ct` carriers at compile time. Sealed so a
/// downstream impl can't reopen the gate.
pub trait NonCt: sealed::Sealed {}

macro_rules! impl_nonct_primitive {
    ($($t:ty),*) => {
        $(
            impl sealed::Sealed for $t {}
            impl NonCt for $t {}
        )*
    };
}

impl_nonct_primitive!(u8, u16, u32, u64, u128, usize);

#[cfg(feature = "fixed-bigint")]
mod fixed_bigint_impl {
    use super::{NonCt, sealed};
    use const_num_traits::Nct;
    use fixed_bigint::{FixedUInt, MachineWord};

    impl<W: MachineWord, const N: usize> sealed::Sealed for FixedUInt<W, N, Nct> {}
    impl<W: MachineWord, const N: usize> NonCt for FixedUInt<W, N, Nct> {}
}

#[cfg(test)]
mod tests {
    use super::NonCt;

    /// Compile-check: primitives impl NonCt. FixedUInt<_, _, Nct> impl
    /// lives behind the `fixed-bigint` feature; the Ct-personality
    /// non-impl is documented in the module docstring and belongs in a
    /// compile-fail fixture.
    #[test]
    fn nonct_impl_check_primitives() {
        fn assert_nonct<T: NonCt>() {}
        assert_nonct::<u8>();
        assert_nonct::<u16>();
        assert_nonct::<u32>();
        assert_nonct::<u64>();
        assert_nonct::<u128>();
        assert_nonct::<usize>();
    }

    #[cfg(feature = "fixed-bigint")]
    #[test]
    fn nonct_impl_check_fixed_bigint() {
        use const_num_traits::Nct;
        use fixed_bigint::FixedUInt;
        fn assert_nonct<T: NonCt>() {}
        assert_nonct::<FixedUInt<u32, 4, Nct>>();
        assert_nonct::<FixedUInt<u32, 8, Nct>>();
    }
}
