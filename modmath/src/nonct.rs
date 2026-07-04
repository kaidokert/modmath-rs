//! Sealed marker for "this type is not the Ct personality" — i.e.
//! safe to use in variable-time algorithms.
//!
//! Closes the gap in the `_pr` schoolbook family: those entries don't
//! bound on `Rem` / `Div` (their precondition is that the caller
//! pre-reduced), so the missing-`Rem`-impl gate that blocks Ct on the
//! non-`_pr` surface doesn't apply. Without this marker, the
//! variable-time bodies (`while b > zero { if b.is_odd() ... ; b >>=
//! 1; }`) would admit Ct operands and leak through the loop count and
//! per-bit branch.
//!
//! `T: NonCt` is the bound that gates those entries. Impl'd for the
//! primitives (conventionally Nct — variable-time at the hardware
//! level on common cores) and for `FixedUInt<W, N, Nct>`. Deliberately
//! **not** impl'd for `FixedUInt<W, N, Ct>`, so attempts to call
//! variable-time entries with a `Ct`-personality carrier fail to
//! compile.
//!
//! ## Why a modmath-local marker
//!
//! Conceptually adjacent to `const_num_traits::Personality`, but with
//! a different consumer surface: the personality trait describes the
//! type's CT shape; `NonCt` gates downstream algorithm choice on it.
//! Both can coexist — the personality is the *attribute*, `NonCt` is
//! the *gate that consumes it*.
//!
//! Lives in modmath rather than cnt because (a) it's a modmath-
//! specific gate, (b) the alternative (adding a `HasPersonality`
//! projection trait to cnt + primitive impls + FB impls) is three
//! crates' worth of coordination for the same enforcement guarantee,
//! and (c) if a second consumer crate ever wants the same vocabulary,
//! we can promote to cnt then. For now, modmath-local.

pub(crate) mod sealed {
    pub trait Sealed {}
}

/// Sealed marker for "type is not Ct-personality and is therefore
/// allowed to flow through variable-time algorithms in modmath."
///
/// Implemented for primitives (`u8`/`u16`/`u32`/`u64`/`u128`/`usize`)
/// and, with the `fixed-bigint` feature enabled (on by default), for
/// `fixed_bigint::FixedUInt<W, N, Nct>`. Not implemented for
/// `FixedUInt<W, N, Ct>` — that's the whole point.
///
/// Used as a bound on the variable-time `_pr` schoolbook entries in
/// `add.rs`, `sub.rs`, `mul.rs`, `exp.rs` and on the variable-time
/// Montgomery `_pr` family in `basic_mont.rs`. The bound is a sealed
/// trait so downstreams can't bypass it by impl'ing `NonCt` for a Ct
/// carrier.
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

    // Deliberately NO impl for FixedUInt<W, N, Ct>. Calling
    // variable-time algorithms with Ct-personality carriers must be a
    // compile error.
}

#[cfg(test)]
mod tests {
    use super::NonCt;
    use const_num_traits::{Ct, Nct};
    use fixed_bigint::FixedUInt;

    /// Compile-check: primitives and Nct-personality FixedUInt impl NonCt.
    #[test]
    fn nonct_impl_check() {
        fn assert_nonct<T: NonCt>() {}
        assert_nonct::<u8>();
        assert_nonct::<u16>();
        assert_nonct::<u32>();
        assert_nonct::<u64>();
        assert_nonct::<u128>();
        assert_nonct::<usize>();
        assert_nonct::<FixedUInt<u32, 4, Nct>>();
        assert_nonct::<FixedUInt<u32, 8, Nct>>();
        // FixedUInt<_, _, Ct> does NOT impl NonCt — uncommenting the
        // line below must fail to compile.
        // assert_nonct::<FixedUInt<u32, 4, Ct>>();
        let _ = core::marker::PhantomData::<Ct>;
    }
}
