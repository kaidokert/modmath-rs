//! CIOS (Coarsely Integrated Operand Scanning) Montgomery multiplication.
//!
//! Interleaves multiplication and reduction in a single pass, using
//! 2N²+N limb multiplies instead of 3N² for separate wide-mul + REDC.
//!
//! The algorithm is generic over [`CiosRowOps`] (in the leaf
//! `modmath-cios` crate), which exposes infallible limb access and the
//! two CIOS row kernels. modmath never touches raw limb arrays.

use const_num_traits::BorrowingSub;
use const_num_traits::ops::ct::CtIsZero;
use const_num_traits::ops::overflowing::OverflowingAdd;
use const_num_traits::{One, WrappingMul, Zero};
use modmath_cios::CiosRowOps;
use subtle::{Choice, ConditionallySelectable};

/// CIOS Montgomery multiplication — variable-time: `a * b * R⁻¹ mod modulus`.
///
/// Both `a` and `b` must be in Montgomery form and in `[0, modulus)`.
/// `n_prime_0` is the lowest word of `−N⁻¹ mod R`.
///
/// The algorithm fuses multiplication and REDC into a single double-loop,
/// saving ~25 % of limb multiplies compared to separate wide-mul + REDC.
///
/// Final reduction is a predicted branch. Timing leaks operand magnitude.
/// Use [`cios_montgomery_mul_ct`] in CT-sensitive paths (private-key
/// operations, secret scalar multiplication). For verify / public-key
/// paths the branched version is significantly faster — see `wide_redc`
/// docs for the cost breakdown.
pub fn cios_montgomery_mul<T>(a: &T, b: &T, modulus: &T, n_prime_0: T::Word) -> T
where
    T: CiosRowOps + PartialOrd + BorrowingSub<Output = T> + core::ops::Sub<Output = T> + Copy,
    T::Word: const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::WrappingMul<Output = T::Word>
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T::Word>
        + core::ops::Add<Output = T::Word>
        + core::ops::Mul<Output = T::Word>,
{
    debug_assert!(a < modulus, "CIOS input a must be in [0, modulus)");
    debug_assert!(b < modulus, "CIOS input b must be in [0, modulus)");

    let n = a.word_count();
    let zero = <T::Word as Zero>::zero();
    let one = <T::Word as One>::one();
    let mut acc = T::default();
    let mut acc_hi = zero;
    let mut acc_hi2 = zero;

    for i in 0..n {
        let ai = a.word(i);

        // Phase 1: acc += a[i] * b
        let carry = T::mul_acc_row(ai, b, &mut acc, zero);
        let (sum, overflow) = acc_hi.overflowing_add(carry);
        acc_hi = sum;
        if overflow {
            acc_hi2 = acc_hi2 + one;
        }

        // Compute reduction factor: m = acc[0] * n_prime_0 (mod word)
        let m = acc.word(0).wrapping_mul(n_prime_0);

        // Phase 2: [acc, acc_hi] = ([acc, acc_hi] + m * modulus) >> word_bits
        let new_overflow = T::mul_acc_shift_row(m, modulus, &mut acc, acc_hi);
        debug_assert!(
            new_overflow == zero || new_overflow == one,
            "mul_acc_shift_row must return 0 or 1"
        );
        acc_hi = acc_hi2 + new_overflow;
        acc_hi2 = zero;
    }

    // Final conditional subtraction
    if acc_hi > zero || acc >= *modulus {
        let (result, _) = <T as BorrowingSub>::borrowing_sub(acc, *modulus, false);
        acc = result;
    }
    acc
}

/// CIOS Montgomery multiplication — constant-time finalize.
///
/// Same algorithm as [`cios_montgomery_mul`] but performs the final
/// conditional subtraction branchlessly via `subtle::ConditionallySelectable`,
/// removing the operand-magnitude side-channel.
///
/// `acc_hi` is a `T::Word`; its `CtIsZero` provides the "high half
/// nonzero" test. The "low half ≥ modulus" test is free — the borrow
/// flag from `borrowing_sub` already encodes `acc < modulus`.
pub fn cios_montgomery_mul_ct<T>(a: &T, b: &T, modulus: &T, n_prime_0: T::Word) -> T
where
    T: CiosRowOps + BorrowingSub<Output = T> + subtle::ConditionallySelectable + Copy,
    T::Word: const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::WrappingMul<Output = T::Word>
        + const_num_traits::CtIsZero
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T::Word>
        + core::ops::Add<Output = T::Word>
        + core::ops::Mul<Output = T::Word>
        + subtle::ConditionallySelectable,
{
    let n = a.word_count();
    let zero = <T::Word as Zero>::zero();
    let one = <T::Word as One>::one();
    let mut acc = T::default();
    let mut acc_hi = zero;
    let mut acc_hi2 = zero;

    for i in 0..n {
        let ai = a.word(i);

        let carry = T::mul_acc_row(ai, b, &mut acc, zero);
        let (sum, overflow) = acc_hi.overflowing_add(carry);
        acc_hi = sum;
        let overflow_word = T::Word::conditional_select(&zero, &one, Choice::from(overflow as u8));
        acc_hi2 = acc_hi2 + overflow_word;

        let m = acc.word(0).wrapping_mul(n_prime_0);

        let new_overflow = T::mul_acc_shift_row(m, modulus, &mut acc, acc_hi);
        debug_assert!(
            new_overflow == zero || new_overflow == one,
            "mul_acc_shift_row must return 0 or 1"
        );
        acc_hi = acc_hi2 + new_overflow;
        acc_hi2 = zero;
    }

    // Final reduction: the borrow flag from BorrowingSub already
    // encodes `acc < modulus` — no separate compare needed.
    let (sub_result, borrow) = <T as BorrowingSub>::borrowing_sub(acc, *modulus, false);
    let acc_hi_nonzero = !acc_hi.ct_is_zero();
    let needs_sub = acc_hi_nonzero | !Choice::from(borrow as u8);
    T::conditional_select(&acc, &sub_result, needs_sub)
}

/// Convenience trait: types that support variable-time CIOS Montgomery multiplication.
///
/// Implemented automatically for any `T: MulAccOps + PartialOrd + BorrowingSub`
/// with the required `Word` bounds. Higher-level code (e.g. `MontgomeryCtx`)
/// can bound on this trait instead of requiring knowledge of `MulAccOps::Word`
/// or the CIOS call signature. For CT-sensitive paths, see [`CiosMontMulCt`].
pub trait CiosMontMul:
    CiosRowOps + PartialOrd + BorrowingSub + core::ops::Sub<Output = Self> + Copy
{
    fn cios_mont_mul(a: &Self, b: &Self, modulus: &Self, n_prime: &Self) -> Self;
}

impl<T> CiosMontMul for T
where
    T: CiosRowOps + PartialOrd + BorrowingSub<Output = T> + core::ops::Sub<Output = T> + Copy,
    T::Word: const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::WrappingMul<Output = T::Word>
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T::Word>
        + core::ops::Add<Output = T::Word>
        + core::ops::Mul<Output = T::Word>,
{
    #[inline]
    fn cios_mont_mul(a: &Self, b: &Self, modulus: &Self, n_prime: &Self) -> Self {
        cios_montgomery_mul(a, b, modulus, n_prime.word(0))
    }
}

/// Constant-time analog of [`CiosMontMul`].
///
/// Implemented automatically for any
/// `T: MulAccOps + BorrowingSub + ConditionallySelectable` with the
/// required `Word` bounds (including `ConstantTimeEq`). Calls
/// [`cios_montgomery_mul_ct`] under the hood.
pub trait CiosMontMulCt:
    CiosRowOps + BorrowingSub<Output = Self> + subtle::ConditionallySelectable + Copy
{
    fn cios_mont_mul_ct(a: &Self, b: &Self, modulus: &Self, n_prime: &Self) -> Self;
}

impl<T> CiosMontMulCt for T
where
    T: CiosRowOps + BorrowingSub<Output = T> + subtle::ConditionallySelectable + Copy,
    T::Word: const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::WrappingMul<Output = T::Word>
        + const_num_traits::CtIsZero
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T::Word>
        + core::ops::Add<Output = T::Word>
        + core::ops::Mul<Output = T::Word>
        + subtle::ConditionallySelectable,
{
    #[inline]
    fn cios_mont_mul_ct(a: &Self, b: &Self, modulus: &Self, n_prime: &Self) -> Self {
        cios_montgomery_mul_ct(a, b, modulus, n_prime.word(0))
    }
}

#[cfg(test)]
mod smoke_tests {
    use super::*;
    use crate::montgomery::{compute_n_prime_newton, type_bit_width};

    /// New CIOS via `CiosRowOps` must match the existing wide-REDC
    /// `wide_montgomery_mul` reference for every `(a, b)` modulo a small
    /// prime. Generated per primitive — validates the trait surface
    /// end-to-end on the built-in single-word impls alone.
    macro_rules! cios_smoke_test {
        ($name:ident, $t:ty) => {
            #[test]
            fn $name() {
                let modulus: $t = 13;
                let n_prime = compute_n_prime_newton(modulus, type_bit_width::<$t>());
                for a in 0..modulus {
                    for b in 0..modulus {
                        let direct = crate::montgomery::basic_mont::wide_montgomery_mul(
                            a, b, modulus, n_prime,
                        );
                        let via_cios = cios_montgomery_mul(&a, &b, &modulus, n_prime);
                        assert_eq!(direct, via_cios, "CIOS mismatch at a={a}, b={b}");
                    }
                }
            }
        };
    }

    cios_smoke_test!(cios_u8_matches_wide_montgomery_mul_small, u8);
    cios_smoke_test!(cios_u16_matches_wide_montgomery_mul_small, u16);
    cios_smoke_test!(cios_u32_matches_wide_montgomery_mul_small, u32);
    cios_smoke_test!(cios_u64_matches_wide_montgomery_mul_small, u64);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::montgomery::basic_mont::{
        compute_n_prime_newton, compute_r_mod_n, compute_r2_mod_n, type_bit_width,
        wide_montgomery_mul, wide_redc,
    };
    use const_num_traits::Ct;
    use fixed_bigint::FixedUInt;

    /// Compile-time check: CiosMontMul must be usable as a generic bound
    /// without requiring the consumer to constrain `T::Word`.
    fn _assert_generic_bound<T: CiosMontMul>() {}

    /// Verify CIOS matches wide_montgomery_mul for u8 FixedUInt.
    #[test]
    fn test_cios_vs_wide_redc_u8() {
        type U16 = FixedUInt<u8, 2>;

        let modulus = U16::from(13u16);
        let w = type_bit_width::<U16>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

        // Test all pairs in [0, 13)
        for a_val in 0u16..13 {
            for b_val in 0u16..13 {
                let a = U16::from(a_val);
                let b = U16::from(b_val);

                // Convert to Montgomery form
                let (lo, hi) = crate::WideMul::wide_mul(&a, &r2_mod_n);
                let a_m = wide_redc(lo, hi, modulus, n_prime);
                let (lo, hi) = crate::WideMul::wide_mul(&b, &r2_mod_n);
                let b_m = wide_redc(lo, hi, modulus, n_prime);

                // Multiply via wide-REDC (reference)
                let expected = wide_montgomery_mul(a_m, b_m, modulus, n_prime);

                // Multiply via CIOS
                let got = cios_montgomery_mul(&a_m, &b_m, &modulus, n_prime.word(0));

                assert_eq!(
                    got, expected,
                    "CIOS mismatch for {a_val}*{b_val} mod 13: got {got:?}, expected {expected:?}"
                );
            }
        }
    }

    /// Verify CIOS with a larger type (u32 words, 4 limbs = 128-bit).
    #[test]
    fn test_cios_u32x4() {
        type U128 = FixedUInt<u32, 4>;

        let modulus = !U128::from(0u64) - U128::from(58u64); // 2^128 - 59 (odd)
        let w = type_bit_width::<U128>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

        let test_vals = [0u64, 1, 2, 42, 0xDEAD_BEEF, 0xCAFE_BABE];
        for &a_val in &test_vals {
            for &b_val in &test_vals {
                let a = U128::from(a_val);
                let b = U128::from(b_val);

                let (lo, hi) = crate::WideMul::wide_mul(&a, &r2_mod_n);
                let a_m = wide_redc(lo, hi, modulus, n_prime);
                let (lo, hi) = crate::WideMul::wide_mul(&b, &r2_mod_n);
                let b_m = wide_redc(lo, hi, modulus, n_prime);

                let expected = wide_montgomery_mul(a_m, b_m, modulus, n_prime);
                let got = cios_montgomery_mul(&a_m, &b_m, &modulus, n_prime.word(0));

                assert_eq!(got, expected, "CIOS mismatch for {a_val:#x}*{b_val:#x}");
            }
        }
    }

    /// Verify CIOS round-trip: to_mont → CIOS multiply → from_mont = (a*b) mod N.
    #[test]
    fn test_cios_roundtrip() {
        type U128 = FixedUInt<u32, 4>;

        let modulus = !U128::from(0u64) - U128::from(58u64);
        let w = type_bit_width::<U128>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

        let a = U128::from(0xDEAD_BEEF_u64);
        let b = U128::from(0xCAFE_BABE_u64);

        // Convert to Montgomery form
        let (lo, hi) = crate::WideMul::wide_mul(&a, &r2_mod_n);
        let a_m = wide_redc(lo, hi, modulus, n_prime);
        let (lo, hi) = crate::WideMul::wide_mul(&b, &r2_mod_n);
        let b_m = wide_redc(lo, hi, modulus, n_prime);

        // CIOS multiply in Montgomery domain
        let result_m = cios_montgomery_mul(&a_m, &b_m, &modulus, n_prime.word(0));

        // Convert back
        let result = wide_redc(result_m, U128::from(0u64), modulus, n_prime);

        // Compare with basic_mod_mul_pr (inputs are tiny vs modulus, so naturally pre-reduced).
        let expected = crate::mul::basic_mod_mul_pr(a, b, modulus);
        assert_eq!(result, expected);
    }

    /// Test CiosMontMul convenience trait.
    #[test]
    fn test_cios_mont_mul_trait() {
        type U128 = FixedUInt<u32, 4>;

        let modulus = !U128::from(0u64) - U128::from(58u64);
        let w = type_bit_width::<U128>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

        let a = U128::from(42u64);
        let b = U128::from(69u64);
        let (lo, hi) = crate::WideMul::wide_mul(&a, &r2_mod_n);
        let a_m = wide_redc(lo, hi, modulus, n_prime);
        let (lo, hi) = crate::WideMul::wide_mul(&b, &r2_mod_n);
        let b_m = wide_redc(lo, hi, modulus, n_prime);

        // Use the trait method
        let result = CiosMontMul::cios_mont_mul(&a_m, &b_m, &modulus, &n_prime);
        let expected = wide_montgomery_mul(a_m, b_m, modulus, n_prime);
        assert_eq!(result, expected);
    }

    /// CT variant must produce the same output as the NCT path for all inputs.
    ///
    /// The precompute (Newton's iteration, r/r² mod N) runs on the Nct side
    /// for speed; the to-Montgomery REDC also runs in Nct. Only the CT CIOS
    /// call itself requires Ct-typed inputs (its `ConditionallySelectable`
    /// bound resolves only there). Convert via `.into()` for the CT call and
    /// `.forget_ct()` to bring the result back for equality.
    #[test]
    fn test_cios_ct_matches_nct_u8() {
        type U16 = FixedUInt<u8, 2>;
        type U16Ct = FixedUInt<u8, 2, Ct>;

        let modulus = U16::from(13u16);
        let w = type_bit_width::<U16>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);
        let modulus_ct: U16Ct = modulus.into();

        for a_val in 0u16..13 {
            for b_val in 0u16..13 {
                let a = U16::from(a_val);
                let b = U16::from(b_val);
                let (lo, hi) = crate::WideMul::wide_mul(&a, &r2_mod_n);
                let a_m = wide_redc(lo, hi, modulus, n_prime);
                let (lo, hi) = crate::WideMul::wide_mul(&b, &r2_mod_n);
                let b_m = wide_redc(lo, hi, modulus, n_prime);

                let n_prime_0 = n_prime.word(0);
                let nct = cios_montgomery_mul(&a_m, &b_m, &modulus, n_prime_0);
                let a_m_ct: U16Ct = a_m.into();
                let b_m_ct: U16Ct = b_m.into();
                let ct = cios_montgomery_mul_ct(&a_m_ct, &b_m_ct, &modulus_ct, n_prime_0);
                assert_eq!(
                    nct,
                    ct.forget_ct(),
                    "CIOS CT mismatch for {a_val}*{b_val} mod 13"
                );
            }
        }
    }

    /// CT variant matches NCT at a larger type (128-bit modulus).
    #[test]
    fn test_cios_ct_matches_nct_u32x4() {
        type U128 = FixedUInt<u32, 4>;
        type U128Ct = FixedUInt<u32, 4, Ct>;

        let modulus = !U128::from(0u64) - U128::from(58u64);
        let w = type_bit_width::<U128>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);
        let modulus_ct: U128Ct = modulus.into();

        let test_vals = [0u64, 1, 2, 42, 0xDEAD_BEEF, 0xCAFE_BABE];
        for &a_val in &test_vals {
            for &b_val in &test_vals {
                let a = U128::from(a_val);
                let b = U128::from(b_val);
                let (lo, hi) = crate::WideMul::wide_mul(&a, &r2_mod_n);
                let a_m = wide_redc(lo, hi, modulus, n_prime);
                let (lo, hi) = crate::WideMul::wide_mul(&b, &r2_mod_n);
                let b_m = wide_redc(lo, hi, modulus, n_prime);

                let n_prime_0 = n_prime.word(0);
                let nct = cios_montgomery_mul(&a_m, &b_m, &modulus, n_prime_0);
                let a_m_ct: U128Ct = a_m.into();
                let b_m_ct: U128Ct = b_m.into();
                let ct = cios_montgomery_mul_ct(&a_m_ct, &b_m_ct, &modulus_ct, n_prime_0);
                assert_eq!(
                    nct,
                    ct.forget_ct(),
                    "CIOS CT mismatch for {a_val:#x}*{b_val:#x}"
                );
            }
        }
    }

    /// Trait CiosMontMulCt is usable as a bound without constraining T::Word.
    fn _assert_ct_trait_bound<T: CiosMontMulCt>() {}
}
