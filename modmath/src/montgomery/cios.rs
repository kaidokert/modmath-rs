//! CIOS (Coarsely Integrated Operand Scanning) Montgomery multiplication.
//!
//! Interleaves multiplication and reduction in a single pass, using
//! 2N²+N limb multiplies instead of 3N² for separate wide-mul + REDC.
//!
//! The algorithm operates entirely through the [`MulAccOps`] trait,
//! never touching raw limb arrays.

use fixed_bigint::MulAccOps;
use fixed_bigint::const_numtraits::ConstBorrowingSub;
use num_traits::ops::overflowing::OverflowingAdd;
use num_traits::{One, WrappingMul, Zero};

/// CIOS Montgomery multiplication: `a * b * R⁻¹ mod modulus`.
///
/// Both `a` and `b` must be in Montgomery form and in `[0, modulus)`.
/// `n_prime_0` is the lowest word of `−N⁻¹ mod R`.
///
/// The algorithm fuses multiplication and REDC into a single double-loop,
/// saving ~25 % of limb multiplies compared to separate wide-mul + REDC.
pub fn cios_montgomery_mul<T: MulAccOps + PartialOrd + ConstBorrowingSub>(
    a: &T,
    b: &T,
    modulus: &T,
    n_prime_0: T::Word,
) -> Option<T>
where
    T::Word: num_traits::Zero
        + num_traits::One
        + num_traits::WrappingMul
        + num_traits::ops::overflowing::OverflowingAdd
        + core::ops::Add<Output = T::Word>,
{
    debug_assert!(a < modulus, "CIOS input a must be in [0, modulus)");
    debug_assert!(b < modulus, "CIOS input b must be in [0, modulus)");

    let n = T::word_count();
    let mut acc = T::default();
    let mut acc_hi = <T::Word as Zero>::zero();
    let mut acc_hi2 = <T::Word as Zero>::zero();

    for i in 0..n {
        let ai = a.get_word(i)?;

        // Phase 1: acc += a[i] * b
        let carry = T::mul_acc_row(ai, b, &mut acc, <T::Word as Zero>::zero());
        let (sum, overflow) = acc_hi.overflowing_add(&carry);
        acc_hi = sum;
        if overflow {
            acc_hi2 = acc_hi2 + <T::Word as One>::one();
        }

        // Compute reduction factor: m = acc[0] * n_prime_0 (mod word)
        let m = acc.get_word(0)?.wrapping_mul(&n_prime_0);

        // Phase 2: [acc, acc_hi] = ([acc, acc_hi] + m * modulus) >> word_bits
        let new_overflow = T::mul_acc_shift_row(m, modulus, &mut acc, acc_hi);
        // Safety: acc_hi2 ∈ {0,1} (reset each iteration, incremented at most once)
        // and new_overflow ∈ {0,1} (bool→word from mul_acc_shift_row), so max sum = 2.
        debug_assert!(
            new_overflow == <T::Word as Zero>::zero() || new_overflow == <T::Word as One>::one(),
            "mul_acc_shift_row must return 0 or 1"
        );
        acc_hi = acc_hi2 + new_overflow;
        acc_hi2 = <T::Word as Zero>::zero();
    }

    // Final conditional subtraction
    if acc_hi > <T::Word as Zero>::zero() || acc >= *modulus {
        let (result, _) = <T as ConstBorrowingSub>::borrowing_sub(acc, *modulus, false);
        acc = result;
    }
    Some(acc)
}

/// Convenience trait: types that support CIOS Montgomery multiplication.
///
/// Implemented automatically for any `T: MulAccOps + PartialOrd + ConstBorrowingSub`
/// with the required `Word` bounds.  Higher-level code (e.g. `MontgomeryCtx`)
/// can bound on this trait instead of requiring knowledge of `MulAccOps::Word`
/// or the CIOS call signature.
pub trait CiosMontMul: MulAccOps + PartialOrd + ConstBorrowingSub {
    fn cios_mont_mul(a: &Self, b: &Self, modulus: &Self, n_prime: &Self) -> Option<Self>;
}

impl<T: MulAccOps + PartialOrd + ConstBorrowingSub> CiosMontMul for T
where
    T::Word: num_traits::Zero
        + num_traits::One
        + num_traits::WrappingMul
        + num_traits::ops::overflowing::OverflowingAdd
        + core::ops::Add<Output = T::Word>,
{
    #[inline]
    fn cios_mont_mul(a: &Self, b: &Self, modulus: &Self, n_prime: &Self) -> Option<Self> {
        cios_montgomery_mul(a, b, modulus, n_prime.get_word(0)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::montgomery::basic_mont::{
        compute_n_prime_newton, compute_r_mod_n, compute_r2_mod_n, type_bit_width,
        wide_montgomery_mul, wide_redc,
    };

    /// Compile-time check: CiosMontMul must be usable as a generic bound
    /// without requiring the consumer to constrain `T::Word`.
    fn _assert_generic_bound<T: CiosMontMul>() {}

    /// Verify CIOS matches wide_montgomery_mul for u8 FixedUInt.
    #[test]
    fn test_cios_vs_wide_redc_u8() {
        use fixed_bigint::FixedUInt;
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
                let got = cios_montgomery_mul(&a_m, &b_m, &modulus, n_prime.get_word(0).unwrap())
                    .unwrap();

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
        use fixed_bigint::FixedUInt;
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
                let got = cios_montgomery_mul(&a_m, &b_m, &modulus, n_prime.get_word(0).unwrap())
                    .unwrap();

                assert_eq!(got, expected, "CIOS mismatch for {a_val:#x}*{b_val:#x}");
            }
        }
    }

    /// Verify CIOS round-trip: to_mont → CIOS multiply → from_mont = (a*b) mod N.
    #[test]
    fn test_cios_roundtrip() {
        use fixed_bigint::FixedUInt;
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
        let result_m =
            cios_montgomery_mul(&a_m, &b_m, &modulus, n_prime.get_word(0).unwrap()).unwrap();

        // Convert back
        let result = wide_redc(result_m, U128::from(0u64), modulus, n_prime);

        // Compare with basic_mod_mul
        let expected = crate::mul::basic_mod_mul(a, b, modulus);
        assert_eq!(result, expected);
    }

    /// Test CiosMontMul convenience trait.
    #[test]
    fn test_cios_mont_mul_trait() {
        use fixed_bigint::FixedUInt;
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
        let result = CiosMontMul::cios_mont_mul(&a_m, &b_m, &modulus, &n_prime).unwrap();
        let expected = wide_montgomery_mul(a_m, b_m, modulus, n_prime);
        assert_eq!(result, expected);
    }
}
