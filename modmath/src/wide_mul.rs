/// Widening multiplication producing a double-width result.
///
/// When the `wide-mul` feature is enabled, types implementing
/// `fixed_bigint::patch_num_traits::WideningMul` use an optimized path,
/// while primitive integers use the fallback bit-by-bit algorithm.
///
/// Without the feature, all types use the generic bit-by-bit
/// algorithm using OverflowingAdd and Parity.
pub trait WideMul {
    fn wide_mul(&self, rhs: &Self) -> (Self, Self)
    where
        Self: Sized;
}

/// Bit-by-bit widening multiplication fallback.
///
/// Uses double-and-add algorithm that works with any type supporting
/// the required traits. Uses `is_zero()` instead of comparison to
/// minimize trait bounds.
///
/// # Overflow invariants
///
/// The algorithm maintains a double-width accumulator (result_lo, result_hi) and
/// a double-width shift register (shift_lo, shift_hi). The product of two W-bit
/// values fits in 2W bits, so:
/// - The shift register holds at most `rhs * 2^k` where k is the current bit position,
///   which is bounded by `rhs * 2^(W-1) < 2^(2W)`.
/// - The result accumulator holds at most `lhs * rhs < 2^(2W)`.
///
/// Therefore, high-half additions that overflow would indicate a bug in the algorithm,
/// not in the inputs. We assert this in debug builds.
#[cfg(not(feature = "wide-mul"))]
fn wide_mul_fallback<T>(lhs: T, rhs: T) -> (T, T)
where
    T: num_traits::ops::overflowing::OverflowingAdd
        + crate::parity::Parity
        + core::ops::Shr<usize, Output = T>
        + num_traits::Zero
        + num_traits::One
        + Copy,
{
    let one = T::one();
    let mut result_lo = T::zero();
    let mut result_hi = T::zero();
    let mut shift_lo = rhs;
    let mut shift_hi = T::zero();
    let mut a_rem = lhs;

    while !a_rem.is_zero() {
        if a_rem.is_odd() {
            let (new_lo, carry) = result_lo.overflowing_add(&shift_lo);
            result_lo = new_lo;
            let (new_hi, overflow1) = result_hi.overflowing_add(&shift_hi);
            result_hi = new_hi;
            if carry {
                let (new_hi, overflow2) = result_hi.overflowing_add(&one);
                result_hi = new_hi;
                // High-half overflow means result > 2^(2W), which is impossible
                // for the product of two W-bit values.
                debug_assert!(!overflow2, "wide_mul: result high-half overflow");
            } else {
                debug_assert!(!overflow1, "wide_mul: result high-half overflow");
            }
        }
        // double the shift register
        let (new_shift_lo, carry) = shift_lo.overflowing_add(&shift_lo);
        shift_lo = new_shift_lo;
        let (new_shift_hi, overflow1) = shift_hi.overflowing_add(&shift_hi);
        shift_hi = new_shift_hi;
        if carry {
            let (new_shift_hi, overflow2) = shift_hi.overflowing_add(&one);
            shift_hi = new_shift_hi;
            // Shift register overflow means shift > 2^(2W), impossible for rhs * 2^k
            // where k < W.
            debug_assert!(!overflow2, "wide_mul: shift high-half overflow");
        } else {
            debug_assert!(!overflow1, "wide_mul: shift high-half overflow");
        }
        a_rem = a_rem >> 1;
    }

    (result_lo, result_hi)
}

// --- Without wide-mul: blanket impl using fallback for all types ---
#[cfg(not(feature = "wide-mul"))]
impl<T> WideMul for T
where
    T: num_traits::ops::overflowing::OverflowingAdd
        + crate::parity::Parity
        + core::ops::Shr<usize, Output = T>
        + num_traits::Zero
        + num_traits::One
        + Copy,
{
    fn wide_mul(&self, rhs: &Self) -> (Self, Self) {
        wide_mul_fallback(*self, *rhs)
    }
}

// --- With wide-mul: optimized impl for WideningMul types ---
#[cfg(feature = "wide-mul")]
impl<T> WideMul for T
where
    T: Copy + fixed_bigint::patch_num_traits::WideningMul<Output = T>,
{
    fn wide_mul(&self, rhs: &Self) -> (Self, Self) {
        fixed_bigint::patch_num_traits::WideningMul::widening_mul(*self, *rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_u8_basics() {
        assert_eq!(0u8.wide_mul(&0), (0, 0));
        assert_eq!(1u8.wide_mul(&1), (1, 0));
        assert_eq!(3u8.wide_mul(&7), (21, 0));
        // 200 * 200 = 40000 = 0x9C40 -> (0x40, 0x9C) = (64, 156)
        assert_eq!(200u8.wide_mul(&200), (64, 156));
        // 255 * 255 = 65025 = 0xFE01 -> (0x01, 0xFE) = (1, 254)
        assert_eq!(255u8.wide_mul(&255), (1, 254));
    }

    #[test]
    fn test_u16_basics() {
        assert_eq!(0u16.wide_mul(&0), (0, 0));
        assert_eq!(1000u16.wide_mul(&1000), (16960, 15));
        assert_eq!(u16::MAX.wide_mul(&u16::MAX), (1, u16::MAX - 1));
    }

    #[test]
    fn test_u32_basics() {
        assert_eq!(0u32.wide_mul(&0), (0, 0));
        assert_eq!(u32::MAX.wide_mul(&u32::MAX), (1, u32::MAX - 1));
        assert_eq!(u32::MAX.wide_mul(&2), (u32::MAX - 1, 1));
    }

    #[test]
    fn test_u64_basics() {
        assert_eq!(0u64.wide_mul(&0), (0, 0));
        assert_eq!(u64::MAX.wide_mul(&u64::MAX), (1, u64::MAX - 1));
        assert_eq!(u64::MAX.wide_mul(&1), (u64::MAX, 0));
    }

    #[test]
    fn test_u8_exhaustive() {
        for a in 0u16..=255 {
            for b in 0u16..=255 {
                let expected = a * b;
                let (lo, hi) = (a as u8).wide_mul(&(b as u8));
                let got = (hi as u16) << 8 | (lo as u16);
                assert_eq!(
                    got, expected,
                    "wide_mul({a}, {b}): expected {expected}, got ({lo}, {hi})"
                );
            }
        }
    }

    #[test]
    fn test_commutativity() {
        let pairs: &[(u32, u32)] = &[
            (0, 0),
            (1, u32::MAX),
            (12345, 67890),
            (u32::MAX, u32::MAX),
            (0x8000_0000, 2),
        ];
        for &(a, b) in pairs {
            assert_eq!(a.wide_mul(&b), b.wide_mul(&a), "commutativity: {a} * {b}");
        }
    }

    #[test]
    fn test_fixed_bigint() {
        use fixed_bigint::FixedUInt;
        type U128 = FixedUInt<u32, 4>;

        let a = U128::from(0xDEAD_BEEF_u64);
        let b = U128::from(0xCAFE_BABE_u64);
        let (lo, hi) = a.wide_mul(&b);

        // Verify against u128 arithmetic
        let a128 = 0xDEAD_BEEF_u128;
        let b128 = 0xCAFE_BABE_u128;
        let expected = a128 * b128;
        let expected_lo = expected as u64;
        let expected_hi = (expected >> 64) as u64;

        assert_eq!(lo, U128::from(expected_lo));
        assert_eq!(hi, U128::from(expected_hi));

        // Build max = U128::MAX via !0
        let max = !U128::from(0u64);
        let one = U128::from(1u64);
        let (lo, hi) = max.wide_mul(&max);
        // MAX * MAX = (MAX-1) * 2^128 + 1
        assert_eq!(lo, one);
        assert_eq!(hi, max - one);

        // Zero and identity
        let zero = U128::from(0u64);
        assert_eq!(a.wide_mul(&zero), (zero, zero));
        assert_eq!(a.wide_mul(&one), (a, zero));

        // Verify via U256: (hi << 128) | lo == a256 * b256
        type U256 = FixedUInt<u32, 8>;

        fn to_u256(lo: FixedUInt<u32, 4>, hi: FixedUInt<u32, 4>) -> FixedUInt<u32, 8> {
            let mut buf = [0u8; 32];
            let mut lo_buf = [0u8; 16];
            let mut hi_buf = [0u8; 16];
            lo.to_le_bytes(&mut lo_buf).unwrap();
            hi.to_le_bytes(&mut hi_buf).unwrap();
            buf[..16].copy_from_slice(&lo_buf);
            buf[16..].copy_from_slice(&hi_buf);
            U256::from_le_bytes(&buf)
        }

        let (lo, hi) = a.wide_mul(&b);
        let a256 = U256::from(0xDEAD_BEEF_u64);
        let b256 = U256::from(0xCAFE_BABE_u64);
        assert_eq!(to_u256(lo, hi), a256 * b256);
    }
}
