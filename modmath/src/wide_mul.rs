/// Widening multiplication producing a double-width `(lo, hi)` result.
///
/// Blanket-impl'd for any `T: const_num_traits::CarryingMul<Unsigned = T>`.
pub trait WideMul {
    fn wide_mul(&self, rhs: &Self) -> (Self, Self)
    where
        Self: Sized;
}

impl<T> WideMul for T
where
    T: Copy
        + const_num_traits::Zero
        + core::ops::Mul<Output = T>
        + const_num_traits::CarryingMul<Unsigned = T, Output = T>,
{
    fn wide_mul(&self, rhs: &Self) -> (Self, Self) {
        const_num_traits::CarryingMul::carrying_mul(*self, *rhs, T::zero())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fixed_bigint::FixedUInt;

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
