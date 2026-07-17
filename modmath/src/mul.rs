use crate::parity::Parity;
use const_num_traits::WithPrecision;
#[cfg(feature = "nightly")]
use const_num_traits::{One, OverflowingAdd, OverflowingSub, Zero};

#[cfg(feature = "nightly")]
c0nst::c0nst! {
    /// # Const Modular Multiplication
    /// Const-evaluable version of modular multiplication.
    pub c0nst fn const_mod_mul<T>(a: T, b: T, m: T) -> T
    where
        T: [c0nst] core::cmp::PartialOrd
            + [c0nst] core::cmp::PartialEq
            + Copy
            + [c0nst] Zero
            + [c0nst] One
            + [c0nst] core::ops::BitAnd<Output = T>
            + [c0nst] OverflowingAdd<Output = T>
            + [c0nst] OverflowingSub<Output = T>
            + [c0nst] core::ops::Shr<usize, Output = T>
            + [c0nst] core::ops::Rem<Output = T>,
    {
        let mut a = a % m;
        let mut b = b % m;
        let mut result = T::zero();

        while b > T::zero() {
            if b & T::one() == T::one() {
                let (sum, overflow) = result.overflowing_add(a);
                result = if sum >= m || overflow {
                    sum.overflowing_sub(m).0
                } else {
                    sum
                };
            }

            b = b >> 1;

            if b > T::zero() {
                let (doubled, overflow) = a.overflowing_add(a);
                a = if doubled >= m || overflow {
                    doubled.overflowing_sub(m).0
                } else {
                    doubled
                };
            }
        }

        result
    }
}

/// # Modular Multiplication (Basic)
/// Simple version that operates on values and copies them. Requires
/// `WrappingAdd` and `WrappingSub` traits to be implemented.
pub fn basic_mod_mul<T>(a: T, b: T, m: T) -> T
where
    T: core::cmp::PartialOrd
        + Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + crate::parity::Parity
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Rem<Output = T>
        + crate::NonCt
        + WithPrecision,
{
    basic_mod_mul_pr(a % m, b % m, m)
}

/// # Modular Multiplication (Basic, proven-non-zero modulus). **Infallible.**
///
/// Same body as [`basic_mod_mul`] but the modulus is taken as
/// `T::NonZero` — every `% m` becomes `rem_nonzero(m)` with no
/// divide-by-zero branch. When the carrier's `DivNonZero` impl elides
/// the underlying `Div` zero-check, this **deletes the `panic_fmt`
/// symbol** that the runtime modulus-zero check would otherwise
/// synthesize into the linked binary on embedded targets.
///
/// Sibling `*_nz` siblings in [`basic_mod_add_nz`](crate::add::basic_mod_add_nz),
/// [`basic_mod_sub_nz`](crate::sub::basic_mod_sub_nz), and
/// [`basic_mod_exp_nz`](crate::exp::basic_mod_exp_nz) cover the
/// remaining schoolbook surface; the
/// [`constrained_mod_*_nz`](crate::constrained::nonzero) and
/// [`strict_mod_*_nz`](crate::strict::nonzero) families mirror the same
/// adapter across the by-reference and reference-only bound profiles.
///
/// Routing note: the `_nz` surface is Nct-aligned because the
/// underlying `_pr` bodies are operand-magnitude-branching. Secret
/// (Ct-personality) moduli go through
/// [`Field::try_new_odd_ct`](crate::Field::try_new_odd_ct) and
/// the Montgomery path (which doesn't divide), not through `_nz`.
pub fn basic_mod_mul_nz<T>(a: T, b: T, m: T::NonZero) -> T
where
    T: core::cmp::PartialOrd
        + Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + crate::parity::Parity
        + const_num_traits::HasNonZero
        + const_num_traits::DivNonZero<Output = T>
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + crate::NonCt
        + WithPrecision,
{
    let m_raw = T::nonzero_get(m);
    basic_mod_mul_pr(a.rem_nonzero(m), b.rem_nonzero(m), m_raw)
}

/// # Modular Multiplication (Basic, pre-reduced)
/// Precondition: `a < m` and `b < m`. No `Rem` bound.
///
/// Note: this only removes the division side-channel from the signature.
/// The double-and-add loop still leaks `bit_length(b)`; for constant-time
/// multiplication, use the Montgomery wide-REDC path.
pub(crate) fn basic_mod_mul_pr<T>(a: T, mut b: T, m: T) -> T
where
    T: core::cmp::PartialOrd
        + Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + crate::parity::Parity
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + crate::NonCt
        + WithPrecision,
{
    let m1 = m;
    // Seed the doubling operand and accumulator at the modulus width so
    // `a + a` / `result + a` overflow at bit W, not the narrow operand's.
    let mut a = a.widen_to_precision_of(&m);
    let mut result = T::zero_with_precision_of(&m);

    while b > T::zero() {
        if b.is_odd() {
            // Inline basic_mod_add but skip the redundant modulo on 'a' since we know it's reduced
            let sum = result.wrapping_add(a);
            result = if sum >= m1 || sum < result {
                sum.wrapping_sub(m1)
            } else {
                sum
            };
        }

        b = b >> 1;

        if b > T::zero() {
            let sum = a.wrapping_add(a);
            a = if sum >= m1 || sum < a {
                sum.wrapping_sub(m1)
            } else {
                sum
            };
        }
    }

    result
}

/// # Modular Multiplication (Constrained)
/// Version that works with references, requires `WrappingAdd` and
/// `WrappingSub` traits to be implemented.
pub fn constrained_mod_mul<T>(mut a: T, b: &T, m: &T) -> T
where
    T: Clone
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialOrd
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + crate::NonCt
        + WithPrecision,
    for<'a> T: core::ops::RemAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>,
    for<'a> &'a T: crate::parity::Parity,
{
    a.rem_assign(m);
    let b = b % m;
    constrained_mod_mul_pr(a, &b, m)
}

/// # Modular Multiplication (Constrained, proven-non-zero modulus). **Infallible.**
pub fn constrained_mod_mul_nz<T>(a: T, b: &T, m: T::NonZero) -> T
where
    T: Clone
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialOrd
        + const_num_traits::HasNonZero
        + const_num_traits::DivNonZero<Output = T>
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + crate::NonCt
        + WithPrecision,
    for<'a> &'a T: crate::parity::Parity,
{
    let m_raw = T::nonzero_get(m);
    let b_mod = b.clone().rem_nonzero(m);
    constrained_mod_mul_pr(a.rem_nonzero(m), &b_mod, &m_raw)
}

/// # Modular Multiplication (Constrained, pre-reduced)
/// Precondition: `a < *m` and `*b < *m`. No `Rem` family bound. See
/// [`basic_mod_mul_pr`] for the CT caveat.
pub(crate) fn constrained_mod_mul_pr<T>(a: T, b: &T, m: &T) -> T
where
    T: Clone
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialOrd
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + crate::NonCt
        + WithPrecision,
    for<'a> &'a T: crate::parity::Parity,
{
    let mut b = b.clone();
    // Seed the doubling operand and accumulator at the modulus width so
    // `a + a` / `result + a` overflow at bit W, not the narrow operand's.
    let mut a = a.widen_to_precision_of(m);
    let mut result = T::zero_with_precision_of(m);

    while b > T::zero() {
        if (&b).is_odd() {
            // Inline constrained_mod_add but skip the redundant modulo on 'a' since we know it's reduced
            let sum = result.clone().wrapping_add(a.clone());
            result = if &sum >= m || sum < result {
                sum.wrapping_sub(m.clone())
            } else {
                sum
            };
        }

        b = b >> 1;

        if b > T::zero() {
            let sum = a.clone().wrapping_add(a.clone());
            a = if &sum >= m || sum < a {
                sum.wrapping_sub(m.clone())
            } else {
                sum
            };
        }
    }

    result
}

/// # Modular Multiplication (Strict)
/// Most constrained version that works with references. Requires
/// `OverflowingAdd` and `OverflowingSub` traits to be implemented.
/// # Modular Multiplication (Strict, proven-non-zero modulus). **Infallible.**
pub fn strict_mod_mul_nz<T>(a: T, b: &T, m: T::NonZero) -> T
where
    T: Clone
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialOrd
        + const_num_traits::HasNonZero
        + const_num_traits::DivNonZero<Output = T>
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::ops::overflowing::OverflowingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + crate::NonCt
        + WithPrecision,
    for<'a> &'a T: crate::parity::Parity,
{
    let m_raw = T::nonzero_get(m);
    let b_mod = b.clone().rem_nonzero(m);
    strict_mod_mul_pr(a.rem_nonzero(m), &b_mod, &m_raw)
}

/// # Modular Multiplication (Strict)
/// Most constrained version that works with references. Reduces both
/// operands, then multiplies via repeated modular doubling.
pub fn strict_mod_mul<T>(mut a: T, b: &T, m: &T) -> T
where
    T: Clone
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialOrd
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::ops::overflowing::OverflowingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + crate::NonCt
        + WithPrecision,
    for<'a> T: core::ops::RemAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>,
    for<'a> &'a T: crate::parity::Parity,
{
    a.rem_assign(m);
    let b = b % m;
    strict_mod_mul_pr(a, &b, m)
}

/// # Modular Multiplication (Strict, pre-reduced)
/// Precondition: `a < *m` and `*b < *m`. No `Rem` family bound. See
/// [`basic_mod_mul_pr`] for the CT caveat.
pub(crate) fn strict_mod_mul_pr<T>(a: T, b: &T, m: &T) -> T
where
    T: Clone
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialOrd
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::ops::overflowing::OverflowingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + crate::NonCt
        + WithPrecision,
    for<'a> &'a T: crate::parity::Parity,
{
    let mut b = b.clone();
    // Seed the doubling operand and accumulator at the modulus width so
    // `a + a` / `result + a` overflow at bit W, not the narrow operand's.
    let mut a = a.widen_to_precision_of(m);
    let mut result = T::zero_with_precision_of(m);

    while b > T::zero() {
        if (&b).is_odd() {
            let (sum, overflow) = result.overflowing_add(a.clone());
            result = if &sum >= m || overflow {
                sum.overflowing_sub(m.clone()).0
            } else {
                sum
            };
        }

        b = b >> 1;

        if b > T::zero() {
            let (doubled, overflow) = a.clone().overflowing_add(a);
            a = if &doubled >= m || overflow {
                doubled.overflowing_sub(m.clone()).0
            } else {
                doubled
            };
        }
    }

    result
}

#[cfg(test)]
macro_rules! select_mod_mul {
    ($mod_mul:path, $t:ty, by_ref) => {
        fn mod_mul(a: $t, b: &$t, m: &$t) -> $t {
            $mod_mul(a, b, m)
        }
    };
    ($mod_mul:path, $t:ty, by_val) => {
        fn mod_mul(a: $t, b: &$t, m: &$t) -> $t {
            $mod_mul(a, *b, *m)
        }
    };
}

#[cfg(test)]
macro_rules! generate_mod_mul_tests_block_1 {
    ($mod_add:path, $t:ty, $by_ref:tt) => {
        select_mod_mul!($mod_add, $t, $by_ref);

        #[test]
        fn test_basic_cases() {
            assert_eq!(mod_mul(7u8, &13u8, &19u8), 15); // (7 * 13) % 19 = 15
            assert_eq!(mod_mul(6u8, &9u8, &7u8), 5); // (6 * 9) % 7 = 5
            assert_eq!(mod_mul(5u8, &5u8, &11u8), 3); // (5 * 5) % 11 = 3
        }

        #[test]
        fn test_a_is_zero() {
            assert_eq!(mod_mul(0u8, &5u8, &7u8), 0); // (0 * 5) % 7 = 0
            assert_eq!(mod_mul(0u8, &255u8, &19u8), 0); // (0 * 255) % 19 = 0
        }

        #[test]
        fn test_b_is_zero() {
            assert_eq!(mod_mul(5u8, &0u8, &7u8), 0); // (5 * 0) % 7 = 0
            assert_eq!(mod_mul(255u8, &0u8, &19u8), 0); // (255 * 0) % 19 = 0
        }

        #[test]
        fn test_modulus_is_one() {
            assert_eq!(mod_mul(7u8, &13u8, &1u8), 0); // (7 * 13) % 1 = 0
            assert_eq!(mod_mul(255u8, &255u8, &1u8), 0); // (255 * 255) % 1 = 0
        }

        #[test]
        #[should_panic]
        fn test_modulus_is_zero() {
            mod_mul(7u8, &13u8, &0u8); // Undefined behavior, expect panic or error
        }

        #[test]
        fn test_max_values() {
            assert_eq!(mod_mul(255u8, &255u8, &19u8), 7); // (255 * 255) % 19 = 7
            assert_eq!(mod_mul(255u8, &255u8, &255u8), 0); // (255 * 255) % 255 = 0
        }

        #[test]
        fn test_multiplication_by_one() {
            assert_eq!(mod_mul(1u8, &5u8, &7u8), 5 % 7); // (1 * 5) % 7 = 5
            assert_eq!(mod_mul(7u8, &1u8, &19u8), 7 % 19); // (7 * 1) % 19 = 7
        }

        #[test]
        fn test_equal_values() {
            assert_eq!(mod_mul(7u8, &7u8, &19u8), (7 * 7) % 19); // (7 * 7) % 19 = 11
            assert_eq!(mod_mul(13u8, &13u8, &19u8), (13 * 13) % 19); // (13 * 13) % 19 = 17
        }

        #[test]
        fn test_prime_moduli() {
            assert_eq!(mod_mul(7u8, &13u8, &19u8), (7 * 13) % 19); // (7 * 13) % 19 = 15
            assert_eq!(mod_mul(8u8, &9u8, &17u8), (8 * 9) % 17); // (8 * 9) % 17 = 4
            assert_eq!(mod_mul(5u8, &11u8, &23u8), (5 * 11) % 23); // (5 * 11) % 23 = 9
        }

        #[test]
        fn test_large_values_small_modulus() {
            assert_eq!(mod_mul(200u8, &200u8, &7u8), 2); // (200 * 200) % 7 = 2
            assert_eq!(mod_mul(255u8, &255u8, &3u8), 0); // (255 * 255) % 3 = 0
        }

        #[test]
        fn test_small_modulus() {
            assert_eq!(mod_mul(7u8, &8u8, &2u8), (7 * 8) % 2); // (7 * 8) % 2 = 0
            assert_eq!(mod_mul(5u8, &6u8, &4u8), (5 * 6) % 4); // (5 * 6) % 4 = 2
        }

        #[test]
        fn test_powers_of_two_modulus() {
            assert_eq!(mod_mul(7u8, &13u8, &8u8), 3); // (7 * 13) % 8 = 3
            assert_eq!(mod_mul(16u8, &16u8, &16u8), 0); // (16 * 16) % 16 = 0
        }

        #[test]
        fn test_modulus_greater_than_a_or_b() {
            assert_eq!(mod_mul(10u8, &12u8, &20u8), (10 * 12) % 20); // (10 * 12) % 20 = 0
            assert_eq!(mod_mul(15u8, &14u8, &30u8), (15 * 14) % 30); // (15 * 14) % 30 = 0
        }

        #[test]
        fn test_a_or_b_equals_m_minus_1() {
            assert_eq!(mod_mul(18u8, &13u8, &19u8), (18 * 13) % 19); // (18 * 13) % 19 = 6
            assert_eq!(mod_mul(7u8, &16u8, &17u8), (7 * 16) % 17); // (7 * 16) % 17 = 10
        }

        #[test]
        fn test_binary_modulus() {
            assert_eq!(mod_mul(5u8, &6u8, &2u8), (5 * 6) % 2); // (5 * 6) % 2 = 0
        }

        #[test]
        fn test_small_moduli_explicit() {
            assert_eq!(mod_mul(10u8, &9u8, &2u8), (10 * 9) % 2); // (10 * 9) % 2 = 0
            assert_eq!(mod_mul(10u8, &9u8, &3u8), (10 * 9) % 3); // (10 * 9) % 3 = 0
        }

        #[test]
        fn test_a_and_b_equals_m_minus_1() {
            assert_eq!(mod_mul(18u8, &18u8, &19u8), 1); // (18 * 18) % 19 = 1
            assert_eq!(mod_mul(254u8, &254u8, &255u8), 1); // (254 * 254) % 255 = 1
        }

        #[test]
        fn test_a_or_b_equals_modulus() {
            assert_eq!(mod_mul(7u8, &8u8, &8u8), 0); // (7 * 8) % 8 = 0
            assert_eq!(mod_mul(8u8, &8u8, &8u8), 0); // (8 * 8) % 8 = 0
        }

        #[test]
        fn test_large_product_small_modulus() {
            assert_eq!(mod_mul(250u8, &240u8, &13u8), 5); // (250 * 240) % 13 = 5
            assert_eq!(mod_mul(200u8, &200u8, &5u8), 0); // (200 * 200) % 5 = 0
        }
    };
}

#[cfg(test)]
macro_rules! generate_mod_mul_tests_block_2 {
    ($mod_add:path, $t:ty, $by_ref:tt) => {
        select_mod_mul!($mod_add, $t, $by_ref);

        #[test]
        fn test_64bit_large_values() {
            // Maximum possible u64 values
            assert_eq!(mod_mul(u64::MAX, &u64::MAX, &u64::MAX), 0); // (MAX * MAX) % MAX = 0

            // Large prime modulus and values close to MAX
            assert_eq!(
                mod_mul(u64::MAX, &(u64::MAX - 1), &1_000_000_007_u64),
                532600269
            );

            // Case where a and b multiply to a value larger than 64-bit but mod reduces
            assert_eq!(
                mod_mul(
                    12345678901234567890_u64,
                    &9876543210987654321_u64,
                    &1_000_000_007_u64
                ),
                77470638
            );

            // Edge case: modulus just above a and b
            assert_eq!(mod_mul(10_u64, &20_u64, &u64::MAX), (10 * 20) % u64::MAX);
        }

        #[test]
        fn test_64bit_overflows() {
            // Overflow scenario for smaller types but correct for u64
            assert_eq!(mod_mul(2_u64.pow(32), &2_u64.pow(32), &u64::MAX), 1);
        }

        #[test]
        fn test_64bit_specific_patterns() {
            // Test with powers of two
            assert_eq!(mod_mul(2_u64.pow(63), &2_u64, &u64::MAX), 1);
            // Near max and near zero modulus
            assert_eq!(
                mod_mul(u64::MAX - 1, &1_u64, &u64::MAX),
                (u64::MAX - 1) % u64::MAX
            ); // Just below MAX
            // Prime modulus large values
            let large_prime = 18_446_744_073_709_551_557_u64; // Largest 64-bit prime

            assert_eq!(mod_mul(u64::MAX, &(u64::MAX - 1), &large_prime), 3306_u64);

            // Stress with small modulus and large values
            assert_eq!(mod_mul(u64::MAX - 1, &(u64::MAX - 1), &2_u64), 0);
        }
    };
}

#[cfg(test)]
mod strict_mod_mul_tests {
    use super::strict_mod_mul;
    mod u8_tests {
        use super::strict_mod_mul;
        generate_mod_mul_tests_block_1!(strict_mod_mul, u8, by_ref);
    }
    mod u64_tests {
        use super::strict_mod_mul;
        generate_mod_mul_tests_block_2!(strict_mod_mul, u64, by_ref);
    }
}

#[cfg(test)]
mod constrained_mod_mul_tests {
    use super::constrained_mod_mul;
    mod u8_tests {
        use super::constrained_mod_mul;
        generate_mod_mul_tests_block_1!(constrained_mod_mul, u8, by_ref);
    }
    mod u64_tests {
        use super::constrained_mod_mul;
        generate_mod_mul_tests_block_2!(constrained_mod_mul, u64, by_ref);
    }
}

#[cfg(test)]
mod basic_mod_mul_tests {
    use super::basic_mod_mul;
    mod u8_tests {
        use super::basic_mod_mul;
        generate_mod_mul_tests_block_1!(basic_mod_mul, u8, by_val);
    }
    mod u64_tests {
        use super::basic_mod_mul;
        generate_mod_mul_tests_block_2!(basic_mod_mul, u64, by_val);
    }
}

#[cfg(test)]
#[cfg(feature = "nightly")]
const _: () = {
    let result = const_mod_mul(7u8, 13u8, 19u8);
    assert!(result == 15u8);
};

#[cfg(test)]
#[cfg(feature = "nightly")]
mod const_mod_mul_tests {
    use super::const_mod_mul;
    mod u8_tests {
        use super::const_mod_mul;
        generate_mod_mul_tests_block_1!(const_mod_mul, u8, by_val);
    }
    mod u64_tests {
        use super::const_mod_mul;
        generate_mod_mul_tests_block_2!(const_mod_mul, u64, by_val);
    }
}

#[cfg(test)]
macro_rules! mul_test_module {
    (
        $stem:ident,           // Base name (e.g., "bnum")
        $type_path:path,       // Full path to the type
        $(type $type_def:ty = $type_expr:ty;)? // Optional type definition
        strict: $strict:expr,
        constrained: $constrained:expr,
        basic: $basic:expr,
    ) => {
        paste::paste! {
            mod [<$stem _tests>] {     // This creates e.g., mod bnum_tests
                #[allow(unused_imports)]
                use $type_path;
                $( type $type_def = $type_expr; )?

                #[test]
                #[allow(unused_variables)]
                fn test_mod_mul_basic() {
                    let a_val = 5u8;
                    let a = U256::from(a_val);
                    let b = U256::from(10u8);
                    let m = U256::from(20u8);
                    let result = U256::from(10u8);

                    crate::maybe_test!($strict, assert_eq!(super::strict_mod_mul(a, &b, &m), result));
                    let a = U256::from(a_val);
                    crate::maybe_test!($constrained, assert_eq!(super::constrained_mod_mul(a, &b, &m), result));
                    let a = U256::from(a_val);
                    crate::maybe_test!($basic, assert_eq!(super::basic_mod_mul(a, b, m), result));

                    let a_val = 12345u32;
                    let a: U256 = a_val.try_into().unwrap();
                    let b: U256 = 6789u32.try_into().unwrap();
                    let m: U256 = 10000u32.try_into().unwrap();
                    let result: U256 = 205u32.try_into().unwrap();

                    crate::maybe_test!($strict, assert_eq!(super::strict_mod_mul(a, &b, &m), result));
                    let a: U256 = a_val.try_into().unwrap();
                    crate::maybe_test!($constrained, assert_eq!(super::constrained_mod_mul(a, &b, &m), result));
                    let a: U256 = a_val.try_into().unwrap();
                    crate::maybe_test!($basic, assert_eq!(super::basic_mod_mul(a, b, m), result));
                }
            }
        }
    };
}

#[cfg(test)]
mod bnum_mul_tests {
    use super::basic_mod_mul;
    use super::constrained_mod_mul;
    use super::strict_mod_mul;

    //     mul_test_module!(
    //         bnum,
    //         bnum::types::U256,
    //         strict: off, // OverflowingAdd + OverflowingSub is not implemented
    //         constrained: on,
    //         basic: on,
    //     );

    mul_test_module!(
        bnum_patched,
        bnum_patched::types::U256,
        strict: on,
        constrained: on,
        basic: on,
    );

    //     mul_test_module!(
    //         crypto_bigint,
    //         crypto_bigint::U256,
    //         strict: off,  // "Missing OverflowingAdd + OverflowingSub" },
    //         constrained: off, // "Rem<'a> is not implemented for U256" },
    //         basic: on,
    //     );

    mul_test_module!(
        crypto_bigint_patched,
        crypto_bigint_patched::U256,
        strict: on,
        constrained: on,
        basic: on,
    );

    //     mul_test_module!(
    //         num_bigint,
    //         num_bigint::BigUint,
    //         type U256 = num_bigint::BigUint;
    //         strict: off, // OverflowingAdd + OverflowingSub is not implemented
    //         constrained: off, // WrappingAdd + WrappingSub is not implemented
    //         basic: off, // Copy cannot be implemented, heap allocation
    //     );

    // num-bigint `FixedWidthBigUint`: heap carrier, Nct, constrained/strict
    // only (not `Copy`, so `basic: off`).
    mul_test_module!(
        num_bigint_patched,
        num_bigint_patched::FixedWidthBigUint,
        type U256 = num_bigint_patched::FixedWidthBigUint;
        strict: on,
        constrained: on,
        basic: off,
    );

    //     mul_test_module!(
    //         ibig,
    //         ibig::UBig,
    //         type U256 = ibig::UBig;
    //         strict: off, // OverflowingAdd + OverflowingSub is not implemented
    //         constrained: off, // WrappingAdd + WrappingSub is not implemented
    //         basic: off, // Copy cannot be implemented, heap allocation
    //     );

    //     mul_test_module!(
    //         ibig_patched,
    //         ibig_patched::UBig,
    //         type U256 = ibig_patched::UBig;
    //         strict: on,
    //         constrained: on,
    //         basic: off, // Copy cannot be implemented, heap allocation
    //     );

    // Default (Nct-personality) FixedUInt provides Rem; the wrapper
    // functions are usable here. The Ct personality intentionally omits
    // Rem (variable-time on operand magnitudes).
    mul_test_module!(
        fixed_bigint,
        fixed_bigint::FixedUInt,
        type U256 = fixed_bigint::FixedUInt<u32, 4>;
        strict: on,
        constrained: on,
        basic: on,
    );

    mul_test_module!(
        heapless_bigint,
        fixed_bigint::FixedUInt,
        type U256 = fixed_bigint::HeaplessBigInt<u32, 4>;
        strict: on,
        constrained: on,
        basic: on,
    );

    // Operands occupy one byte, modulus spans two: the peasant loop doubles
    // `a` and must overflow at the modulus width, not the operand's. Regression
    // guard for the len(a) < len(m) case the single-word modules can't reach.
    #[test]
    fn heapless_narrow_operand_wide_modulus() {
        use fixed_bigint::HeaplessBigInt;
        type H = HeaplessBigInt<u8, 4>;
        let m: H = 1000u16.into();
        let want: H = 104u16.into(); // 234 * 56 = 13104 ≡ 104 (mod 1000)
        assert_eq!(super::basic_mod_mul(H::from(234u8), H::from(56u8), m), want);
        assert_eq!(
            super::constrained_mod_mul(H::from(234u8), &H::from(56u8), &m),
            want
        );
        assert_eq!(
            super::strict_mod_mul(H::from(234u8), &H::from(56u8), &m),
            want
        );
    }
}

#[cfg(test)]
mod fixed_bigint_pr_tests {
    use super::{basic_mod_mul_pr, constrained_mod_mul_pr, strict_mod_mul_pr};
    type U256 = fixed_bigint::FixedUInt<u32, 4>;

    #[test]
    fn test_mod_mul_basic_pr() {
        let m = U256::from(19u8);
        let a = U256::from(7u8);
        let b = U256::from(13u8);
        let expected = U256::from(15u8); // (7 * 13) mod 19 = 91 mod 19 = 15

        assert_eq!(strict_mod_mul_pr(a, &b, &m), expected);
        let a = U256::from(7u8);
        assert_eq!(constrained_mod_mul_pr(a, &b, &m), expected);
        let a = U256::from(7u8);
        assert_eq!(basic_mod_mul_pr(a, b, m), expected);
    }

    #[test]
    fn test_mod_mul_zero_pr() {
        let m = U256::from(19u8);
        let zero = U256::from(0u8);
        let b = U256::from(13u8);

        assert_eq!(strict_mod_mul_pr(zero, &b, &m), zero);
        assert_eq!(constrained_mod_mul_pr(zero, &b, &m), zero);
        assert_eq!(basic_mod_mul_pr(zero, b, m), zero);
    }
}

#[cfg(test)]
mod nz_tests {
    use super::*;
    use const_num_traits::HasNonZero;

    #[test]
    fn basic_mod_mul_nz_matches_basic_mod_mul() {
        let m: u32 = 97;
        let m_nz = m.into_nonzero().unwrap();
        for a in [0u32, 1, 5, 96, 200, 0xFFFF] {
            for b in [0u32, 1, 50, 96, 300, 0xABCD] {
                assert_eq!(basic_mod_mul_nz(a, b, m_nz), basic_mod_mul(a, b, m));
            }
        }
    }

    #[test]
    fn constrained_mod_mul_nz_matches_constrained_mod_mul() {
        let m: u32 = 97;
        let m_nz = m.into_nonzero().unwrap();
        for a in [0u32, 1, 5, 96, 200, 0xFFFF] {
            for b in [0u32, 1, 50, 96, 300, 0xABCD] {
                assert_eq!(
                    constrained_mod_mul_nz(a, &b, m_nz),
                    constrained_mod_mul(a, &b, &m)
                );
            }
        }
    }

    #[test]
    fn strict_mod_mul_nz_matches_strict_mod_mul() {
        let m: u32 = 97;
        let m_nz = m.into_nonzero().unwrap();
        for a in [0u32, 1, 5, 96, 200, 0xFFFF] {
            for b in [0u32, 1, 50, 96, 300, 0xABCD] {
                assert_eq!(strict_mod_mul_nz(a, &b, m_nz), strict_mod_mul(a, &b, &m));
            }
        }
    }

    /// End-to-end: the `*_nz` surface resolves against fixed-bigint's
    /// `NonZeroFixedUInt`-driven `DivNonZero` impl. Pins that the FB
    /// carrier-side ship (`HasNonZero` + `DivNonZero` for `FixedUInt`)
    /// is consumable through modmath without per-call-site adapter
    /// boilerplate.
    #[test]
    fn basic_mod_mul_nz_with_fixed_bigint() {
        use fixed_bigint::FixedUInt;
        type U128 = FixedUInt<u32, 4>;

        let m = U128::from(97u32);
        let m_nz = m.into_nonzero().expect("97 != 0");

        let a = U128::from(50u32);
        let b = U128::from(50u32);

        let via_nz = basic_mod_mul_nz(a, b, m_nz);
        let via_old = basic_mod_mul(a, b, m);
        assert_eq!(via_nz, via_old);
        // 50 * 50 = 2500; 2500 % 97 = 2500 - 25*97 = 2500 - 2425 = 75
        assert_eq!(via_nz, U128::from(75u32));
    }
}
