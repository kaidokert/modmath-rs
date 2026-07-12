#[cfg(feature = "nightly")]
use const_num_traits::{OverflowingAdd, OverflowingSub};

#[cfg(feature = "nightly")]
c0nst::c0nst! {
    /// # Const Modular Subtraction
    /// Const-evaluable version of modular subtraction.
    pub c0nst fn const_mod_sub<T>(a: T, b: T, m: T) -> T
    where
        T: [c0nst] core::cmp::PartialOrd
            + Copy
            + [c0nst] OverflowingAdd<Output = T>
            + [c0nst] OverflowingSub<Output = T>
            + [c0nst] core::ops::Rem<Output = T>,
    {
        let a = a % m;
        let (diff, overflow) = a.overflowing_sub(b % m);
        if overflow {
            m.overflowing_add(diff).0
        } else {
            diff
        }
    }
}

/// # Modular Subtraction (Basic)
/// Simple version that operates on values and copies them. Requires
/// `WrappingAdd` and `WrappingSub` traits to be implemented.
pub fn basic_mod_sub<T>(a: T, b: T, m: T) -> T
where
    T: core::cmp::PartialOrd
        + Copy
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + core::ops::Rem<Output = T>
        + crate::NonCt,
{
    basic_mod_sub_pr(a % m, b % m, m)
}

/// # Modular Subtraction (Basic, proven-non-zero modulus). **Infallible.**
pub fn basic_mod_sub_nz<T>(a: T, b: T, m: T::NonZero) -> T
where
    T: core::cmp::PartialOrd
        + Copy
        + const_num_traits::HasNonZero
        + const_num_traits::DivNonZero<Output = T>
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + crate::NonCt,
{
    let m_raw = T::nonzero_get(m);
    basic_mod_sub_pr(a.rem_nonzero(m), b.rem_nonzero(m), m_raw)
}

/// # Modular Subtraction (Basic, pre-reduced)
/// Precondition: `a < m` and `b < m`. No `Rem` bound.
pub fn basic_mod_sub_pr<T>(a: T, b: T, m: T) -> T
where
    T: core::cmp::PartialOrd
        + Copy
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + crate::NonCt,
{
    let diff = a.wrapping_sub(b);
    if diff > a {
        // If we wrapped around (underflow)
        diff.wrapping_add(m)
    } else {
        diff
    }
}

/// # Modular Subtraction (Constrained)
/// Version that works with references, requires `WrappingAdd` and
/// `WrappingSub` traits to be implemented.
pub fn constrained_mod_sub<T>(a: T, b: &T, m: &T) -> T
where
    T: core::cmp::PartialOrd
        + Clone
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + crate::NonCt,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>,
{
    let a_mod = &a % m;
    let b_mod = b % m;
    constrained_mod_sub_pr(a_mod, &b_mod, m)
}

/// # Modular Subtraction (Constrained, proven-non-zero modulus). **Infallible.**
pub fn constrained_mod_sub_nz<T>(a: T, b: &T, m: T::NonZero) -> T
where
    T: core::cmp::PartialOrd
        + Clone
        + const_num_traits::HasNonZero
        + const_num_traits::DivNonZero<Output = T>
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + crate::NonCt,
{
    let m_raw = T::nonzero_get(m);
    let b_mod = b.clone().rem_nonzero(m);
    constrained_mod_sub_pr(a.rem_nonzero(m), &b_mod, &m_raw)
}

/// # Modular Subtraction (Constrained, pre-reduced)
/// Precondition: `a < *m` and `*b < *m`. No `Rem` family bound.
pub fn constrained_mod_sub_pr<T>(a: T, b: &T, m: &T) -> T
where
    T: core::cmp::PartialOrd
        + Clone
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + crate::NonCt,
{
    let diff = a.clone().wrapping_sub(b.clone());
    if diff > a {
        // If we wrapped around (underflow)
        diff.wrapping_add(m.clone())
    } else {
        diff
    }
}

/// # Modular Subtraction (Strict)
/// Most constrained version that works with references. Requires
/// `OverflowingAdd` and `OverflowingSub` traits to be implemented.
pub fn strict_mod_sub<T>(mut a: T, b: &T, m: &T) -> T
where
    T: core::cmp::PartialOrd
        + Clone
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::ops::overflowing::OverflowingSub<Output = T>
        + crate::NonCt,
    for<'b> T: core::ops::RemAssign<&'b T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>,
{
    a.rem_assign(m);
    let b_mod = b % m;
    strict_mod_sub_pr(a, &b_mod, m)
}

/// # Modular Subtraction (Strict, proven-non-zero modulus). **Infallible.**
pub fn strict_mod_sub_nz<T>(a: T, b: &T, m: T::NonZero) -> T
where
    T: core::cmp::PartialOrd
        + Clone
        + const_num_traits::HasNonZero
        + const_num_traits::DivNonZero<Output = T>
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::ops::overflowing::OverflowingSub<Output = T>
        + crate::NonCt,
{
    let m_raw = T::nonzero_get(m);
    let b_mod = b.clone().rem_nonzero(m);
    strict_mod_sub_pr(a.rem_nonzero(m), &b_mod, &m_raw)
}

/// # Modular Subtraction (Strict, pre-reduced)
/// Precondition: `a < *m` and `*b < *m`. No `Rem` family bound.
pub fn strict_mod_sub_pr<T>(a: T, b: &T, m: &T) -> T
where
    T: core::cmp::PartialOrd
        + Clone
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::ops::overflowing::OverflowingSub<Output = T>
        + crate::NonCt,
{
    let (diff, overflow) = a.overflowing_sub(b.clone());
    if overflow {
        m.clone().overflowing_add(diff).0
    } else {
        diff
    }
}

#[cfg(test)]
macro_rules! select_mod_sub {
    ($mod_sub:path, $t:ty, by_ref) => {
        fn mod_sub(a: $t, b: &$t, m: &$t) -> $t {
            $mod_sub(a, b, m)
        }
    };
    ($mod_sub:path, $t:ty, by_val) => {
        fn mod_sub(a: $t, b: &$t, m: &$t) -> $t {
            $mod_sub(a, *b, *m)
        }
    };
}

#[cfg(test)]
macro_rules! generate_mod_sub_tests {
    ($mod_sub:path, $t:ty, $by_ref:tt) => {
        select_mod_sub!($mod_sub, $t, $by_ref);

        #[test]
        fn test_mod_sub_basic() {
            assert_eq!(mod_sub(5u8, &10u8, &20u8), 15u8); // 5 - 10 = -5 + 20 = 15
            assert_eq!(mod_sub(7u8, &3u8, &11u8), 4u8); // 7 - 3 = 4
            assert_eq!(mod_sub(0u8, &0u8, &10u8), 0u8); // 0 - 0 = 0
        }

        #[test]
        fn test_mod_sub_res_equals_modulus() {
            assert_eq!(mod_sub(10u8, &10u8, &20u8), 0u8); // 10 - 10 = 0
            assert_eq!(mod_sub(5u8, &5u8, &7u8), 0u8); // 5 - 5 = 0
        }

        #[test]
        fn test_mod_sub_res_less_than_modulus() {
            assert_eq!(mod_sub(15u8, &5u8, &20u8), 10u8); // 15 - 5 = 10
            assert_eq!(mod_sub(8u8, &3u8, &10u8), 5u8); // 8 - 3 = 5
        }

        #[test]
        fn test_mod_sub_with_large_numbers() {
            assert_eq!(mod_sub(255u8, &254u8, &100u8), 1u8); // (255 % 100) - (254 % 100) = 55 - 54 = 1
            assert_eq!(mod_sub(200u8, &100u8, &50u8), 0u8); // (200 % 50) - (100 % 50) = 0 - 0 = 0
        }

        #[test]
        fn test_mod_sub_with_zero() {
            assert_eq!(mod_sub(0u8, &5u8, &10u8), 5u8); // 0 - 5 = -5 + 10 = 5
            assert_eq!(mod_sub(5u8, &0u8, &10u8), 5u8); // 5 - 0 = 5
        }

        #[test]
        fn test_mod_sub_with_max_values() {
            assert_eq!(mod_sub(255u8, &255u8, &100u8), 0u8); // (255 % 100) - (255 % 100) = 55 - 55 = 0
            assert_eq!(mod_sub(255u8, &1u8, &100u8), 54u8); // (255 % 100) - (1 % 100) = 55 - 1 = 54
        }

        #[test]
        fn test_mod_sub_modulus_is_one() {
            assert_eq!(mod_sub(10u8, &20u8, &1u8), 0u8); // Everything mod 1 is 0
            assert_eq!(mod_sub(255u8, &255u8, &1u8), 0u8);
        }

        #[test]
        #[should_panic]
        fn test_mod_sub_modulus_is_zero() {
            mod_sub(10u8, &20u8, &0u8); // Should panic - division by zero
        }

        #[test]
        fn test_mod_sub_operands_exceed_modulus() {
            assert_eq!(mod_sub(100u8, &50u8, &30u8), 20u8); // (100 % 30) - (50 % 30) = 10 - 20 + 30 = 20
            assert_eq!(mod_sub(75u8, &80u8, &20u8), 15u8); // (75 % 20) - (80 % 20) = 15 - 0 = 15
        }

        #[test]
        fn test_mod_sub_edge_cases() {
            assert_eq!(mod_sub(u8::MAX, &u8::MAX, &100u8), 0u8);
            assert_eq!(mod_sub(u8::MAX, &0u8, &100u8), 55u8); // 255 % 100 = 55
            assert_eq!(mod_sub(0u8, &u8::MAX, &100u8), 45u8); // -(255 % 100) + 100 = -55 + 100 = 45
        }

        #[test]
        fn test_mod_sub_result_equals_modulus() {
            assert_eq!(mod_sub(25u8, &5u8, &20u8), 0u8); // (25 % 20) - (5 % 20) = 5 - 5 = 0
            assert_eq!(mod_sub(45u8, &25u8, &20u8), 0u8); // (45 % 20) - (25 % 20) = 5 - 5 = 0
        }

        #[test]
        fn test_mod_sub_overflow() {
            assert_eq!(mod_sub(200u8, &100u8, &50u8), 0u8); // (200 % 50) - (100 % 50) = 0
            assert_eq!(mod_sub(255u8, &254u8, &100u8), 1u8); // Overflow in intermediate calculation
        }
    };
}

#[cfg(test)]
mod strict_mod_sub_tests {
    generate_mod_sub_tests!(super::strict_mod_sub, u8, by_ref);
}

#[cfg(test)]
mod constrained_mod_sub_tests {
    generate_mod_sub_tests!(super::constrained_mod_sub, u8, by_ref);
}

#[cfg(test)]
mod basic_mod_sub_tests {
    generate_mod_sub_tests!(super::basic_mod_sub, u8, by_val);
}

#[cfg(test)]
#[cfg(feature = "nightly")]
const _: () = {
    let result = const_mod_sub(5u8, 10u8, 20u8);
    assert!(result == 15u8);
};

#[cfg(test)]
#[cfg(feature = "nightly")]
mod const_mod_sub_tests {
    generate_mod_sub_tests!(super::const_mod_sub, u8, by_val);
}

#[cfg(test)]
macro_rules! sub_test_module {
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
                fn test_mod_sub_basic() {
                    let a = U256::from(10u8);
                    let b = U256::from(5u8);
                    let m = U256::from(20u8);
                    let result = U256::from(5u8);

                    crate::maybe_test!($strict, assert_eq!(super::strict_mod_sub(a, &b, &m), result));
                    let a = U256::from(10u8);
                    crate::maybe_test!($constrained, assert_eq!(super::constrained_mod_sub(a, &b, &m), result));
                    let a = U256::from(10u8);
                    crate::maybe_test!($basic, assert_eq!(super::basic_mod_sub(a, b, m), result));
                }
            }
        }
    };
}

#[cfg(test)]
mod bnum_sub_tests {
    use super::basic_mod_sub;
    use super::constrained_mod_sub;
    use super::strict_mod_sub;

    //     sub_test_module!(
    //         bnum,
    //         bnum::types::U256,
    //         strict: off, // OverflowingAdd + OverflowingSub is not implemented
    //         constrained: on,
    //         basic: on,
    //     );

    sub_test_module!(
        bnum_patched,
        bnum_patched::types::U256,
        strict: on,
        constrained: on,
        basic: on,
    );

    //     sub_test_module!(
    //         crypto_bigint,
    //         crypto_bigint::U256,
    //         strict: off,  // "Missing OverflowingAdd + OverflowingSub" },
    //         constrained: off, // "Rem<'a> is not implemented for U256" },
    //         basic: on,
    //     );

    sub_test_module!(
        crypto_bigint_patched,
        crypto_bigint_patched::U256,
        strict: on,
        constrained: on,
        basic: on,
    );

    //     sub_test_module!(
    //         num_bigint,
    //         num_bigint::BigUint,
    //         type U256 = num_bigint::BigUint;
    //         strict: off, // OverflowingAdd + OverflowingSub is not implemented
    //         constrained: off, // WrappingAdd + WrappingSub is not implemented
    //         basic: off, // Copy cannot be implemented, heap allocation
    //     );

    //     sub_test_module!(
    //         num_bigint_patched,
    //         num_bigint_patched::BigUint,
    //         type U256 = num_bigint_patched::BigUint;
    //         strict: on,
    //         constrained: on,
    //         basic: off, // Copy cannot be implemented, heap allocation
    //     );

    //     sub_test_module!(
    //         ibig,
    //         ibig::UBig,
    //         type U256 = ibig::UBig;
    //         strict: off, // OverflowingAdd + OverflowingSub is not implemented
    //         constrained: off, // WrappingAdd + WrappingSub is not implemented
    //         basic: off, // Copy cannot be implemented, heap allocation
    //     );

    //     sub_test_module!(
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
    sub_test_module!(
        fixed_bigint,
        fixed_bigint::FixedUInt,
        type U256 = fixed_bigint::FixedUInt<u32, 4>;
        strict: on,
        constrained: on,
        basic: on,
    );

    sub_test_module!(
        heapless_bigint,
        fixed_bigint::FixedUInt,
        type U256 = fixed_bigint::HeaplessBigInt<u32, 4>;
        strict: on,
        constrained: on,
        basic: on,
    );
}

#[cfg(test)]
mod fixed_bigint_pr_tests {
    use super::{basic_mod_sub_pr, constrained_mod_sub_pr, strict_mod_sub_pr};
    type U256 = fixed_bigint::FixedUInt<u32, 4>;

    #[test]
    fn test_mod_sub_basic_pr() {
        let m = U256::from(20u8);
        let a = U256::from(10u8);
        let b = U256::from(5u8);
        let expected = U256::from(5u8);

        assert_eq!(strict_mod_sub_pr(a, &b, &m), expected);
        let a = U256::from(10u8);
        assert_eq!(constrained_mod_sub_pr(a, &b, &m), expected);
        let a = U256::from(10u8);
        assert_eq!(basic_mod_sub_pr(a, b, m), expected);
    }

    #[test]
    fn test_mod_sub_underflow_pr() {
        // a < b: wraps via + m
        let m = U256::from(20u8);
        let a = U256::from(3u8);
        let b = U256::from(7u8);
        let expected = U256::from(16u8); // 3 - 7 = -4 ≡ 16 (mod 20)

        assert_eq!(strict_mod_sub_pr(a, &b, &m), expected);
        let a = U256::from(3u8);
        assert_eq!(constrained_mod_sub_pr(a, &b, &m), expected);
        let a = U256::from(3u8);
        assert_eq!(basic_mod_sub_pr(a, b, m), expected);
    }
}

#[cfg(test)]
mod nz_tests {
    use super::*;
    use const_num_traits::HasNonZero;

    #[test]
    fn basic_mod_sub_nz_matches_basic_mod_sub() {
        let m: u32 = 97;
        let m_nz = m.into_nonzero().unwrap();
        for a in [0u32, 1, 5, 96, 200, u32::MAX] {
            for b in [0u32, 1, 50, 96, 300, u32::MAX - 1] {
                assert_eq!(basic_mod_sub_nz(a, b, m_nz), basic_mod_sub(a, b, m));
            }
        }
    }

    #[test]
    fn constrained_mod_sub_nz_matches_constrained_mod_sub() {
        let m: u32 = 97;
        let m_nz = m.into_nonzero().unwrap();
        for a in [0u32, 1, 5, 96, 200, u32::MAX] {
            for b in [0u32, 1, 50, 96, 300, u32::MAX - 1] {
                assert_eq!(
                    constrained_mod_sub_nz(a, &b, m_nz),
                    constrained_mod_sub(a, &b, &m)
                );
            }
        }
    }

    #[test]
    fn strict_mod_sub_nz_matches_strict_mod_sub() {
        let m: u32 = 97;
        let m_nz = m.into_nonzero().unwrap();
        for a in [0u32, 1, 5, 96, 200, u32::MAX] {
            for b in [0u32, 1, 50, 96, 300, u32::MAX - 1] {
                assert_eq!(strict_mod_sub_nz(a, &b, m_nz), strict_mod_sub(a, &b, &m));
            }
        }
    }
}
