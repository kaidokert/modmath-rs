pub mod safegcd;
mod signed;
use signed::Signed;

/*
Based on this pseudocode:

function inverse(a, n)
    t := 0;     newt := 1
    r := n;     newr := a

    while newr ≠ 0 do
        quotient := r div newr
        (t, newt) := (newt, t − quotient × newt)
        (r, newr) := (newr, r − quotient × newr)

    if r > 1 then
        return "a is not invertible"
    if t < 0 then
        t := t + n

    return t
*/

/// # Modular Inverse (Strict)
/// Most constrained version that works with references. Requires
/// reference-based operations for division and subtraction.
///
/// # Panics
///
/// Panics if `T` is too narrow to hold the extended-GCD intermediates
/// (checked coefficient arithmetic overflows) — a carrier-precondition
/// violation, distinct from the `None` return for a non-invertible input.
pub fn strict_mod_inv<T>(a: T, modulus: &T) -> Option<T>
where
    // Checked mul/add throughout (see `basic_mod_inv`): both the `Signed`
    // coefficients and the raw `T` r-sequence refuse to rely on `*`'s
    // implementation-defined overflow.
    T: const_num_traits::Zero
        + const_num_traits::One
        + PartialEq
        + const_num_traits::CheckedAdd<Output = T>
        + const_num_traits::CheckedMul<Output = T>
        + core::ops::Sub<Output = T>
        + const_num_traits::WithPrecision
        + core::cmp::PartialOrd,
    for<'a> T: core::ops::Sub<&'a T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::cmp::PartialOrd,
    for<'a> &'a T: core::ops::Div<&'a T, Output = T> + core::ops::Sub<T, Output = T>,
{
    // Seed the signed coefficients at the modulus width: on a runtime-width
    // carrier `T::zero()` is width 0, and the signed magnitude arithmetic
    // (canonicalization, magnitude subtraction) misbehaves when a width-0
    // magnitude meets field-width values. `zero_with_precision_of` establishes
    // the ring width; identity on fixed-width carriers.
    let mut t = Signed::new(T::zero_with_precision_of(modulus), false);
    let mut new_t = Signed::new(T::one_with_precision_of(modulus), false);
    // makes a clone of modulus
    let mut r = T::zero_with_precision_of(modulus) + modulus;
    let mut new_r = a;

    while new_r != T::zero() {
        let quotient = &r / &new_r;
        // strict avoids a `Clone` bound: duplicate `quotient` via the same
        // reference-add trick used for `r`/`t`, so each coefficient update can
        // consume an owned copy through the checked multiply.
        let quotient2 = T::zero_with_precision_of(modulus) + &quotient;

        // clone
        let tmp_t = Signed::new(T::zero_with_precision_of(modulus), false) + &new_t;
        new_t = t - new_t * quotient;
        t = tmp_t;

        // clone
        let tmp_r = T::zero_with_precision_of(modulus) + &new_r;
        new_r = r - new_r.checked_mul(quotient2).expect(signed::OVERFLOW_MSG);
        r = tmp_r;
    }

    if r > T::one() {
        return None;
    }

    if t < T::zero().into() {
        t = t + modulus;
    }

    Some(t.into_inner())
}

// Another version with less constraints
// Still uses input references. Requires Clone

/// # Modular Inverse (Constrained)
/// Version that works with references. Requires Clone and
/// reference-based operations.
///
/// # Panics
///
/// Panics if `T` is too narrow to hold the extended-GCD intermediates
/// (checked coefficient arithmetic overflows) — a carrier-precondition
/// violation, distinct from the `None` return for a non-invertible input.
pub fn constrained_mod_inv<T>(a: T, modulus: &T) -> Option<T>
where
    // Checked mul/add throughout (see `basic_mod_inv`): both the `Signed`
    // coefficients and the raw `T` r-sequence refuse to rely on `*`'s
    // implementation-defined overflow.
    T: const_num_traits::Zero
        + const_num_traits::One
        + Clone
        + PartialEq
        + core::cmp::PartialOrd
        + const_num_traits::CheckedAdd<Output = T>
        + core::ops::Sub<Output = T>
        + const_num_traits::CheckedMul<Output = T>
        + const_num_traits::WithPrecision,
    for<'a> T: core::ops::Add<&'a T, Output = T> + core::ops::Sub<&'a T, Output = T>,
    for<'a> &'a T: core::ops::Sub<T, Output = T> + core::ops::Div<&'a T, Output = T>,
{
    // Seed the signed coefficients at the modulus width (as strict_mod_inv):
    // on a runtime-width carrier `T::zero()`/`one()` are minimal-width and the
    // signed magnitude arithmetic corrupts the inverse otherwise.
    let mut t = Signed::new(T::zero_with_precision_of(modulus), false);
    let mut new_t = Signed::new(T::one_with_precision_of(modulus), false);
    let mut r = modulus.clone();
    let mut new_r = a.widen_to_precision_of(modulus);

    while new_r != T::zero() {
        let quotient = &r / &new_r;

        let tmp_t = new_t.clone();
        new_t = t - new_t * quotient.clone();
        t = tmp_t;

        let tmp_r = new_r.clone();
        new_r = r - new_r.checked_mul(quotient).expect(signed::OVERFLOW_MSG);
        r = tmp_r;
    }

    if r > T::one() {
        return None;
    }

    if t < T::zero().into() {
        t = t + modulus;
    }

    Some(t.into_inner())
}

/// # Modular Inverse (Basic)
/// Simple version that operates on values and copies them.
///
/// # Panics
///
/// Panics if `T` is too narrow to hold the extended-GCD intermediates
/// (checked coefficient arithmetic overflows) — a carrier-precondition
/// violation, distinct from the `None` return for a non-invertible input.
pub fn basic_mod_inv<T>(a: T, modulus: T) -> Option<T>
where
    // `Signed<T>` uses checked mul/add, not plain `*`/`+` whose overflow on
    // `T` is implementation-defined. Products fit any carrier sized for the
    // modulus, so this guards no live overflow — it refuses to rely on
    // unspecified behavior, turning a too-narrow carrier into a clear panic
    // rather than a wrong inverse. `Sub` stays plain (smaller from larger).
    T: const_num_traits::Zero
        + const_num_traits::One
        + Copy
        + PartialEq
        + core::ops::Div<Output = T>
        + const_num_traits::CheckedAdd<Output = T>
        + core::ops::Sub<Output = T>
        + const_num_traits::CheckedMul<Output = T>
        + const_num_traits::WithPrecision
        + core::cmp::PartialOrd,
{
    // Seed the signed coefficients at the modulus width (as strict_mod_inv):
    // on a runtime-width carrier `T::zero()`/`one()` are minimal-width and the
    // signed magnitude arithmetic corrupts the inverse otherwise.
    let mut t = Signed::new(T::zero_with_precision_of(&modulus), false);
    let mut new_t = Signed::new(T::one_with_precision_of(&modulus), false);
    let mut r = Signed::new(modulus, false);
    let mut new_r = Signed::new(a.widen_to_precision_of(&modulus), false);

    while new_r != Signed::new(T::zero(), false) {
        let quotient = r / new_r;

        let tmp_t = new_t;
        new_t = t - new_t * quotient;
        t = tmp_t;

        let tmp_r = new_r;
        new_r = r - new_r * quotient;
        r = tmp_r;
    }

    if r > T::one().into() {
        return None;
    }

    if t < T::zero().into() {
        t = t + modulus.into();
    }

    Some(t.into_inner())
}

#[cfg(test)]
macro_rules! select_mod_inv {
    ($mod_inv:path, $t:ty, by_ref) => {
        fn mod_inv(a: $t, modulus: &$t) -> Option<$t> {
            $mod_inv(a, modulus)
        }
    };
    ($mod_inv:path, $t:ty, by_val) => {
        fn mod_inv(a: $t, modulus: &$t) -> Option<$t> {
            $mod_inv(a, *modulus)
        }
    };
}

#[cfg(test)]
macro_rules! generate_mod_inv_tests_block_1 {
    ($mod_inv:path, $t:ty , $by_ref:tt) => {
        select_mod_inv!($mod_inv, $t, $by_ref);

        #[test]
        fn test_mod_inv_1_mod_13() {
            assert_eq!(mod_inv(0u32, &7u32), None);
            assert_eq!(mod_inv(1u32, &7u32).unwrap(), 1);
            assert_eq!(mod_inv(6u32, &8u32), None);
            assert_eq!(mod_inv(1u32, &13u32).unwrap(), 1u32); // 1 * 1 ≡ 1 (mod 13)
            assert_eq!(mod_inv(8u32, &13u32).unwrap(), 5u32); // Check: 8×5=40≡1mod  13.8×5=40≡1mod13.
            assert_eq!(mod_inv(12u32, &13u32).unwrap(), 12u32); // 1 * 1 ≡ 1 (mod 13)
            assert_eq!(mod_inv(14u32, &13u32).unwrap(), 1); // 14 * 10 ≡ 1 (mod 13)
            assert_eq!(mod_inv(15u32, &13u32).unwrap(), 7u32); // 15 * 9 ≡ 1 (mod 13)
            assert_eq!(mod_inv(16u32, &13u32).unwrap(), 9u32); // 16 * 8 ≡ 1 (mod 13)
            assert_eq!(mod_inv(10u32, &17).unwrap(), 12);
        }
    };
}

#[cfg(test)]
mod strict_mod_inv_tests {
    generate_mod_inv_tests_block_1!(super::strict_mod_inv, u32, by_ref);
}

#[cfg(test)]
mod constrained_mod_inv_tests {
    generate_mod_inv_tests_block_1!(super::constrained_mod_inv, u32, by_ref);
}

#[cfg(test)]
mod basic_mod_inv_tests {
    generate_mod_inv_tests_block_1!(super::basic_mod_inv, u32, by_val);
}

#[cfg(test)]
macro_rules! inv_test_module {
    (
        $stem:ident,
        $type_path:path,
        $(type $type_def:ty = $type_expr:ty;)? // Optional type definition
        strict: $strict:expr,
        constrained: $constrained:expr,
        basic: $basic:expr,
    ) => {
        paste::paste! {
            mod [<$stem _tests>] {
                #[allow(unused_imports)]
                use $type_path;
                $( type $type_def = $type_expr; )?

                #[test]
                #[allow(unused_variables)]
                fn test_mod_inv_basic() {
                    let a_val = 5u8;
                    let a = U256::from(a_val);
                    let modulus = U256::from(13u8);
                    let result_val = 8u8;
                    let result = U256::from(result_val);

                    crate::maybe_test!($strict, assert_eq!(super::strict_mod_inv(a, &modulus), Some(result)));
                    let a = U256::from(a_val);
                    let result = U256::from(result_val);
                    crate::maybe_test!($constrained, assert_eq!(super::constrained_mod_inv(a, &modulus), Some(result)));
                    let a = U256::from(a_val);
                    let result = U256::from(result_val);
                    crate::maybe_test!($basic, assert_eq!(super::basic_mod_inv(a, modulus), Some(result)));
                }
            }
        }
    };
}

#[cfg(test)]
mod bnum_inv_tests {
    use super::basic_mod_inv;
    use super::constrained_mod_inv;
    use super::strict_mod_inv;

    //     inv_test_module!(
    //         bnum,
    //         bnum::types::U256,
    //         strict: on,
    //         constrained: on,
    //         basic: on,
    //     );

    inv_test_module!(
        bnum_patched,
        bnum_patched::types::U256,
        strict: on,
        constrained: on,
        basic: on,
    );

    //     inv_test_module!(
    //         crypto_bigint,
    //         crypto_bigint::U256,
    //         strict: off, // &'a Div missing
    //         constrained: off, // &'a Div missing
    //         basic: on,
    //     );

    inv_test_module!(
        crypto_bigint_patched,
        crypto_bigint_patched::U256,
        strict: on,
        constrained: on,
        basic: on,
    );

    //     inv_test_module!(
    //         num_bigint,
    //         num_bigint::BigUint,
    //         type U256 = num_bigint::BigUint;
    //         strict: on,
    //         constrained: on,
    //         basic: off, // Copy is not implemented, heap
    //     );

    // num-bigint `FixedWidthBigUint`: heap carrier, Nct, constrained/strict
    // only (not `Copy`, so `basic: off`).
    inv_test_module!(
        num_bigint_patched,
        num_bigint_patched::FixedWidthBigUint,
        type U256 = num_bigint_patched::FixedWidthBigUint;
        strict: on,
        constrained: on,
        basic: off,
    );

    //     inv_test_module!(
    //         ibig,
    //         ibig::UBig,
    //         type U256 = ibig::UBig;
    //         strict: on,
    //         constrained: on,
    //         basic: off, // Copy is not implemented, heap
    //     );

    //     inv_test_module!(
    //         ibig_patched,
    //         ibig_patched::UBig,
    //         type U256 = ibig_patched::UBig;
    //         strict: on,
    //         constrained: on,
    //         basic: off, // Copy is not implemented, heap
    //     );

    // fixed_bigint inv tests run on the Nct personality (the type alias
    // default). The Ct personality intentionally lacks `Div`/`Rem` because
    // EEA's inner loop count depends on operand magnitudes — schoolbook
    // modular inverse is fundamentally variable-time and so does not
    // belong on a CT-typed FixedUInt. Calling `basic_mod_inv` etc. with a
    // Ct-typed value is a compile error rather than silent NCT execution,
    // which is the safety property the personality typestate exists to
    // enforce.
    inv_test_module!(
        fixed_bigint,
        fixed_bigint::FixedUInt,
        type U256 = fixed_bigint::FixedUInt<u8, 4>;
        strict: on,
        constrained: on,
        basic: on,
    );

    inv_test_module!(
        heapless_bigint,
        fixed_bigint::FixedUInt,
        type U256 = fixed_bigint::HeaplessBigInt<u8, 4>;
        strict: on,
        constrained: on,
        basic: on,
    );

    // Regression for the width-0 seed bug in the EEA inverses: on a runtime-width
    // carrier (len < CAP) `zero()`/`one()`-seeded coefficients corrupted the
    // inverse for a large fraction of residues. All three flavors now seed at the
    // modulus width; cross-check every residue against the u32 oracle.
    #[test]
    fn mod_inv_all_flavors_runtime_width_full_range() {
        use fixed_bigint::HeaplessBigInt;
        type H = HeaplessBigInt<u8, 4>;
        for m in [13u8, 17, 251] {
            let hm = H::from(m);
            for a in 1..m {
                let ha = H::from(a);
                let want = super::basic_mod_inv(a as u32, m as u32).map(|w| H::from(w as u8));
                assert_eq!(super::basic_mod_inv(ha, hm), want, "basic inv({a}) mod {m}");
                assert_eq!(
                    super::constrained_mod_inv(ha, &hm),
                    want,
                    "constrained inv({a}) mod {m}"
                );
                assert_eq!(
                    super::strict_mod_inv(ha, &hm),
                    want,
                    "strict inv({a}) mod {m}"
                );
            }
        }
    }

    // Multi-limb guard for num-bigint's FixedWidthBigUint (const-7 minimized its
    // width handling): cross-check inv + schoolbook mul against a u128 oracle at
    // a >1-limb odd modulus. The single-limb macro tests can't see a multi-limb
    // width regression; this can.
    #[test]
    fn num_bigint_multilimb_matches_u128_oracle() {
        use num_bigint_patched::{BigUint, FixedWidthBigUint as FW};
        const N: usize = 4; // >= 128 bits even for u32 digits; holds m + EEA headroom
        let m_u128: u128 = 1_180_591_620_717_411_303_451; // ~2^70, odd, multi-limb
        let fw = |v: u128| FW::new(BigUint::from(v), N);
        let m = fw(m_u128);
        for a in [2u128, 3, 5, 65_537, 999_999_937, m_u128 - 1] {
            for (label, gotv, wantv) in [
                (
                    "strict",
                    super::strict_mod_inv(fw(a), &m),
                    super::strict_mod_inv(a, &m_u128),
                ),
                (
                    "constrained",
                    super::constrained_mod_inv(fw(a), &m),
                    super::constrained_mod_inv(a, &m_u128),
                ),
            ] {
                match (gotv, wantv) {
                    (Some(g), Some(w)) => assert_eq!(g, fw(w), "{label} inv({a})"),
                    (None, None) => {}
                    (g, w) => {
                        panic!("{label} inv({a}): FW.is_some={} u128={w:?}", g.is_some())
                    }
                }
            }
        }
        // Multi-limb schoolbook mul against the u128 oracle.
        for (a, b) in [
            (3u128, 5u128),
            (123_456_789, 987_654_321),
            (m_u128 - 1, m_u128 - 2),
        ] {
            let got = crate::mul::constrained_mod_mul(fw(a), &fw(b), &m);
            let want = crate::mul::constrained_mod_mul(a, &b, &m_u128);
            assert_eq!(got, fw(want), "mul({a},{b})");
        }
    }
}
