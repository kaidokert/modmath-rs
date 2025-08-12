// Constrained Montgomery arithmetic functions
// These work with references to avoid unnecessary copies, following the pattern from exp.rs

use super::basic_mont::NPrimeMethod;
use crate::inv::constrained_mod_inv;

/// Compute N' using trial search method - O(R) complexity (constrained version)
/// Finds N' such that modulus * N' ≡ -1 (mod R)
/// Returns None if N' cannot be found (should not happen for valid Montgomery parameters)
fn compute_n_prime_trial_search_constrained<T>(modulus: &T, r: &T) -> Option<T>
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T> + core::ops::Mul<&'a T, Output = T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>,
{
    // We need to find N' where modulus * N' ≡ R - 1 (mod R)
    let target = r.clone().wrapping_sub(&T::one()); // This is -1 mod R

    // Simple trial search for N'
    let mut n_prime = T::one();
    loop {
        if (modulus.clone() * &n_prime) % r == target {
            return Some(n_prime);
        }
        n_prime = n_prime.wrapping_add(&T::one());

        // Safety check to avoid infinite loop
        if &n_prime >= r {
            return None; // Could not find N' - should not happen for valid inputs
        }
    }
}

/// Compute N' using Extended Euclidean Algorithm - O(log R) complexity (constrained version)
/// Finds N' such that modulus * N' ≡ -1 (mod R)
/// Returns None if modular inverse cannot be found
fn compute_n_prime_extended_euclidean_constrained<T>(modulus: &T, r: &T) -> Option<T>
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub,
    for<'a> T: core::ops::Add<&'a T, Output = T> + core::ops::Sub<&'a T, Output = T>,
    for<'a> &'a T: core::ops::Sub<T, Output = T> + core::ops::Div<&'a T, Output = T>,
{
    // We need to solve: modulus * N' ≡ -1 (mod R)
    // This is equivalent to: N' ≡ -modulus^(-1) (mod R)

    // Use constrained_mod_inv to find modulus^(-1) mod R
    if let Some(modulus_inv) = constrained_mod_inv(modulus.clone(), r) {
        // N' = -modulus^(-1) mod R = R - modulus^(-1) mod R
        if modulus_inv == T::zero() {
            Some(r.clone().wrapping_sub(&T::one())) // Handle edge case where inverse is 0
        } else {
            Some(r.clone().wrapping_sub(&modulus_inv))
        }
    } else {
        None // Could not find modular inverse - gcd(modulus, R) should be 1 for valid Montgomery parameters
    }
}

/// Compute N' using Hensel's lifting - O(log R) complexity, optimized for R = 2^k (constrained version)
/// Finds N' such that modulus * N' ≡ -1 (mod R)
/// Returns None if Hensel's lifting fails to produce correct N'
fn compute_n_prime_hensels_lifting_constrained<T>(modulus: &T, r: &T, r_bits: usize) -> Option<T>
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shl<usize, Output = T>
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T> + core::ops::Mul<&'a T, Output = T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T> + core::ops::BitAnd<&'a T, Output = T>,
{
    // Hensel's lifting requires modulus to be odd (prerequisite for Montgomery arithmetic)
    debug_assert!(
        modulus & &T::one() == T::one(),
        "Hensel's lifting requires an odd modulus for Montgomery arithmetic"
    );

    // Hensel's lifting for N' computation when R = 2^k
    let mut n_prime = T::one();

    // Lift from 2^1 to 2^r_bits using Newton's method
    for k in 2..=r_bits {
        let target_mod = T::one() << k; // 2^k
        let temp_prod = modulus.clone() * &n_prime;
        let temp_sum = temp_prod.wrapping_add(&T::one());
        let check_val = &temp_sum % &target_mod;

        if check_val != T::zero() {
            let prev_power = T::one() << (k - 1); // 2^(k-1)
            if check_val == prev_power {
                n_prime = n_prime.wrapping_add(&prev_power);
            }
        }
    }

    // Final check
    let final_check = (modulus.clone() * &n_prime) % r;
    let target = r.clone().wrapping_sub(&T::one()); // -1 mod R

    if final_check != target {
        None // Hensel lifting failed to produce correct N'
    } else {
        Some(n_prime)
    }
}

/// Montgomery parameter computation (Constrained)
/// Computes R, R^(-1) mod N, N', and R bit length for Montgomery arithmetic
/// Returns None if N' computation fails or R^(-1) mod N cannot be found
pub fn constrained_compute_montgomery_params_with_method<T>(
    modulus: &T,
    method: NPrimeMethod,
) -> Option<(T, T, T, usize)>
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shl<usize, Output = T>
        + core::ops::Sub<Output = T>
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::Add<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::RemAssign<&'a T>,
    for<'a> &'a T: core::ops::Sub<T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Rem<&'a T, Output = T>
        + core::ops::BitAnd<&'a T, Output = T>,
{
    // Step 1: Find R = 2^k where R > modulus
    let mut r = T::one();
    let mut r_bits = 0usize;

    while &r <= modulus {
        r = r << 1; // r *= 2
        r_bits += 1;
    }

    // Step 2: Compute R^(-1) mod modulus
    let r_inv = constrained_mod_inv(r.clone(), modulus)?;

    // Step 3: Compute N' such that N * N' ≡ -1 (mod R) using selected method
    let n_prime = match method {
        NPrimeMethod::TrialSearch => compute_n_prime_trial_search_constrained(modulus, &r)?,
        NPrimeMethod::ExtendedEuclidean => {
            compute_n_prime_extended_euclidean_constrained(modulus, &r)?
        }
        NPrimeMethod::HenselsLifting => {
            compute_n_prime_hensels_lifting_constrained(modulus, &r, r_bits)?
        }
    };

    Some((r, r_inv, n_prime, r_bits))
}

/// Montgomery parameter computation (Constrained) with default method
/// Computes R, R^(-1) mod N, N', and R bit length for Montgomery arithmetic
/// Returns None if parameter computation fails
pub fn constrained_compute_montgomery_params<T>(modulus: &T) -> Option<(T, T, T, usize)>
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shl<usize, Output = T>
        + core::ops::Sub<Output = T>
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::Add<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::RemAssign<&'a T>,
    for<'a> &'a T: core::ops::Sub<T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Rem<&'a T, Output = T>
        + core::ops::BitAnd<&'a T, Output = T>,
{
    constrained_compute_montgomery_params_with_method(modulus, NPrimeMethod::default())
}

/// Convert to Montgomery form (Constrained): a -> (a * R) mod N
pub fn constrained_to_montgomery<T>(a: T, modulus: &T, r: &T) -> T
where
    T: num_traits::Zero
        + num_traits::One
        + PartialOrd
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shr<usize, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T> + core::ops::BitAnd<Output = T>,
{
    crate::mul::constrained_mod_mul(a, r, modulus)
}

/// Convert from Montgomery form (Constrained): (a * R) -> a mod N
/// Uses Montgomery reduction algorithm
pub fn constrained_from_montgomery<T>(a_mont: T, modulus: &T, n_prime: &T, r_bits: usize) -> T
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Shr<usize, Output = T>
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::Mul<&'a T, Output = T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>,
{
    // Montgomery reduction algorithm:
    let r = T::one() << r_bits; // R = 2^r_bits

    // Step 1: m = (a_mont * N') mod R
    let m = (a_mont.clone() * n_prime) % &r;

    // Step 2: t = (a_mont + m * N) / R
    let temp_prod = m.clone() * modulus;
    let temp_sum = a_mont.wrapping_add(&temp_prod);
    let t = temp_sum >> r_bits; // Divide by R = 2^r_bits

    // Step 3: Final reduction
    if &t >= modulus {
        t.wrapping_sub(modulus)
    } else {
        t
    }
}

/// Montgomery multiplication (Constrained): (a * R) * (b * R) -> (a * b * R) mod N
pub fn constrained_montgomery_mul<T>(
    a_mont: &T,
    b_mont: &T,
    modulus: &T,
    n_prime: &T,
    r_bits: usize,
) -> T
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Shr<usize, Output = T>
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T> + core::ops::Mul<&'a T, Output = T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T> + core::ops::BitAnd<Output = T>,
{
    // Step 1: Regular modular multiplication in Montgomery domain
    let product = crate::mul::constrained_mod_mul(a_mont.clone(), b_mont, modulus);

    // Step 2: Apply Montgomery reduction to get result in Montgomery form
    constrained_from_montgomery(product, modulus, n_prime, r_bits)
}

/// Complete Montgomery modular multiplication with method selection (Constrained): A * B mod N
/// Returns None if Montgomery parameter computation fails
pub fn constrained_montgomery_mod_mul_with_method<T>(
    a: T,
    b: &T,
    modulus: &T,
    method: NPrimeMethod,
) -> Option<T>
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shl<usize, Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Sub<Output = T>
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::Add<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::RemAssign<&'a T>,
    for<'a> &'a T: core::ops::Sub<T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Rem<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    let (r, _r_inv, n_prime, r_bits) =
        constrained_compute_montgomery_params_with_method(modulus, method)?;
    let a_mont = constrained_to_montgomery(a, modulus, &r);
    let b_mont = constrained_to_montgomery(b.clone(), modulus, &r);
    let result_mont = constrained_montgomery_mul(&a_mont, &b_mont, modulus, &n_prime, r_bits);
    Some(constrained_from_montgomery(
        result_mont,
        modulus,
        &n_prime,
        r_bits,
    ))
}

/// Complete Montgomery modular multiplication (Constrained): A * B mod N
/// Returns None if Montgomery parameter computation fails
pub fn constrained_montgomery_mod_mul<T>(a: T, b: &T, modulus: &T) -> Option<T>
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shl<usize, Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Sub<Output = T>
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::Add<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::RemAssign<&'a T>,
    for<'a> &'a T: core::ops::Sub<T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Rem<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    constrained_montgomery_mod_mul_with_method(a, b, modulus, NPrimeMethod::default())
}

/// Montgomery-based modular exponentiation with method selection (Constrained): base^exponent mod modulus
/// Uses Montgomery arithmetic for efficient repeated multiplication
/// Returns None if Montgomery parameter computation fails
pub fn constrained_montgomery_mod_exp_with_method<T>(
    mut base: T,
    exponent: &T,
    modulus: &T,
    method: NPrimeMethod,
) -> Option<T>
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shl<usize, Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::ShrAssign<usize>
        + core::ops::Sub<Output = T>
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>,
    for<'a> &'a T: core::ops::Sub<T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Rem<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    // Compute Montgomery parameters using specified method
    let (r, _r_inv, n_prime, r_bits) =
        constrained_compute_montgomery_params_with_method(modulus, method)?;

    // Reduce base and convert to Montgomery form
    base.rem_assign(modulus);
    base = constrained_to_montgomery(base, modulus, &r);

    // Montgomery form of 1 (the initial result)
    let mut result = constrained_to_montgomery(T::one(), modulus, &r);

    // Copy exponent for manipulation
    let mut exp = exponent.clone();
    let two = T::one().clone().wrapping_add(&T::one());

    // Binary exponentiation using Montgomery multiplication
    while exp > T::zero() {
        // If exponent is odd, multiply result by current base power
        if &exp % &two == T::one() {
            result = constrained_montgomery_mul(&result, &base, modulus, &n_prime, r_bits);
        }

        // Square the base for next iteration
        exp >>= 1;
        base = constrained_montgomery_mul(&base, &base, modulus, &n_prime, r_bits);
    }

    // Convert result back from Montgomery form
    Some(constrained_from_montgomery(
        result, modulus, &n_prime, r_bits,
    ))
}

/// Montgomery-based modular exponentiation (Constrained): base^exponent mod modulus
/// Uses Montgomery arithmetic for efficient repeated multiplication
/// Returns None if Montgomery parameter computation fails
pub fn constrained_montgomery_mod_exp<T>(base: T, exponent: &T, modulus: &T) -> Option<T>
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shl<usize, Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::ShrAssign<usize>
        + core::ops::Sub<Output = T>
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>,
    for<'a> &'a T: core::ops::Sub<T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Rem<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    constrained_montgomery_mod_exp_with_method(base, exponent, modulus, NPrimeMethod::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constrained_compute_n_prime_trial_search_failure() {
        // Test case where no N' can be found - this happens when gcd(modulus, R) != 1
        // Example: N = 4, R = 8 (both even, so gcd(4, 8) = 4 != 1)
        let modulus = 4u32;
        let r = 8u32;
        let result = compute_n_prime_trial_search_constrained(&modulus, &r);
        assert!(
            result.is_none(),
            "Should return None for invalid modulus-R pair"
        );
    }

    #[test]
    fn test_constrained_compute_montgomery_params_failure() {
        // Test Montgomery parameter computation failure with even modulus
        let even_modulus = 4u32;
        let result = constrained_compute_montgomery_params(&even_modulus);
        assert!(result.is_none(), "Should return None for even modulus");
    }

    #[test]
    fn test_constrained_compute_montgomery_params_failure_with_method() {
        // Test all N' computation methods with invalid inputs
        let invalid_modulus = 4u32; // Even modulus

        // Trial search should fail
        let trial_result = constrained_compute_montgomery_params_with_method(
            &invalid_modulus,
            NPrimeMethod::TrialSearch,
        );
        assert!(
            trial_result.is_none(),
            "Trial search should fail with even modulus"
        );

        // Extended Euclidean should fail
        let euclidean_result = constrained_compute_montgomery_params_with_method(
            &invalid_modulus,
            NPrimeMethod::ExtendedEuclidean,
        );
        assert!(
            euclidean_result.is_none(),
            "Extended Euclidean should fail with even modulus"
        );

        // Hensel's lifting should fail
        let hensels_result = constrained_compute_montgomery_params_with_method(
            &invalid_modulus,
            NPrimeMethod::HenselsLifting,
        );
        assert!(
            hensels_result.is_none(),
            "Hensel's lifting should fail with even modulus"
        );
    }

    #[test]
    fn test_constrained_montgomery_mod_mul_parameter_failure() {
        // Test that montgomery_mod_mul returns None when parameter computation fails
        let invalid_modulus = 4u32;
        let a = 2u32;
        let b = 3u32;

        let result = constrained_montgomery_mod_mul(a, &b, &invalid_modulus);
        assert!(
            result.is_none(),
            "Montgomery mod_mul should return None for invalid modulus"
        );
    }

    #[test]
    fn test_constrained_montgomery_mod_exp_parameter_failure() {
        // Test that montgomery_mod_exp returns None when parameter computation fails
        let invalid_modulus = 4u32;
        let base = 2u32;
        let exponent = 3u32;

        let result = constrained_montgomery_mod_exp(base, &exponent, &invalid_modulus);
        assert!(
            result.is_none(),
            "Montgomery mod_exp should return None for invalid modulus"
        );
    }

    #[test]
    fn test_constrained_montgomery_reduction_final_subtraction() {
        // Test to trigger t >= modulus branch in constrained_from_montgomery
        let modulus = 15u32;
        let (r, _r_inv, n_prime, r_bits) = constrained_compute_montgomery_params(&modulus).unwrap();

        // Test with maximum value to potentially trigger final subtraction
        let high_value = 14u32;
        let mont_high = constrained_to_montgomery(high_value, &modulus, &r);
        let result = constrained_from_montgomery(mont_high, &modulus, &n_prime, r_bits);
        assert_eq!(result, high_value);

        // Test with another high value
        let mont_13 = constrained_to_montgomery(13u32, &modulus, &r);
        let result_13 = constrained_from_montgomery(mont_13, &modulus, &n_prime, r_bits);
        assert_eq!(result_13, 13u32);
    }

    #[test]
    fn test_constrained_hensel_lifting_branches() {
        // Test Hensel's lifting with moduli that may trigger different conditional paths
        let test_moduli = [9u32, 15u32, 21u32, 35u32, 45u32]; // Various composite odd moduli

        for &modulus in &test_moduli {
            let hensels_result = constrained_compute_montgomery_params_with_method(
                &modulus,
                crate::montgomery::NPrimeMethod::HenselsLifting,
            );

            assert!(
                hensels_result.is_some(),
                "Hensel's lifting should work for modulus {}",
                modulus
            );

            // Verify the result is mathematically correct
            if let Some((r, _r_inv, n_prime, _r_bits)) = &hensels_result {
                let check = (modulus * n_prime.clone()) % r.clone();
                let expected = r.clone() - 1; // Should equal R - 1 (which is -1 mod R)
                assert_eq!(
                    check, expected,
                    "N' verification failed for modulus {} with Hensel's lifting",
                    modulus
                );
            }
        }
    }

    #[test]
    fn test_constrained_multiplication_stress() {
        // Test multiplication with values designed to stress different code paths
        let modulus = 33u32; // 33 = 3 * 11, composite modulus
        let (r, _r_inv, n_prime, r_bits) = constrained_compute_montgomery_params(&modulus).unwrap();

        // Test with values that may cause intermediate results needing reduction
        let test_pairs = [(31u32, 32u32), (29u32, 30u32), (25u32, 27u32)];

        for (a, b) in test_pairs.iter() {
            let a_mont = constrained_to_montgomery(*a, &modulus, &r);
            let b_mont = constrained_to_montgomery(*b, &modulus, &r);

            // This may hit different branches in Montgomery multiplication
            let result_mont =
                constrained_montgomery_mul(&a_mont, &b_mont, &modulus, &n_prime, r_bits);
            let result = constrained_from_montgomery(result_mont, &modulus, &n_prime, r_bits);

            let expected = (a * b) % modulus;
            assert_eq!(result, expected, "Failed for {} * {} mod {}", a, b, modulus);
        }
    }

    #[test]
    fn test_constrained_exponentiation_conditional_branches() {
        // Test exponentiation with specific values to hit different loop branches
        let modulus = 19u32; // Prime modulus

        // Test with base near modulus and various exponent patterns
        let test_cases = [
            (18u32, 2u32),  // Base near modulus, small exponent
            (17u32, 7u32),  // Various combinations
            (15u32, 31u32), // Larger exponent with specific bit pattern
            (2u32, 127u32), // Small base, large exponent
        ];

        for (base, exponent) in test_cases.iter() {
            let result = constrained_montgomery_mod_exp(*base, exponent, &modulus).unwrap();
            let expected = crate::exp::constrained_mod_exp(*base, exponent, &modulus);
            assert_eq!(
                result, expected,
                "Failed for {}^{} mod {}",
                base, exponent, modulus
            );
        }
    }

    #[test]
    fn test_constrained_extended_euclidean_edge_cases() {
        // Test Extended Euclidean N' computation with various moduli
        let edge_moduli = [7u32, 9u32, 25u32, 49u32, 121u32]; // Powers of primes and products

        for &modulus in &edge_moduli {
            let euclidean_result = constrained_compute_montgomery_params_with_method(
                &modulus,
                crate::montgomery::NPrimeMethod::ExtendedEuclidean,
            );

            assert!(
                euclidean_result.is_some(),
                "Extended Euclidean should work for modulus {}",
                modulus
            );

            // Cross-validate with trial search
            let trial_result = constrained_compute_montgomery_params_with_method(
                &modulus,
                crate::montgomery::NPrimeMethod::TrialSearch,
            );

            assert_eq!(
                euclidean_result, trial_result,
                "Extended Euclidean vs Trial Search mismatch for modulus {}",
                modulus
            );
        }
    }
}
