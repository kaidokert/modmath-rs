// Basic Montgomery arithmetic functions
// These require Copy trait but have minimal constraints

use crate::inv::basic_mod_inv;

/// Methods for computing N' in Montgomery parameter computation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NPrimeMethod {
    /// Trial search - O(R) complexity, simple but slow for large R
    TrialSearch,
    /// Extended Euclidean Algorithm - O(log R) complexity
    #[default]
    ExtendedEuclidean,
    /// Hensel's lifting - O(log R) complexity, optimized for R = 2^k
    HenselsLifting,
}

/// Compute N' using trial search method - O(R) complexity
/// Finds N' such that modulus * N' ≡ -1 (mod R)
/// Returns None if N' cannot be found
fn compute_n_prime_trial_search<T>(modulus: T, r: T) -> Option<T>
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Add<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Mul<Output = T>
        + core::ops::Rem<Output = T>,
{
    // We need to find N' where modulus * N' ≡ R - 1 (mod R)
    let target = r - T::one(); // This is -1 mod R

    // Simple trial search for N'
    // TODO: Replace with Extended Euclidean Algorithm for O(log R) complexity instead of O(R)
    // Current implementation is fine for small numbers but inefficient for large moduli
    let mut n_prime = T::one();
    loop {
        if (modulus * n_prime) % r == target {
            return Some(n_prime);
        }
        n_prime = n_prime + T::one();

        // Safety check to avoid infinite loop
        if n_prime >= r {
            return None; // Could not find N' - should not happen for valid inputs
        }
    }
}

/// Compute N' using Extended Euclidean Algorithm - O(log R) complexity
/// Finds N' such that modulus * N' ≡ -1 (mod R)
/// This is equivalent to N' ≡ -modulus^(-1) (mod R)
/// Returns None if modular inverse cannot be found
fn compute_n_prime_extended_euclidean<T>(modulus: T, r: T) -> Option<T>
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Add<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Mul<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::Div<Output = T>,
{
    // We need to solve: modulus * N' ≡ -1 (mod R)
    // This is equivalent to: modulus * N' ≡ R - 1 (mod R)
    // So: N' ≡ (R - 1) * modulus^(-1) (mod R)
    // Or: N' ≡ -modulus^(-1) (mod R)

    // Use basic_mod_inv to find modulus^(-1) mod R
    if let Some(modulus_inv) = basic_mod_inv(modulus, r) {
        // N' = -modulus^(-1) mod R = R - modulus^(-1) mod R
        if modulus_inv == T::zero() {
            Some(r - T::one()) // Handle edge case where inverse is 0
        } else {
            Some(r - modulus_inv)
        }
    } else {
        None // Could not find modular inverse - gcd(modulus, R) should be 1 for valid Montgomery parameters
    }
}

/// Compute N' using Hensel's lifting - O(log R) complexity, optimized for R = 2^k
/// Finds N' such that modulus * N' ≡ -1 (mod R)
/// Uses Newton's method to iteratively lift from small powers to full R
/// Returns None if Hensel's lifting fails
fn compute_n_prime_hensels_lifting<T>(modulus: T, r: T, r_bits: usize) -> Option<T>
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Add<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Mul<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::Shl<usize, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    // Hensel's lifting for N' computation when R = 2^k
    // Start with base case: find N' such that modulus * N' ≡ -1 (mod 2)
    // Then iteratively lift to larger powers of 2

    // Base case: modulus * N' ≡ -1 ≡ 1 (mod 2)
    // Since modulus is odd (required for Montgomery), modulus ≡ 1 (mod 2)
    // So we need N' ≡ 1 (mod 2), hence N' starts as 1
    let mut n_prime = T::one();

    // Lift from 2^1 to 2^r_bits using Newton's method
    for k in 2..=r_bits {
        // We have: modulus * n_prime ≡ -1 (mod 2^(k-1))
        // We want: modulus * n_prime_new ≡ -1 (mod 2^k)

        // Newton iteration: x_new = x - f(x)/f'(x)
        // Where f(x) = modulus * x + 1
        // f'(x) = modulus
        // So: x_new = x - (modulus * x + 1) / modulus
        //     x_new = x - x - 1/modulus  (but we work mod powers of 2)

        let target_mod = T::one() << k; // 2^k
        let check_val = (modulus * n_prime + T::one()) % target_mod;

        if check_val != T::zero() {
            // Need to adjust n_prime
            // If modulus * n_prime + 1 = t * 2^(k-1) for odd t, add 2^(k-1) to n_prime
            let prev_power = T::one() << (k - 1); // 2^(k-1)

            if check_val == prev_power {
                n_prime = n_prime + prev_power;
            }
        }
    }

    // Final check and adjustment to ensure modulus * N' ≡ -1 (mod R)
    let final_check = (modulus * n_prime) % r;
    let target = r - T::one(); // -1 mod R

    if final_check != target {
        // This shouldn't happen with correct Hensel lifting, but safety check
        None // Hensel lifting failed to produce correct N'
    } else {
        Some(n_prime)
    }
}

/// Montgomery parameter computation (Basic)
/// Computes R, R^(-1) mod N, N', and R bit length for Montgomery arithmetic
/// Returns None if parameter computation fails
pub fn basic_compute_montgomery_params_with_method<T>(
    modulus: T,
    method: NPrimeMethod,
) -> Option<(T, T, T, usize)>
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Div<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Mul<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::Add<Output = T>
        + core::ops::BitAnd<Output = T>,
{
    // Step 1: Find R = 2^k where R > modulus
    let mut r = T::one();
    let mut r_bits = 0usize;

    while r <= modulus {
        r = r << 1; // r *= 2
        r_bits += 1;
    }

    // Step 2: Compute R^(-1) mod modulus
    let r_inv = basic_mod_inv(r, modulus)?;

    // Step 3: Compute N' such that N * N' ≡ -1 (mod R) using selected method
    let n_prime = match method {
        NPrimeMethod::TrialSearch => compute_n_prime_trial_search(modulus, r)?,
        NPrimeMethod::ExtendedEuclidean => compute_n_prime_extended_euclidean(modulus, r)?,
        NPrimeMethod::HenselsLifting => compute_n_prime_hensels_lifting(modulus, r, r_bits)?,
    };

    Some((r, r_inv, n_prime, r_bits))
}

/// Montgomery parameter computation (Basic) with default method
/// Computes R, R^(-1) mod N, N', and R bit length for Montgomery arithmetic using trial search
/// Returns None if parameter computation fails
pub fn basic_compute_montgomery_params<T>(modulus: T) -> Option<(T, T, T, usize)>
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Div<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Mul<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::Add<Output = T>
        + core::ops::BitAnd<Output = T>,
{
    basic_compute_montgomery_params_with_method(modulus, NPrimeMethod::default())
}

/// Convert to Montgomery form (Basic): a -> (a * R) mod N
pub fn basic_to_montgomery<T>(a: T, modulus: T, r: T) -> T
where
    T: core::cmp::PartialOrd
        + Copy
        + num_traits::Zero
        + num_traits::One
        + core::ops::BitAnd<Output = T>
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shr<usize, Output = T>
        + core::ops::Rem<Output = T>,
{
    crate::mul::basic_mod_mul(a, r, modulus)
}

/// Convert from Montgomery form (Basic): (a * R) -> a mod N
/// Uses Montgomery reduction algorithm
pub fn basic_from_montgomery<T>(a_mont: T, modulus: T, n_prime: T, r_bits: usize) -> T
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialOrd
        + core::ops::Mul<Output = T>
        + core::ops::Add<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Shl<usize, Output = T>,
{
    // Montgomery reduction algorithm:
    // Input: a_mont (Montgomery form), N (modulus), N', r_bits
    // 1. R = 2^r_bits
    // 2. m = (a_mont * N') mod R
    // 3. t = (a_mont + m * N) / R
    // 4. if t >= N then return t - N else return t

    let r = T::one() << r_bits; // R = 2^r_bits

    // Step 1: m = (a_mont * N') mod R
    let m = (a_mont * n_prime) % r;

    // Step 2: t = (a_mont + m * N) / R
    let t = (a_mont + m * modulus) >> r_bits; // Divide by R = 2^r_bits

    // Step 3: Final reduction
    if t >= modulus { t - modulus } else { t }
}

/// Montgomery multiplication (Basic): (a * R) * (b * R) -> (a * b * R) mod N
pub fn basic_montgomery_mul<T>(a_mont: T, b_mont: T, modulus: T, n_prime: T, r_bits: usize) -> T
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialOrd
        + core::ops::Mul<Output = T>
        + core::ops::Add<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Shl<usize, Output = T>
        + core::ops::BitAnd<Output = T>
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub,
{
    // Montgomery multiplication algorithm:
    // Input: a_mont, b_mont (both in Montgomery form), modulus N, N', r_bits
    // 1. Compute product = a_mont * b_mont (mod N)
    // 2. Apply Montgomery reduction to get (a * b * R) mod N

    // Step 1: Regular modular multiplication in Montgomery domain
    let product = crate::mul::basic_mod_mul(a_mont, b_mont, modulus);

    // Step 2: Apply Montgomery reduction to get result in Montgomery form
    basic_from_montgomery(product, modulus, n_prime, r_bits)
}

/// Complete Montgomery modular multiplication (Basic): A * B mod N
/// Returns None if Montgomery parameter computation fails
pub fn basic_montgomery_mod_mul<T>(a: T, b: T, modulus: T) -> Option<T>
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Div<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Mul<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::BitAnd<Output = T>
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shr<usize, Output = T>,
{
    let (r, _r_inv, n_prime, r_bits) = basic_compute_montgomery_params(modulus)?;
    let a_mont = basic_to_montgomery(a, modulus, r);
    let b_mont = basic_to_montgomery(b, modulus, r);
    let result_mont = basic_montgomery_mul(a_mont, b_mont, modulus, n_prime, r_bits);
    Some(basic_from_montgomery(result_mont, modulus, n_prime, r_bits))
}

/// Montgomery-based modular exponentiation (Basic): base^exponent mod modulus
/// Uses Montgomery arithmetic for efficient repeated multiplication
/// Returns None if Montgomery parameter computation fails
pub fn basic_montgomery_mod_exp<T>(mut base: T, exponent: T, modulus: T) -> Option<T>
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Div<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Mul<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::BitAnd<Output = T>
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shr<usize, Output = T>
        + core::ops::ShrAssign<usize>,
{
    // Compute Montgomery parameters
    let (r, _r_inv, n_prime, r_bits) = basic_compute_montgomery_params(modulus)?;

    // Convert base to Montgomery form
    base = basic_to_montgomery(base % modulus, modulus, r); // Reduce base first

    // Montgomery form of 1 (the initial result)
    let mut result = basic_to_montgomery(T::one(), modulus, r);

    // Copy exponent for manipulation
    let mut exp = exponent;

    // Binary exponentiation using Montgomery multiplication
    while exp > T::zero() {
        // If exponent is odd, multiply result by current base power
        if exp & T::one() == T::one() {
            result = basic_montgomery_mul(result, base, modulus, n_prime, r_bits);
        }

        // Square the base for next iteration
        exp >>= 1;
        if exp > T::zero() {
            base = basic_montgomery_mul(base, base, modulus, n_prime, r_bits);
        }
    }

    // Convert result back from Montgomery form
    Some(basic_from_montgomery(result, modulus, n_prime, r_bits))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_compute_n_prime_trial_search_failure() {
        // Test case where no N' can be found - this happens when gcd(modulus, R) != 1
        // Example: N = 4, R = 8 (both even, so gcd(4, 8) = 4 != 1)
        let modulus = 4u32;
        let r = 8u32;
        let result = compute_n_prime_trial_search(modulus, r);
        assert!(
            result.is_none(),
            "Should return None for invalid modulus-R pair"
        );
    }

    #[test]
    fn test_basic_compute_montgomery_params_failure() {
        // Test Montgomery parameter computation failure with even modulus
        let even_modulus = 4u32;
        let result = basic_compute_montgomery_params(even_modulus);
        assert!(result.is_none(), "Should return None for even modulus");
    }

    #[test]
    fn test_basic_compute_montgomery_params_failure_with_method() {
        // Test all N' computation methods with invalid inputs
        let invalid_modulus = 4u32; // Even modulus

        // Trial search should fail
        let trial_result =
            basic_compute_montgomery_params_with_method(invalid_modulus, NPrimeMethod::TrialSearch);
        assert!(
            trial_result.is_none(),
            "Trial search should fail with even modulus"
        );

        // Extended Euclidean should fail
        let euclidean_result = basic_compute_montgomery_params_with_method(
            invalid_modulus,
            NPrimeMethod::ExtendedEuclidean,
        );
        assert!(
            euclidean_result.is_none(),
            "Extended Euclidean should fail with even modulus"
        );

        // Hensel's lifting should fail
        let hensels_result = basic_compute_montgomery_params_with_method(
            invalid_modulus,
            NPrimeMethod::HenselsLifting,
        );
        assert!(
            hensels_result.is_none(),
            "Hensel's lifting should fail with even modulus"
        );
    }

    #[test]
    fn test_basic_montgomery_mod_mul_parameter_failure() {
        // Test that montgomery_mod_mul returns None when parameter computation fails
        let invalid_modulus = 4u32;
        let a = 2u32;
        let b = 3u32;

        let result = basic_montgomery_mod_mul(a, b, invalid_modulus);
        assert!(
            result.is_none(),
            "Montgomery mod_mul should return None for invalid modulus"
        );
    }

    #[test]
    fn test_basic_montgomery_mod_exp_parameter_failure() {
        // Test that montgomery_mod_exp returns None when parameter computation fails
        let invalid_modulus = 4u32;
        let base = 2u32;
        let exponent = 3u32;

        let result = basic_montgomery_mod_exp(base, exponent, invalid_modulus);
        assert!(
            result.is_none(),
            "Montgomery mod_exp should return None for invalid modulus"
        );
    }

    #[test]
    fn test_basic_montgomery_reduction_final_subtraction() {
        // Test to trigger t >= modulus branch in basic_from_montgomery
        let modulus = 15u32;
        let (r, _r_inv, n_prime, r_bits) = basic_compute_montgomery_params(modulus).unwrap();

        // Use values designed to need final subtraction in Montgomery reduction
        let high_value = 14u32; // Near maximum for this modulus
        let mont_high = basic_to_montgomery(high_value, modulus, r);
        let result = basic_from_montgomery(mont_high, modulus, n_prime, r_bits);
        assert_eq!(result, high_value);

        // Test with maximum value - 1 to stress the >= check
        let mont_max = basic_to_montgomery(modulus - 1, modulus, r);
        let result_max = basic_from_montgomery(mont_max, modulus, n_prime, r_bits);
        assert_eq!(result_max, modulus - 1);
    }

    #[test]
    fn test_basic_multiplication_edge_cases() {
        // Test Montgomery multiplication with values that stress the algorithm
        let modulus = 21u32; // Composite modulus
        let (r, _r_inv, n_prime, r_bits) = basic_compute_montgomery_params(modulus).unwrap();

        // Test with values that may trigger intermediate results >= modulus
        let a = 20u32; // Near maximum
        let b = 19u32; // Near maximum

        let a_mont = basic_to_montgomery(a, modulus, r);
        let b_mont = basic_to_montgomery(b, modulus, r);

        // This may hit different code paths in Montgomery multiplication
        let result_mont = basic_montgomery_mul(a_mont, b_mont, modulus, n_prime, r_bits);
        let result = basic_from_montgomery(result_mont, modulus, n_prime, r_bits);

        let expected = (a * b) % modulus;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_basic_n_prime_computation_edge_cases() {
        // Test N' computation with moduli that may hit different branches
        let test_moduli = [9u32, 15u32, 21u32, 25u32, 27u32]; // Various composite odd numbers

        for &modulus in &test_moduli {
            // Test all N' computation methods
            let trial_result =
                basic_compute_montgomery_params_with_method(modulus, NPrimeMethod::TrialSearch);
            let euclidean_result = basic_compute_montgomery_params_with_method(
                modulus,
                NPrimeMethod::ExtendedEuclidean,
            );
            let hensels_result =
                basic_compute_montgomery_params_with_method(modulus, NPrimeMethod::HenselsLifting);

            // All should succeed for valid odd moduli
            assert!(
                trial_result.is_some(),
                "Trial search failed for modulus {}",
                modulus
            );
            assert!(
                euclidean_result.is_some(),
                "Extended Euclidean failed for modulus {}",
                modulus
            );
            assert!(
                hensels_result.is_some(),
                "Hensel's lifting failed for modulus {}",
                modulus
            );

            // And all should produce same results
            assert_eq!(
                trial_result, euclidean_result,
                "Methods disagree for modulus {}",
                modulus
            );
            assert_eq!(
                trial_result, hensels_result,
                "Methods disagree for modulus {}",
                modulus
            );
        }
    }

    #[test]
    fn test_basic_exponentiation_loop_branches() {
        // Test binary exponentiation with values that hit different loop conditions
        let modulus = 17u32; // Prime modulus

        // Test exponentiation that exercises specific loop branches
        let base = 16u32; // Base = modulus - 1
        let exponent = 15u32; // Exponent with specific bit pattern

        let result = basic_montgomery_mod_exp(base, exponent, modulus).unwrap();
        let expected = crate::exp::basic_mod_exp(base, exponent, modulus);
        assert_eq!(result, expected);

        // Test with power of 2 exponent to hit specific branches
        let exp_pow2 = 16u32; // 2^4, will exercise specific bit patterns
        let result_pow2 = basic_montgomery_mod_exp(3u32, exp_pow2, modulus).unwrap();
        let expected_pow2 = crate::exp::basic_mod_exp(3u32, exp_pow2, modulus);
        assert_eq!(result_pow2, expected_pow2);

        // Test with odd exponent
        let exp_odd = 255u32; // Large odd exponent
        let result_odd = basic_montgomery_mod_exp(2u32, exp_odd, modulus).unwrap();
        let expected_odd = crate::exp::basic_mod_exp(2u32, exp_odd, modulus);
        assert_eq!(result_odd, expected_odd);
    }

    #[test]
    fn test_basic_extended_euclidean_none_case() {
        // Test the case where basic_mod_inv returns None in compute_n_prime_extended_euclidean
        // This happens when gcd(modulus, r) > 1, i.e., when modulus and r are not coprime

        // Since r is always a power of 2 in Montgomery arithmetic, we need an even modulus
        // to make gcd(modulus, r) > 1
        let even_modulus = 6u32; // Even modulus
        let r = 8u32; // Power of 2

        // gcd(6, 8) = 2 > 1, so basic_mod_inv(6, 8) should return None
        assert!(
            crate::inv::basic_mod_inv(even_modulus, r).is_none(),
            "basic_mod_inv should return None for non-coprime inputs"
        );

        // This should trigger the None path in compute_n_prime_extended_euclidean
        let result = compute_n_prime_extended_euclidean(even_modulus, r);
        assert!(
            result.is_none(),
            "Should return None when basic_mod_inv fails"
        );

        // Test with other even moduli to ensure the None path is consistently hit
        let test_cases = [
            (4u32, 8u32),   // gcd(4, 8) = 4
            (6u32, 12u32),  // gcd(6, 12) = 6
            (10u32, 16u32), // gcd(10, 16) = 2
            (12u32, 8u32),  // gcd(12, 8) = 4
        ];

        for (modulus, r) in test_cases.iter() {
            // Verify basic_mod_inv returns None for these non-coprime pairs
            assert!(
                crate::inv::basic_mod_inv(*modulus, *r).is_none(),
                "basic_mod_inv({}, {}) should return None",
                modulus,
                r
            );

            // Verify compute_n_prime_extended_euclidean returns None
            let result = compute_n_prime_extended_euclidean(*modulus, *r);
            assert!(
                result.is_none(),
                "compute_n_prime_extended_euclidean({}, {}) should return None",
                modulus,
                r
            );
        }
    }
}
