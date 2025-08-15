/// Compute N' using trial search method - O(R) complexity (Strict)
/// Finds N' such that modulus * N' ≡ -1 (mod R)
/// Uses reference-based operations to minimize copying of large integers
/// Returns None if N' cannot be found
pub fn strict_compute_n_prime_trial_search<T>(modulus: &T, r: &T) -> Option<T>
where
    T: num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub,
    for<'a> T: core::ops::RemAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>,
{
    // We need to find N' where modulus * N' ≡ R - 1 (mod R)
    let target = r - &T::one(); // This is -1 mod R

    // Simple trial search for N'
    // TODO: Replace with Extended Euclidean Algorithm for O(log R) complexity instead of O(R)
    // Current implementation is fine for small numbers but inefficient for large moduli
    let mut n_prime = T::one();
    loop {
        // Check if (modulus * n_prime) % r == target
        // Use references to avoid copying large integers
        let product = modulus * &n_prime;
        let remainder = &product % r;
        if remainder == target {
            return Some(n_prime);
        }

        // Increment n_prime: n_prime = n_prime + 1
        // Use overflowing_add for strict arithmetic
        let (incremented, _overflow) = n_prime.overflowing_add(&T::one());
        n_prime = incremented;

        // Safety check to avoid infinite loop
        if &n_prime >= r {
            return None; // Could not find N' - should not happen for valid inputs
        }
    }
}

/// Compute N' using Extended Euclidean Algorithm - O(log R) complexity (strict version)
/// Finds N' such that modulus * N' ≡ -1 (mod R)
/// This is equivalent to N' ≡ -modulus^(-1) (mod R)
/// Returns None if modular inverse cannot be found
fn strict_compute_n_prime_extended_euclidean<T>(modulus: &T, r: &T) -> Option<T>
where
    T: num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub,
    for<'a> T: core::ops::RemAssign<&'a T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::AddAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<T, Output = T>
        + core::ops::Add<&'a T, Output = T>,
{
    // We need to solve: modulus * N' ≡ -1 (mod R)
    // This is equivalent to: N' ≡ -modulus^(-1) (mod R)

    // Use strict_mod_inv to find modulus^(-1) mod R
    let modulus_clone = &T::zero() + modulus; // Clone using references
    if let Some(modulus_inv) = crate::inv::strict_mod_inv(modulus_clone, r) {
        // N' = -modulus^(-1) mod R = R - modulus^(-1) mod R
        if modulus_inv == T::zero() {
            Some(r - &T::one()) // Handle edge case where inverse is 0
        } else {
            Some(r - &modulus_inv)
        }
    } else {
        None // Could not find modular inverse - gcd(modulus, R) should be 1 for valid Montgomery parameters
    }
}

/// Compute N' using Hensel's lifting - O(log R) complexity, optimized for R = 2^k (strict version)
/// Finds N' such that modulus * N' ≡ -1 (mod R)
/// Uses Hensel lifting, 1-bit per step, from 2^1 to 2^r_bits
/// Returns None if Hensel's lifting fails
fn strict_compute_n_prime_hensels_lifting<T>(modulus: &T, r: &T, r_bits: usize) -> Option<T>
where
    T: num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T>,
    for<'a> &'a T: core::ops::Add<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Rem<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    // Hensel's lifting for N' computation when R = 2^k
    // Start with base case: find N' such that modulus * N' ≡ -1 (mod 2)
    // Then iteratively lift to larger powers of 2

    // Base case: modulus * N' ≡ -1 ≡ 1 (mod 2)
    // Since modulus is odd (required for Montgomery), modulus ≡ 1 (mod 2)
    // So we need N' ≡ 1 (mod 2), hence N' starts as 1
    let mut n_prime = T::one();

    // Lift from 2^1 to 2^r_bits using Hensel lifting, 1-bit per step
    for k in 2..=r_bits {
        // We have: modulus * n_prime ≡ -1 (mod 2^(k-1))
        // We want: modulus * n_prime_new ≡ -1 (mod 2^k)

        let target_mod = T::one() << k; // 2^k
        let temp_prod = modulus * &n_prime;
        let (temp_sum, _overflow) = temp_prod.overflowing_add(&T::one());
        let check_val = &temp_sum % &target_mod;

        if check_val != T::zero() {
            // Need to adjust n_prime
            // If modulus * n_prime + 1 = t * 2^(k-1) for odd t, add 2^(k-1) to n_prime
            let prev_power = T::one() << (k - 1); // 2^(k-1)

            if check_val == prev_power {
                let (adjusted, _overflow) = n_prime.overflowing_add(&prev_power);
                n_prime = adjusted;
            }
        }
    }

    // Final check and adjustment to ensure modulus * N' ≡ -1 (mod R)
    let final_check = (modulus * &n_prime) % r;
    let target = r - &T::one(); // -1 mod R

    if final_check != target {
        // This shouldn't happen with correct Hensel lifting, but safety check
        None // Hensel lifting failed to produce correct N'
    } else {
        // Canonicalize N' to [0, R) range
        let canonical_n_prime = &n_prime % r;
        Some(canonical_n_prime)
    }
}

/// Montgomery parameter computation (Strict)
/// Computes R, R^(-1) mod N, N', and R bit length for Montgomery arithmetic using default method
/// Uses reference-based operations to minimize copying of large integers
pub fn strict_compute_montgomery_params<T>(modulus: &T) -> Option<(T, T, T, usize)>
where
    T: num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::AddAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    strict_compute_montgomery_params_with_method(
        modulus,
        crate::montgomery::NPrimeMethod::default(),
    )
}

/// Montgomery parameter computation with method selection (Strict)
/// Computes R, R^(-1) mod N, N', and R bit length for Montgomery arithmetic using specified method
/// Uses reference-based operations to minimize copying of large integers
/// Returns None if parameter computation fails
pub fn strict_compute_montgomery_params_with_method<T>(
    modulus: &T,
    method: crate::montgomery::NPrimeMethod,
) -> Option<(T, T, T, usize)>
where
    T: num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::AddAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    // Step 1: Find R = 2^k where R > modulus
    let mut r = T::one();
    let mut r_bits = 0usize;

    while &r <= modulus {
        r = r << 1; // r *= 2
        r_bits += 1;
    }

    // Create a clone using reference-based addition to avoid moving r
    let r_clone = &r + &T::zero(); // Clone r using references

    // Step 2: Compute R^(-1) mod modulus using strict_mod_inv
    let r_inv = crate::inv::strict_mod_inv(r_clone, modulus)?;

    // Step 3: Compute N' such that N * N' ≡ -1 (mod R) using selected method
    let n_prime = match method {
        crate::montgomery::NPrimeMethod::TrialSearch => {
            strict_compute_n_prime_trial_search(modulus, &r)?
        }
        crate::montgomery::NPrimeMethod::ExtendedEuclidean => {
            strict_compute_n_prime_extended_euclidean(modulus, &r)?
        }
        crate::montgomery::NPrimeMethod::HenselsLifting => {
            strict_compute_n_prime_hensels_lifting(modulus, &r, r_bits)?
        }
    };

    Some((r, r_inv, n_prime, r_bits))
}

/// Montgomery multiplication (Strict): (a * R) * (b * R) -> (a * b * R) mod N
/// Multiplies two values already in Montgomery form and returns result in Montgomery form
pub fn strict_montgomery_mul<T>(a_mont: T, b_mont: &T, modulus: &T, n_prime: &T, r_bits: usize) -> T
where
    T: num_traits::Zero
        + num_traits::One
        + PartialOrd
        + core::ops::Sub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Shl<usize, Output = T>
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub,
    for<'a> T: core::ops::Mul<&'a T, Output = T>
        + core::ops::RemAssign<&'a T>
        + core::ops::Sub<&'a T, Output = T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    // Step 1: Raw multiplication (no modular reduction yet)
    let product = a_mont * b_mont;
    // Step 2: Apply Montgomery reduction to get result in Montgomery form
    strict_from_montgomery(product, modulus, n_prime, r_bits)
}

/// Convert from Montgomery form (Strict): (a * R) -> a mod N
/// Uses Montgomery reduction algorithm with reference-based operations
/// to minimize copying of large integers
pub fn strict_from_montgomery<T>(a_mont: T, modulus: &T, n_prime: &T, r_bits: usize) -> T
where
    T: num_traits::Zero
        + num_traits::One
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Sub<Output = T>
        + num_traits::ops::overflowing::OverflowingAdd,
    for<'a> T: core::ops::Mul<&'a T, Output = T>
        + core::ops::RemAssign<&'a T>
        + core::ops::Sub<&'a T, Output = T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::BitAnd<&'a T, Output = T>,
{
    // Montgomery reduction algorithm:
    // Input: a_mont (Montgomery form), N (modulus), N', r_bits
    // 1. R = 2^r_bits, mask = R - 1 = (1 << r_bits) - 1
    // 2. m = ((a_mont & mask) * N') & mask  [only low bits, no expensive modulo!]
    // 3. t = (a_mont + m * N) >> r_bits     [bit shift, no division!]
    // 4. if t >= N then return t - N else return t

    let mask = (T::one() << r_bits) - T::one(); // mask = 2^r_bits - 1

    // Step 1: m = ((a_mont & mask) * N') & mask
    // Only use low r_bits of a_mont, then mask result to low r_bits
    let a_low = &a_mont & &mask;
    let product = &a_low * n_prime;
    let m = &product & &mask;

    // Step 2: t = (a_mont + m * N) >> r_bits
    // Use overflowing_add for strict arithmetic
    let m_times_n = &m * modulus;
    let (sum, _overflow) = a_mont.overflowing_add(&m_times_n);
    let t = sum >> r_bits; // Divide by R = 2^r_bits using bit shift

    // Step 3: Final reduction
    if &t >= modulus { t - modulus } else { t }
}

/// Convert to Montgomery form (Strict): a -> (a * R) mod N
/// Uses reference-based operations to minimize copying of large integers
pub fn strict_to_montgomery<T>(a: T, modulus: &T, r: &T) -> T
where
    T: PartialOrd
        + num_traits::Zero
        + num_traits::One
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub
        + core::ops::Shr<usize, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T> + core::ops::BitAnd<Output = T>,
{
    crate::mul::strict_mod_mul(a, r, modulus)
}

/// Complete Montgomery modular multiplication with method selection (Strict): A * B mod N
/// Uses reference-based operations throughout to minimize copying of large integers
/// Returns None if Montgomery parameter computation fails
pub fn strict_montgomery_mod_mul_with_method<T>(
    a: T,
    b: &T,
    modulus: &T,
    method: crate::montgomery::NPrimeMethod,
) -> Option<T>
where
    T: num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::AddAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    // Step 1: Compute Montgomery parameters using specified method
    let (r, _r_inv, n_prime, r_bits) =
        strict_compute_montgomery_params_with_method(modulus, method)?;

    // Step 2: Convert inputs to Montgomery form
    let a_mont = strict_to_montgomery(a, modulus, &r);
    let b_clone = b + &T::zero(); // Clone b to pass by value
    let b_mont = strict_to_montgomery(b_clone, modulus, &r);

    // Step 3: Montgomery multiplication in Montgomery domain
    // a_mont * b_mont = (a*R) * (b*R) = a*b*R² mod N
    let product = crate::mul::strict_mod_mul(a_mont, &b_mont, modulus);

    // Apply Montgomery reduction once to get a*b*R mod N (still in Montgomery form)
    let result_mont = strict_from_montgomery(product, modulus, &n_prime, r_bits);

    // Step 4: Convert final result back from Montgomery form to get (a * b) mod N
    Some(strict_from_montgomery(
        result_mont,
        modulus,
        &n_prime,
        r_bits,
    ))
}

/// Complete Montgomery modular multiplication (Strict): A * B mod N
/// Uses reference-based operations throughout to minimize copying of large integers
/// Returns None if Montgomery parameter computation fails
pub fn strict_montgomery_mod_mul<T>(a: T, b: &T, modulus: &T) -> Option<T>
where
    T: num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub
        + for<'a> core::ops::Rem<&'a T, Output = T>,
    for<'a> T: core::ops::RemAssign<&'a T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::AddAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    strict_montgomery_mod_mul_with_method(a, b, modulus, crate::montgomery::NPrimeMethod::default())
}

/// Montgomery-based modular exponentiation with method selection (Strict): base^exponent mod modulus
/// Uses Montgomery arithmetic for efficient repeated multiplication with reference-based operations
/// to minimize copying of large integers
/// Returns None if Montgomery parameter computation fails
pub fn strict_montgomery_mod_exp_with_method<T>(
    mut base: T,
    exponent: &T,
    modulus: &T,
    method: crate::montgomery::NPrimeMethod,
) -> Option<T>
where
    T: num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::ShrAssign<usize>
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub,
    for<'a> T: core::ops::RemAssign<&'a T>
        + core::ops::Rem<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::AddAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    // Step 1: Compute Montgomery parameters using specified method
    let (r, _r_inv, n_prime, r_bits) =
        strict_compute_montgomery_params_with_method(modulus, method)?;

    // Step 2: Reduce base and convert to Montgomery form
    base.rem_assign(modulus); // Reduce base first to ensure base < modulus
    base = strict_to_montgomery(base, modulus, &r);

    // Step 3: Initialize result to Montgomery form of 1
    let mut result = strict_to_montgomery(T::one(), modulus, &r);

    // Step 4: Clone exponent for manipulation to avoid moving the reference
    let mut exp = exponent + &T::zero(); // Clone exponent

    // Step 5: Binary exponentiation staying in Montgomery form
    // Use Montgomery multiplication to stay in Montgomery domain throughout
    while exp > T::zero() {
        // If exponent is odd, multiply result by current base power
        if &exp & &T::one() == T::one() {
            result = strict_montgomery_mul(result, &base, modulus, &n_prime, r_bits);
        }

        // Square the base for next iteration
        exp >>= 1;
        if exp > T::zero() {
            let base_clone = &base + &T::zero(); // Clone base using references
            base = strict_montgomery_mul(base, &base_clone, modulus, &n_prime, r_bits);
        }
    }

    // Step 6: Convert result back from Montgomery form
    Some(strict_from_montgomery(result, modulus, &n_prime, r_bits))
}

/// Montgomery-based modular exponentiation (Strict): base^exponent mod modulus
/// Uses Montgomery arithmetic for efficient repeated multiplication with reference-based operations
/// to minimize copying of large integers
/// Returns None if Montgomery parameter computation fails
pub fn strict_montgomery_mod_exp<T>(base: T, exponent: &T, modulus: &T) -> Option<T>
where
    T: num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::ShrAssign<usize>
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub,
    for<'a> T: core::ops::RemAssign<&'a T>
        + core::ops::Rem<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::AddAssign<&'a T>,
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>
        + core::ops::Mul<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::ops::Div<&'a T, Output = T>
        + core::ops::Sub<T, Output = T>
        + core::ops::Add<&'a T, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    strict_montgomery_mod_exp_with_method(
        base,
        exponent,
        modulus,
        crate::montgomery::NPrimeMethod::default(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strict_compute_n_prime_trial_search() {
        let modulus = 13u32;
        let r = 16u32;
        let n_prime = strict_compute_n_prime_trial_search(&modulus, &r).unwrap();
        assert_eq!(n_prime, 11);
    }

    #[test]
    fn test_strict_with_fixed_bigint() {
        type U256 = fixed_bigint::FixedUInt<u32, 4>;
        let modulus = U256::from(13u32);
        let r = U256::from(16u32);
        let n_prime = strict_compute_n_prime_trial_search(&modulus, &r).unwrap();
        assert_eq!(n_prime, U256::from(11u32));
    }

    #[test]
    fn test_strict_compute_montgomery_params() {
        // Test with our documented example: N = 13
        // Expected: R = 16, R^(-1) = 9, N' = 11, r_bits = 4
        let modulus = 13u32;
        let (r, r_inv, n_prime, r_bits) = strict_compute_montgomery_params(&modulus).unwrap();

        assert_eq!(r, 16);
        assert_eq!(r_inv, 9);
        assert_eq!(n_prime, 11);
        assert_eq!(r_bits, 4);

        // Verify the mathematical properties
        // 1. R * R^(-1) ≡ 1 (mod N)
        assert_eq!((r * r_inv) % 13, 1);

        // 2. N * N' ≡ -1 (mod R) which means N * N' ≡ R - 1 (mod R)
        assert_eq!((13 * n_prime) % r, r - 1);

        // 3. R should be > N and a power of 2
        assert!(r > 13);
        assert_eq!(r, 1u32 << r_bits);
    }

    #[test]
    fn test_strict_compute_montgomery_params_with_method() {
        // Test that the parametrized version produces identical results
        let modulus = 13u32;
        let default_result = strict_compute_montgomery_params(&modulus).unwrap();
        let explicit_trial_result = strict_compute_montgomery_params_with_method(
            &modulus,
            crate::montgomery::NPrimeMethod::TrialSearch,
        )
        .unwrap();

        // Both should produce identical results since TrialSearch is the default
        assert_eq!(default_result, explicit_trial_result);

        // Verify the explicit method call produces correct values
        let (r, r_inv, n_prime, r_bits) = explicit_trial_result;
        assert_eq!(r, 16);
        assert_eq!(r_inv, 9);
        assert_eq!(n_prime, 11);
        assert_eq!(r_bits, 4);

        // Test that all methods produce the same result (for now they all use trial search)
        let euclidean_result = strict_compute_montgomery_params_with_method(
            &modulus,
            crate::montgomery::NPrimeMethod::ExtendedEuclidean,
        )
        .unwrap();
        let hensels_result = strict_compute_montgomery_params_with_method(
            &modulus,
            crate::montgomery::NPrimeMethod::HenselsLifting,
        )
        .unwrap();

        assert_eq!(default_result, euclidean_result);
        assert_eq!(default_result, hensels_result);
    }

    #[test]
    fn test_strict_montgomery_params_with_fixed_bigint() {
        type U256 = fixed_bigint::FixedUInt<u32, 4>;
        let modulus = U256::from(13u32);
        let (r, r_inv, n_prime, r_bits) = strict_compute_montgomery_params(&modulus).unwrap();

        assert_eq!(r, U256::from(16u32));
        assert_eq!(r_inv, U256::from(9u32));
        assert_eq!(n_prime, U256::from(11u32));
        assert_eq!(r_bits, 4);

        // Verify mathematical properties
        let thirteen = U256::from(13u32);
        let one = U256::from(1u32);

        // R * R^(-1) ≡ 1 (mod N)
        assert_eq!((r * r_inv) % thirteen, one);

        // N * N' ≡ -1 (mod R) which means N * N' ≡ R - 1 (mod R)
        assert_eq!((thirteen * n_prime) % r, r - one);

        // R should be > N and a power of 2
        assert!(r > thirteen);
        assert_eq!(r, one << r_bits);
    }

    #[test]
    fn test_strict_from_montgomery() {
        // Test with our documented example: N = 13
        let modulus = 13u32;
        let (_r, _r_inv, n_prime, r_bits) = strict_compute_montgomery_params(&modulus).unwrap();

        // Test conversion from Montgomery form back to normal form
        // We know from the basic tests that 7 -> Montgomery form is 8
        // So 8 -> normal form should be 7
        let mont_form = 8u32; // This is 7 in Montgomery form
        let normal_form = strict_from_montgomery(mont_form, &modulus, &n_prime, r_bits);
        assert_eq!(normal_form, 7u32);

        // Test with 5 -> Montgomery (2) -> back to normal (should be 5)
        let mont_5 = 2u32; // This is 5 in Montgomery form
        let normal_5 = strict_from_montgomery(mont_5, &modulus, &n_prime, r_bits);
        assert_eq!(normal_5, 5u32);

        // Test edge cases
        let mont_0 = 0u32; // 0 in Montgomery form is still 0
        let normal_0 = strict_from_montgomery(mont_0, &modulus, &n_prime, r_bits);
        assert_eq!(normal_0, 0u32);
    }

    #[test]
    fn test_strict_from_montgomery_with_fixed_bigint() {
        type U256 = fixed_bigint::FixedUInt<u32, 4>;
        let modulus = U256::from(13u32);
        let (_r, _r_inv, n_prime, r_bits) = strict_compute_montgomery_params(&modulus).unwrap();

        // Test conversion from Montgomery form back to normal form
        let mont_form = U256::from(8u32); // This is 7 in Montgomery form
        let normal_form = strict_from_montgomery(mont_form, &modulus, &n_prime, r_bits);
        assert_eq!(normal_form, U256::from(7u32));

        let mont_5 = U256::from(2u32); // This is 5 in Montgomery form
        let normal_5 = strict_from_montgomery(mont_5, &modulus, &n_prime, r_bits);
        assert_eq!(normal_5, U256::from(5u32));

        // Test edge case
        let mont_0 = U256::from(0u32);
        let normal_0 = strict_from_montgomery(mont_0, &modulus, &n_prime, r_bits);
        assert_eq!(normal_0, U256::from(0u32));
    }

    #[test]
    fn test_strict_round_trip_conversion() {
        // Test round-trip: normal -> Montgomery -> normal using strict functions
        // We'll need to implement strict_to_montgomery first to complete this test
        // For now, we can verify that from_montgomery works with known Montgomery values

        let modulus = 13u32;
        let (r, _r_inv, n_prime, r_bits) = strict_compute_montgomery_params(&modulus).unwrap();

        // Test all values 0..13 using basic_to_montgomery to get Montgomery form,
        // then strict_from_montgomery to get back to normal
        for i in 0u32..13u32 {
            // Use basic function to convert to Montgomery form
            let mont = crate::montgomery::basic_mont::basic_to_montgomery(i, modulus, r);
            // Use strict function to convert back
            let back = strict_from_montgomery(mont, &modulus, &n_prime, r_bits);
            assert_eq!(
                back, i,
                "Round-trip failed for {}: {} -> {} -> {}",
                i, i, mont, back
            );
        }
    }

    #[test]
    fn test_strict_to_montgomery() {
        // Test with our documented example: N = 13, R = 16
        let modulus = 13u32;
        let (r, _r_inv, _n_prime, _r_bits) = strict_compute_montgomery_params(&modulus).unwrap();

        // From EXAMPLE1_COMPUTE_PARAM.md:
        // 7 -> Montgomery: 7 * 16 mod 13 = 112 mod 13 = 8
        // 5 -> Montgomery: 5 * 16 mod 13 = 80 mod 13 = 2
        assert_eq!(strict_to_montgomery(7u32, &modulus, &r), 8u32);
        assert_eq!(strict_to_montgomery(5u32, &modulus, &r), 2u32);

        // Test edge cases
        assert_eq!(strict_to_montgomery(0u32, &modulus, &r), 0u32); // 0 * R mod N = 0
        assert_eq!(strict_to_montgomery(1u32, &modulus, &r), 3u32); // 1 * 16 mod 13 = 3
    }

    #[test]
    fn test_strict_to_montgomery_with_fixed_bigint() {
        type U256 = fixed_bigint::FixedUInt<u32, 4>;
        let modulus = U256::from(13u32);
        let (r, _r_inv, _n_prime, _r_bits) = strict_compute_montgomery_params(&modulus).unwrap();

        // Test conversions
        assert_eq!(
            strict_to_montgomery(U256::from(7u32), &modulus, &r),
            U256::from(8u32)
        );
        assert_eq!(
            strict_to_montgomery(U256::from(5u32), &modulus, &r),
            U256::from(2u32)
        );
        assert_eq!(
            strict_to_montgomery(U256::from(0u32), &modulus, &r),
            U256::from(0u32)
        );
        assert_eq!(
            strict_to_montgomery(U256::from(1u32), &modulus, &r),
            U256::from(3u32)
        );
    }

    #[test]
    fn test_strict_montgomery_mod_mul() {
        // Test the complete strict Montgomery workflow end-to-end
        let modulus = 13u32;

        // Test: 7 * 5 mod 13 = 9
        let a = 7u32;
        let b = 5u32;
        let result = strict_montgomery_mod_mul(a, &b, &modulus).unwrap();
        assert_eq!(result, 9u32);

        // Verify against regular modular multiplication
        let regular_result = crate::mul::strict_mod_mul(a, &b, &modulus);
        assert_eq!(result, regular_result);

        // Test more cases to ensure correctness
        let test_cases = [
            (0u32, 5u32, 0u32),   // 0 * 5 mod 13 = 0
            (7u32, 0u32, 0u32),   // 7 * 0 mod 13 = 0
            (1u32, 7u32, 7u32),   // 1 * 7 mod 13 = 7
            (7u32, 1u32, 7u32),   // 7 * 1 mod 13 = 7
            (12u32, 12u32, 1u32), // 12 * 12 mod 13 = 144 mod 13 = 1
            (6u32, 9u32, 2u32),   // 6 * 9 mod 13 = 54 mod 13 = 2
        ];

        for (a, b, expected) in test_cases.iter() {
            let result = strict_montgomery_mod_mul(*a, b, &modulus).unwrap();
            assert_eq!(result, *expected, "Failed for {} * {} mod 13", a, b);
        }
    }

    #[test]
    fn test_strict_montgomery_mod_mul_with_fixed_bigint() {
        type U256 = fixed_bigint::FixedUInt<u32, 4>;
        let modulus = U256::from(13u32);

        // Test: 7 * 5 mod 13 = 9
        let a = U256::from(7u32);
        let b = U256::from(5u32);
        let result = strict_montgomery_mod_mul(a, &b, &modulus).unwrap();
        assert_eq!(result, U256::from(9u32));

        // Verify against regular modular multiplication
        let regular_result = crate::mul::strict_mod_mul(U256::from(7u32), &b, &modulus);
        assert_eq!(result, regular_result);
    }

    #[test]
    fn test_strict_complete_round_trip() {
        // Test the complete workflow: to_montgomery -> from_montgomery
        let modulus = 13u32;
        let (r, _r_inv, n_prime, r_bits) = strict_compute_montgomery_params(&modulus).unwrap();

        // Test all values 0..13 for complete round-trip
        for i in 0u32..13u32 {
            let mont = strict_to_montgomery(i, &modulus, &r);
            let back = strict_from_montgomery(mont, &modulus, &n_prime, r_bits);
            assert_eq!(
                back, i,
                "Complete round-trip failed for {}: {} -> {} -> {}",
                i, i, mont, back
            );
        }
    }

    #[test]
    fn test_strict_montgomery_mod_exp() {
        // Test Montgomery-based exponentiation against regular exponentiation
        let modulus = 13u32;

        // Test: 7^5 mod 13 = 16807 mod 13 = 11
        let base = 7u32;
        let exponent = 5u32;
        let montgomery_result = strict_montgomery_mod_exp(base, &exponent, &modulus).unwrap();
        let regular_result = crate::exp::strict_mod_exp(base, &exponent, &modulus);
        assert_eq!(montgomery_result, regular_result);
        assert_eq!(montgomery_result, 11u32);

        // Test edge cases
        assert_eq!(
            strict_montgomery_mod_exp(0u32, &5u32, &modulus).unwrap(),
            0u32
        ); // 0^5 = 0
        assert_eq!(
            strict_montgomery_mod_exp(7u32, &0u32, &modulus).unwrap(),
            1u32
        ); // 7^0 = 1
        assert_eq!(
            strict_montgomery_mod_exp(1u32, &100u32, &modulus).unwrap(),
            1u32
        ); // 1^100 = 1
        assert_eq!(
            strict_montgomery_mod_exp(7u32, &1u32, &modulus).unwrap(),
            7u32
        ); // 7^1 = 7

        // Test more cases
        let test_cases = [
            (2u32, 10u32, 13u32, 10u32), // 2^10 mod 13 = 1024 mod 13 = 10
            (3u32, 4u32, 13u32, 3u32),   // 3^4 mod 13 = 81 mod 13 = 3
            (12u32, 2u32, 13u32, 1u32),  // 12^2 mod 13 = 144 mod 13 = 1
            (5u32, 3u32, 13u32, 8u32),   // 5^3 mod 13 = 125 mod 13 = 8
        ];

        for (base, exp, mod_val, expected) in test_cases.iter() {
            let result = strict_montgomery_mod_exp(*base, exp, mod_val).unwrap();
            assert_eq!(
                result, *expected,
                "Failed for {}^{} mod {}",
                base, exp, mod_val
            );
        }
    }

    #[test]
    fn test_strict_montgomery_mod_exp_with_fixed_bigint() {
        type U256 = fixed_bigint::FixedUInt<u32, 4>;
        let modulus = U256::from(13u32);

        // Test: 7^5 mod 13 = 11
        let base = U256::from(7u32);
        let exponent = U256::from(5u32);
        let montgomery_result = strict_montgomery_mod_exp(base, &exponent, &modulus).unwrap();
        let regular_result = crate::exp::strict_mod_exp(U256::from(7u32), &exponent, &modulus);
        assert_eq!(montgomery_result, regular_result);
        assert_eq!(montgomery_result, U256::from(11u32));

        // Test edge cases
        assert_eq!(
            strict_montgomery_mod_exp(U256::from(0u32), &U256::from(5u32), &modulus).unwrap(),
            U256::from(0u32)
        ); // 0^5 = 0

        assert_eq!(
            strict_montgomery_mod_exp(U256::from(7u32), &U256::from(0u32), &modulus).unwrap(),
            U256::from(1u32)
        ); // 7^0 = 1

        assert_eq!(
            strict_montgomery_mod_exp(U256::from(1u32), &U256::from(100u32), &modulus).unwrap(),
            U256::from(1u32)
        ); // 1^100 = 1
    }

    #[test]
    fn test_strict_montgomery_mod_exp_comprehensive() {
        // Comprehensive test: verify Montgomery exp matches regular exp for all small values
        let modulus = 13u32;

        for base in 0u32..13u32 {
            for exponent in 0u32..10u32 {
                let montgomery_result =
                    strict_montgomery_mod_exp(base, &exponent, &modulus).unwrap();
                let regular_result = crate::exp::strict_mod_exp(base, &exponent, &modulus);
                assert_eq!(
                    montgomery_result, regular_result,
                    "Montgomery vs regular exp mismatch: {}^{} mod 13: {} != {}",
                    base, exponent, montgomery_result, regular_result
                );
            }
        }
    }

    #[test]
    fn test_strict_montgomery_mod_exp_large_exponents() {
        // Test with larger exponents to verify efficiency benefits would apply
        let modulus = 13u32;

        assert_eq!(
            strict_montgomery_mod_exp(2u32, &100u32, &modulus).unwrap(),
            crate::exp::strict_mod_exp(2u32, &100u32, &modulus)
        );

        assert_eq!(
            strict_montgomery_mod_exp(3u32, &1000u32, &modulus).unwrap(),
            crate::exp::strict_mod_exp(3u32, &1000u32, &modulus)
        );

        // Test with very large exponent
        assert_eq!(
            strict_montgomery_mod_exp(7u32, &999999u32, &modulus).unwrap(),
            crate::exp::strict_mod_exp(7u32, &999999u32, &modulus)
        );
    }

    #[test]
    fn test_strict_compute_n_prime_trial_search_failure() {
        // Test case where no N' can be found - this happens when gcd(modulus, R) != 1
        // For Montgomery arithmetic, we need gcd(N, R) = 1, but let's create a scenario
        // where this fails by using an even modulus with R = power of 2

        // Example: N = 4, R = 8 (both even, so gcd(4, 8) = 4 != 1)
        // There should be no N' such that 4 * N' ≡ -1 (mod 8) since 4 and 8 share factors
        let modulus = 4u32;
        let r = 8u32;
        let result = strict_compute_n_prime_trial_search(&modulus, &r);
        assert!(
            result.is_none(),
            "Should return None for invalid modulus-R pair"
        );
    }

    #[test]
    fn test_strict_compute_montgomery_params_failure_even_modulus() {
        // Test Montgomery parameter computation failure with even modulus
        // Montgomery arithmetic requires odd modulus, but the parameter computation
        // should handle this gracefully and return None
        let even_modulus = 4u32;
        let result = strict_compute_montgomery_params(&even_modulus);
        assert!(result.is_none(), "Should return None for even modulus");
    }

    #[test]
    fn test_strict_compute_montgomery_params_failure_with_method() {
        // Test all N' computation methods with invalid inputs
        let invalid_modulus = 4u32; // Even modulus

        // Trial search should fail
        let trial_result = strict_compute_montgomery_params_with_method(
            &invalid_modulus,
            crate::montgomery::NPrimeMethod::TrialSearch,
        );
        assert!(
            trial_result.is_none(),
            "Trial search should fail with even modulus"
        );

        // Extended Euclidean should fail
        let euclidean_result = strict_compute_montgomery_params_with_method(
            &invalid_modulus,
            crate::montgomery::NPrimeMethod::ExtendedEuclidean,
        );
        assert!(
            euclidean_result.is_none(),
            "Extended Euclidean should fail with even modulus"
        );

        // Hensel's lifting should fail
        let hensels_result = strict_compute_montgomery_params_with_method(
            &invalid_modulus,
            crate::montgomery::NPrimeMethod::HenselsLifting,
        );
        assert!(
            hensels_result.is_none(),
            "Hensel's lifting should fail with even modulus"
        );
    }

    #[test]
    fn test_strict_montgomery_mod_mul_parameter_failure() {
        // Test that montgomery_mod_mul returns None when parameter computation fails
        let invalid_modulus = 4u32;
        let a = 2u32;
        let b = 3u32;

        let result = strict_montgomery_mod_mul(a, &b, &invalid_modulus);
        assert!(
            result.is_none(),
            "Montgomery mod_mul should return None for invalid modulus"
        );
    }

    #[test]
    fn test_strict_montgomery_mod_exp_parameter_failure() {
        // Test that montgomery_mod_exp returns None when parameter computation fails
        let invalid_modulus = 4u32;
        let base = 2u32;
        let exponent = 3u32;

        let result = strict_montgomery_mod_exp(base, &exponent, &invalid_modulus);
        assert!(
            result.is_none(),
            "Montgomery mod_exp should return None for invalid modulus"
        );
    }

    #[test]
    fn test_montgomery_reduction_final_subtraction() {
        // Test case designed to trigger t >= modulus in Montgomery reduction
        // We need values where the intermediate result requires final subtraction

        // Using specific values that cause Montgomery reduction to need final subtraction
        let modulus = 15u32; // N = 15
        let (r, _r_inv, n_prime, r_bits) = strict_compute_montgomery_params(&modulus).unwrap();

        // Create Montgomery forms that when processed will result in t >= modulus
        // This requires careful construction to hit the >= modulus branch in from_montgomery
        let a_mont = strict_to_montgomery(14u32, &modulus, &r); // Near maximum value

        // The from_montgomery reduction should trigger the t >= modulus branch
        let result = strict_from_montgomery(a_mont, &modulus, &n_prime, r_bits);
        assert_eq!(result, 14u32);

        // Another case with maximum value
        let max_mont = strict_to_montgomery(14u32, &modulus, &r);
        let result_max = strict_from_montgomery(max_mont, &modulus, &n_prime, r_bits);
        assert_eq!(result_max, 14u32);
    }

    #[test]
    fn test_hensel_lifting_conditional_branches() {
        // Test Hensel's lifting with values that trigger conditional branches
        let modulus = 15u32; // Composite number that's harder for some algorithms

        // Test with different N' computation methods to hit different branches
        let hensels_result = strict_compute_montgomery_params_with_method(
            &modulus,
            crate::montgomery::NPrimeMethod::HenselsLifting,
        );

        // Should still work but may hit different code paths
        assert!(
            hensels_result.is_some(),
            "Hensel's lifting should work with modulus 15"
        );

        // Try with a larger odd composite to stress test the algorithm
        let large_modulus = 255u32; // 255 = 3 * 5 * 17, odd composite
        let result_large = strict_compute_montgomery_params_with_method(
            &large_modulus,
            crate::montgomery::NPrimeMethod::HenselsLifting,
        );
        assert!(
            result_large.is_some(),
            "Should handle larger composite moduli"
        );
    }

    #[test]
    fn test_specific_montgomery_edge_cases() {
        // Test specific mathematical edge cases that may hit uncovered branches

        // Case 1: Modulus where R^(-1) computation is at edge cases
        let edge_modulus = 9u32; // 9 = 3^2, may hit different inv computation paths
        let params = strict_compute_montgomery_params(&edge_modulus);
        assert!(params.is_some(), "Should handle modulus 9");

        if let Some((r, _r_inv, n_prime, r_bits)) = params {
            // Test Montgomery operations with edge values
            let edge_val = 8u32; // Maximum value for this modulus
            let mont_edge = strict_to_montgomery(edge_val, &edge_modulus, &r);
            let back = strict_from_montgomery(mont_edge, &edge_modulus, &n_prime, r_bits);
            assert_eq!(back, edge_val);

            // Test multiplication that might hit t >= modulus branch
            let a = 7u32;
            let b = 8u32;
            let result = strict_montgomery_mod_mul(a, &b, &edge_modulus).unwrap();
            let expected = (a * b) % edge_modulus;
            assert_eq!(result, expected);
        }

        // Case 2: Values designed to stress Montgomery multiplication
        let modulus2 = 21u32; // 21 = 3 * 7
        if let Some((_r, _r_inv, _n_prime, _r_bits)) = strict_compute_montgomery_params(&modulus2) {
            // This Montgomery multiplication may hit the final reduction branch
            let result = strict_montgomery_mod_mul(20u32, &19u32, &modulus2).unwrap();

            let expected = (20u32 * 19u32) % modulus2;
            assert_eq!(result, expected);
        }
    }

    #[test]
    fn test_montgomery_exponentiation_edge_paths() {
        // Test exponentiation with values designed to hit specific loop branches
        let modulus = 11u32; // Prime modulus

        // Test with exponent that exercises different loop paths
        let base = 10u32; // Base near modulus
        let exponent = 10u32; // Exponent that will exercise binary exponentiation loops

        let result = strict_montgomery_mod_exp(base, &exponent, &modulus).unwrap();
        let expected = crate::exp::strict_mod_exp(base, &exponent, &modulus);
        assert_eq!(result, expected);

        // Test with specific values that may hit the exp > T::zero() conditional
        let small_exp = 1u32;
        let result_small = strict_montgomery_mod_exp(base, &small_exp, &modulus).unwrap();
        assert_eq!(result_small, base);

        // Test with larger exponent to exercise more iterations
        let large_exp = 100u32;
        let result_large = strict_montgomery_mod_exp(2u32, &large_exp, &modulus).unwrap();
        let expected_large = crate::exp::strict_mod_exp(2u32, &large_exp, &modulus);
        assert_eq!(result_large, expected_large);
    }
}
