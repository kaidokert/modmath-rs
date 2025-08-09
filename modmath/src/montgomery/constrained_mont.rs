// Constrained Montgomery arithmetic functions
// These work with references to avoid unnecessary copies, following the pattern from exp.rs

use super::basic_mont::NPrimeMethod;
use crate::inv::constrained_mod_inv;

/// Compute N' using trial search method - O(R) complexity (constrained version)
/// Finds N' such that modulus * N' ≡ -1 (mod R)
fn compute_n_prime_trial_search_constrained<T>(modulus: &T, r: &T) -> T
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
            return n_prime;
        }
        n_prime = n_prime.wrapping_add(&T::one());

        // Safety check to avoid infinite loop
        if &n_prime >= r {
            panic!("Could not find N' - should not happen for valid inputs");
        }
    }
}

/// Compute N' using Extended Euclidean Algorithm - O(log R) complexity (constrained version)
/// Finds N' such that modulus * N' ≡ -1 (mod R)
fn compute_n_prime_extended_euclidean_constrained<T>(modulus: &T, r: &T) -> T
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
            r.clone().wrapping_sub(&T::one()) // Handle edge case where inverse is 0
        } else {
            r.clone().wrapping_sub(&modulus_inv)
        }
    } else {
        panic!("Could not find modular inverse - gcd(modulus, R) should be 1 for valid Montgomery parameters");
    }
}

/// Compute N' using Hensel's lifting - O(log R) complexity, optimized for R = 2^k (constrained version)
/// Finds N' such that modulus * N' ≡ -1 (mod R)
fn compute_n_prime_hensels_lifting_constrained<T>(modulus: &T, r: &T, r_bits: usize) -> T
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
    for<'a> &'a T: core::ops::Rem<&'a T, Output = T>,
{
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
        panic!("Hensel lifting failed to produce correct N'");
    }

    n_prime
}

/// Montgomery parameter computation (Constrained)
/// Computes R, R^(-1) mod N, N', and R bit length for Montgomery arithmetic
pub fn constrained_compute_montgomery_params_with_method<T>(
    modulus: &T,
    method: NPrimeMethod,
) -> (T, T, T, usize)
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
        + core::ops::Rem<&'a T, Output = T>,
{
    // Step 1: Find R = 2^k where R > modulus
    let mut r = T::one();
    let mut r_bits = 0usize;

    while &r <= modulus {
        r = r << 1; // r *= 2
        r_bits += 1;
    }

    // Step 2: Compute R^(-1) mod modulus
    let r_inv =
        constrained_mod_inv(r.clone(), modulus).expect("R should always be invertible mod N");

    // Step 3: Compute N' such that N * N' ≡ -1 (mod R) using selected method
    let n_prime = match method {
        NPrimeMethod::TrialSearch => compute_n_prime_trial_search_constrained(modulus, &r),
        NPrimeMethod::ExtendedEuclidean => {
            compute_n_prime_extended_euclidean_constrained(modulus, &r)
        }
        NPrimeMethod::HenselsLifting => {
            compute_n_prime_hensels_lifting_constrained(modulus, &r, r_bits)
        }
    };

    (r, r_inv, n_prime, r_bits)
}

/// Montgomery parameter computation (Constrained) with default method
/// Computes R, R^(-1) mod N, N', and R bit length for Montgomery arithmetic
pub fn constrained_compute_montgomery_params<T>(modulus: &T) -> (T, T, T, usize)
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
        + core::ops::Rem<&'a T, Output = T>,
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

/// Complete Montgomery modular multiplication (Constrained): A * B mod N
pub fn constrained_montgomery_mod_mul<T>(a: T, b: &T, modulus: &T) -> T
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
    let (r, _r_inv, n_prime, r_bits) = constrained_compute_montgomery_params(modulus);
    let a_mont = constrained_to_montgomery(a, modulus, &r);
    let b_mont = constrained_to_montgomery(b.clone(), modulus, &r);
    let result_mont = constrained_montgomery_mul(&a_mont, &b_mont, modulus, &n_prime, r_bits);
    constrained_from_montgomery(result_mont, modulus, &n_prime, r_bits)
}

/// Montgomery-based modular exponentiation (Constrained): base^exponent mod modulus
/// Uses Montgomery arithmetic for efficient repeated multiplication
pub fn constrained_montgomery_mod_exp<T>(mut base: T, exponent: &T, modulus: &T) -> T
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
    // Compute Montgomery parameters
    let (r, _r_inv, n_prime, r_bits) = constrained_compute_montgomery_params(modulus);

    // Reduce base and convert to Montgomery form
    base.rem_assign(modulus);
    base = constrained_to_montgomery(base, modulus, &r);

    // Montgomery form of 1 (the initial result)
    let mut result = constrained_to_montgomery(T::one(), modulus, &r);

    // Copy exponent for manipulation
    let mut exp = exponent.clone();
    let two = T::one().wrapping_add(&T::one());

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
    constrained_from_montgomery(result, modulus, &n_prime, r_bits)
}
