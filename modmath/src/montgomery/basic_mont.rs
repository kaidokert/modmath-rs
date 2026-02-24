// Basic Montgomery arithmetic functions
// These require Copy trait but have minimal constraints

use crate::inv::basic_mod_inv;
use crate::parity::Parity;
use crate::wide_mul::WideMul;

/// Methods for computing N' in Montgomery parameter computation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NPrimeMethod {
    /// Trial search - O(R) complexity, simple but slow for large R.
    /// Only works when R fits in the type (R > N, smallest power of 2).
    TrialSearch,
    /// Extended Euclidean Algorithm - O(log R) complexity.
    /// Only works when R fits in the type (R > N, smallest power of 2).
    ExtendedEuclidean,
    /// Hensel's lifting - O(log R) complexity, optimized for R = 2^k.
    /// Only works when R fits in the type (R > N, smallest power of 2).
    HenselsLifting,
    /// Newton's method - O(log W) complexity, works with R = 2^W (full type width).
    /// Uses wrapping arithmetic so R never needs to be represented explicitly.
    /// This is the method used internally by the wide-REDC path; other methods
    /// are accepted for API compatibility but all use Newton internally.
    #[default]
    Newton,
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
        // N' is already in [0, R) after Hensel lifting (starts at 1, accumulates powers of 2)
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
        NPrimeMethod::Newton => {
            // In this legacy param path (R > N, smallest power of 2),
            // delegate to ExtendedEuclidean which computes -N^{-1} mod R
            // using the same mathematical identity.
            compute_n_prime_extended_euclidean(modulus, r)?
        }
    };

    Some((r, r_inv, n_prime, r_bits))
}

/// Montgomery parameter computation (Basic) with default method
/// Computes R, R^(-1) mod N, N', and R bit length for Montgomery arithmetic
/// Uses the default NPrimeMethod (Newton).
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
        + num_traits::ops::wrapping::WrappingAdd
        + num_traits::ops::wrapping::WrappingSub
        + core::ops::Shr<usize, Output = T>
        + core::ops::Rem<Output = T>
        + crate::parity::Parity,
{
    crate::mul::basic_mod_mul(a, r, modulus)
}

/// Convert from Montgomery form (Basic): (a * R) -> a mod N
/// Uses Montgomery reduction algorithm with legacy R > N semantics.
///
/// **Warning**: This function can overflow for large moduli where m * N exceeds
/// the type width. For overflow-free reduction, use [`wide_from_montgomery`]
/// which uses wide-REDC with R = 2^W.
pub fn basic_from_montgomery<T>(a_mont: T, modulus: T, n_prime: T, r_bits: usize) -> T
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialOrd
        + core::ops::Mul<Output = T>
        + core::ops::Add<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Shl<usize, Output = T>
        + core::ops::BitAnd<Output = T>,
{
    // Montgomery reduction algorithm:
    // Input: a_mont (Montgomery form), N (modulus), N', r_bits
    // 1. mask = 2^r_bits - 1
    // 2. m = ((a_mont & mask) * N') & mask  [only low bits, no expensive modulo!]
    // 3. t = (a_mont + m * N) >> r_bits     [bit shift, no division!]
    // 4. if t >= N then return t - N else return t

    // Fast path for R=1 (r_bits == 0): Montgomery reduction simplifies to conditional subtraction
    if r_bits == 0 {
        return if a_mont >= modulus {
            a_mont - modulus
        } else {
            a_mont
        };
    }

    let mask = (T::one() << r_bits) - T::one(); // mask = 2^r_bits - 1

    // Step 1: m = ((a_mont & mask) * N') & mask
    let m = ((a_mont & mask) * n_prime) & mask;

    // Step 2: t = (a_mont + m * N) >> r_bits
    // WARNING: m * N can overflow for large moduli (m < R, N < R, so
    // m*N can reach R^2). Use wide_from_montgomery for overflow-free reduction.
    let t = (a_mont + m * modulus) >> r_bits;

    // Step 3: Final reduction
    if t >= modulus { t - modulus } else { t }
}

/// Convert from Montgomery form using wide-REDC: (a * R) -> a mod N
///
/// Uses R = 2^W (full type width) semantics with overflow-free reduction.
/// The N' parameter must be computed for R = 2^W (e.g., via Newton's method).
pub fn wide_from_montgomery<T>(a_mont: T, modulus: T, n_prime: T) -> T
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialOrd
        + WideMul
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::WrappingMul
        + num_traits::WrappingSub,
{
    wide_redc(a_mont, T::zero(), modulus, n_prime)
}

/// Montgomery multiplication (Basic): (a * R) * (b * R) -> (a * b * R) mod N
///
/// **Warning**: This building-block function uses the legacy reduction path which
/// can overflow for large moduli (see `basic_from_montgomery` warning). For
/// overflow-free multiplication, use [`basic_montgomery_mod_mul`] which uses
/// wide-REDC internally.
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
        + num_traits::ops::wrapping::WrappingSub
        + crate::parity::Parity,
{
    // Montgomery multiplication algorithm:
    // Input: a_mont, b_mont (both in Montgomery form), modulus N, N', r_bits
    // 1. Compute product = a_mont * b_mont (mod N)
    // 2. Apply Montgomery reduction to get (a * b * R) mod N

    // TODO(Phase 1): Replace mod_mul + from_montgomery with proper Montgomery
    // reduction using WideningMul. Current mod_mul is O(k) double-and-add which
    // defeats the performance purpose of Montgomery, and from_montgomery's
    // m * N intermediate can overflow at key sizes. See ROADMAP.md Phase 1.
    let product = crate::mul::basic_mod_mul(a_mont, b_mont, modulus);
    basic_from_montgomery(product, modulus, n_prime, r_bits)
}

// ---------------------------------------------------------------------------
// Wide-REDC private helpers (R = 2^W, full type width)
// ---------------------------------------------------------------------------

/// Bit width of type T (e.g. 32 for u32, 128 for FixedUInt<u32,4>).
const fn type_bit_width<T>() -> usize {
    core::mem::size_of::<T>() * 8
}

/// Modular doubling: (val + val) mod modulus, handling overflow.
fn mod_double<T>(val: T, modulus: T) -> T
where
    T: Copy + PartialOrd + num_traits::ops::overflowing::OverflowingAdd + num_traits::WrappingSub,
{
    let (doubled, overflow) = val.overflowing_add(&val);
    if overflow || doubled >= modulus {
        doubled.wrapping_sub(&modulus)
    } else {
        doubled
    }
}

/// Newton's method for N' = -N^{-1} mod 2^W.
///
/// Invariant: after each iteration, `modulus * x ≡ 1 (mod 2^precision)`.
/// We return `0 - x` so that `modulus * N' ≡ -1 (mod 2^W)`.
fn compute_n_prime_newton<T>(modulus: T, w: usize) -> T
where
    T: Copy
        + num_traits::One
        + num_traits::Zero
        + num_traits::WrappingMul
        + num_traits::WrappingSub
        + num_traits::WrappingAdd,
{
    let two = T::one().wrapping_add(&T::one());
    let mut x = T::one(); // modulus * 1 ≡ 1 (mod 2) for odd modulus
    let mut precision = 1usize;
    while precision < w {
        // x = x * (2 - modulus * x)   mod 2^(2*precision)
        x = x.wrapping_mul(&two.wrapping_sub(&modulus.wrapping_mul(&x)));
        precision *= 2;
    }
    // N' = -x mod 2^W  (wrapping_sub from 0 gives two's complement negation)
    T::zero().wrapping_sub(&x)
}

/// Compute (val * 2^w) mod N via w modular doublings.
fn mod_exp2<T>(val: T, modulus: T, w: usize) -> T
where
    T: Copy
        + PartialEq
        + PartialOrd
        + num_traits::Zero
        + num_traits::One
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::WrappingSub,
{
    // For modulus == 1, any value mod 1 == 0
    if modulus == T::one() {
        return T::zero();
    }
    let mut result = val;
    for _ in 0..w {
        result = mod_double(result, modulus);
    }
    result
}

/// Compute R mod N = 2^W mod N via W modular doublings starting from 1.
fn compute_r_mod_n<T>(modulus: T, w: usize) -> T
where
    T: Copy
        + PartialEq
        + PartialOrd
        + num_traits::Zero
        + num_traits::One
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::WrappingSub,
{
    mod_exp2(T::one(), modulus, w)
}

/// Compute R^2 mod N via W more modular doublings from (R mod N).
fn compute_r2_mod_n<T>(r_mod_n: T, modulus: T, w: usize) -> T
where
    T: Copy
        + PartialEq
        + PartialOrd
        + num_traits::Zero
        + num_traits::One
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::WrappingSub,
{
    mod_exp2(r_mod_n, modulus, w)
}

/// Accumulate the high half with carries from low-half addition.
///
/// Given the low-half carry `carry1`, adds it to `result` and returns the
/// combined extra-bit flag (true if there was overflow from any addition).
///
/// Invariants maintained:
/// - `result` is in range [0, modulus + 1] on entry (sum of two values < modulus plus carry)
/// - Returns (result + carry1, extra_bit) where extra_bit indicates overflow
fn accumulate_high_half_carry<T>(result: T, carry1: bool, carry2: bool) -> (T, bool)
where
    T: Copy + num_traits::One + num_traits::ops::overflowing::OverflowingAdd,
{
    if carry1 {
        let (r2, carry3) = result.overflowing_add(&T::one());
        (r2, carry2 || carry3)
    } else {
        (result, carry2)
    }
}

/// REDC on a double-width input (t_lo, t_hi).
///
/// Computes  (t_lo + t_hi * 2^W) * R^{-1}  mod N.
///
/// Invariants:
/// - Inputs t_lo, t_hi represent a value T < N * R (product of two values < N)
/// - modulus is odd and non-zero
/// - n_prime = -N^{-1} mod 2^W
fn wide_redc<T>(t_lo: T, t_hi: T, modulus: T, n_prime: T) -> T
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialOrd
        + WideMul
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::WrappingMul
        + num_traits::WrappingSub,
{
    // m = t_lo * N'  (mod 2^W -- wrapping mul gives that for free)
    let m = t_lo.wrapping_mul(&n_prime);

    // (m_lo, m_hi) = m * modulus  (full double-width product)
    let (m_lo, m_hi) = m.wide_mul(&modulus);

    // low half:  t_lo + m_lo  -- the low W bits cancel by construction,
    // we only need the carry.
    let (_discard_lo, carry1) = t_lo.overflowing_add(&m_lo);

    // high half:  t_hi + m_hi + carry1
    let (result, carry2) = t_hi.overflowing_add(&m_hi);
    let (result, extra_bit) = accumulate_high_half_carry(result, carry1, carry2);

    if extra_bit || result >= modulus {
        result.wrapping_sub(&modulus)
    } else {
        result
    }
}

// ---------------------------------------------------------------------------
// Wide-REDC Montgomery multiply helper
// ---------------------------------------------------------------------------

/// Montgomery multiplication using wide REDC: REDC(a_mont * b_mont).
fn wide_montgomery_mul<T>(a_mont: T, b_mont: T, modulus: T, n_prime: T) -> T
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialOrd
        + WideMul
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::WrappingMul
        + num_traits::WrappingSub,
{
    let (lo, hi) = a_mont.wide_mul(&b_mont);
    wide_redc(lo, hi, modulus, n_prime)
}

// ---------------------------------------------------------------------------
// Public API: mod_mul / mod_exp using wide REDC
// ---------------------------------------------------------------------------

/// Complete Montgomery modular multiplication with method selection (Basic): A * B mod N
///
/// Uses wide REDC (R = 2^W) for correct overflow-free Montgomery reduction.
/// Note: The `method` parameter is accepted for API compatibility but all methods
/// use Newton internally since it's designed for R = 2^W.
/// Returns None if modulus is even or zero.
pub fn basic_montgomery_mod_mul_with_method<T>(
    a: T,
    b: T,
    modulus: T,
    _method: NPrimeMethod,
) -> Option<T>
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::WrappingMul
        + num_traits::WrappingAdd
        + num_traits::WrappingSub
        + Parity
        + core::ops::Rem<Output = T>,
{
    // All methods use Newton internally, so delegate to the simple version
    basic_montgomery_mod_mul(a, b, modulus)
}

/// Reduce value modulo modulus.
///
/// **Note**: This function assumes unsigned types. For signed types, the `%`
/// operator may return negative values, which would violate the `[0, modulus)`
/// range required by Montgomery arithmetic.
fn reduce_mod<T>(val: T, modulus: T) -> T
where
    T: Copy + core::ops::Rem<Output = T>,
{
    val % modulus
}

/// Complete Montgomery modular multiplication (Basic): A * B mod N
///
/// Uses wide REDC (R = 2^W) with Newton's method for N'.
/// Inputs are reduced modulo N before conversion to Montgomery form.
/// Returns None if modulus is even or zero.
pub fn basic_montgomery_mod_mul<T>(a: T, b: T, modulus: T) -> Option<T>
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::WrappingMul
        + num_traits::WrappingAdd
        + num_traits::WrappingSub
        + Parity
        + core::ops::Rem<Output = T>,
{
    if modulus == T::zero() || modulus.is_even() {
        return None;
    }
    let w = type_bit_width::<T>();
    let n_prime = compute_n_prime_newton(modulus, w);
    let r_mod_n = compute_r_mod_n(modulus, w);
    let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

    // Reduce inputs to [0, modulus) before Montgomery conversion
    let a_red = reduce_mod(a, modulus);
    let b_red = reduce_mod(b, modulus);

    // Convert to Montgomery form via REDC(a * R^2)
    let (lo, hi) = a_red.wide_mul(&r2_mod_n);
    let a_m = wide_redc(lo, hi, modulus, n_prime);
    let (lo, hi) = b_red.wide_mul(&r2_mod_n);
    let b_m = wide_redc(lo, hi, modulus, n_prime);

    // Multiply in Montgomery domain
    let r_m = wide_montgomery_mul(a_m, b_m, modulus, n_prime);

    // Convert back: REDC(r_m, 0)
    Some(wide_redc(r_m, T::zero(), modulus, n_prime))
}

/// Montgomery-based modular exponentiation with method selection (Basic): base^exponent mod modulus
///
/// Uses wide REDC (R = 2^W) for correct overflow-free Montgomery reduction.
/// Note: The `method` parameter is accepted for API compatibility but all methods
/// use Newton internally since it's designed for R = 2^W.
/// Returns None if modulus is even or zero.
pub fn basic_montgomery_mod_exp_with_method<T>(
    base: T,
    exponent: T,
    modulus: T,
    _method: NPrimeMethod,
) -> Option<T>
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::WrappingMul
        + num_traits::WrappingAdd
        + num_traits::WrappingSub
        + Parity
        + core::ops::Rem<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::ShrAssign<usize>,
{
    // All methods use Newton internally, so delegate to the simple version
    basic_montgomery_mod_exp(base, exponent, modulus)
}

/// Montgomery-based modular exponentiation (Basic): base^exponent mod modulus
///
/// Uses wide REDC (R = 2^W) with Newton's method for N'.
/// The base is reduced modulo N before conversion to Montgomery form.
/// Returns None if modulus is even or zero.
pub fn basic_montgomery_mod_exp<T>(base: T, exponent: T, modulus: T) -> Option<T>
where
    T: Copy
        + num_traits::Zero
        + num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::WrappingMul
        + num_traits::WrappingAdd
        + num_traits::WrappingSub
        + Parity
        + core::ops::Rem<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::ShrAssign<usize>,
{
    if modulus == T::zero() || modulus.is_even() {
        return None;
    }
    let w = type_bit_width::<T>();
    let n_prime = compute_n_prime_newton(modulus, w);
    let r_mod_n = compute_r_mod_n(modulus, w);
    let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

    // 1 in Montgomery form = R mod N (since REDC(1 * R²) = 1 * R² * R⁻¹ = R mod N)
    let one_mont = r_mod_n;

    // Reduce base to [0, modulus) before Montgomery conversion
    let base_red = reduce_mod(base, modulus);
    let (lo, hi) = base_red.wide_mul(&r2_mod_n);
    let mut base_mont = wide_redc(lo, hi, modulus, n_prime);
    let mut result = one_mont;
    let mut exp = exponent;

    while exp > T::zero() {
        if exp.is_odd() {
            result = wide_montgomery_mul(result, base_mont, modulus, n_prime);
        }
        exp >>= 1;
        if exp > T::zero() {
            base_mont = wide_montgomery_mul(base_mont, base_mont, modulus, n_prime);
        }
    }

    Some(wide_redc(result, T::zero(), modulus, n_prime))
}

#[cfg(test)]
mod tests {
    use super::*;

    // -- Old basic_mont tests (param computation, N' methods, etc.) ----------

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

    // -- Wide REDC helper tests (moved from wide_mont.rs) --------------------

    #[test]
    fn test_type_bit_width() {
        assert_eq!(type_bit_width::<u8>(), 8);
        assert_eq!(type_bit_width::<u16>(), 16);
        assert_eq!(type_bit_width::<u32>(), 32);
        assert_eq!(type_bit_width::<u64>(), 64);
    }

    #[test]
    fn test_mod_double() {
        // Normal case: 5+5=10, 10 < 13 -> 10
        assert_eq!(mod_double(5u8, 13), 10);
        // Needs reduction: 8+8=16 >= 13 -> 3
        assert_eq!(mod_double(8u8, 13), 3);
        // Overflow case: 200+200 wraps in u8, result should still be correct
        // 400 mod 201 = 400 - 201 = 199
        // In u8: 200+200 = 144 (wrap), overflow=true -> 144 - 201 wraps to 199
        assert_eq!(mod_double(200u8, 201), 199);
        // Edge: 0 doubled is 0
        assert_eq!(mod_double(0u8, 13), 0);
    }

    #[test]
    fn test_compute_n_prime_newton() {
        // For several small odd moduli, verify N * N' == -1 (mod 2^W)
        let test_moduli: &[u8] = &[3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 255];
        for &n in test_moduli {
            let np = compute_n_prime_newton(n, 8);
            // N * N' should wrap to 0xFF (which is -1 mod 256)
            let product = n.wrapping_mul(np);
            assert_eq!(
                product, 0xFF,
                "n={n}: n*n_prime={product:#04x}, expected 0xFF"
            );
        }

        // Also check u32
        let np32 = compute_n_prime_newton(13u32, 32);
        assert_eq!(13u32.wrapping_mul(np32), u32::MAX);
    }

    #[test]
    fn test_compute_r_mod_n() {
        // R = 2^8 = 256 for u8
        // 256 mod 13 = 256 - 19*13 = 256-247 = 9
        assert_eq!(compute_r_mod_n(13u8, 8), 9);
        // 256 mod 255 = 1
        assert_eq!(compute_r_mod_n(255u8, 8), 1);
        // 256 mod 3 = 1
        assert_eq!(compute_r_mod_n(3u8, 8), 1);

        // R = 2^32 for u32: 2^32 mod 13 = 4294967296 mod 13 = 9
        assert_eq!(compute_r_mod_n(13u32, 32), 9);
    }

    #[test]
    fn test_wide_redc() {
        // With u8, R=256, N=13
        // REDC(35, 0) should give 35 * R^{-1} mod 13
        // R^{-1} mod 13: R=256, 256 mod 13 = 9, inv(9,13)=3 (9*3=27==1 mod 13)
        // So REDC(35,0) = 35*3 mod 13 = 105 mod 13 = 1
        let n: u8 = 13;
        let n_prime = compute_n_prime_newton(n, 8);
        let result = wide_redc(35u8, 0u8, n, n_prime);
        assert_eq!(result, 1);
    }

    // -- Round-trip tests ----------------------------------------------------

    #[test]
    fn test_wide_montgomery_roundtrip() {
        let modulus = 13u8;
        let w = type_bit_width::<u8>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);
        for a in 0u8..13 {
            let (lo, hi) = a.wide_mul(&r2_mod_n);
            let a_m = wide_redc(lo, hi, modulus, n_prime);
            let back = wide_redc(a_m, 0u8, modulus, n_prime);
            assert_eq!(back, a, "roundtrip failed for a={a}");
        }
    }

    #[test]
    fn test_wide_montgomery_roundtrip_u32() {
        let modulus = 0xFFFF_FFF1u32; // large odd u32
        let w = type_bit_width::<u32>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);
        let test_vals = [0u32, 1, 2, 100, modulus - 1, modulus - 2, 0x7FFF_FFFF];
        for &a in &test_vals {
            let (lo, hi) = a.wide_mul(&r2_mod_n);
            let a_m = wide_redc(lo, hi, modulus, n_prime);
            let back = wide_redc(a_m, 0u32, modulus, n_prime);
            assert_eq!(back, a, "roundtrip failed for a={a:#x}");
        }
    }

    // -- Mod mul tests -------------------------------------------------------

    #[test]
    fn test_wide_mod_mul_u8_exhaustive() {
        let modulus = 13u8;
        for a in 0u8..13 {
            for b in 0u8..13 {
                let expected = ((a as u16 * b as u16) % 13) as u8;
                let got = basic_montgomery_mod_mul(a, b, modulus).unwrap();
                assert_eq!(
                    got, expected,
                    "{a}*{b} mod 13: got {got}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn test_wide_mod_mul_u32() {
        // Compare against basic_mod_mul for a range of values with a large modulus
        let modulus = 0xFFFF_FFF1u32;
        let vals = [0u32, 1, 2, 7, 1000, 0x7FFF_FFFF, modulus - 1, modulus - 2];
        for &a in &vals {
            for &b in &vals {
                let expected = crate::mul::basic_mod_mul(a, b, modulus);
                let got = basic_montgomery_mod_mul(a, b, modulus).unwrap();
                assert_eq!(
                    got, expected,
                    "{a:#x}*{b:#x} mod {modulus:#x}: got {got:#x}, expected {expected:#x}"
                );
            }
        }
    }

    // -- Mod exp tests -------------------------------------------------------

    #[test]
    fn test_wide_mod_exp_small() {
        let modulus = 13u8;
        for base in 0u8..13 {
            for exp in 0u8..20 {
                let expected = crate::exp::basic_mod_exp(base, exp, modulus);
                let got = basic_montgomery_mod_exp(base, exp, modulus).unwrap();
                assert_eq!(
                    got, expected,
                    "{base}^{exp} mod 13: got {got}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn test_wide_mod_exp_u32() {
        let modulus = 0xFFFF_FFF1u32;
        // 2^100 mod modulus
        let expected = crate::exp::basic_mod_exp(2u32, 100, modulus);
        let got = basic_montgomery_mod_exp(2, 100, modulus).unwrap();
        assert_eq!(got, expected);

        // 3^1000 mod modulus
        let expected = crate::exp::basic_mod_exp(3u32, 1000, modulus);
        let got = basic_montgomery_mod_exp(3, 1000, modulus).unwrap();
        assert_eq!(got, expected);

        // (modulus-1)^(modulus-2) mod modulus  (Fermat inverse if prime-ish)
        let expected = crate::exp::basic_mod_exp(modulus - 1, modulus - 2, modulus);
        let got = basic_montgomery_mod_exp(modulus - 1, modulus - 2, modulus).unwrap();
        assert_eq!(got, expected);
    }

    #[test]
    fn test_wide_mod_exp_large_u64() {
        // Values that would overflow with the old mod_mul-based approach
        let modulus = 0xFFFF_FFFF_FFFF_FFC5u64; // large odd u64
        let base = 0xDEAD_BEEF_CAFE_BABEu64;
        let exp = 0x1234_5678u64;
        let expected = crate::exp::basic_mod_exp(base, exp, modulus);
        let got = basic_montgomery_mod_exp(base, exp, modulus).unwrap();
        assert_eq!(got, expected);
    }

    // -- Error paths ---------------------------------------------------------

    #[test]
    fn test_wide_params_even_modulus() {
        assert!(basic_montgomery_mod_mul(2u32, 3u32, 4u32).is_none());
        assert!(basic_montgomery_mod_mul(2u32, 3u32, 0u32).is_none());
    }

    #[test]
    fn test_wide_mod_mul_even_modulus() {
        assert!(basic_montgomery_mod_mul(2u32, 3u32, 4u32).is_none());
    }

    #[test]
    fn test_wide_mod_exp_even_modulus() {
        assert!(basic_montgomery_mod_exp(2u32, 3u32, 4u32).is_none());
    }

    // -- FixedUInt test ------------------------------------------------------

    #[test]
    fn test_wide_fixed_bigint() {
        use fixed_bigint::FixedUInt;
        type U128 = FixedUInt<u32, 4>;

        // Pick a 128-bit-ish odd modulus close to type max
        let modulus = !U128::from(0u64) - U128::from(58u64); // 2^128 - 59 (odd)

        // Round-trip test via wide helpers directly
        let w = type_bit_width::<U128>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

        let a = U128::from(0xDEAD_BEEF_u64);
        let (lo, hi) = a.wide_mul(&r2_mod_n);
        let a_m = wide_redc(lo, hi, modulus, n_prime);
        let back = wide_redc(a_m, U128::from(0u64), modulus, n_prime);
        assert_eq!(back, a);

        // Multiplication test: compare with basic_mod_mul
        let b = U128::from(0xCAFE_BABE_u64);
        let expected = crate::mul::basic_mod_mul(a, b, modulus);
        let got = basic_montgomery_mod_mul(a, b, modulus).unwrap();
        assert_eq!(got, expected);

        // Exponentiation test
        let base = U128::from(42u64);
        let exp = U128::from(1000u64);
        let expected = crate::exp::basic_mod_exp(base, exp, modulus);
        let got = basic_montgomery_mod_exp(base, exp, modulus).unwrap();
        assert_eq!(got, expected);
    }

    // -- Input reduction tests -----------------------------------------------

    #[test]
    fn test_input_reduction_mod_mul() {
        // Test that inputs >= modulus are handled correctly via reduction
        let modulus = 13u32;

        // Inputs larger than modulus should be reduced before Montgomery conversion
        let a = 27u32; // 27 mod 13 = 1
        let b = 39u32; // 39 mod 13 = 0
        let got = basic_montgomery_mod_mul(a, b, modulus).unwrap();
        assert_eq!(got, 0, "27 * 39 mod 13 should be 0");

        // Another case: both inputs need reduction
        let a = 100u32; // 100 mod 13 = 9
        let b = 200u32; // 200 mod 13 = 5
        let expected = (9 * 5) % 13; // 45 mod 13 = 6
        let got = basic_montgomery_mod_mul(a, b, modulus).unwrap();
        assert_eq!(got, expected);

        // Edge case: input equals modulus (should reduce to 0)
        let got = basic_montgomery_mod_mul(13u32, 5u32, modulus).unwrap();
        assert_eq!(got, 0, "modulus * 5 mod modulus should be 0");
    }

    #[test]
    fn test_input_reduction_mod_exp() {
        // Test that base >= modulus is handled correctly
        let modulus = 13u32;

        // Base larger than modulus
        let base = 27u32; // 27 mod 13 = 1
        let exp = 100u32;
        let expected = 1u32; // 1^100 = 1
        let got = basic_montgomery_mod_exp(base, exp, modulus).unwrap();
        assert_eq!(got, expected);

        // Another case
        let base = 100u32; // 100 mod 13 = 9
        let exp = 3u32;
        let expected = crate::exp::basic_mod_exp(9u32, 3, modulus); // 729 mod 13 = 1
        let got = basic_montgomery_mod_exp(base, exp, modulus).unwrap();
        assert_eq!(got, expected);
    }

    #[test]
    fn test_wide_from_montgomery() {
        // Test the wide-REDC based from_montgomery function
        let modulus = 13u8;
        let w = type_bit_width::<u8>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

        // Convert to Montgomery form and back for various values
        for a in 0u8..13 {
            let (lo, hi) = a.wide_mul(&r2_mod_n);
            let a_m = wide_redc(lo, hi, modulus, n_prime);
            let back = wide_from_montgomery(a_m, modulus, n_prime);
            assert_eq!(back, a, "wide_from_montgomery roundtrip failed for {a}");
        }
    }
}
