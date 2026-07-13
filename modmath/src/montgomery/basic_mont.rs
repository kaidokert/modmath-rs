// Basic Montgomery arithmetic functions
// These require Copy trait but have minimal constraints

use crate::inv::basic_mod_inv;
use crate::parity::Parity;
use crate::wide_mul::WideMul;
use const_num_traits::Odd;
use subtle::Choice;

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
    /// The default. The wide-REDC path uses Newton unconditionally (it's the
    /// only N' method that fits R = 2^W). On the R>N path (where R is the
    /// smallest power of 2 above N rather than the full type width),
    /// `*_compute_montgomery_params_with_method` maps `Newton` to
    /// `ExtendedEuclidean` since Newton's R = 2^W assumption doesn't hold there.
    #[default]
    Newton,
}

/// Compute N' using trial search method - O(R) complexity
/// Finds N' such that modulus * N' ≡ -1 (mod R)
/// Returns None if N' cannot be found
fn compute_n_prime_trial_search<T>(modulus: T, r: T) -> Option<T>
where
    T: Copy
        + const_num_traits::One
        + core::ops::Add<Output = T>
        + const_num_traits::CheckedMul<Output = T>
        + PartialEq
        + PartialOrd
        + core::ops::Sub<Output = T>
        + core::ops::Rem<Output = T>,
{
    // We need to find N' where modulus * N' ≡ R - 1 (mod R)
    let target = r - T::one(); // This is -1 mod R

    // Simple O(R) trial search for N'. Adequate for small moduli; the
    // ExtendedEuclidean and HenselsLifting NPrimeMethod variants offer
    // O(log R) alternatives for larger moduli.
    let mut n_prime = T::one();
    loop {
        if modulus
            .checked_mul(n_prime)
            .expect(crate::montgomery::OVERFLOW_MSG)
            % r
            == target
        {
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
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::CheckedAdd<Output = T>
        + const_num_traits::CheckedMul<Output = T>
        + PartialEq
        + PartialOrd
        + core::ops::Sub<Output = T>
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
        + const_num_traits::Zero
        + core::ops::Add<Output = T>
        + const_num_traits::CheckedMul<Output = T>
        + const_num_traits::One
        + PartialEq
        + core::ops::Sub<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::Shl<usize, Output = T>,
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
        let check_val = (modulus
            .checked_mul(n_prime)
            .expect(crate::montgomery::OVERFLOW_MSG)
            + T::one())
            % target_mod;

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
    let final_check = modulus
        .checked_mul(n_prime)
        .expect(crate::montgomery::OVERFLOW_MSG)
        % r;
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
        + const_num_traits::Zero
        + core::ops::Mul<Output = T>
        + const_num_traits::One
        + const_num_traits::CheckedAdd<Output = T>
        + const_num_traits::CheckedMul<Output = T>
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Div<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::Add<Output = T>,
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
            // In this R > N param path (smallest power of 2 above N),
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
        + const_num_traits::Zero
        + core::ops::Mul<Output = T>
        + const_num_traits::One
        + const_num_traits::CheckedAdd<Output = T>
        + const_num_traits::CheckedMul<Output = T>
        + PartialEq
        + PartialOrd
        + core::ops::Shl<usize, Output = T>
        + core::ops::Div<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::Add<Output = T>,
{
    basic_compute_montgomery_params_with_method(modulus, NPrimeMethod::default())
}

/// Convert to Montgomery form (Basic): a -> (a * R) mod N
pub fn basic_to_montgomery<T>(a: T, modulus: T, r: T) -> T
where
    T: core::cmp::PartialOrd
        + Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Rem<Output = T>
        + crate::parity::Parity
        + crate::NonCt,
{
    crate::mul::basic_mod_mul(a, r, modulus)
}

/// Convert to Montgomery form (Basic, pre-reduced): a -> (a * R) mod N
/// Precondition: `a < modulus` and `r < modulus`. No `Rem` bound.
pub fn basic_to_montgomery_pr<T>(a: T, modulus: T, r: T) -> T
where
    T: core::cmp::PartialOrd
        + Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + crate::parity::Parity
        + crate::NonCt,
{
    crate::mul::basic_mod_mul_pr(a, r, modulus)
}

/// Convert from Montgomery form (Basic): (a * R) -> a mod N
/// Uses Montgomery reduction algorithm with R > N semantics.
///
/// **Warning**: This function can overflow for large moduli where m * N exceeds
/// the type width. For overflow-free reduction, use [`wide_redc`] (call as
/// `wide_redc(a_mont, T::zero(), modulus, n_prime)`) which uses wide-REDC
/// with R = 2^W.
pub fn basic_from_montgomery<T>(a_mont: T, modulus: T, n_prime: T, r_bits: usize) -> T
where
    T: Copy
        + const_num_traits::One
        + core::ops::Add<Output = T>
        + const_num_traits::CheckedMul<Output = T>
        + PartialOrd
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
    let m = (a_mont & mask)
        .checked_mul(n_prime)
        .expect(crate::montgomery::OVERFLOW_MSG)
        & mask;

    // Step 2: t = (a_mont + m * N) >> r_bits
    // m * N reaches up to R² (m < R, N < R); checked_mul turns a
    // carrier-too-narrow product into a panic rather than a wrapped, wrong
    // reduction. Use wide_redc(a_mont, T::zero(), modulus, n_prime) for
    // overflow-free reduction.
    let t = (a_mont
        + m.checked_mul(modulus)
            .expect(crate::montgomery::OVERFLOW_MSG))
        >> r_bits;

    // Step 3: Final reduction
    if t >= modulus { t - modulus } else { t }
}

/// Montgomery multiplication (Basic): (a * R) * (b * R) -> (a * b * R) mod N
///
/// **Warning**: This building-block function uses the R > N reduction path which
/// can overflow for large moduli (see [`basic_from_montgomery`] warning). For
/// overflow-free multiplication, use [`crate::basic::montgomery::mod_mul`]
/// which uses wide-REDC internally.
pub fn basic_montgomery_mul<T>(a_mont: T, b_mont: T, modulus: T, n_prime: T, r_bits: usize) -> T
where
    T: Copy
        + const_num_traits::Zero
        + core::ops::Add<Output = T>
        + const_num_traits::CheckedMul<Output = T>
        + const_num_traits::One
        + PartialOrd
        + core::ops::Sub<Output = T>
        + core::ops::Rem<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Shl<usize, Output = T>
        + core::ops::BitAnd<Output = T>
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + crate::parity::Parity
        + crate::NonCt,
{
    // Montgomery multiplication algorithm:
    // Input: a_mont, b_mont (both in Montgomery form), modulus N, N', r_bits
    // 1. Compute product = a_mont * b_mont (mod N)
    // 2. Apply Montgomery reduction to get (a * b * R) mod N

    // Note: this R>N path uses double-and-add mod_mul (O(k)), which
    // defeats Montgomery's perf purpose; the m*N intermediate in
    // from_montgomery can also overflow at key sizes. For overflow-free
    // multiplication, route through wide-REDC / CIOS instead.
    let product = crate::mul::basic_mod_mul(a_mont, b_mont, modulus);
    basic_from_montgomery(product, modulus, n_prime, r_bits)
}

/// Montgomery multiplication (Basic, pre-reduced)
/// Precondition: `a_mont < modulus` and `b_mont < modulus`. No `Rem` bound.
pub fn basic_montgomery_mul_pr<T>(a_mont: T, b_mont: T, modulus: T, n_prime: T, r_bits: usize) -> T
where
    T: Copy
        + const_num_traits::Zero
        + core::ops::Add<Output = T>
        + const_num_traits::CheckedMul<Output = T>
        + const_num_traits::One
        + PartialOrd
        + core::ops::Sub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Shl<usize, Output = T>
        + core::ops::BitAnd<Output = T>
        + const_num_traits::ops::wrapping::WrappingAdd<Output = T>
        + const_num_traits::ops::wrapping::WrappingSub<Output = T>
        + crate::parity::Parity
        + crate::NonCt,
{
    let product = crate::mul::basic_mod_mul_pr(a_mont, b_mont, modulus);
    basic_from_montgomery(product, modulus, n_prime, r_bits)
}

// ---------------------------------------------------------------------------
// Wide-REDC private helpers (R = 2^W, full type width)
// ---------------------------------------------------------------------------

/// Modular doubling: (val + val) mod modulus, handling overflow.
///
/// **Variable-time.** Branches on the `overflow || doubled >= modulus`
/// magnitude check; used by [`compute_r_mod_n`] / [`compute_r2_mod_n`]
/// for the Nct precompute path where the modulus is public. For the Ct
/// path (secret modulus), see [`mod_double_ct`].
fn mod_double<T>(val: T, modulus: T) -> T
where
    T: Copy
        + PartialOrd
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>,
{
    let (doubled, overflow) = val.overflowing_add(val);
    if overflow || doubled >= modulus {
        doubled.wrapping_sub(modulus)
    } else {
        doubled
    }
}

/// CT modular doubling: `(val + val) mod modulus`, no value-dependent
/// branches. Always computes both the post-sub and pre-sub candidates,
/// selects via `subtle::Choice`.
///
/// Used by [`compute_r_mod_n_ct`] / [`compute_r2_mod_n_ct`] for the Ct
/// precompute path called from [`Field::new_odd_ct`](crate::Field::new_odd_ct).
fn mod_double_ct<T>(val: T, modulus: T) -> T
where
    T: Copy
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + subtle::ConditionallySelectable
        + subtle::ConstantTimeLess,
{
    let doubled = val.wrapping_add(val);
    // Wraparound ⇔ sum < operand. subtle's ct_lt is pure mask
    // arithmetic on every ISA; overflowing_add's carry flag lowers to
    // an equality-branch chain on targets without a flags register
    // (riscv32).
    let overflow = doubled.ct_lt(&val);
    let sub_result = doubled.wrapping_sub(modulus);
    let doubled_ge_modulus = !doubled.ct_lt(&modulus);
    let needs_sub = overflow | doubled_ge_modulus;
    T::conditional_select(&doubled, &sub_result, needs_sub)
}

/// Newton's method for N' = -N^{-1} mod 2^W.
///
/// Invariant: after each iteration, `modulus * x ≡ 1 (mod 2^precision)`.
/// We return `0 - x` so that `modulus * N' ≡ -1 (mod 2^W)`.
pub fn compute_n_prime_newton<T>(modulus: T, w: usize) -> T
where
    T: Copy
        + const_num_traits::One
        + const_num_traits::Zero
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + const_num_traits::WrappingAdd<Output = T>,
{
    let two = T::one().wrapping_add(T::one());
    let mut x = T::one(); // modulus * 1 ≡ 1 (mod 2) for odd modulus
    let mut precision = 1usize;
    while precision < w {
        // x = x * (2 - modulus * x)   mod 2^(2*precision)
        x = x.wrapping_mul(two.wrapping_sub(modulus.wrapping_mul(x)));
        precision *= 2;
    }
    // N' = -x mod 2^W  (wrapping_sub from 0 gives two's complement negation)
    T::zero().wrapping_sub(x)
}

/// Compute (val * 2^w) mod N via w modular doublings.
///
/// **Variable-time** (delegates to [`mod_double`]). Used for the Nct
/// precompute path; for the Ct path see [`mod_exp2_ct`].
fn mod_exp2<T>(val: T, modulus: T, w: usize) -> T
where
    T: Copy
        + PartialEq
        + PartialOrd
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>,
{
    // For modulus == 1, any value mod 1 == 0
    if modulus == T::one() {
        return T::zero();
    }
    // Widen the start value to the modulus width so mod_double's
    // overflow flag fires at bit W, not at the narrow value's width.
    let mut result = modulus.wrapping_sub(modulus).overflowing_add(val).0;
    for _ in 0..w {
        result = mod_double(result, modulus);
    }
    result
}

/// CT modular doubling iterated `w` times: `val · 2^w mod modulus`.
/// Constant-time over `val` and `modulus`. The `modulus == 1` edge
/// case (every value reduces to 0) is handled branchlessly via a
/// final `conditional_select` rather than the variable-time early
/// return in [`mod_exp2`].
fn mod_exp2_ct<T>(val: T, modulus: T, w: usize) -> T
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + subtle::ConditionallySelectable
        + subtle::ConstantTimeEq
        + subtle::ConstantTimeLess,
{
    // Widen the start value to the modulus width so mod_double_ct's
    // overflow flag fires at bit W, not at the narrow value's width.
    let mut result = modulus.wrapping_sub(modulus).wrapping_add(val);
    for _ in 0..w {
        result = mod_double_ct(result, modulus);
    }
    let m_is_one = modulus.ct_eq(&T::one());
    T::conditional_select(&result, &T::zero(), m_is_one)
}

/// Compute R mod N = 2^W mod N via W modular doublings starting from 1.
///
/// **Variable-time.** Used by [`Field::new_odd`](crate::Field::new_odd) for the Nct precompute
/// path (public modulus); the CT sibling [`compute_r_mod_n_ct`] is
/// used by [`Field::new_odd_ct`](crate::Field::new_odd_ct) when the modulus is secret.
pub fn compute_r_mod_n<T>(modulus: T, w: usize) -> T
where
    T: Copy
        + PartialEq
        + PartialOrd
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>,
{
    mod_exp2(T::one(), modulus, w)
}

/// CT version of [`compute_r_mod_n`]: `2^W mod modulus` with no
/// value-dependent branches. Used by [`Field::new_odd_ct`](crate::Field::new_odd_ct) for the
/// secret-modulus precompute path.
pub fn compute_r_mod_n_ct<T>(modulus: T, w: usize) -> T
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + subtle::ConditionallySelectable
        + subtle::ConstantTimeLess,
{
    mod_exp2_ct(T::one(), modulus, w)
}

/// Compute R^2 mod N via W more modular doublings from (R mod N).
///
/// **Variable-time.** See [`compute_r2_mod_n_ct`] for the CT sibling.
pub fn compute_r2_mod_n<T>(r_mod_n: T, modulus: T, w: usize) -> T
where
    T: Copy
        + PartialEq
        + PartialOrd
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>,
{
    mod_exp2(r_mod_n, modulus, w)
}

/// CT version of [`compute_r2_mod_n`]: `R² mod modulus` via W modular
/// doublings from `R mod N`. No value-dependent branches.
pub fn compute_r2_mod_n_ct<T>(r_mod_n: T, modulus: T, w: usize) -> T
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + subtle::ConditionallySelectable
        + subtle::ConstantTimeLess,
{
    mod_exp2_ct(r_mod_n, modulus, w)
}

/// Accumulate the high half with carries from low-half addition.
///
/// Given the low-half carry `carry1`, adds it to `result` and returns the
/// combined extra-bit flag (true if there was overflow from any addition).
///
/// Invariants maintained:
/// - `result` is in range [0, modulus + 1] on entry (sum of two values < modulus plus carry)
/// - Returns (result + carry1, extra_bit) where extra_bit indicates overflow
///
/// This is the variable-time path: branches on `carry1` and on `carry3`
/// via `||` short-circuit. Used by the NCT REDC functions. For the CT
/// REDC path, see `accumulate_high_half_carry_ct`.
fn accumulate_high_half_carry<T>(result: T, carry1: bool, carry2: bool) -> (T, bool)
where
    T: const_num_traits::One + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>,
{
    if carry1 {
        let (r2, carry3) = result.overflowing_add(T::one());
        (r2, carry2 || carry3)
    } else {
        (result, carry2)
    }
}

/// Constant-time analog of [`accumulate_high_half_carry`].
///
/// Same semantics — adds `carry1` to `result` and returns the combined
/// extra-bit flag — but eliminates the value-dependent branches. The
/// addition is always performed; the result is selected branchlessly via
/// `subtle::ConditionallySelectable`. Carries stay `Choice` end to end:
/// bitwise `&` / `|` never short-circuit, and no `bool` materializes for
/// the optimizer to re-branch on.
///
/// Called by the CT REDC functions (`wide_redc_ct`, `strict_wide_redc_ct`).
fn accumulate_high_half_carry_ct<T>(result: T, carry1: Choice, carry2: Choice) -> (T, Choice)
where
    T: const_num_traits::One
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::CtIsZero
        + subtle::ConditionallySelectable,
{
    // Always compute the addition; branchlessly choose whether to keep
    // it. `+1` wraps ⇔ the sum is zero — no overflow flag needed.
    let r2 = result.wrapping_add(T::one());
    let carry3 = r2.ct_is_zero();
    let chosen = T::conditional_select(&result, &r2, carry1);
    (chosen, carry2 | (carry1 & carry3))
}

/// REDC on a double-width input (t_lo, t_hi) — variable-time.
///
/// Computes  (t_lo + t_hi * 2^W) * R^{-1}  mod N.
///
/// Final reduction is a predicted branch. Timing leaks operand magnitude.
/// Use [`wide_redc_ct`] in constant-time-sensitive contexts (signing keys,
/// scalar multiplication on secret data). For verify/public-key paths the
/// branched version is significantly faster — on Cortex-M3 the branchless
/// finalize is ~6× slower because the cond-sub branch is highly predictable
/// while the branchless mask construction does limb-wise work every call.
///
/// Invariants:
/// - Inputs t_lo, t_hi represent a value T < N * R (product of two values < N)
/// - modulus is odd and non-zero
/// - n_prime = -N^{-1} mod 2^W
pub fn wide_redc<T>(t_lo: T, t_hi: T, modulus: T, n_prime: T) -> T
where
    T: Copy
        + const_num_traits::One
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingSub<Output = T>,
{
    // m = t_lo * N'  (mod 2^W -- wrapping mul gives that for free)
    let m = t_lo.wrapping_mul(n_prime);

    // (m_lo, m_hi) = m * modulus  (full double-width product)
    let (m_lo, m_hi) = m.wide_mul(&modulus);

    // low half:  t_lo + m_lo  -- the low W bits cancel by construction,
    // we only need the carry.
    let (_discard_lo, carry1) = t_lo.overflowing_add(m_lo);

    // high half:  t_hi + m_hi + carry1
    let (result, carry2) = t_hi.overflowing_add(m_hi);
    let (result, extra_bit) = accumulate_high_half_carry(result, carry1, carry2);

    if extra_bit || result >= modulus {
        result.wrapping_sub(modulus)
    } else {
        result
    }
}

/// REDC on a double-width input — variable-time, reference-based inputs.
///
/// Same algorithm as [`wide_redc`] but takes all operands by reference,
/// avoiding the per-call value copy. For `Copy` types the compiler
/// register-passes either way; the by-ref form matters for non-`Copy`
/// bigint backends where each value pass would clone.
///
/// Takes operands by reference; `T: Copy` is still required (the body
/// deref-copies into the word-level kernel).
pub fn strict_wide_redc<T>(t_lo: &T, t_hi: &T, modulus: &T, n_prime: &T) -> T
where
    T: Copy
        + const_num_traits::One
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingSub<Output = T>,
{
    let m = (*t_lo).wrapping_mul(*n_prime);
    let (m_lo, m_hi) = m.wide_mul(modulus);
    let (_discard_lo, carry1) = (*t_lo).overflowing_add(m_lo);
    let (result, carry2) = (*t_hi).overflowing_add(m_hi);
    let (result, extra_bit) = accumulate_high_half_carry(result, carry1, carry2);

    if extra_bit || &result >= modulus {
        result.wrapping_sub(*modulus)
    } else {
        result
    }
}

/// REDC on a double-width input (t_lo, t_hi) — constant-time finalize.
///
/// Same algorithm as [`wide_redc`] but performs the final reduction step
/// branchlessly via `subtle::ConditionallySelectable` and
/// `subtle::ConstantTimeLess`, removing the operand-magnitude side-channel
/// on the `result >= modulus` comparison.
///
/// Use this in CT-sensitive paths (private-key operations, secret scalar
/// multiplication). On platforms with branch prediction the branchless
/// finalize is a few percent slower than [`wide_redc`]; on simple cores
/// (Cortex-M0/M3) the gap is much larger (~5–6×). Choose at the call site.
pub fn wide_redc_ct<T>(t_lo: T, t_hi: T, modulus: T, n_prime: T) -> T
where
    T: Copy
        + const_num_traits::One
        + WideMul
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::CtIsZero
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + subtle::ConditionallySelectable
        + subtle::ConstantTimeLess,
{
    let m = t_lo.wrapping_mul(n_prime);
    let (m_lo, m_hi) = m.wide_mul(&modulus);
    // Carries via ct_lt (wraparound ⇔ sum < operand): overflowing_add's
    // flag lowers to an equality-branch chain on targets without a
    // flags register (riscv32).
    let sum_lo = t_lo.wrapping_add(m_lo);
    let carry1 = sum_lo.ct_lt(&t_lo);
    let result = t_hi.wrapping_add(m_hi);
    let carry2 = result.ct_lt(&t_hi);
    let (result, extra_bit) = accumulate_high_half_carry_ct(result, carry1, carry2);

    // Branchless final reduction: needs_sub = extra_bit | !(result < modulus)
    let sub_result = result.wrapping_sub(modulus);
    let result_lt_modulus = result.ct_lt(&modulus);
    let needs_sub = extra_bit | !result_lt_modulus;
    T::conditional_select(&result, &sub_result, needs_sub)
}

/// REDC on a double-width input — constant-time finalize, reference-based inputs.
///
/// Same algorithm as [`wide_redc_ct`] but takes all operands by reference,
/// avoiding the per-call value copy. `T: Copy` is still required — see
/// [`strict_wide_redc`].
pub fn strict_wide_redc_ct<T>(t_lo: &T, t_hi: &T, modulus: &T, n_prime: &T) -> T
where
    T: Copy
        + const_num_traits::One
        + WideMul
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::CtIsZero
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + subtle::ConditionallySelectable
        + subtle::ConstantTimeLess,
{
    let m = (*t_lo).wrapping_mul(*n_prime);
    let (m_lo, m_hi) = m.wide_mul(modulus);
    // Carry idiom rationale in `wide_redc_ct`.
    let sum_lo = (*t_lo).wrapping_add(m_lo);
    let carry1 = sum_lo.ct_lt(t_lo);
    let result = (*t_hi).wrapping_add(m_hi);
    let carry2 = result.ct_lt(t_hi);
    let (result, extra_bit) = accumulate_high_half_carry_ct(result, carry1, carry2);

    let sub_result = result.wrapping_sub(*modulus);
    let result_lt_modulus = result.ct_lt(modulus);
    let needs_sub = extra_bit | !result_lt_modulus;
    T::conditional_select(&result, &sub_result, needs_sub)
}

// ---------------------------------------------------------------------------
// Wide-REDC Montgomery multiply helper
// ---------------------------------------------------------------------------

/// Montgomery multiplication using wide REDC: REDC(a_mont * b_mont).
pub fn wide_montgomery_mul<T>(a_mont: T, b_mont: T, modulus: T, n_prime: T) -> T
where
    T: Copy
        + const_num_traits::One
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingSub<Output = T>,
{
    let (lo, hi) = a_mont.wide_mul(&b_mont);
    wide_redc(lo, hi, modulus, n_prime)
}

/// Montgomery multiplication using wide REDC — constant-time finalize.
///
/// Same shape as [`wide_montgomery_mul`] but routes through [`wide_redc_ct`].
/// Use in CT-sensitive paths.
pub fn wide_montgomery_mul_ct<T>(a_mont: T, b_mont: T, modulus: T, n_prime: T) -> T
where
    T: Copy
        + const_num_traits::One
        + WideMul
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::CtIsZero
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + subtle::ConditionallySelectable
        + subtle::ConstantTimeLess,
{
    let (lo, hi) = a_mont.wide_mul(&b_mont);
    wide_redc_ct(lo, hi, modulus, n_prime)
}

/// Montgomery multiplication using wide REDC — reference-based inputs.
///
/// Same shape as [`wide_montgomery_mul`] but takes operands by reference,
/// avoiding per-call copies on non-Copy bigint backends. Delegates to
/// [`strict_wide_redc`] for the reduction.
pub fn strict_wide_montgomery_mul<T>(a_mont: &T, b_mont: &T, modulus: &T, n_prime: &T) -> T
where
    T: Copy
        + const_num_traits::One
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingSub<Output = T>,
{
    let (lo, hi) = a_mont.wide_mul(b_mont);
    strict_wide_redc(&lo, &hi, modulus, n_prime)
}

/// Montgomery multiplication using wide REDC — constant-time finalize,
/// reference-based inputs.
///
/// Same shape as [`wide_montgomery_mul_ct`] but takes operands by reference.
/// Delegates to [`strict_wide_redc_ct`].
pub fn strict_wide_montgomery_mul_ct<T>(a_mont: &T, b_mont: &T, modulus: &T, n_prime: &T) -> T
where
    T: Copy
        + const_num_traits::One
        + WideMul
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::CtIsZero
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + subtle::ConditionallySelectable
        + subtle::ConstantTimeLess,
{
    let (lo, hi) = a_mont.wide_mul(b_mont);
    strict_wide_redc_ct(&lo, &hi, modulus, n_prime)
}

/// Wide multiply-accumulate: `(acc_lo, acc_hi) += a * b`.
///
/// Produces `a * b` as a double-width pair and folds it into the
/// caller-held accumulator. Result is **not reduced** — call
/// [`wide_redc`] once at the end of accumulation.
///
/// # Bound
///
/// Caller must keep the accumulator within `2·WIDTH(T)`. With Mont-form
/// `a, b ∈ [0, q)`, summing `N` products requires `N ≤ R/q` where
/// `R = 2^WIDTH(T)`. Out-of-bound use silently wraps the high word.
pub fn wide_montgomery_mul_acc<T>(acc_lo: T, acc_hi: T, a: T, b: T) -> (T, T)
where
    T: Copy
        + const_num_traits::One
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>,
{
    let (m_lo, m_hi) = a.wide_mul(&b);
    let (new_lo, carry1) = acc_lo.overflowing_add(m_lo);
    let (sum_hi, _) = acc_hi.overflowing_add(m_hi);
    // Carry-out is discarded under the documented `N ≤ R/q` bound.
    let new_hi = if carry1 {
        let (r, _) = sum_hi.overflowing_add(T::one());
        r
    } else {
        sum_hi
    };
    (new_lo, new_hi)
}

/// Wide multiply-accumulate — constant-time carry propagation.
///
/// Same shape as [`wide_montgomery_mul_acc`] but the carry select is
/// branchless (`ct_lt`-derived `Choice` into `conditional_select`), so
/// the per-term cost has no value-dependent branches. Use in
/// CT-sensitive paths; pair with [`wide_redc_ct`] for the final
/// reduction.
///
/// Same bound as [`wide_montgomery_mul_acc`] — caller-verified `N ≤ R/q`.
pub fn wide_montgomery_mul_acc_ct<T>(acc_lo: T, acc_hi: T, a: T, b: T) -> (T, T)
where
    T: Copy
        + const_num_traits::One
        + WideMul
        + const_num_traits::WrappingAdd<Output = T>
        + subtle::ConstantTimeLess
        + subtle::ConditionallySelectable,
{
    let (m_lo, m_hi) = a.wide_mul(&b);
    // Carry idiom rationale in `wide_redc_ct`; the two adds whose
    // flags were discarded are plain wrapping adds.
    let new_lo = acc_lo.wrapping_add(m_lo);
    let carry1 = new_lo.ct_lt(&acc_lo);
    let sum_hi = acc_hi.wrapping_add(m_hi);
    let r = sum_hi.wrapping_add(T::one());
    let new_hi = T::conditional_select(&sum_hi, &r, carry1);
    (new_lo, new_hi)
}

/// Wide multiply-accumulate — reference-based inputs.
///
/// Same shape as [`wide_montgomery_mul_acc`] but takes operands by
/// reference, for non-`Copy` bigint backends. Pair with
/// [`strict_wide_redc`] for the final reduction.
pub fn strict_wide_montgomery_mul_acc<T>(acc_lo: &T, acc_hi: &T, a: &T, b: &T) -> (T, T)
where
    T: Copy
        + const_num_traits::One
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>,
{
    let (m_lo, m_hi) = a.wide_mul(b);
    let (new_lo, carry1) = (*acc_lo).overflowing_add(m_lo);
    let (sum_hi, _) = (*acc_hi).overflowing_add(m_hi);
    let new_hi = if carry1 {
        let (r, _) = sum_hi.overflowing_add(T::one());
        r
    } else {
        sum_hi
    };
    (new_lo, new_hi)
}

/// Wide multiply-accumulate — constant-time carry, reference-based inputs.
///
/// Same shape as [`wide_montgomery_mul_acc_ct`] but reference-based.
/// Pair with [`strict_wide_redc_ct`] for the final reduction.
pub fn strict_wide_montgomery_mul_acc_ct<T>(acc_lo: &T, acc_hi: &T, a: &T, b: &T) -> (T, T)
where
    T: Copy
        + const_num_traits::One
        + WideMul
        + const_num_traits::WrappingAdd<Output = T>
        + subtle::ConstantTimeLess
        + subtle::ConditionallySelectable,
{
    let (m_lo, m_hi) = a.wide_mul(b);
    // Carry idiom rationale in `wide_redc_ct`.
    let new_lo = (*acc_lo).wrapping_add(m_lo);
    let carry1 = new_lo.ct_lt(acc_lo);
    let sum_hi = (*acc_hi).wrapping_add(m_hi);
    let r = sum_hi.wrapping_add(T::one());
    let new_hi = T::conditional_select(&sum_hi, &r, carry1);
    (new_lo, new_hi)
}

// ---------------------------------------------------------------------------
// Public API: mod_mul / mod_exp using wide REDC
// ---------------------------------------------------------------------------

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
pub fn basic_montgomery_mod_mul_odd<T>(a: T, b: T, modulus: Odd<T>) -> T
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + core::ops::Rem<Output = T>
        + const_num_traits::BitsPrecision,
{
    let m = modulus.get();
    basic_montgomery_mod_mul_pr_odd(reduce_mod(a, m), reduce_mod(b, m), modulus)
}

/// Returns None if modulus is even or zero. Thin wrapper around
/// [`basic_montgomery_mod_mul_odd`].
pub fn basic_montgomery_mod_mul<T>(a: T, b: T, modulus: T) -> Option<T>
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + Parity
        + core::ops::Rem<Output = T>
        + const_num_traits::BitsPrecision,
{
    Odd::new(modulus).map(|m| basic_montgomery_mod_mul_odd(a, b, m))
}

/// Complete Montgomery modular multiplication (Basic, pre-reduced,
/// proven-odd modulus): A * B mod N.
///
/// **Infallible.** The `Odd<T>` typestate carries the "modulus is odd and
/// nonzero" precondition; no internal `Option` plumbing or runtime parity
/// check. Use [`basic_montgomery_mod_mul_pr`] if the proof has to be done
/// at runtime.
///
/// Precondition (unchanged from the `Option`-returning sibling): `a < modulus`
/// and `b < modulus`.
pub fn basic_montgomery_mod_mul_pr_odd<T>(a: T, b: T, modulus: Odd<T>) -> T
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + const_num_traits::BitsPrecision,
{
    let modulus = modulus.get();
    let w = modulus.bits_precision() as usize;
    let n_prime = compute_n_prime_newton(modulus, w);
    let r_mod_n = compute_r_mod_n(modulus, w);
    let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

    // Convert to Montgomery form via REDC(a * R^2)
    let (lo, hi) = a.wide_mul(&r2_mod_n);
    let a_m = wide_redc(lo, hi, modulus, n_prime);
    let (lo, hi) = b.wide_mul(&r2_mod_n);
    let b_m = wide_redc(lo, hi, modulus, n_prime);

    // Multiply in Montgomery domain
    let r_m = wide_montgomery_mul(a_m, b_m, modulus, n_prime);

    // Convert back: REDC(r_m, 0)
    wide_redc(r_m, T::zero(), modulus, n_prime)
}

/// Complete Montgomery modular multiplication (Basic, pre-reduced): A * B mod N
///
/// Precondition: `a < modulus` and `b < modulus`. No `Rem` bound. Returns
/// None only if modulus is even or zero. Thin wrapper around
/// [`basic_montgomery_mod_mul_pr_odd`] that performs the parity proof at
/// runtime — prefer the `_odd` form to keep the panic path out of the
/// linked binary.
pub fn basic_montgomery_mod_mul_pr<T>(a: T, b: T, modulus: T) -> Option<T>
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + Parity
        + const_num_traits::BitsPrecision,
{
    Odd::new(modulus).map(|m| basic_montgomery_mod_mul_pr_odd(a, b, m))
}

/// Montgomery-based modular exponentiation (Basic, proven-odd modulus):
/// base^exponent mod modulus. **Infallible.**
pub fn basic_montgomery_mod_exp_odd<T>(base: T, exponent: T, modulus: Odd<T>) -> T
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + Parity
        + core::ops::Rem<Output = T>
        + core::ops::ShrAssign<usize>
        + const_num_traits::BitsPrecision,
{
    let m = modulus.get();
    basic_montgomery_mod_exp_pr_odd(reduce_mod(base, m), exponent, modulus)
}

/// Montgomery-based modular exponentiation (Basic): base^exponent mod modulus
///
/// Uses wide REDC (R = 2^W) with Newton's method for N'.
/// The base is reduced modulo N before conversion to Montgomery form.
/// Returns None if modulus is even or zero. Thin wrapper around
/// [`basic_montgomery_mod_exp_odd`].
pub fn basic_montgomery_mod_exp<T>(base: T, exponent: T, modulus: T) -> Option<T>
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + Parity
        + core::ops::Rem<Output = T>
        + core::ops::ShrAssign<usize>
        + const_num_traits::BitsPrecision,
{
    Odd::new(modulus).map(|m| basic_montgomery_mod_exp_odd(base, exponent, m))
}

/// Complete Montgomery modular exponentiation (Basic, pre-reduced,
/// proven-odd modulus). **Infallible.** Precondition: `base < modulus`.
pub fn basic_montgomery_mod_exp_pr_odd<T>(base: T, exponent: T, modulus: Odd<T>) -> T
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + Parity
        + core::ops::ShrAssign<usize>
        + const_num_traits::BitsPrecision,
{
    let modulus = modulus.get();
    let w = modulus.bits_precision() as usize;
    let n_prime = compute_n_prime_newton(modulus, w);
    let r_mod_n = compute_r_mod_n(modulus, w);
    let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

    // 1 in Montgomery form = R mod N (since REDC(1 * R²) = 1 * R² * R⁻¹ = R mod N)
    let one_mont = r_mod_n;

    let (lo, hi) = base.wide_mul(&r2_mod_n);
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

    wide_redc(result, T::zero(), modulus, n_prime)
}

/// Complete Montgomery modular exponentiation (Basic, pre-reduced): base^exponent mod modulus
///
/// Precondition: `base < modulus`. No `Rem` bound. Returns None only if
/// modulus is even or zero. Thin wrapper around
/// [`basic_montgomery_mod_exp_pr_odd`].
pub fn basic_montgomery_mod_exp_pr<T>(base: T, exponent: T, modulus: T) -> Option<T>
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + Parity
        + core::ops::ShrAssign<usize>
        + const_num_traits::BitsPrecision,
{
    Odd::new(modulus).map(|m| basic_montgomery_mod_exp_pr_odd(base, exponent, m))
}

/// Complete Montgomery modular exponentiation (Basic, CT, pre-reduced):
/// base^exponent mod modulus, constant-time over `exponent` (and `base`).
///
/// Precondition: `base < modulus`. No `Rem` bound.
///
/// Implements a fixed-iteration constant-time square-and-multiply:
///   - Iterates over **every** bit position of the exponent type (not just
///     significant bits), so the loop count does not leak `bit_length(exp)`.
///   - Performs the squaring and the conditional multiply on **every**
///     iteration, with `subtle::conditional_select` choosing whether to keep
///     the multiplied result based on the current exponent bit — so timing
///     and memory access patterns are independent of the exponent bit
///     pattern.
///   - All Montgomery reductions route through [`wide_redc_ct`].
///
/// Precomputation (`compute_n_prime_newton`, `compute_r_mod_n`,
/// `compute_r2_mod_n`) operates only on the modulus, which is public; using
/// the NCT compute_* helpers there is intentional and does not leak any
/// secret.
///
/// Complete Montgomery modular exponentiation (Basic, CT, pre-reduced,
/// proven-odd modulus). **Infallible.** Precondition: `base < modulus`.
pub fn basic_montgomery_mod_exp_pr_odd_ct<T>(base: T, exponent: T, modulus: Odd<T>) -> T
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + const_num_traits::CtIsZero
        + core::ops::Shr<usize, Output = T>
        + core::ops::BitAnd<Output = T>
        + subtle::ConditionallySelectable
        + subtle::ConstantTimeLess
        + const_num_traits::BitsPrecision,
{
    let modulus = modulus.get();
    let w = modulus.bits_precision() as usize;
    let n_prime = compute_n_prime_newton(modulus, w);
    let r_mod_n = compute_r_mod_n(modulus, w);
    let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

    // 1 in Montgomery form = R mod N
    let one_mont = r_mod_n;

    // Convert base to Montgomery form via REDC(base * R²)
    let base_mont = wide_montgomery_mul_ct(base, r2_mod_n, modulus, n_prime);

    let mut result = one_mont;

    // Constant-time square-and-multiply, MSB to LSB. Iteration count is the
    // exponent width (public shape); equals w on a fixed carrier, tracks the
    // exponent independently of the modulus R width on a runtime-len carrier.
    let exp_bits = exponent.bits_precision() as usize;
    let one = T::one();
    for i in (0..exp_bits).rev() {
        // Always square
        result = wide_montgomery_mul_ct(result, result, modulus, n_prime);

        // Always compute the conditional product
        let multiplied = wide_montgomery_mul_ct(result, base_mont, modulus, n_prime);

        // Extract bit i of exponent (i is public; the shift amount is the
        // loop index, not derived from exp). Bit value is secret.
        let bit_t = (exponent >> i) & one;
        let choice = !bit_t.ct_is_zero();
        result = T::conditional_select(&result, &multiplied, choice);
    }

    // Convert back from Montgomery form: REDC(result, 0)
    wide_redc_ct(result, T::zero(), modulus, n_prime)
}

/// Returns None if modulus is even or zero. Thin wrapper around
/// [`basic_montgomery_mod_exp_pr_odd_ct`].
pub fn basic_montgomery_mod_exp_pr_ct<T>(base: T, exponent: T, modulus: T) -> Option<T>
where
    T: Copy
        + const_num_traits::Zero
        + const_num_traits::One
        + PartialEq
        + PartialOrd
        + WideMul
        + const_num_traits::WrappingMul<Output = T>
        + const_num_traits::WrappingAdd<Output = T>
        + const_num_traits::ops::overflowing::OverflowingAdd<Output = T>
        + const_num_traits::WrappingSub<Output = T>
        + const_num_traits::CtIsZero
        + Parity
        + core::ops::Shr<usize, Output = T>
        + core::ops::BitAnd<Output = T>
        + subtle::ConditionallySelectable
        + subtle::ConstantTimeLess
        + const_num_traits::BitsPrecision,
{
    Odd::new(modulus).map(|m| basic_montgomery_mod_exp_pr_odd_ct(base, exponent, m))
}

#[cfg(test)]
mod tests {
    use super::*;
    use const_num_traits::BitsPrecision;
    use const_num_traits::Ct;
    use fixed_bigint::FixedUInt;

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
    fn odd_surface_matches_option_surface_mul() {
        // The `_odd` (infallible) and `Option`-returning entry points must
        // produce identical results when the modulus is in fact odd. This
        // pins the contract that the wrapper is purely a parity-proof
        // adapter — no behavioural divergence.
        let m: u32 = 97;
        let modulus_odd = Odd::new(m).expect("97 is odd");
        for a in [0u32, 1, 2, 42, 50, 96] {
            for b in [0u32, 1, 13, 49, 95] {
                let via_odd_pr = basic_montgomery_mod_mul_pr_odd(a % m, b % m, modulus_odd);
                let via_opt_pr = basic_montgomery_mod_mul_pr(a % m, b % m, m).unwrap();
                assert_eq!(via_odd_pr, via_opt_pr, "_pr divergence at ({a}, {b})");

                let via_odd = basic_montgomery_mod_mul_odd(a, b, modulus_odd);
                let via_opt = basic_montgomery_mod_mul(a, b, m).unwrap();
                assert_eq!(via_odd, via_opt, "mul divergence at ({a}, {b})");
                assert_eq!(via_odd, (a * b) % m);
            }
        }
    }

    #[test]
    fn odd_surface_matches_option_surface_exp() {
        let m: u32 = 97;
        let modulus_odd = Odd::new(m).expect("97 is odd");
        for base in [0u32, 1, 2, 5, 96] {
            for exp in [0u32, 1, 5, 96] {
                let via_odd_pr = basic_montgomery_mod_exp_pr_odd(base % m, exp, modulus_odd);
                let via_opt_pr = basic_montgomery_mod_exp_pr(base % m, exp, m).unwrap();
                assert_eq!(via_odd_pr, via_opt_pr, "_pr divergence at ({base}, {exp})");

                let via_odd = basic_montgomery_mod_exp_odd(base, exp, modulus_odd);
                let via_opt = basic_montgomery_mod_exp(base, exp, m).unwrap();
                assert_eq!(via_odd, via_opt, "exp divergence at ({base}, {exp})");

                let via_ct = basic_montgomery_mod_exp_pr_odd_ct(base % m, exp, modulus_odd);
                assert_eq!(via_ct, via_odd, "ct/nct divergence at ({base}, {exp})");
            }
        }
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
        let w = modulus.bits_precision() as usize;
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
        let w = modulus.bits_precision() as usize;
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

    /// Sibling of [`test_wide_fixed_bigint_pr`] for the Rem-bound wrappers
    /// (`basic_mod_mul`, `basic_mod_exp`). Default (Nct-personality)
    /// FixedUInt provides Rem and so these wrappers are reachable here.
    /// A Ct-typed FixedUInt would fail to satisfy the Rem bound — which is
    /// the correct outcome, since reduce-then-multiply is variable-time.
    #[test]
    fn test_wide_fixed_bigint() {
        type U128 = FixedUInt<u32, 4>;

        let modulus = !U128::from(0u64) - U128::from(58u64); // 2^128 - 59 (odd)

        let a = U128::from(0xDEAD_BEEF_u64);
        let b = U128::from(0xCAFE_BABE_u64);
        let expected = crate::mul::basic_mod_mul(a, b, modulus);
        let got = basic_montgomery_mod_mul(a, b, modulus).unwrap();
        assert_eq!(got, expected);

        let base = U128::from(42u64);
        let exp = U128::from(1000u64);
        let expected = crate::exp::basic_mod_exp(base, exp, modulus);
        let got = basic_montgomery_mod_exp(base, exp, modulus).unwrap();
        assert_eq!(got, expected);
    }

    #[test]
    fn test_wide_fixed_bigint_pr() {
        type U128 = FixedUInt<u32, 4>;

        // Pick a 128-bit-ish odd modulus close to type max
        let modulus = !U128::from(0u64) - U128::from(58u64); // 2^128 - 59 (odd)

        // Round-trip test via wide helpers directly
        let w = modulus.bits_precision() as usize;
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

        // All test inputs are tiny (< 2^32) and modulus is ~2^128, so the
        // pre-reduced precondition is naturally satisfied.
        let a = U128::from(0xDEAD_BEEF_u64);
        let (lo, hi) = a.wide_mul(&r2_mod_n);
        let a_m = wide_redc(lo, hi, modulus, n_prime);
        let back = wide_redc(a_m, U128::from(0u64), modulus, n_prime);
        assert_eq!(back, a);

        // Multiplication: compare _pr Montgomery against _pr double-and-add.
        let b = U128::from(0xCAFE_BABE_u64);
        let expected = crate::mul::basic_mod_mul_pr(a, b, modulus);
        let got = basic_montgomery_mod_mul_pr(a, b, modulus).unwrap();
        assert_eq!(got, expected);

        // Exponentiation
        let base = U128::from(42u64);
        let exp = U128::from(1000u64);
        let expected = crate::exp::basic_mod_exp_pr(base, exp, modulus);
        let got = basic_montgomery_mod_exp_pr(base, exp, modulus).unwrap();
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

    /// wide_redc_ct must produce the same output as wide_redc for all inputs.
    /// Exhaustive over a small u8 modulus to cover both the "needs subtraction"
    /// and "doesn't need subtraction" branches plus the extra_bit case.
    #[test]
    fn test_wide_redc_ct_matches_nct_u8() {
        let modulus = 13u8;
        for t_lo in 0u8..=255 {
            for t_hi in 0u8..modulus {
                // Pick an arbitrary odd n_prime — exact value doesn't matter for
                // comparing the two implementations; they share inputs.
                let n_prime = 11u8;
                let nct = wide_redc(t_lo, t_hi, modulus, n_prime);
                let ct = wide_redc_ct(t_lo, t_hi, modulus, n_prime);
                assert_eq!(nct, ct, "wide_redc_ct mismatch at t_lo={t_lo} t_hi={t_hi}");
            }
        }
    }

    /// All four `strict_wide_*` functions must produce identical outputs to
    /// their by-value siblings on Copy types. Same input grid as the
    /// `wide_redc_ct` equivalence test; local `basic_*` aliases keep the
    /// comparison readable.
    #[test]
    fn test_strict_wide_matches_basic_u8() {
        // Local aliases for "the by-value (basic-flavor) version" purely so
        // the asserts read symmetrically.
        use super::{
            wide_montgomery_mul as basic_wide_montgomery_mul,
            wide_montgomery_mul_ct as basic_wide_montgomery_mul_ct, wide_redc as basic_wide_redc,
            wide_redc_ct as basic_wide_redc_ct,
        };

        let modulus = 13u8;
        let n_prime = 11u8;
        for t_lo in 0u8..=255 {
            for t_hi in 0u8..modulus {
                assert_eq!(
                    basic_wide_redc(t_lo, t_hi, modulus, n_prime),
                    strict_wide_redc(&t_lo, &t_hi, &modulus, &n_prime),
                    "strict_wide_redc mismatch at t_lo={t_lo} t_hi={t_hi}"
                );
                assert_eq!(
                    basic_wide_redc_ct(t_lo, t_hi, modulus, n_prime),
                    strict_wide_redc_ct(&t_lo, &t_hi, &modulus, &n_prime),
                    "strict_wide_redc_ct mismatch at t_lo={t_lo} t_hi={t_hi}"
                );
            }
        }

        // mul wrappers exercise a different code path (wide_mul + redc).
        for a in 0u8..modulus {
            for b in 0u8..modulus {
                assert_eq!(
                    basic_wide_montgomery_mul(a, b, modulus, n_prime),
                    strict_wide_montgomery_mul(&a, &b, &modulus, &n_prime),
                    "strict_wide_montgomery_mul mismatch at a={a} b={b}"
                );
                assert_eq!(
                    basic_wide_montgomery_mul_ct(a, b, modulus, n_prime),
                    strict_wide_montgomery_mul_ct(&a, &b, &modulus, &n_prime),
                    "strict_wide_montgomery_mul_ct mismatch at a={a} b={b}"
                );
            }
        }
    }

    /// `mul_acc` with a zero accumulator + `redc` must produce the same
    /// result as the existing fused `wide_montgomery_mul`. Identity check
    /// on the new primitive.
    #[test]
    fn test_wide_mul_acc_identity_u8() {
        let modulus = 13u8;
        let n_prime = compute_n_prime_newton(modulus, 8);
        for a in 0u8..modulus {
            for b in 0u8..modulus {
                let direct = wide_montgomery_mul(a, b, modulus, n_prime);
                let (lo, hi) = wide_montgomery_mul_acc(0u8, 0u8, a, b);
                let via_acc = wide_redc(lo, hi, modulus, n_prime);
                assert_eq!(direct, via_acc, "identity mismatch at a={a} b={b}");
            }
        }
    }

    /// Σ a_i * b_i computed via `mul_acc` per term + one `redc` at the
    /// end must equal the direct residue-domain sum of products. Proves
    /// the accumulation semantics — the core promise of the primitive.
    #[test]
    fn test_wide_mul_acc_dot_product_u8() {
        let modulus = 13u8;
        let n_prime = compute_n_prime_newton(modulus, 8);
        let r_mod_n = compute_r_mod_n(modulus, 8);
        let r2 = compute_r2_mod_n(r_mod_n, modulus, 8);
        let to_mont = |x: u8| {
            let (lo, hi) = x.wide_mul(&r2);
            wide_redc(lo, hi, modulus, n_prime)
        };
        let pairs: &[(u8, u8)] = &[(2, 3), (5, 7), (11, 4), (1, 12)];

        let (mut acc_lo, mut acc_hi) = (0u8, 0u8);
        for &(a, b) in pairs {
            let (l, h) = wide_montgomery_mul_acc(acc_lo, acc_hi, to_mont(a), to_mont(b));
            acc_lo = l;
            acc_hi = h;
        }
        let result_mont = wide_redc(acc_lo, acc_hi, modulus, n_prime);
        let result = wide_redc(result_mont, 0u8, modulus, n_prime);

        let expected = pairs
            .iter()
            .fold(0u8, |acc, &(a, b)| (acc + (a * b) % modulus) % modulus);
        assert_eq!(result, expected);
    }

    /// CT and NCT `mul_acc` variants must produce bit-identical output
    /// over a representative input grid.
    #[test]
    fn test_wide_mul_acc_ct_matches_nct_u8() {
        let modulus = 13u8;
        for a in 0u8..modulus {
            for b in 0u8..modulus {
                for acc_lo in [0u8, 1, 50, 100, 200, 255] {
                    for acc_hi in [0u8, 1, 5, modulus - 1] {
                        let nct = wide_montgomery_mul_acc(acc_lo, acc_hi, a, b);
                        let ct = wide_montgomery_mul_acc_ct(acc_lo, acc_hi, a, b);
                        assert_eq!(
                            nct, ct,
                            "ct/nct mismatch at lo={acc_lo} hi={acc_hi} a={a} b={b}"
                        );
                    }
                }
            }
        }
    }

    /// Reference-based `strict_*_mul_acc` siblings must match the
    /// by-value forms on `Copy` types.
    #[test]
    fn test_strict_wide_mul_acc_matches_basic_u8() {
        let modulus = 13u8;
        for a in 0u8..modulus {
            for b in 0u8..modulus {
                for acc_lo in [0u8, 1, 50, 100, 255] {
                    for acc_hi in [0u8, 1, 5, modulus - 1] {
                        assert_eq!(
                            wide_montgomery_mul_acc(acc_lo, acc_hi, a, b),
                            strict_wide_montgomery_mul_acc(&acc_lo, &acc_hi, &a, &b),
                            "strict mul_acc mismatch at lo={acc_lo} hi={acc_hi} a={a} b={b}"
                        );
                        assert_eq!(
                            wide_montgomery_mul_acc_ct(acc_lo, acc_hi, a, b),
                            strict_wide_montgomery_mul_acc_ct(&acc_lo, &acc_hi, &a, &b),
                            "strict mul_acc_ct mismatch at lo={acc_lo} hi={acc_hi} a={a} b={b}"
                        );
                    }
                }
            }
        }
    }

    /// ML-KEM-shaped dot product at `q=3329`, `T=u32`. Exercises the
    /// actual use-case: 4 Mont products summed via `mul_acc`, single
    /// `redc` at the end.
    #[test]
    fn test_wide_mul_acc_dot_product_u32_mlkem_q() {
        let modulus = 3329u32;
        let w = modulus.bits_precision() as usize;
        let n_prime = compute_n_prime_newton(modulus, w);
        let r2 = compute_r2_mod_n(compute_r_mod_n(modulus, w), modulus, w);
        let to_mont = |x: u32| {
            let (lo, hi) = x.wide_mul(&r2);
            wide_redc(lo, hi, modulus, n_prime)
        };
        let pairs: &[(u32, u32)] = &[(123, 456), (789, 1011), (2222, 3000), (1, 3328)];

        let (mut acc_lo, mut acc_hi) = (0u32, 0u32);
        for &(a, b) in pairs {
            let (l, h) = wide_montgomery_mul_acc(acc_lo, acc_hi, to_mont(a), to_mont(b));
            acc_lo = l;
            acc_hi = h;
        }
        let result_mont = wide_redc(acc_lo, acc_hi, modulus, n_prime);
        let result = wide_redc(result_mont, 0u32, modulus, n_prime);

        let expected = pairs.iter().fold(0u64, |acc, &(a, b)| {
            acc + (a as u64 * b as u64) % modulus as u64
        }) % modulus as u64;
        assert_eq!(result as u64, expected);
    }

    /// wide_redc_ct matches at FixedUInt sizes.
    ///
    /// Under fixed-bigint's personality typestate, `wide_redc_ct`'s
    /// `ConditionallySelectable`/`ConstantTimeLess` bounds resolve only for
    /// Ct-typed FixedUInts; the Nct-typed precompute values are bridged via
    /// the free `.into()` upgrade, and the CT result is brought back via
    /// `forget_ct()` for cross-personality equality.
    #[test]
    fn test_wide_redc_ct_matches_nct_fixed() {
        type U128 = FixedUInt<u32, 4>;
        type U128Ct = FixedUInt<u32, 4, Ct>;

        let modulus = !U128::from(0u64) - U128::from(58u64);
        let w = modulus.bits_precision() as usize;
        let n_prime = compute_n_prime_newton(modulus, w);
        let modulus_ct: U128Ct = modulus.into();
        let n_prime_ct: U128Ct = n_prime.into();

        // A mix: small values, values near modulus, and full-width-ish values
        let test_vals = [
            U128::from(0u64),
            U128::from(1u64),
            U128::from(0xDEAD_BEEF_u64),
            modulus - U128::from(1u64),
            !U128::from(0u64), // 2^128 - 1 (well above modulus)
        ];
        for &t_lo in &test_vals {
            for &t_hi in &test_vals {
                let nct = wide_redc(t_lo, t_hi, modulus, n_prime);
                let t_lo_ct: U128Ct = t_lo.into();
                let t_hi_ct: U128Ct = t_hi.into();
                let ct = wide_redc_ct(t_lo_ct, t_hi_ct, modulus_ct, n_prime_ct);
                assert_eq!(
                    nct,
                    ct.forget_ct(),
                    "wide_redc_ct mismatch at t_lo={t_lo:?} t_hi={t_hi:?}"
                );
            }
        }
    }

    /// basic_montgomery_mod_exp_pr_ct must produce the same output as
    /// the NCT version when the caller pre-reduces base.
    #[test]
    fn test_basic_montgomery_mod_exp_pr_ct_matches_nct_u32_with_external_reduce() {
        let modulus = 13u32;
        for base in 0u32..modulus {
            for exp in 0u32..32 {
                let nct = basic_montgomery_mod_exp(base, exp, modulus).unwrap();
                // External pre-reduction: caller's responsibility on
                // the CT path. Inputs are already in [0, modulus) for
                // this test, so the `% modulus` is a no-op observable.
                let ct = basic_montgomery_mod_exp_pr_ct(base % modulus, exp, modulus).unwrap();
                assert_eq!(nct, ct, "mod_exp_pr_ct mismatch at base={base} exp={exp}");
            }
        }
    }

    /// basic_montgomery_mod_exp_pr_ct must produce the same output as the NCT
    /// _pr version for pre-reduced inputs.
    #[test]
    fn test_basic_montgomery_mod_exp_pr_ct_matches_nct_u32() {
        let modulus = 13u32;
        for base in 0u32..modulus {
            for exp in 0u32..32 {
                let nct = basic_montgomery_mod_exp_pr(base, exp, modulus).unwrap();
                let ct = basic_montgomery_mod_exp_pr_ct(base, exp, modulus).unwrap();
                assert_eq!(nct, ct, "mod_exp_pr_ct mismatch at base={base} exp={exp}");
            }
        }
    }

    /// CT exp matches NCT exp at FixedUInt sizes.
    ///
    /// See [`test_wide_redc_ct_matches_nct_fixed`] for the personality-bridge
    /// rationale — Ct-typed inputs are required for `_ct` functions.
    #[test]
    fn test_basic_montgomery_mod_exp_pr_ct_matches_nct_fixed() {
        type U128 = FixedUInt<u32, 4>;
        type U128Ct = FixedUInt<u32, 4, Ct>;

        let modulus = !U128::from(0u64) - U128::from(58u64); // 2^128 - 59 (odd prime)
        let modulus_ct: U128Ct = modulus.into();
        let bases = [U128::from(2u64), U128::from(0xDEAD_BEEF_u64)];
        let exps = [
            U128::from(0u64),
            U128::from(1u64),
            U128::from(7u64),
            U128::from(0xCAFE_BABE_u64),
        ];
        for &base in &bases {
            for &exp in &exps {
                let nct = basic_montgomery_mod_exp_pr(base, exp, modulus).unwrap();
                let base_ct: U128Ct = base.into();
                let exp_ct: U128Ct = exp.into();
                let ct = basic_montgomery_mod_exp_pr_ct(base_ct, exp_ct, modulus_ct).unwrap();
                assert_eq!(
                    nct,
                    ct.forget_ct(),
                    "mod_exp_pr_ct mismatch at base={base:?} exp={exp:?}"
                );
            }
        }
    }

    /// Equivalence of the NCT and CT high-half accumulators. Exhaustive over
    /// `result` (u8), both carry flags. Confirms `accumulate_high_half_carry_ct`
    /// produces the same `(result, extra_bit)` pair as `accumulate_high_half_carry`
    /// for every input — the swap inside `wide_redc_ct` / `strict_wide_redc_ct`
    /// is purely a side-channel hardening, not a semantic change.
    #[test]
    fn test_accumulate_high_half_carry_ct_matches_nct_u8() {
        for result in 0u8..=255 {
            for &carry1 in &[false, true] {
                for &carry2 in &[false, true] {
                    let nct = accumulate_high_half_carry(result, carry1, carry2);
                    let ct = accumulate_high_half_carry_ct(
                        result,
                        Choice::from(carry1 as u8),
                        Choice::from(carry2 as u8),
                    );
                    let ct = (ct.0, bool::from(ct.1));
                    assert_eq!(
                        nct, ct,
                        "mismatch at result={result} carry1={carry1} carry2={carry2}"
                    );
                }
            }
        }
    }
}
