// Clippy's `clone_on_copy` / `op_ref` lints misfire on generic code that
// uses `.clone()` or `&x + &T::zero()` as portable "produce owned T"
// idioms across trait-bound flavors where T may or may not be Copy.
#![allow(clippy::clone_on_copy, clippy::op_ref)]

//! Bernstein-Yang constant-time modular inverse via "safegcd" divsteps.
//!
//! Computes `a⁻¹ mod n` in constant time over the value being inverted,
//! for **any modulus** — composite or prime. Used by RSA private-key
//! blinding (where `n = p·q` is composite, so Fermat's little theorem
//! doesn't apply) and anywhere else CT inversion of arbitrary residues
//! is needed.
//!
//! Reference: Bernstein & Yang, *"Fast constant-time gcd computation
//! and modular inversion"*, IACR ePrint 2019/266
//! (<https://gcd.cr.yp.to/>).
//!
//! ## Algorithm shape
//!
//! Operates on:
//! - `f, g`: signed two's-complement values stored in unsigned `T`.
//!   `f` starts as `modulus`, `g` starts as `value`. They evolve via
//!   divsteps until `g == 0` and `f == ±1` (iff `gcd == 1`).
//! - `d, e`: modular coefficients kept in `[0, modulus)` throughout,
//!   such that `d * value ≡ f mod modulus` and `e * value ≡ g mod
//!   modulus`. At convergence, `d` is `±value⁻¹ mod modulus`.
//! - `delta`: small signed integer that tracks the divsteps state
//!   machine. Stored as `i64`.
//!
//! ## Why this is "signed" without naming a signed bigint type
//!
//! The algorithm conceptually operates on signed values `f`, `g` that
//! range over `(-modulus, modulus)`. **The implementation never types
//! a signed bigint.** Every "signed" op is two's-complement-on-unsigned:
//!
//! - signed `f + g`  ≡ `T::wrapping_add(f, g)`         on unsigned `T`
//! - signed `-f`     ≡ `T::wrapping_neg(f)`            on unsigned `T`
//! - signed `f < 0`  ≡ MSB of unsigned `f` is set
//! - arithmetic shr  ≡ `(x >> 1) | sign_extend_mask`   on unsigned `T`
//!
//! ## Bound precondition
//!
//! `modulus` must fit in `T` with at least one bit of headroom — i.e.
//! `2 * modulus` does not overflow `T`. The algorithm maintains
//! `|d|, |e| < modulus` strictly, so intermediate sums `d + e` are
//! bounded by `2 * modulus`. Practical types (e.g.
//! `FixedUInt<u32, 64>` for a 2048-bit RSA modulus) leave plenty of
//! headroom.
//!
//! ## Implementation note
//!
//! This is a **per-step (unbatched) implementation**: each divstep
//! does multi-precision arithmetic on the full `T`. Asymptotically
//! `O(n²)` for an `n`-bit modulus. The batched ("jumpdivstep")
//! variant brings this to `O(n²/W)` where `W` is the limb width — a
//! follow-up optimization. The unbatched version is **correct** and
//! suitable for cases where the latency budget allows it (RSA blinding
//! at 2048 bits is dominated by the main exponentiation anyway).

use const_num_traits::{CtIsZero, CtParity, One, WrappingAdd, WrappingSub, Zero};
use modmath_cios::CiosRowOps;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, ConstantTimeLess, CtOption};

/// Total number of divsteps for an `n`-bit operand. Per Theorem 11.2
/// of the Bernstein-Yang paper, for the *integer* case with δ₀=1 (our
/// initialization at [`safegcd_inv_ct`]) the proven sufficient count
/// is `⌊(49·n + 80)/17⌋` for `n < 46` and `⌊(49·n + 57)/17⌋` for
/// `n ≥ 46`. Under-iterating causes worst-case coprime pairs to exit
/// divsteps with `g ≠ 0`, silently masking valid inverses to `None`.
/// (The `(45·n + 64)/19` bound is the *polynomial* case from Section
/// 7 — not applicable here.)
pub(crate) const fn divsteps_total(modulus_bits: usize) -> usize {
    if modulus_bits < 46 {
        (49 * modulus_bits + 80) / 17
    } else {
        (49 * modulus_bits + 57) / 17
    }
}

/// Returns `Choice::from(1)` iff the low bit of `value` is set.
///
/// Delegates to cnt's [`CtParity`] on the low limb. The CT contract
/// for "is this primitive odd" lives upstream; we compose by
/// extracting `word(0)` (which is a primitive once the limb level is
/// reached).
#[inline]
fn ct_low_bit<T>(value: &T) -> Choice
where
    T: CiosRowOps,
    T::Word: CtParity,
{
    value.word(0).ct_is_odd()
}

/// Returns `Choice::from(1)` iff the most-significant bit of `value`'s
/// bit pattern is set (i.e. value is "negative" in two's-complement
/// interpretation).
///
/// Exposed at crate visibility so `Field::inv_safegcd_ct` can fold the
/// "modulus has one bit of headroom" precondition into its CtOption
/// instead of relying on documentation alone.
#[inline]
pub(crate) fn ct_msb_set<T>(value: &T) -> Choice
where
    T: CiosRowOps,
    T::Word: core::ops::BitAnd<Output = T::Word>
        + core::ops::Shl<usize, Output = T::Word>
        + One
        + CtIsZero,
{
    let n = value.word_count();
    let word_bits = core::mem::size_of::<T::Word>() * 8;
    let msb_mask = T::Word::one() << (word_bits - 1);
    let masked = value.word(n - 1) & msb_mask;
    // Delegates to cnt's CtIsZero rather than ct_eq(&T::Word::zero()).
    // Identical semantics; the dedicated trait is cnt-tested upstream.
    !masked.ct_is_zero()
}

/// Constant-time `delta > 0` for the i64 state variable.
#[inline]
fn ct_i64_positive(delta: i64) -> Choice {
    let delta_u = delta as u64;
    // nonzero = ((x | -x) >> 63) — top bit set iff x != 0
    let nonzero_top = (delta_u | delta_u.wrapping_neg()) >> 63;
    // sign_bit_clear = 1 iff x's high bit is 0 (i.e. x >= 0 as i64)
    let sign_bit_clear = (!delta_u) >> 63;
    Choice::from(((nonzero_top & sign_bit_clear) & 1) as u8)
}

/// Arithmetic right shift by 1 (sign-extending). For T interpreted as
/// two's-complement, returns `value / 2` with floor rounding for
/// negative values.
#[inline]
fn arithmetic_shr_one<T>(value: &T) -> T
where
    T: CiosRowOps
        + Clone
        + ConditionallySelectable
        + One
        + Zero
        + core::ops::Shr<usize, Output = T>
        + core::ops::Shl<usize, Output = T>
        + core::ops::BitOr<Output = T>,
    T::Word: core::ops::BitAnd<Output = T::Word>
        + core::ops::Shl<usize, Output = T::Word>
        + One
        + CtIsZero,
{
    let logical = value.clone() >> 1;
    let msb_set = ct_msb_set(value);
    let n_bits = value.word_count() * core::mem::size_of::<T::Word>() * 8;
    let top_bit_mask = T::one() << (n_bits - 1);
    let with_sign_ext = logical.clone() | top_bit_mask;
    T::conditional_select(&logical, &with_sign_ext, msb_set)
}

/// Reduce `sum` to `[0, m)` when `sum` is already in `[0, 2·m)`.
/// One conditional subtract suffices.
#[inline]
fn reduce_lt_2m_ct<T>(sum: T, m: &T) -> T
where
    T: Clone + ConditionallySelectable + WrappingSub<Output = T> + ConstantTimeLess,
{
    let sum_minus_m = sum.clone().wrapping_sub(m.clone());
    let need_sub = !sum.ct_lt(m);
    T::conditional_select(&sum, &sum_minus_m, need_sub)
}

/// Modular add: `(a + b) mod m`, assuming `a, b ∈ [0, m)` so the sum
/// is in `[0, 2·m)`. Precondition: `2·m` fits in `T` (one bit of
/// headroom over `m`).
#[inline]
fn add_mod_ct<T>(a: &T, b: &T, m: &T) -> T
where
    T: Clone
        + ConditionallySelectable
        + WrappingAdd<Output = T>
        + WrappingSub<Output = T>
        + ConstantTimeLess,
{
    let sum = a.clone().wrapping_add(b.clone());
    reduce_lt_2m_ct(sum, m)
}

/// Modular negate: `(m - x) mod m`, preserving the `[0, m)` invariant.
/// Returns 0 when `x == 0` (since `-0 mod m = 0`).
///
/// Delegates to cnt's [`CtIsZero`] for the masked zero check rather
/// than hand-rolling `ct_eq(&T::zero())`. Identical semantically; the
/// dedicated trait is cnt-tested.
#[inline]
fn neg_mod_ct<T>(x: &T, m: &T) -> T
where
    T: Clone + Zero + ConditionallySelectable + WrappingSub<Output = T> + CtIsZero,
{
    let zero = T::zero();
    let x_is_zero = x.ct_is_zero();
    let neg = m.clone().wrapping_sub(x.clone());
    T::conditional_select(&neg, &zero, x_is_zero)
}

/// Modular halving: returns `y ∈ [0, m)` such that `2·y ≡ x mod m`.
/// Precondition: `x ∈ [0, m)`, `m` is odd.
///
/// For odd `m`, the inverse of 2 mod `m` is `(m + 1) / 2`; computing
/// `x / 2 mod m` via the standard branch-free trick:
/// - if `x` even: result is `x >> 1`.
/// - if `x` odd: result is `(x + m) >> 1` (since `x + m` is even when
///   `m` is odd, and the sum is in `[m, 2m)` which fits with one bit
///   of headroom).
#[inline]
fn half_mod_ct<T>(x: &T, m: &T) -> T
where
    T: CiosRowOps
        + Clone
        + ConditionallySelectable
        + WrappingAdd<Output = T>
        + core::ops::Shr<usize, Output = T>,
    T::Word: CtParity,
{
    let x_odd = ct_low_bit(x);
    let x_plus_m = x.clone().wrapping_add(m.clone());
    let candidate_odd = x_plus_m >> 1;
    let candidate_even = x.clone() >> 1;
    T::conditional_select(&candidate_even, &candidate_odd, x_odd)
}

/// Compute `value⁻¹ mod modulus` in constant time over `value`. Works
/// for any modulus (composite or prime), provided the modulus is odd
/// and `2 * modulus` fits in `T`.
///
/// Returns `CtOption::Some(inv)` when `gcd(value, modulus) == 1`, or
/// `CtOption::None` masked when no inverse exists. The failure path
/// timing is independent of input magnitudes.
///
/// The `modulus is odd` and `modulus > 1` preconditions are folded
/// into the returned mask: passing an even modulus or `1` yields
/// `CtOption::None`, matching the behavior for "no inverse exists" —
/// the divsteps loop still runs the full step count in constant time.
///
/// # Preconditions
///
/// - `value < modulus` (caller should reduce first)
/// - `2 * modulus` fits in `T` without overflow (i.e. `T` is at least
///   one bit wider than `modulus`)
pub fn safegcd_inv_ct<T>(value: &T, modulus: &T) -> CtOption<T>
where
    T: CiosRowOps
        + Clone
        + ConditionallySelectable
        + ConstantTimeEq
        + ConstantTimeLess
        + CtIsZero
        + Zero
        + One
        + WrappingAdd<Output = T>
        + WrappingSub<Output = T>
        + core::ops::Shr<usize, Output = T>
        + core::ops::Shl<usize, Output = T>
        + core::ops::BitOr<Output = T>,
    T::Word: Copy
        + ConditionallySelectable
        + ConstantTimeEq
        + CtIsZero
        + CtParity
        + One
        + Zero
        + core::ops::BitAnd<Output = T::Word>
        + core::ops::Shl<usize, Output = T::Word>,
{
    let n_bits = value.word_count() * core::mem::size_of::<T::Word>() * 8;
    let total_steps = divsteps_total(n_bits);

    let mut f = modulus.clone();
    let mut g = value.clone();
    let mut d = T::zero();
    let mut e = T::one();
    let mut delta: i64 = 1;

    for _ in 0..total_steps {
        let delta_pos = ct_i64_positive(delta);
        let g_odd = ct_low_bit(&g);
        let swap = delta_pos & g_odd;

        // swap_choice: (f, g) ← (g, -f); (d, e) ← (e, -d mod m); delta ← -delta
        let new_f_if_swap = g.clone();
        let new_g_if_swap = f.clone().wrapping_neg_two_complement();
        let new_d_if_swap = e.clone();
        let new_e_if_swap = neg_mod_ct(&d, modulus);

        f = T::conditional_select(&f, &new_f_if_swap, swap);
        g = T::conditional_select(&g, &new_g_if_swap, swap);
        d = T::conditional_select(&d, &new_d_if_swap, swap);
        e = T::conditional_select(&e, &new_e_if_swap, swap);

        let neg_delta = (delta as u64).wrapping_neg() as i64;
        delta = i64::conditional_select(&delta, &neg_delta, swap);
        delta = delta.wrapping_add(1);

        // g_odd was determined before swap; after swap+negate the parity
        // of g equals that of the original g (both old g and -old_f are
        // odd when the swap fired, since old f is odd by loop invariant).
        // In the non-swap case g_odd tracks current g's parity directly.
        let g_odd_now = g_odd;
        let to_add_to_g = T::conditional_select(&T::zero(), &f, g_odd_now);
        g = g.wrapping_add(to_add_to_g);

        let add_to_e = add_mod_ct(&e, &d, modulus);
        e = T::conditional_select(&e, &add_to_e, g_odd_now);

        g = arithmetic_shr_one(&g);
        e = half_mod_ct(&e, modulus);
    }

    // After total_steps divsteps, g == 0 and f == ±gcd. For the
    // invertible case (gcd == 1), f is either +1 (bit pattern T::one())
    // or -1 (bit pattern T::one().wrapping_neg()).
    let one = T::one();
    let neg_one_pattern = one.clone().wrapping_neg_two_complement();

    let f_is_one = f.ct_eq(&one);
    let f_is_neg_one = f.ct_eq(&neg_one_pattern);
    let has_inverse = f_is_one | f_is_neg_one;

    // Result: d if f == 1; (m - d) = -d mod m if f == -1.
    let neg_d = neg_mod_ct(&d, modulus);
    let result = T::conditional_select(&d, &neg_d, f_is_neg_one);

    // Fold in the modulus preconditions. `half_mod_ct`'s (x + m) >> 1
    // trick is invalid unless m is odd; modulus = 1 collapses residues
    // to a single class. Divsteps still ran (constant time), the
    // resulting `d`/`e` are garbage in those cases, and the mask
    // discards them.
    let modulus_is_odd = ct_low_bit(modulus);
    let modulus_is_one = modulus.ct_eq(&one);
    let modulus_ok = modulus_is_odd & !modulus_is_one;

    CtOption::new(result, has_inverse & modulus_ok)
}

// Inline WrappingNeg-style two's-complement negate. const-num-traits'
// WrappingNeg has an `Output` associated type that doesn't necessarily
// equal T for generic T; this trait gives a clean `T -> T` shape as
// `0 - x` via wrapping_sub.
trait WrappingNegT: Sized {
    fn wrapping_neg_two_complement(self) -> Self;
}

impl<T> WrappingNegT for T
where
    T: Zero + WrappingSub<Output = T>,
{
    #[inline]
    fn wrapping_neg_two_complement(self) -> Self {
        T::zero().wrapping_sub(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::inv::basic_mod_inv;

    #[test]
    fn divsteps_total_matches_paper_bound() {
        assert!(divsteps_total(256) >= (45usize * 256 + 64).div_ceil(19));
        assert!(divsteps_total(2048) >= (45usize * 2048 + 64).div_ceil(19));
        assert!(divsteps_total(4096) >= (45usize * 4096 + 64).div_ceil(19));
    }

    /// Exhaustive cross-check against `basic_mod_inv` over a small
    /// odd-prime modulus. T = u32 gives plenty of headroom for the
    /// `2 * modulus` precondition.
    #[test]
    fn matches_basic_mod_inv_mod_small_primes_u32() {
        for &m in &[3u32, 5, 7, 11, 13, 17, 19, 23, 29, 31, 97, 251] {
            for v in 1..m {
                let got = safegcd_inv_ct::<u32>(&v, &m);
                let want = basic_mod_inv(v, m);
                match (got.into_option(), want) {
                    (Some(g), Some(w)) => {
                        assert_eq!(
                            g, w,
                            "value={v} modulus={m}: safegcd gave {g}, basic_mod_inv gave {w}"
                        );
                        assert_eq!((g as u64 * v as u64) % m as u64, 1, "inv check failed");
                    }
                    (Some(_), None) | (None, Some(_)) => {
                        panic!("disagreement on value={v} modulus={m}");
                    }
                    (None, None) => {}
                }
            }
        }
    }

    /// Composite modulus — the load-bearing case for RSA blinding,
    /// where Fermat's little theorem doesn't apply but safegcd still
    /// computes the inverse correctly.
    #[test]
    fn handles_composite_modulus_u32() {
        // n = 3 * 5 = 15. Values coprime to 15: 1, 2, 4, 7, 8, 11, 13, 14.
        let m: u32 = 15;
        let coprime_values = [1u32, 2, 4, 7, 8, 11, 13, 14];
        for &v in &coprime_values {
            let got = safegcd_inv_ct::<u32>(&v, &m).into_option();
            let want = basic_mod_inv(v, m);
            assert_eq!(got, want, "value={v} modulus=15: composite case mismatch");
            if let Some(inv) = got {
                assert_eq!((inv as u64 * v as u64) % 15, 1);
            }
        }
        for &v in &[3u32, 5, 6, 9, 10, 12] {
            assert!(
                safegcd_inv_ct::<u32>(&v, &m).into_option().is_none(),
                "value={v} modulus=15: expected None (not coprime)"
            );
        }
    }

    /// FixedUInt end-to-end test. Confirms the algorithm works for
    /// multi-limb T as it does for primitive T. This is the
    /// load-bearing case — RSA uses FixedUInt, not primitives.
    #[test]
    fn fixed_bigint_smoke_test() {
        use const_num_traits::Ct;
        use fixed_bigint::FixedUInt;
        type U64Ct = FixedUInt<u32, 2, Ct>;
        type U64Nct = FixedUInt<u32, 2>;

        let m_raw: u64 = 0x7FFF_FFFF_FFFF_FFE7; // 2^63 - 25, prime
        let m = U64Ct::from(m_raw);
        let m_nct = U64Nct::from(m_raw);
        let test_vals = [1u64, 2, 7, 42, 0xDEAD_BEEF];
        for &v_raw in &test_vals {
            let v = U64Ct::from(v_raw);
            let v_nct = U64Nct::from(v_raw);
            let got = safegcd_inv_ct(&v, &m).into_option();
            assert!(got.is_some(), "expected inverse to exist for v={v_raw}");
            let inv = got.unwrap();
            // Cross-check by converting Ct → Nct and running basic_mod_inv
            let baseline = basic_mod_inv(v_nct, m_nct).expect("baseline inv");
            let inv_nct: U64Nct = inv.forget_ct();
            assert_eq!(
                inv_nct, baseline,
                "FixedUInt v={v_raw}: mismatch with basic_mod_inv"
            );
        }
    }

    /// u64 composite modulus. Coprime values must invert; 0xDEAD_BEEF
    /// (shares factor 11 with n = 0x100000707000031) must return
    /// CtOption::None.
    #[test]
    fn u64_composite_coprime() {
        let n: u64 = 0x1_0000_0007 * 0x100_0007;
        for v in [1u64, 2, 3, 0xCAFE_BABE, 0xFEED_FACE] {
            let got = safegcd_inv_ct::<u64>(&v, &n).into_option();
            assert!(got.is_some(), "u64 case: v={v} mod {n:#x} expected Some");
            let inv = got.unwrap();
            let prod = (inv as u128 * v as u128) % n as u128;
            assert_eq!(prod, 1, "u64 v={v}: inv*v mod n != 1");
        }
        // Non-coprime: 0xDEAD_BEEF shares factor 11 with n
        assert!(
            safegcd_inv_ct::<u64>(&0xDEAD_BEEF, &n)
                .into_option()
                .is_none()
        );
    }

    /// Larger primitive: u64 modulus, sanity check on a handful of pairs.
    /// Uses a modulus < 2^63 to satisfy the `2 * modulus fits in T`
    /// precondition.
    #[test]
    fn u64_smoke_test() {
        let m: u64 = 0x7FFF_FFFF_FFFF_FFE7; // 2^63 - 25, a prime
        let test_vals = [1u64, 2, 7, 0xDEAD_BEEF, 0xCAFE_BABE];
        for &v in &test_vals {
            let got = safegcd_inv_ct::<u64>(&v, &m).into_option();
            let want = basic_mod_inv(v, m);
            assert_eq!(got, want, "value={v} modulus={m}: u64 case mismatch");
            if let Some(inv) = got {
                let prod = (inv as u128 * v as u128) % m as u128;
                assert_eq!(prod, 1, "u64 inv * value mod m != 1");
            }
        }
    }

    /// Even modulus and modulus = 1 violate the algorithm's precondition
    /// (`half_mod_ct`'s trick assumes odd m; modulus = 1 has a single
    /// residue class). The mask folds those into CtOption::None.
    #[test]
    fn modulus_precondition_masks() {
        for m in [2u32, 4, 6, 100, 0xFFFF_FFFE] {
            for v in [1u32, 3, 7, 0xCAFE_BABE] {
                assert!(
                    safegcd_inv_ct::<u32>(&v, &m).into_option().is_none(),
                    "even modulus {m:#x} with value {v:#x}: expected None"
                );
            }
        }
        for v in [0u32, 1, 7] {
            assert!(
                safegcd_inv_ct::<u32>(&v, &1u32).into_option().is_none(),
                "modulus = 1 with value {v}: expected None"
            );
        }
    }
}
