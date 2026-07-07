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
//! - `f, g`: signed two's-complement values stored in [`SignedExt`]
//!   (unsigned `T` plus one extension word). `f` starts as `modulus`,
//!   `g` starts as `value`. They evolve via divsteps until `g == 0`
//!   and `f == ±1` (iff `gcd == 1`).
//! - `d, e`: modular coefficients kept in `[0, modulus)` throughout,
//!   such that `d * value ≡ f mod modulus` and `e * value ≡ g mod
//!   modulus`. At convergence, `d` is `±value⁻¹ mod modulus`.
//! - `delta`: small signed integer that tracks the divsteps state
//!   machine. Stored as `i64`.
//!
//! ## Why this is "signed" without naming a signed bigint type
//!
//! The algorithm conceptually operates on signed values `f`, `g` that
//! range over `(-modulus, modulus)`, transiently `(-2·modulus,
//! 2·modulus)`. **The implementation never types a signed bigint.**
//! Every "signed" op is two's complement over [`SignedExt`] — the
//! unsigned `T` bit pattern plus one extension word for sign and
//! overflow:
//!
//! - signed `f + g`  ≡ wrapping add on `T`, carry (detected as
//!   `sum < a`) propagated into the extension word
//! - signed `-f`     ≡ `0 - f` with the borrow propagated likewise
//! - signed `f < 0`  ≡ MSB of the extension word is set
//! - arithmetic shr  ≡ `lo >> 1` with the extension word's low bit
//!   shifted into the top of `lo`
//!
//! ## Full-width carriers
//!
//! The modulus may occupy the full width of `T` — no headroom bit is
//! required, so a 2048-bit RSA modulus works in an exact 2048-bit
//! carrier. The signed intermediates can exceed `T` by up to two
//! bits; [`SignedExt`]'s extension word absorbs them without the bit
//! width ever being known at compile time. The modular coefficients
//! `d`, `e` stay in `[0, modulus)` and their update formulas are
//! written so no intermediate exceeds the carrier (see [`add_mod_ct`]
//! and [`half_mod_ct`]).
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

/// Two's-complement signed value spanning `T` plus one extension word.
///
/// Represents `hi · 2^k + lo` where `k` is `T`'s bit width, `lo` is
/// the unsigned low part and `hi` (a two's-complement word) carries
/// sign and overflow. The divsteps intermediates range over
/// `(-2·modulus, 2·modulus)` — up to two bits beyond what `T` holds
/// when the modulus occupies the full carrier width — so `hi` only
/// ever takes values in `{-2, -1, 0, 1}`.
#[derive(Clone)]
struct SignedExt<T> {
    lo: T,
    hi: u64,
}

#[inline]
fn se_select<T: ConditionallySelectable>(
    a: &SignedExt<T>,
    b: &SignedExt<T>,
    choice: Choice,
) -> SignedExt<T> {
    SignedExt {
        lo: T::conditional_select(&a.lo, &b.lo, choice),
        hi: u64::conditional_select(&a.hi, &b.hi, choice),
    }
}

/// Signed add over the extended representation. The carry out of the
/// low limb is detected as `sum < a` (unsigned wrap) and propagated
/// into the extension word.
#[inline]
fn se_add<T>(a: &SignedExt<T>, b: &SignedExt<T>) -> SignedExt<T>
where
    T: Clone + WrappingAdd<Output = T> + ConstantTimeLess,
{
    let lo = a.lo.clone().wrapping_add(b.lo.clone());
    let carry = lo.ct_lt(&a.lo);
    let hi =
        a.hi.wrapping_add(b.hi)
            .wrapping_add(carry.unwrap_u8() as u64);
    SignedExt { lo, hi }
}

/// Signed negate (`0 - a`) over the extended representation. The
/// borrow out of the low limb is 1 unless `lo == 0`.
#[inline]
fn se_neg<T>(a: &SignedExt<T>) -> SignedExt<T>
where
    T: Clone + Zero + WrappingSub<Output = T> + CtIsZero,
{
    let lo = T::zero().wrapping_sub(a.lo.clone());
    let borrow = !a.lo.ct_is_zero();
    let hi = 0u64
        .wrapping_sub(a.hi)
        .wrapping_sub(borrow.unwrap_u8() as u64);
    SignedExt { lo, hi }
}

/// Arithmetic right shift by 1 (sign-extending) over the extended
/// representation: the extension word's low bit shifts into the top
/// of `lo`, and the extension word itself shifts arithmetically.
#[inline]
fn se_shr1<T>(a: &SignedExt<T>, n_bits: usize) -> SignedExt<T>
where
    T: Clone
        + ConditionallySelectable
        + One
        + core::ops::Shr<usize, Output = T>
        + core::ops::Shl<usize, Output = T>
        + core::ops::BitOr<Output = T>,
{
    let logical = a.lo.clone() >> 1;
    let top_bit_mask = T::one() << (n_bits - 1);
    let with_hi_bit = logical.clone() | top_bit_mask;
    let hi_low_bit = Choice::from((a.hi & 1) as u8);
    SignedExt {
        lo: T::conditional_select(&logical, &with_hi_bit, hi_low_bit),
        hi: ((a.hi as i64) >> 1) as u64,
    }
}

/// Modular add: `(a + b) mod m`, assuming `a, b ∈ [0, m)`. Needs no
/// carrier headroom: `a + b ≥ m` iff `a ≥ m - b`, so select between
/// `a + b` (exact when the sum stays below `m`) and `a - (m - b)`
/// (exact otherwise) — every kept intermediate stays below `2^k`.
#[inline]
fn add_mod_ct<T>(a: &T, b: &T, m: &T) -> T
where
    T: Clone
        + ConditionallySelectable
        + WrappingAdd<Output = T>
        + WrappingSub<Output = T>
        + ConstantTimeLess,
{
    let m_minus_b = m.clone().wrapping_sub(b.clone());
    let need_sub = !a.ct_lt(&m_minus_b);
    let sum = a.clone().wrapping_add(b.clone());
    let reduced = a.clone().wrapping_sub(m_minus_b);
    T::conditional_select(&sum, &reduced, need_sub)
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
/// Via the standard branch-free trick:
/// - if `x` even: result is `x >> 1`.
/// - if `x` odd: result is `(x + m) / 2`, computed term-wise as
///   `(x >> 1) + (m >> 1) + 1` (exact since both are odd) so no
///   intermediate exceeds `m` — needs no carrier headroom over the
///   `x + m` sum.
#[inline]
fn half_mod_ct<T>(x: &T, m: &T) -> T
where
    T: CiosRowOps
        + Clone
        + ConditionallySelectable
        + One
        + WrappingAdd<Output = T>
        + core::ops::Shr<usize, Output = T>,
    T::Word: CtParity,
{
    let x_odd = ct_low_bit(x);
    let candidate_even = x.clone() >> 1;
    let candidate_odd = candidate_even
        .clone()
        .wrapping_add(m.clone() >> 1)
        .wrapping_add(T::one());
    T::conditional_select(&candidate_even, &candidate_odd, x_odd)
}

/// Compute `value⁻¹ mod modulus` in constant time over `value`. Works
/// for any odd modulus (composite or prime), including one that
/// occupies the full width of `T` — no carrier headroom is required.
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
    T::Word: CtParity,
{
    let n_bits = value.word_count() * core::mem::size_of::<T::Word>() * 8;
    let total_steps = divsteps_total(n_bits);

    let mut f = SignedExt {
        lo: modulus.clone(),
        hi: 0,
    };
    let mut g = SignedExt {
        lo: value.clone(),
        hi: 0,
    };
    let mut d = T::zero();
    let mut e = T::one();
    let mut delta: i64 = 1;

    for _ in 0..total_steps {
        let delta_pos = ct_i64_positive(delta);
        let g_odd = ct_low_bit(&g.lo);
        let swap = delta_pos & g_odd;

        // swap_choice: (f, g) ← (g, -f); (d, e) ← (e, -d mod m); delta ← -delta
        let new_f = se_select(&f, &g, swap);
        let neg_f = se_neg(&f);
        g = se_select(&g, &neg_f, swap);
        f = new_f;

        let new_d = T::conditional_select(&d, &e, swap);
        let neg_d_mod = neg_mod_ct(&d, modulus);
        e = T::conditional_select(&e, &neg_d_mod, swap);
        d = new_d;

        let neg_delta = (delta as u64).wrapping_neg() as i64;
        delta = i64::conditional_select(&delta, &neg_delta, swap);
        delta = delta.wrapping_add(1);

        // g_odd was determined before swap; after swap+negate the parity
        // of g equals that of the original g (both old g and -old_f are
        // odd when the swap fired, since old f is odd by loop invariant).
        // In the non-swap case g_odd tracks current g's parity directly.
        let g_odd_now = g_odd;
        let zero_ext = SignedExt {
            lo: T::zero(),
            hi: 0,
        };
        let to_add_to_g = se_select(&zero_ext, &f, g_odd_now);
        g = se_add(&g, &to_add_to_g);

        let add_to_e = add_mod_ct(&e, &d, modulus);
        e = T::conditional_select(&e, &add_to_e, g_odd_now);

        g = se_shr1(&g, n_bits);
        e = half_mod_ct(&e, modulus);
    }

    // After total_steps divsteps, g == 0 and f == ±gcd. For the
    // invertible case (gcd == 1), f is either +1 (lo == 1, hi == 0)
    // or -1 (lo all-ones, hi all-ones).
    let one = T::one();
    let neg_one_lo = T::zero().wrapping_sub(T::one());

    let f_is_one = f.lo.ct_eq(&one) & f.hi.ct_eq(&0u64);
    let f_is_neg_one = f.lo.ct_eq(&neg_one_lo) & f.hi.ct_eq(&u64::MAX);
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
    /// odd-prime modulus.
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

    /// Exhaustive sweep of the entire u8 carrier: every odd modulus
    /// (including all full-width, MSB-set ones) against every value.
    /// The baseline runs in u32 where `basic_mod_inv`'s signed
    /// intermediates have room.
    #[test]
    fn exhaustive_u8_all_odd_moduli() {
        for m in (3u16..=255).step_by(2) {
            for v in 1..m {
                let got = safegcd_inv_ct::<u8>(&(v as u8), &(m as u8)).into_option();
                let want = basic_mod_inv(v as u32, m as u32);
                match (got, want) {
                    (Some(g), Some(w)) => {
                        assert_eq!(g as u32, w, "value={v} modulus={m}: mismatch");
                    }
                    (None, None) => {}
                    _ => panic!("value={v} modulus={m}: Some/None disagreement"),
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

    /// Full-width u32 modulus (MSB set): 2^32 - 5, the largest 32-bit
    /// prime. Exercises the no-headroom path where the modulus
    /// occupies every bit of the carrier.
    #[test]
    fn u32_full_width_modulus() {
        let m: u32 = 0xFFFF_FFFB; // 2^32 - 5, prime
        for v in [1u32, 2, 7, 0xDEAD_BEEF, 0xFFFF_FFFA] {
            let got = safegcd_inv_ct::<u32>(&v, &m).into_option();
            let inv = got.unwrap_or_else(|| panic!("expected inverse for v={v:#x}"));
            assert_eq!(
                (inv as u64 * v as u64) % m as u64,
                1,
                "full-width u32 v={v:#x}: inv*v mod m != 1"
            );
        }
        // m - 1 ≡ -1 is its own inverse.
        assert_eq!(
            safegcd_inv_ct::<u32>(&(m - 1), &m).into_option(),
            Some(m - 1)
        );
    }

    /// Full-width u64 prime modulus: 2^64 - 59, the largest 64-bit
    /// prime.
    #[test]
    fn u64_full_width_prime() {
        let m: u64 = 0xFFFF_FFFF_FFFF_FFC5; // 2^64 - 59, prime
        for v in [1u64, 2, 7, 0xDEAD_BEEF, 0xCAFE_BABE_F00D_D00D, m - 1] {
            let got = safegcd_inv_ct::<u64>(&v, &m).into_option();
            let inv = got.unwrap_or_else(|| panic!("expected inverse for v={v:#x}"));
            let prod = (inv as u128 * v as u128) % m as u128;
            assert_eq!(prod, 1, "full-width u64 v={v:#x}: inv*v mod m != 1");
        }
    }

    /// Full-width u64 composite modulus — the RSA-blinding shape at
    /// carrier-exact width: n = p·q with p = 2^32 - 5, q = 2^32 - 17
    /// (both prime), so n has its MSB set. Coprime values invert;
    /// multiples of p return None.
    #[test]
    fn u64_full_width_composite() {
        let p: u64 = 0xFFFF_FFFB; // 2^32 - 5
        let q: u64 = 0xFFFF_FFEF; // 2^32 - 17
        let n = p * q;
        assert!(n >> 63 == 1, "test setup: n must have MSB set");
        for v in [1u64, 2, 3, 0xCAFE_BABE, 0xFEED_FACE_DEAD_BEEF % n] {
            let got = safegcd_inv_ct::<u64>(&v, &n).into_option();
            let inv = got.unwrap_or_else(|| panic!("expected inverse for v={v:#x}"));
            let prod = (inv as u128 * v as u128) % n as u128;
            assert_eq!(prod, 1, "full-width composite v={v:#x}: inv*v mod n != 1");
        }
        assert!(
            safegcd_inv_ct::<u64>(&p, &n).into_option().is_none(),
            "factor p must not be invertible mod p*q"
        );
    }

    /// Full-width multi-limb carrier: FixedUInt<u32, 2> holding the
    /// full 64-bit prime 2^64 - 59. Cross-checked against the u64
    /// primitive path (validated independently above).
    #[test]
    fn fixed_bigint_full_width_modulus() {
        use const_num_traits::Ct;
        use fixed_bigint::FixedUInt;
        type U64Ct = FixedUInt<u32, 2, Ct>;
        type U64Nct = FixedUInt<u32, 2>;

        let m_raw: u64 = 0xFFFF_FFFF_FFFF_FFC5; // 2^64 - 59, prime
        let m = U64Ct::from(m_raw);
        for v_raw in [1u64, 2, 7, 42, 0xDEAD_BEEF, m_raw - 1] {
            let v = U64Ct::from(v_raw);
            let got = safegcd_inv_ct(&v, &m).into_option();
            let inv = got.unwrap_or_else(|| panic!("expected inverse for v={v_raw:#x}"));
            let expected = safegcd_inv_ct::<u64>(&v_raw, &m_raw)
                .into_option()
                .expect("u64 baseline");
            let inv_nct: U64Nct = inv.forget_ct();
            assert_eq!(
                inv_nct,
                U64Nct::from(expected),
                "FixedUInt full-width v={v_raw:#x}: mismatch with u64 baseline"
            );
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
