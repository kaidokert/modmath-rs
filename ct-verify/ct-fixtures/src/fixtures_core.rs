//! Free-function CT primitives: wide-REDC family, CIOS, Montgomery
//! modexp. For the mul/REDC entries every input is secret from the
//! driver's point of view — `n_prime` is derived from the modulus, so
//! it carries the same secrecy as the modulus itself. The modexp
//! entries document the modulus as public (their precompute is
//! intentionally NCT), so the driver taints only base and exponent
//! there; secret-modulus exponentiation is covered on the
//! `Field::try_new_odd_ct` path in `fixtures_field`.

use crate::bb;
use const_num_traits::Ct;
use fixed_bigint::FixedUInt;
use modmath::CiosMontMulCt;
use modmath::basic::montgomery::ct::pre_reduced as mont_ct_pr;
use modmath::basic::montgomery::wide::ct as wide_ct;

type Fb256 = FixedUInt<u32, 8, Ct>;

#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__wide_mont_mul__u64__N1(
    a: *const u64,
    b: *const u64,
    m: *const u64,
    np: *const u64,
    out: *mut u64,
) {
    let a = bb(unsafe { *a });
    let b = bb(unsafe { *b });
    let m = bb(unsafe { *m }) | 1;
    let np = bb(unsafe { *np });
    let r = wide_ct::mul(a, b, m, np);
    unsafe { *out = bb(r) }
}

#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__wide_redc__u64__N1(
    t_lo: *const u64,
    t_hi: *const u64,
    m: *const u64,
    np: *const u64,
    out: *mut u64,
) {
    let t_lo = bb(unsafe { *t_lo });
    let t_hi = bb(unsafe { *t_hi });
    let m = bb(unsafe { *m }) | 1;
    let np = bb(unsafe { *np });
    let r = wide_ct::redc(t_lo, t_hi, m, np);
    unsafe { *out = bb(r) }
}

#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__cios_mont_mul__u64__N1(
    a: *const u64,
    b: *const u64,
    m: *const u64,
    np: *const u64,
    out: *mut u64,
) {
    let a = bb(unsafe { *a });
    let b = bb(unsafe { *b });
    let m = bb(unsafe { *m }) | 1;
    let np = bb(unsafe { *np });
    let r = CiosMontMulCt::cios_mont_mul_ct(&a, &b, &m, &np);
    unsafe { *out = bb(r) }
}

#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__cios_mont_mul__fb32__N8(
    a: *const [u32; 8],
    b: *const [u32; 8],
    m: *const [u32; 8],
    np: *const [u32; 8],
    out: *mut [u32; 8],
) {
    let a = Fb256::from(bb(unsafe { *a }));
    let b = Fb256::from(bb(unsafe { *b }));
    let mut mw = bb(unsafe { *m });
    mw[0] |= 1;
    let m = Fb256::from(mw);
    let np = Fb256::from(bb(unsafe { *np }));
    let r = CiosMontMulCt::cios_mont_mul_ct(&a, &b, &m, &np);
    unsafe { *out = bb(*r.words()) }
}

#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__mod_exp_pr_odd__u64__N1(
    base: *const u64,
    e: *const u64,
    m: *const u64,
    out: *mut u64,
) {
    let base = bb(unsafe { *base });
    let e = bb(unsafe { *e });
    let m = bb(unsafe { *m }) | 1;
    // SAFETY: low bit forced above; the proof is honest.
    let m = unsafe { modmath::Odd::new_unchecked(m) };
    let r = mont_ct_pr::mod_exp_odd(base, e, m);
    unsafe { *out = bb(r) }
}

#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__mod_exp_pr_odd__fb32__N8(
    base: *const [u32; 8],
    e: *const [u32; 8],
    m: *const [u32; 8],
    out: *mut [u32; 8],
) {
    let base = Fb256::from(bb(unsafe { *base }));
    let e = Fb256::from(bb(unsafe { *e }));
    let mut mw = bb(unsafe { *m });
    mw[0] |= 1;
    // SAFETY: low bit forced above; the proof is honest.
    let m = unsafe { modmath::Odd::new_unchecked(Fb256::from(mw)) };
    let r = mont_ct_pr::mod_exp_odd(base, e, m);
    unsafe { *out = bb(*r.words()) }
}
