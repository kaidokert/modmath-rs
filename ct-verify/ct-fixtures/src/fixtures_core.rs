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

// Modexp moved to the shipped `Field<T, Ct>::exp` surface — see
// `ct_fix__field_exp__*` in `fixtures_field`. The retired
// `montgomery::ct::pre_reduced` free-function is no longer part of the
// public surface, so there is nothing to taint here.
