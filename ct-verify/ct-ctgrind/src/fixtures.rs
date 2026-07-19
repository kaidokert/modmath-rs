//! Taint wrappers, one per `ct-fixtures` symbol. Each mirrors its
//! fixture's ABI one-for-one; extending the fixture set means adding
//! the symbol there and the matching wrapper here.
//!
//! Zeroed buffers are deliberate: Valgrind's V-bits carry the "secret"
//! marking, not the values, and zero inputs are the ones proven not to
//! reach any panic path (see the fixture-side docs on forced-odd
//! moduli).

use core::hint::black_box;
use krabi_caliper::ctgrind_fixture;
use krabi_caliper::host::ctgrind::{taint_val, untaint_val};

// ============================================================================
// Core primitives — u64 carrier.
// ============================================================================

macro_rules! quad_u64 {
    ($name:ident) => {
        unsafe extern "C" {
            fn $name(a: *const u64, b: *const u64, m: *const u64, np: *const u64, out: *mut u64);
        }
        ctgrind_fixture!($name, {
            let (a, b, m, np) = (0u64, 0u64, 0u64, 0u64);
            let mut out = 0u64;
            taint_val(&a);
            taint_val(&b);
            taint_val(&m);
            taint_val(&np);
            unsafe { $name(&a, &b, &m, &np, &mut out) }
            untaint_val(&out);
            let _ = black_box(out);
        });
    };
}

quad_u64!(ct_fix__wide_mont_mul__u64__N1);
quad_u64!(ct_fix__wide_redc__u64__N1);
quad_u64!(ct_fix__cios_mont_mul__u64__N1);

macro_rules! tri_u64 {
    ($name:ident) => {
        unsafe extern "C" {
            fn $name(a: *const u64, b: *const u64, c: *const u64, out: *mut u64);
        }
        ctgrind_fixture!($name, {
            let (a, b, c) = (0u64, 0u64, 0u64);
            let mut out = 0u64;
            taint_val(&a);
            taint_val(&b);
            taint_val(&c);
            unsafe { $name(&a, &b, &c, &mut out) }
            untaint_val(&out);
            let _ = black_box(out);
        });
    };
}

// The modexp taint model mirrors the documented contract: CT over
// `base` and `exponent`, modulus **public** (its rustdoc says the NCT
// precompute helpers are intentional — R/R² derivation branches on
// modulus magnitude). Tainting the modulus here flags that precompute,
// which the contract explicitly permits; secret-modulus exponentiation
// is covered by the `Field::try_new_odd_ct`-based fixtures instead.
unsafe extern "C" {
    fn ct_fix__field_exp__u64__N1(base: *const u64, e: *const u64, m: *const u64, out: *mut u64);
}
ctgrind_fixture!(ct_fix__field_exp__u64__N1, {
    let (base, e, m) = (0u64, 0u64, 0u64);
    let mut out = 0u64;
    taint_val(&base);
    taint_val(&e);
    unsafe { ct_fix__field_exp__u64__N1(&base, &e, &m, &mut out) }
    untaint_val(&out);
    let _ = black_box(out);
});

unsafe extern "C" {
    fn ct_fix__field_inv_safegcd__u64__N1(m: *const u64, v: *const u64) -> u8;
}
ctgrind_fixture!(ct_fix__field_inv_safegcd__u64__N1, {
    let (m, v) = (0u64, 0u64);
    taint_val(&m);
    taint_val(&v);
    let r = unsafe { ct_fix__field_inv_safegcd__u64__N1(&m, &v) };
    untaint_val(&r);
    let _ = black_box(r);
});

// ============================================================================
// Core primitives — FixedUInt<u32, 8, Ct> carrier.
// ============================================================================

type W8 = [u32; 8];

unsafe extern "C" {
    fn ct_fix__cios_mont_mul__fb32__N8(
        a: *const W8,
        b: *const W8,
        m: *const W8,
        np: *const W8,
        out: *mut W8,
    );
}
ctgrind_fixture!(ct_fix__cios_mont_mul__fb32__N8, {
    let (a, b, m, np) = ([0u32; 8], [0u32; 8], [0u32; 8], [0u32; 8]);
    let mut out = [0u32; 8];
    taint_val(&a);
    taint_val(&b);
    taint_val(&m);
    taint_val(&np);
    unsafe { ct_fix__cios_mont_mul__fb32__N8(&a, &b, &m, &np, &mut out) }
    untaint_val(&out);
    let _ = black_box(out);
});

unsafe extern "C" {
    fn ct_fix__field_exp__fb32__N8(base: *const W8, e: *const W8, m: *const W8, out: *mut W8);
}
// Modulus untainted — same public-modulus contract as the u64 variant
// above.
ctgrind_fixture!(ct_fix__field_exp__fb32__N8, {
    let (base, e, m) = ([0u32; 8], [0u32; 8], [0u32; 8]);
    let mut out = [0u32; 8];
    taint_val(&base);
    taint_val(&e);
    unsafe { ct_fix__field_exp__fb32__N8(&base, &e, &m, &mut out) }
    untaint_val(&out);
    let _ = black_box(out);
});

// ============================================================================
// Field typestate flows.
// ============================================================================

unsafe extern "C" {
    fn ct_fix__field_inv_safegcd__fb32__N8(m: *const W8, v: *const W8) -> u8;
}
ctgrind_fixture!(ct_fix__field_inv_safegcd__fb32__N8, {
    let (m, v) = ([0u32; 8], [0u32; 8]);
    taint_val(&m);
    taint_val(&v);
    let r = unsafe { ct_fix__field_inv_safegcd__fb32__N8(&m, &v) };
    untaint_val(&r);
    let _ = black_box(r);
});

unsafe extern "C" {
    fn ct_fix__field_blind_path__fb32__N8(
        m: *const W8,
        x: *const W8,
        e: *const W8,
        out: *mut W8,
    ) -> u8;
}
ctgrind_fixture!(ct_fix__field_blind_path__fb32__N8, {
    let (m, x, e) = ([0u32; 8], [0u32; 8], [0u32; 8]);
    let mut out = [0u32; 8];
    taint_val(&m);
    taint_val(&x);
    taint_val(&e);
    let r = unsafe { ct_fix__field_blind_path__fb32__N8(&m, &x, &e, &mut out) };
    untaint_val(&out);
    untaint_val(&r);
    let _ = black_box(out);
    let _ = black_box(r);
});

unsafe extern "C" {
    fn ct_fix__field_cswap_eq__fb32__N8(
        m: *const W8,
        a: *const W8,
        b: *const W8,
        choice: *const u8,
        out: *mut W8,
    ) -> u8;
}
ctgrind_fixture!(ct_fix__field_cswap_eq__fb32__N8, {
    let (m, a, b) = ([0u32; 8], [0u32; 8], [0u32; 8]);
    let choice = 0u8;
    let mut out = [0u32; 8];
    taint_val(&m);
    taint_val(&a);
    taint_val(&b);
    taint_val(&choice);
    let r = unsafe { ct_fix__field_cswap_eq__fb32__N8(&m, &a, &b, &choice, &mut out) };
    untaint_val(&out);
    untaint_val(&r);
    let _ = black_box(out);
    let _ = black_box(r);
});

// ============================================================================
// Asymmetric fixtures — only the secret slot is tainted; the public
// operands are constants on the symbol side.
// ============================================================================

unsafe extern "C" {
    fn ct_fix__ASYM__field_inv_safegcd_public_m__u64__N1(v: *const u64) -> u8;
}
ctgrind_fixture!(ct_fix__ASYM__field_inv_safegcd_public_m__u64__N1, {
    let v = 0u64;
    taint_val(&v);
    let r = unsafe { ct_fix__ASYM__field_inv_safegcd_public_m__u64__N1(&v) };
    untaint_val(&r);
    let _ = black_box(r);
});

unsafe extern "C" {
    fn ct_fix__ASYM__field_exp_secret_e__fb32__N8(e: *const W8, out: *mut W8);
}
ctgrind_fixture!(ct_fix__ASYM__field_exp_secret_e__fb32__N8, {
    let e = [0u32; 8];
    let mut out = [0u32; 8];
    taint_val(&e);
    unsafe { ct_fix__ASYM__field_exp_secret_e__fb32__N8(&e, &mut out) }
    untaint_val(&out);
    let _ = black_box(out);
});

unsafe extern "C" {
    fn ct_fix__ASYM__field_cswap_choice__fb32__N8(choice: *const u8, out: *mut W8);
}
ctgrind_fixture!(ct_fix__ASYM__field_cswap_choice__fb32__N8, {
    let choice = 0u8;
    let mut out = [0u32; 8];
    taint_val(&choice);
    unsafe { ct_fix__ASYM__field_cswap_choice__fb32__N8(&choice, &mut out) }
    untaint_val(&out);
    let _ = black_box(out);
});

// ============================================================================
// Negative controls.
// ============================================================================

unsafe extern "C" {
    fn nct_fix__neg__eea_inv__u64__N1(a: *const u64, m: *const u64, out: *mut u64);
}
ctgrind_fixture!(nct_fix__neg__eea_inv__u64__N1, {
    let (a, m) = (0u64, 0u64);
    let mut out = 0u64;
    taint_val(&a);
    taint_val(&m);
    unsafe { nct_fix__neg__eea_inv__u64__N1(&a, &m, &mut out) }
    untaint_val(&out);
    let _ = black_box(out);
});

tri_u64!(nct_fix__neg__schoolbook_exp__u64__N1);

unsafe extern "C" {
    fn nct_fix__neg__table_lookup__u64__N1(a: *const u64, out: *mut u64);
}
ctgrind_fixture!(nct_fix__neg__table_lookup__u64__N1, {
    let a = 0u64;
    let mut out = 0u64;
    taint_val(&a);
    unsafe { nct_fix__neg__table_lookup__u64__N1(&a, &mut out) }
    untaint_val(&out);
    let _ = black_box(out);
});
