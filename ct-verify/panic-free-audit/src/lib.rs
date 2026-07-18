//! Linker-DCE audit for modmath's CT surface.
//!
//! Each `#[no_mangle] pub extern "C"` symbol exercises one CT entry
//! point the way a deployed consumer would — no `unwrap`, no `expect`,
//! `Odd::new_unchecked` at the trust boundary, `CtOption` results
//! observed through `black_box` instead of extracted. After
//! cross-building with the workspace release profile, `check.sh`
//! asserts the resulting archive contains no `core::panicking`
//! machinery: for CT code a reachable panic is both a DoS edge and a
//! timing oracle (the panic formatting path's cost depends on the
//! values being formatted).
//!
//! `black_box` at every boundary so the optimizer can't fold inputs
//! through and DCE the body wholesale (same hygiene as ct-fixtures).

// no_std + the local #[panic_handler] only under the `panic-handler`
// feature (the cross-built audit shape, enabled by check.sh) — host-side
// workspace builds (clippy, tests) link std, which supplies its own.
#![cfg_attr(feature = "panic-handler", no_std)]
#![allow(non_snake_case)]
#![allow(clippy::not_unsafe_ptr_arg_deref)]

use core::hint::black_box;

use const_num_traits::Ct;
use fixed_bigint::FixedUInt;
use modmath::basic::montgomery::wide::ct as wide_ct;
use modmath::{CiosMontMulCt, Field};

type Fb256 = FixedUInt<u32, 8, Ct>;

#[unsafe(no_mangle)]
pub extern "C" fn panic_audit__wide_mont_mul_ct__u64(
    a: u64,
    b: u64,
    m: u64,
    np: u64,
    out: *mut u64,
) {
    let r = wide_ct::mul(black_box(a), black_box(b), black_box(m) | 1, black_box(np));
    unsafe { *out = black_box(r) }
}

#[unsafe(no_mangle)]
pub extern "C" fn panic_audit__wide_redc_ct__u64(
    t_lo: u64,
    t_hi: u64,
    m: u64,
    np: u64,
    out: *mut u64,
) {
    let r = wide_ct::redc(
        black_box(t_lo),
        black_box(t_hi),
        black_box(m) | 1,
        black_box(np),
    );
    unsafe { *out = black_box(r) }
}

#[unsafe(no_mangle)]
pub extern "C" fn panic_audit__cios_mont_mul_ct__fb32(
    a: *const [u32; 8],
    b: *const [u32; 8],
    m: *const [u32; 8],
    np: *const [u32; 8],
    out: *mut [u32; 8],
) {
    let a = Fb256::from(black_box(unsafe { *a }));
    let b = Fb256::from(black_box(unsafe { *b }));
    let mut mw = black_box(unsafe { *m });
    mw[0] |= 1;
    let m = Fb256::from(mw);
    let np = Fb256::from(black_box(unsafe { *np }));
    let r = CiosMontMulCt::cios_mont_mul_ct(&a, &b, &m, &np);
    unsafe { *out = black_box(*r.words()) }
}

#[unsafe(no_mangle)]
pub extern "C" fn panic_audit__field_exp__fb32(
    base: *const [u32; 8],
    e: *const [u32; 8],
    m: *const [u32; 8],
    out: *mut [u32; 8],
) {
    let base = Fb256::from(black_box(unsafe { *base }));
    let e = Fb256::from(black_box(unsafe { *e }));
    let mut mw = black_box(unsafe { *m });
    mw[0] |= 1;
    // SAFETY: low bit forced above. Infallible `new_odd_ct` (no `.unwrap()`)
    // keeps the panic-symbol audit clean.
    let proof = unsafe { modmath::Odd::new_unchecked(Fb256::from(mw)) };
    let f = Field::<Fb256, Ct>::new_odd_ct(proof);
    let r = f.exp(&f.reduce(&base), &e);
    unsafe { *out = black_box(*f.into_raw(&r).words()) }
}

/// The constructor a consumer calls with a secret modulus. The returned
/// `CtOption<Field>` is observed whole — extracting it is the caller's
/// panic site, not the library's.
#[unsafe(no_mangle)]
pub extern "C" fn panic_audit__try_new_odd_ct__fb32(m: *const [u32; 8]) -> u8 {
    let mut mw = black_box(unsafe { *m });
    mw[0] |= 1;
    let f = Field::<Fb256, Ct>::try_new_odd_ct(Fb256::from(mw));
    let is_some = f.is_some().unwrap_u8();
    let _ = black_box(f);
    black_box(is_some)
}

/// The full blinding-shaped pipeline with zero extraction points:
/// infallible construction via the `Odd` proof, reduce → exp →
/// inv_safegcd_ct (observed, not unwrapped) → into_raw, plus the
/// ladder-support cswap / ct_eq.
#[unsafe(no_mangle)]
pub extern "C" fn panic_audit__field_ct_pipeline__fb32(
    m: *const [u32; 8],
    x: *const [u32; 8],
    e: *const [u32; 8],
    choice: u8,
    out: *mut [u32; 8],
) -> u8 {
    let mut mw = black_box(unsafe { *m });
    mw[0] |= 1;
    let x = Fb256::from(black_box(unsafe { *x }));
    let e = Fb256::from(black_box(unsafe { *e }));
    // SAFETY: low bit forced above.
    let proof = unsafe { modmath::Odd::new_unchecked(Fb256::from(mw)) };
    let f = Field::<Fb256, Ct>::new_odd_ct(proof);

    let mut r = f.reduce(&x);
    let mut powed = f.exp(&r, &e);
    let inv = f.inv_safegcd_ct(&powed);
    let is_some = inv.is_some().unwrap_u8();
    let _ = black_box(inv);
    modmath::Residue::cswap(
        subtle::Choice::from(black_box(choice) & 1),
        &mut r,
        &mut powed,
    );
    let eq = r.ct_eq(&powed);
    unsafe { *out = black_box(*f.into_raw(&r).words()) }
    black_box(is_some ^ eq.unwrap_u8())
}

#[cfg(feature = "panic-handler")]
#[panic_handler]
fn panic(_: &core::panic::PanicInfo) -> ! {
    loop {}
}
