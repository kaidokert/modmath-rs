//! `Field<T, Ct>` typestate fixtures — the surface downstream consumers
//! (rsa_heapless blinding) actually call. `safegcd_inv_ct` has no public
//! free-function path, so its coverage routes through
//! `Field::inv_safegcd_ct`, which is also the deployment call site.
//!
//! `try_new_odd_ct(..).unwrap()` inside a fixture is deliberate: the
//! parity `Choice` derives from the modulus's low bit, which the fixture
//! defines via `|= 1`, so the unwrap branches on a *defined* bit —
//! Valgrind-clean — and never actually panics.
//!
//! `inv_safegcd_ct`'s `CtOption<Residue>` can't be unwrapped here: its
//! `is_some` mask derives from fully-tainted gcd state, so any branch on
//! it is a real violation (and a potential panic). The fixtures observe
//! the whole `CtOption` through `bb` — keeping the divsteps loop alive
//! against DCE — and export only the mask byte.

use crate::bb;
use const_num_traits::Ct;
use fixed_bigint::FixedUInt;
use modmath::Field;

type Fb256 = FixedUInt<u32, 8, Ct>;

/// CT modular exponentiation through the shipped `Field<T, Ct>::exp`
/// ladder. Modulus public (forced odd via `|= 1`, so `try_new_odd_ct`'s
/// `CtOption` unwrap branches on a fixture-defined bit); base and
/// exponent secret.
#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__field_exp__u64__N1(
    base: *const u64,
    e: *const u64,
    m: *const u64,
    out: *mut u64,
) {
    let base = bb(unsafe { *base });
    let e = bb(unsafe { *e });
    let m = bb(unsafe { *m }) | 1;
    let f = Field::<u64, Ct>::try_new_odd_ct(m).unwrap();
    let r = f.exp(&f.reduce(&base), &e);
    unsafe { *out = bb(f.into_raw(&r)) }
}

#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__field_exp__fb32__N8(
    base: *const [u32; 8],
    e: *const [u32; 8],
    m: *const [u32; 8],
    out: *mut [u32; 8],
) {
    let base = Fb256::from(bb(unsafe { *base }));
    let e = Fb256::from(bb(unsafe { *e }));
    let mut mw = bb(unsafe { *m });
    mw[0] |= 1;
    let f = Field::<Fb256, Ct>::try_new_odd_ct(Fb256::from(mw)).unwrap();
    let r = f.exp(&f.reduce(&base), &e);
    unsafe { *out = bb(*f.into_raw(&r).words()) }
}

#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__field_inv_safegcd__u64__N1(m: *const u64, v: *const u64) -> u8 {
    let m = bb(unsafe { *m }) | 1;
    let v = bb(unsafe { *v });
    let f = Field::<u64, Ct>::try_new_odd_ct(m).unwrap();
    let r = f.reduce(&v);
    let inv = bb(f.inv_safegcd_ct(&r));
    bb(inv.is_some().unwrap_u8())
}

#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__field_inv_safegcd__fb32__N8(
    m: *const [u32; 8],
    v: *const [u32; 8],
) -> u8 {
    let mut mw = bb(unsafe { *m });
    mw[0] |= 1;
    let v = Fb256::from(bb(unsafe { *v }));
    let f = Field::<Fb256, Ct>::try_new_odd_ct(Fb256::from(mw)).unwrap();
    let r = f.reduce(&v);
    let inv = bb(f.inv_safegcd_ct(&r));
    bb(inv.is_some().unwrap_u8())
}

/// The RSA-blinding-shaped flow: secret modulus precompute → reduce →
/// CT ladder exp → CT inverse → back to raw. One fixture, every stage.
#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__field_blind_path__fb32__N8(
    m: *const [u32; 8],
    x: *const [u32; 8],
    e: *const [u32; 8],
    out: *mut [u32; 8],
) -> u8 {
    let mut mw = bb(unsafe { *m });
    mw[0] |= 1;
    let x = Fb256::from(bb(unsafe { *x }));
    let e = Fb256::from(bb(unsafe { *e }));
    let f = Field::<Fb256, Ct>::try_new_odd_ct(Fb256::from(mw)).unwrap();
    let r = f.reduce(&x);
    let powed = f.exp(&r, &e);
    let inv = bb(f.inv_safegcd_ct(&powed));
    let raw = f.into_raw(&powed);
    unsafe { *out = bb(*raw.words()) }
    bb(inv.is_some().unwrap_u8())
}

/// Ladder-support primitives: `Residue::cswap` on a tainted choice plus
/// `Residue::ct_eq` on the swapped values.
#[unsafe(no_mangle)]
pub extern "C" fn ct_fix__field_cswap_eq__fb32__N8(
    m: *const [u32; 8],
    a: *const [u32; 8],
    b: *const [u32; 8],
    choice: *const u8,
    out: *mut [u32; 8],
) -> u8 {
    let mut mw = bb(unsafe { *m });
    mw[0] |= 1;
    let a = Fb256::from(bb(unsafe { *a }));
    let b = Fb256::from(bb(unsafe { *b }));
    let choice = subtle::Choice::from(bb(unsafe { *choice }));
    let f = Field::<Fb256, Ct>::try_new_odd_ct(Fb256::from(mw)).unwrap();
    let mut ra = f.reduce(&a);
    let mut rb = f.reduce(&b);
    modmath::Residue::cswap(choice, &mut ra, &mut rb);
    let eq = ra.ct_eq(&rb);
    unsafe { *out = bb(*f.into_raw(&ra).words()) }
    bb(eq.unwrap_u8())
}
