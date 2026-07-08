//! Asymmetric fixtures: the public operand is an in-body constant, so
//! LLVM is free to const-propagate it against the secret input. This is
//! the configuration where the XOR-select idiom historically gets
//! rewritten into a conditional move whose flag comes from the secret —
//! symmetric-taint fixtures never expose it because with both operands
//! opaque there is nothing to propagate.

use crate::bb;
use const_num_traits::Ct;
use fixed_bigint::FixedUInt;
use modmath::Field;

type Fb256 = FixedUInt<u32, 8, Ct>;

/// 2^64 − 59, prime. Value is irrelevant to the gate — what matters is
/// that it is a compile-time constant with the low bit set.
const M64: u64 = 0xFFFF_FFFF_FFFF_FFC5;

/// 2^256 − 59 as little-endian u32 limbs. Odd; that is all the fixture
/// needs.
const M256: [u32; 8] = [
    0xFFFF_FFC5,
    0xFFFF_FFFF,
    0xFFFF_FFFF,
    0xFFFF_FFFF,
    0xFFFF_FFFF,
    0xFFFF_FFFF,
    0xFFFF_FFFF,
    0xFFFF_FFFF,
];

/// Public modulus, secret value — the plain-RSA blinding shape (modulus
/// is the public N; the blinding factor is secret).
#[no_mangle]
pub extern "C" fn ct_fix__ASYM__field_inv_safegcd_public_m__u64__N1(v: *const u64) -> u8 {
    let v = bb(unsafe { *v });
    let f = Field::<u64, Ct>::try_new_odd_ct(M64).unwrap();
    let r = f.reduce(&v);
    let inv = bb(f.inv_safegcd_ct(&r));
    bb(inv.is_some().unwrap_u8())
}

/// Public modulus and base, secret exponent — the classic ladder-leak
/// shape (fixed-base exponentiation with a private exponent).
#[no_mangle]
pub extern "C" fn ct_fix__ASYM__field_exp_secret_e__fb32__N8(
    e: *const [u32; 8],
    out: *mut [u32; 8],
) {
    let e = Fb256::from(bb(unsafe { *e }));
    let f = Field::<Fb256, Ct>::try_new_odd_ct(Fb256::from(M256)).unwrap();
    let base = f.reduce(&Fb256::from(3u8));
    let powed = f.exp(&base, &e);
    unsafe { *out = bb(*f.into_raw(&powed).words()) }
}

/// Constant residues, secret swap flag. Exercises the
/// `conditional_swap` path with everything except the flag visible to
/// the optimizer.
#[no_mangle]
pub extern "C" fn ct_fix__ASYM__field_cswap_choice__fb32__N8(
    choice: *const u8,
    out: *mut [u32; 8],
) {
    let choice = subtle::Choice::from(bb(unsafe { *choice }));
    let f = Field::<Fb256, Ct>::try_new_odd_ct(Fb256::from(M256)).unwrap();
    let mut ra = f.reduce(&Fb256::from(2u8));
    let mut rb = f.reduce(&Fb256::from(5u8));
    modmath::Residue::cswap(choice, &mut ra, &mut rb);
    let a_out = f.into_raw(&ra);
    let b_out = bb(f.into_raw(&rb));
    unsafe { *out = bb(*a_out.words()) }
    let _ = bb(*b_out.words());
}
