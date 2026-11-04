//! Negative controls — Nct primitives and hand-written leaks that MUST
//! trip the taint gate. If any of these stops tripping, the harness is
//! broken (compiler folded the body away, taint marks not applied, or
//! Valgrind not actually running), not the library fixed.
//!
//! All of them leak through mechanisms cmov-rewrites can't neutralize:
//! data-dependent *loop* conditions (memcheck must flag the branch) and
//! a secret-derived *address* (memcheck flags the access). A plain
//! `if a < b` canary would be unreliable here — LLVM freely rewrites it
//! to a conditional move, which memcheck only propagates through rather
//! than reporting.

use crate::bb;

/// Extended-Euclid inverse: the loop count and the per-iteration
/// quotients depend on operand magnitudes. Value 0 for `a` skips the
/// loop body at runtime (so the internal division never executes on a
/// zero divisor), but the loop-condition branch itself already reads
/// tainted data.
#[unsafe(no_mangle)]
pub extern "C" fn nct_fix__neg__eea_inv__u64__N1(a: *const u64, m: *const u64, out: *mut u64) {
    let a = bb(unsafe { *a });
    let m = bb(unsafe { *m }) | 1;
    // In-wrapper data-dependent guard: the callee's leaky loops may sit
    // in a non-inlined helper symbol, but the asm-grep gate asserts on
    // THIS symbol's body, so it must contain a visible branch on data.
    if a == u64::MAX {
        unsafe { *out = 1 }
        return;
    }
    let r = modmath::basic::inv(a, m).unwrap_or(0);
    unsafe { *out = bb(r) }
}

/// Schoolbook square-and-multiply: loop length and per-bit work depend
/// on the exponent.
#[unsafe(no_mangle)]
pub extern "C" fn nct_fix__neg__schoolbook_exp__u64__N1(
    base: *const u64,
    e: *const u64,
    m: *const u64,
    out: *mut u64,
) {
    let base = bb(unsafe { *base });
    let e = bb(unsafe { *e });
    let m = bb(unsafe { *m }) | 1;
    // Same in-wrapper guard rationale as `nct_fix__neg__eea_inv`.
    if e == u64::MAX {
        unsafe { *out = 1 }
        return;
    }
    let r = modmath::basic::exp(base, e, m);
    unsafe { *out = bb(r) }
}

/// Secret-indexed table load — the class of leak asm-grep can never see
/// (no branch involved) and the reason the taint layer exists.
#[unsafe(no_mangle)]
pub extern "C" fn nct_fix__neg__table_lookup__u64__N1(a: *const u64, out: *mut u64) {
    static TABLE: [u64; 16] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
    let a = bb(unsafe { *a });
    let r = TABLE[(a & 0xF) as usize];
    unsafe { *out = bb(r) }
}

/// Index of the first zero limb. The `i < 4` exit is a counted loop the
/// asm-grep gate accepts; the `a[i] != 0` exit is decided by loaded data, so
/// this trips only if counted-loop acceptance rejects data-decided loops.
#[unsafe(no_mangle)]
pub extern "C" fn nct_fix__neg__value_loop__u32__N4(a: *const [u32; 4]) -> u32 {
    let a = bb(unsafe { *a });
    let mut i = 0;
    while i < 4 && a[i] != 0 {
        i += 1;
    }
    bb(i as u32)
}

/// `while d < n { d <<= 1 }`, the normalization step of binary division and
/// gcd. The doubling counter has the shape of a public ladder, so this trips
/// only if the gate checks where the counter's value came from.
#[unsafe(no_mangle)]
pub extern "C" fn nct_fix__neg__doubling_normalize__u32__N1(d: *const u32, n: *const u32) -> u32 {
    let mut d = bb(unsafe { *d }) | 1;
    let n = bb(unsafe { *n });
    let mut shifts = 0u32;
    while d < n && d >> 31 == 0 {
        d <<= 1;
        shifts += 1;
    }
    bb(shifts)
}
