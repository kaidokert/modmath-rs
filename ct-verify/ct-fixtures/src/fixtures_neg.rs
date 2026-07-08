//! Negative controls — Nct primitives and hand-written leaks that MUST
//! trip the taint gate. If any of these stops tripping, the harness is
//! broken (compiler folded the body away, taint marks not applied, or
//! Valgrind not actually running), not the library fixed.
//!
//! All three leak through mechanisms cmov-rewrites can't neutralize:
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
#[no_mangle]
pub extern "C" fn nct_fix__neg__eea_inv__u64__N1(a: *const u64, m: *const u64, out: *mut u64) {
    let a = bb(unsafe { *a });
    let m = bb(unsafe { *m }) | 1;
    let r = modmath::basic::inv(a, m).unwrap_or(0);
    unsafe { *out = bb(r) }
}

/// Schoolbook square-and-multiply: branches on each exponent bit.
#[no_mangle]
pub extern "C" fn nct_fix__neg__schoolbook_exp__u64__N1(
    base: *const u64,
    e: *const u64,
    m: *const u64,
    out: *mut u64,
) {
    let base = bb(unsafe { *base });
    let e = bb(unsafe { *e });
    let m = bb(unsafe { *m }) | 1;
    let r = modmath::basic::exp(base, e, m);
    unsafe { *out = bb(r) }
}

/// Secret-indexed table load — the class of leak asm-grep can never see
/// (no branch involved) and the reason the taint layer exists.
#[no_mangle]
pub extern "C" fn nct_fix__neg__table_lookup__u64__N1(a: *const u64, out: *mut u64) {
    static TABLE: [u64; 16] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
    let a = bb(unsafe { *a });
    let r = TABLE[(a & 0xF) as usize];
    unsafe { *out = bb(r) }
}
