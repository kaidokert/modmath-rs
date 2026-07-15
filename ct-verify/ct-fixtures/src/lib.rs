//! Inspection fixtures for the modmath CT-verify harness.
//!
//! Each fixture is a `#[no_mangle] pub extern "C"` symbol exercising one
//! CT entry point at one concrete carrier. The taint driver (`ct-ctgrind`)
//! calls each symbol with its secret inputs tagged `MAKE_MEM_UNDEFINED`
//! under Valgrind; a future asm-grep driver disassembles the same symbols
//! per embedded target.
//!
//! # The `black_box` discipline
//!
//! `#[no_mangle] pub extern "C"` prevents inlining of the *symbol* but
//! says nothing about the body: without `core::hint::black_box` at both
//! ends the optimizer folds pointer loads through, sees the result is
//! never observed, and DCEs the computation — a false pass. Every input
//! goes through `bb` immediately after the load, every output goes
//! through `bb` before the store. Reviewers must reject any fixture
//! that skips this.
//!
//! The surface is small and heterogeneous — field construction, mixed
//! carrier/word types, multi-step flows — so fixtures are written
//! longhand rather than macro-generated; the discipline above is the
//! contract.
//!
//! # Naming
//!
//! `ct_fix__<op>__<carrier>__N<limbs>` — carrier `u64` is the primitive,
//! `fb32` is `FixedUInt<u32, N, Ct>`. `ct_fix__ASYM__*` fixtures keep
//! their public operands as in-body constants so LLVM can const-propagate
//! them — the configuration that historically triggers the XOR-select →
//! conditional-move rewrite that symmetric-taint fixtures miss.
//! `nct_fix__neg__*` are negative controls that MUST trip the gate.
//!
//! # Input values vs. taint
//!
//! The driver passes all-zero buffers; Valgrind's V-bits, not the values,
//! carry the "secret" marking. Where an algorithm needs an odd modulus to
//! avoid a panic path, the fixture ORs the low limb with 1 — a bitwise op
//! on partially-undefined data leaves that one bit defined, so the parity
//! check inside (`Odd::new_ct` / `try_new_odd_ct`) branches on a defined
//! bit and the rest of the modulus stays tainted.

#![cfg_attr(feature = "no-std", no_std)]
#![allow(non_snake_case)]
// The raw-pointer ABI is the harness contract; the only callers are the
// generated ct-ctgrind wrappers, which always pass valid buffers. Marking
// each fixture `unsafe fn` would only move the unsafe block into a
// wrapper with no additional reviewable invariant.
#![allow(clippy::not_unsafe_ptr_arg_deref)]

mod fixtures_asym;
mod fixtures_core;
mod fixtures_field;
mod fixtures_neg;

/// No-op Rust-visible anchor. ct-ctgrind calls this once so rustc links
/// the rlib at all — without a Rust-side reference the crate is dropped
/// from the link line and the `extern "C"` fixture symbols come back
/// undefined.
pub fn link_anchor() {}

// Only under the `panic-handler` feature (cross-built staticlib): std
// consumers (ct-ctgrind) supply their own and the two would collide.
#[cfg(feature = "panic-handler")]
#[panic_handler]
fn panic(_: &core::panic::PanicInfo) -> ! {
    loop {}
}

/// The load-bearing opacifier. See the module docs.
#[inline(always)]
pub(crate) fn bb<T>(v: T) -> T {
    core::hint::black_box(v)
}
