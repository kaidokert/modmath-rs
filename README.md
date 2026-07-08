### Modular math functions

[![Build](https://github.com/kaidokert/modmath-rs/actions/workflows/rust.yaml/badge.svg)](https://github.com/kaidokert/modmath-rs/actions/workflows/rust.yaml)
[![crate](https://img.shields.io/crates/v/modmath.svg)](https://crates.io/crates/modmath)
[![documentation](https://docs.rs/modmath/badge.svg)](https://docs.rs/modmath/)
[![Coverage Status](https://coveralls.io/repos/github/kaidokert/modmath-rs/badge.svg?branch=main)](https://coveralls.io/github/kaidokert/modmath-rs?branch=main)


Yet another mod math implementation, but written for _traits_. Functions
are constrained by `core::ops::` and [`const-num-traits`](https://crates.io/crates/const-num-traits)
traits, so any integer backend implementing them works — no concrete
bigint type is named anywhere.

Implements:
- Unsigned modular addition and subtraction
- Unsigned modular multiplication
- Unsigned modular exponentiation
- Unsigned modular inverse (extended Euclidean, and constant-time
  Bernstein–Yang for odd moduli)
- Unsigned modular Montgomery multiply (wide-REDC and CIOS)
- Unsigned modular Montgomery exponentiation

Two ways in:
- Free functions under the `basic` / `constrained` / `strict` modules —
  the same algorithms over three operand shapes (`Copy` by-value,
  `Clone` mixed, reference-based).
- `Field<T, P>` — precomputed Montgomery parameters with lifetime-branded
  residues. The `P` typestate selects the personality: `Nct` for
  variable-time verify paths, `Ct` for constant-time sign/decrypt paths,
  where operations finalize branchlessly and variable-time entry points
  are rejected at compile time.

The constant-time surface is checked in CI by the harnesses under
`ct-verify/`: Valgrind taint analysis on x86_64/aarch64, per-target
disassembly scanning on Cortex-M and RISC-V, and a panic-free linker
audit. See `ct-verify/README.md` for what each layer does and doesn't
cover.

`no_std`, no-alloc. The generic code trades speed for portability; if
you need maximum performance, use a dedicated bigint crate's own
modular arithmetic.
