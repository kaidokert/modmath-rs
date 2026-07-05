# modmath-cios

Row-op trait surface for CIOS Montgomery multiplication.

`CiosRowOps` is the minimal interface a multi-limb integer type exposes
so a CIOS (Coarsely Integrated Operand Scanning) Montgomery
multiplication loop can drive it: limb access plus the two fused row
kernels. The algorithm itself lives with the consumer — see
[`modmath`](https://crates.io/crates/modmath) — while bigint backends
implement the trait on their own type.

Single-word degenerate impls for `u8`/`u16`/`u32`/`u64` are provided,
treating each primitive as a 1-limb value.

`#![no_std]`, zero dependencies.

## License

Apache-2.0
