# Modmath CYCCNT fixtures

Runs the existing `ct-fixtures` ABI entry points on the J-Trace STM32F407VG
with concrete, contract-valid operands. Secret-modulus fixtures vary the odd
modulus and its matching Montgomery `n_prime`; public-modulus fixtures keep
the modulus fixed and vary only the documented secret operands.

Measurements use RTT, interrupt-free DWT `CYCCNT`, equal warm-up, and balanced
ABBA trial order. Positive fixtures have an absolute 32-cycle combined-spread
limit. The EEA inverse and schoolbook exponentiation ABI fixtures are negative
controls and must produce disjoint A/B ranges.

The carrier uses `krabi-caliper::PairedSuite` for DWT sampling, exact-path
warmups, policy evaluation, versioned reporting, diagnostics, and totals. It
emits lossless `EM_*` schema 1 records plus legacy `CT_*` records during host
tooling migration.

Run the complete hardware gate from this directory with the shared Rust
campaign runner:

```sh
cargo krabi-caliper run modmath-jtrace-f407
```

The declarative profile in `krabi-caliper.toml` owns the release build,
explicit SWD/J-Trace selection, RTT completion, deadline, ELF retention, raw
logs, and JSON/Markdown results below
`target/krabi-caliper/modmath-jtrace-f407/`. The retained metadata includes
the exact `cargo build` and `probe-rs run` commands for direct reproduction.

This timing-regression layer complements, rather than replaces, ctgrind and
the cross-target conditional-branch audit.

Initial STM32F407/J-Trace calibration passed all 13 positive fixtures and both
timing-negative controls (`passed:15 failed:0`). Positive combined spreads
were 0–10 cycles. EEA inverse and schoolbook exponentiation separated by 4,219
and 193,970 cycles respectively. The secret-indexed table lookup was same-cycle
on this target, as expected for an address-only leak; ctgrind remains the
authoritative detector for that class.
