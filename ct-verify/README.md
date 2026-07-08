# ct-verify

Runtime taint verification for modmath's constant-time surface,
instantiating the pattern from fixed-bigint's `ct-verify/` harness.
This directory is the runtime-taint layer; a cross-target asm-grep
layer comes later.

Two members:

- **`ct-fixtures`** — one `#[no_mangle] pub extern "C"` symbol per
  (CT entry point, carrier) pair, `core::hint::black_box` at both ends
  so the optimizer can't fold the body away. Carriers: `u64` and
  `FixedUInt<u32, 8, Ct>` (256-bit). `ct_fix__ASYM__*` fixtures keep
  public operands as in-body constants — the configuration where LLVM
  historically rewrites the XOR-select idiom into a secret-flag
  conditional move. `nct_fix__neg__*` fixtures are negative controls
  (EEA inverse, schoolbook exp, secret-indexed table load) that must
  trip the gate, or the harness itself is broken.
- **`ct-ctgrind`** — runs each fixture under Valgrind memcheck with its
  secret inputs tagged `MAKE_MEM_UNDEFINED` via crabgrind. Valgrind
  flags any conditional jump or memory access that depends on the
  tainted bytes — including secret-derived loads that no branch-level
  inspection can see. The driver reads Valgrind's error counter between
  fixtures to attribute violations, and decides pass/fail itself
  (positives must not trip, negatives must).

Which inputs are tainted mirrors each entry's documented secrecy
contract, not a blanket "everything is secret": the wide-mul/REDC and
`Field::try_new_odd_ct` paths taint the modulus (RSA-CRT: p and q are
secret), while `basic_montgomery_mod_exp_pr_*_ct` documents its
modulus as public (NCT precompute is intentional there) and gets only
base and exponent tainted. When adding a fixture, read the entry's
rustdoc first and encode exactly the secrecy it promises.

## Running locally

On a Linux host:

```
cargo build --release -p ct-ctgrind
valgrind --tool=memcheck --error-exitcode=0 -q target/release/ct-ctgrind
```

On macOS (no Valgrind), use the Docker image — see
`ct-ctgrind/Dockerfile` for build/run commands. On Apple Silicon the
native arm64 image matches the CI aarch64 row; the amd64 image works
at baseline but SIGILLs inside Valgrind's translation when the binary
is built with `-C target-feature=+lzcnt,+bmi1` (an emulation artifact
— the CI x86_64 row runs those flags on real hardware).

CI: `.github/workflows/ct-ctgrind.yml`, x86_64 + aarch64 rows, both
hard-fail.

## Extending

Add the symbol in `ct-fixtures` (black_box discipline is
non-negotiable — see the crate docs) and the matching taint wrapper in
`ct-ctgrind/src/fixtures.rs`; the two mirror each other one-for-one.
All inputs cross the ABI by pointer so taint marks can't be
const-propagated away; keep it that way.
