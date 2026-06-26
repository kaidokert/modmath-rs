# HOW TO FIX: tracking `const-num-traits` from `modmath`

Two independent issues bite when chasing the upstream fork. Both surface as
build failures; neither needs source migration inside `modmath`.

## TL;DR

| Issue | Where to fix | Effort |
|---|---|---|
| 1. Cargo dep `package = "num-traits"` rename is stale | `modmath/Cargo.toml` | one-liner |
| 2. Nightly: `OverflowingAdd` / `OverflowingSub` missing at `const_num_traits` crate root | `../num-traits/src/lib.rs` (upstream) | one-liner upstream |

`modmath/src/**` source code does **not** need the
`impl c0nst Trait` → `c0nst impl Trait` migration that bit `fixed-bigint`
on nightly ≥ 2026-06-19. Modmath only uses the `c0nst` macro at
**function** and **bound** position (`pub c0nst fn …`, `T: [c0nst] Trait`),
neither of which moved in
[rust-lang/rust#158009](https://github.com/rust-lang/rust/pull/158009).
The keyword shift was strictly on `impl`-block heads.

## Issue 1: Cargo dep rename is stale

### Symptom

```text
$ cargo +nightly build --features nightly
error: no matching package named `num-traits` found
location searched: /opt/m/rust/bigint/num-traits
required by package `modmath v0.3.1 (/opt/m/rust/bigint/modmath-rs/modmath)`
```

### Cause

The path-dep dep stanza renamed the package via `package = "num-traits"`:

```toml
const-num-traits = { package = "num-traits", path = "/opt/m/rust/bigint/num-traits", default-features = false }
```

This was correct when the fork's `Cargo.toml` published itself as
`name = "num-traits"`. The fork has since renamed its own published
name to `const-num-traits` (matches the dep alias), so the rename is
not just unnecessary, it's actively wrong — Cargo looks for a package
literally named `num-traits` at the path and finds `const-num-traits`
instead.

### Fix

Drop the `package =` field:

```toml
const-num-traits = { path = "/opt/m/rust/bigint/num-traits", default-features = false }
```

Verify:

```bash
cargo build                                 # stable — unchanged
cargo build --features nightly              # plain check; nightly issue may still bite (see #2)
```

## Issue 2: `OverflowingAdd`/`OverflowingSub` missing at the crate root

### Symptom

```text
$ cargo +nightly build --features nightly
error[E0432]: unresolved imports `const_num_traits::OverflowingAdd`,
                                  `const_num_traits::OverflowingSub`
 --> modmath/src/add.rs:2:24
  | use const_num_traits::{OverflowingAdd, OverflowingSub};
  |                        ^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^ no `OverflowingAdd` in the root
```

…repeated for `add.rs`, `exp.rs`, `mul.rs`, `sub.rs` (8 errors, 2 per file).

### Cause

This is an **upstream gap** in `const-num-traits`, not a `modmath` bug.
`crate::ops::overflowing` defines `OverflowingAdd` / `OverflowingSub`
correctly, but `src/lib.rs` doesn't re-export them at the crate root
(unlike, say, `CarryingAdd`, `CheckedAdd`, `WrappingAdd`, which all
*are* re-exported). The mixed-signed flavors (`OverflowingAddSigned`,
`OverflowingAddUnsigned`, `OverflowingSubSigned`,
`OverflowingSubUnsigned`) from `ops::mixed` *are* re-exported — only
the plain `OverflowingAdd`/`Sub` from `ops::overflowing` are missing.

Modmath compiles fine on stable because the four affected files gate
their `c0nst!` block (where these imports live) on
`#[cfg(feature = "nightly")]`. The plain (non-`c0nst`) bodies in the
same files use `num_traits::OverflowingAdd` from stock `num-traits 0.2`,
which is unaffected.

### Fix (upstream, in `../num-traits/src/lib.rs`)

Add the re-export next to the other `ops::` re-exports — somewhere
near the `crate::ops::wrapping` line:

```rust
pub use crate::ops::overflowing::{
    OverflowingAdd, OverflowingDiv, OverflowingMul, OverflowingNeg,
    OverflowingRem, OverflowingShl, OverflowingShr, OverflowingSub,
};
```

(Adjust the trait list to what `ops/overflowing.rs` actually exports.)

Once that lands in the fork, modmath's nightly build resolves with no
modmath-side edits — the existing `use const_num_traits::{OverflowingAdd,
OverflowingSub};` in `modmath/src/{add,exp,mul,sub}.rs` finds the
re-exports at the root.

### Workaround (modmath-side, if upstream can't be touched)

Import from the module path directly:

```rust
use const_num_traits::ops::overflowing::{OverflowingAdd, OverflowingSub};
```

This works today against the unmodified fork. It's a workaround, not a
fix — once the upstream re-export lands, prefer the root-level import
for consistency with the rest of modmath.

## Why no `impl c0nst Trait` migration is needed in modmath

For context: `fixed-bigint` had to mass-migrate `impl<...> c0nst Trait`
→ `c0nst impl<...> Trait` across ~200 sites because of
[rust-lang/rust#158009](https://github.com/rust-lang/rust/pull/158009)
("Reject `impl const Trait` since the right syntax is `const impl Trait`
now", merged 2026-06-18, in `nightly-2026-06-19+`). The doc for that
migration lives at
[`/opt/m/rust/bigint/num-traits/HOW_TO_FIX.md`](../num-traits/HOW_TO_FIX.md).

Modmath isn't affected because it never writes `impl c0nst Trait`
heads. All four `c0nst::c0nst!` blocks in `modmath/src/` are
*function* declarations:

```rust
c0nst::c0nst! {
    pub c0nst fn const_mod_add<T>(a: T, b: T, m: T) -> T
    where
        T: [c0nst] core::cmp::PartialOrd
            + [c0nst] OverflowingAdd
            + [c0nst] OverflowingSub
            + …
```

`pub c0nst fn` and `[c0nst]` bound syntax are **untouched** by #158009 —
only the `impl` keyword position moved. If modmath ever grows `impl
c0nst Trait for FooT` blocks (e.g. for a custom typestate), it'll need
the same `impl c0nst` → `c0nst impl` mechanical pass; until then, the
syntax migration is a no-op here.

## Verification

After applying both fixes:

```bash
cargo build                                  # stable, default features
cargo build --features nightly               # nightly, c0nst surface
cargo test --features nightly                # full test suite

# expected counts on a clean tree (modulo whatever local edits):
#   stable:  ~459 tests
#   nightly: at least the same; nightly enables more c0nst-gated tests
```

## Reference fixes

- `fixed-bigint` migration commit (`c5d1221`,
  `experiment/external-const-num-traits` branch): full
  `impl c0nst Trait` → `c0nst impl Trait` pass, 200+ sites,
  for comparison. The `[A-Z]`-anchored regex misses path-prefixed
  trait names (`core::cmp::Eq`, `crate::*`) — see the gap-fix commit
  message for details.
- `const-num-traits` migration commit (upstream, `800be14`,
  `squashy` branch): "Migrate impl c0nst → c0nst impl
  (rust-lang/rust#158009); drop nightly pin".
- The procedure: [`../num-traits/HOW_TO_FIX.md`](../num-traits/HOW_TO_FIX.md).
