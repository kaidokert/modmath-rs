# Constant-time hardening plan — CI verification for modmath

> **Status (2026-07):** phase 1 shipped as `ct-verify/` + the
> `ct-ctgrind.yml` workflow — the taint layer landed *first*, inverting
> the step order below, because it dodges both open design questions
> (loop policy, private-helper exposure: taint flows through public
> entries and Valgrind stack traces localize for free) and is the
> loop-tolerant authority anyway. Steps 2–4 (asm-grep matrix) are the
> next wave. One taint-model lesson from phase 1: mirror each entry's
> *documented* secrecy contract, not a blanket "everything is secret" —
> `basic_montgomery_mod_exp_pr_*_ct` documents its modulus as public
> (NCT precompute intentional), so only base/exponent are tainted
> there, while the `Field::try_new_odd_ct` path taints the modulus.

modmath's CT story today is enforced at exactly one layer: the
`Personality` typestate blocks variable-time entries from Ct carriers
at compile time. Nothing verifies that the code emitted for the
`_ct` surface is *actually* branch-free, on any target. fixed-bigint
has built that verification stack over the last months
(`fixed-bigint-rs/ct-verify/`, documented in
`fixed-bigint-rs/notes/CT_VERIFY.md`), and its downstream-consumer
section names modmath explicitly as the next crate to replicate the
pattern. This document is that replication plan, adapted for the one
structural way modmath differs from fixed-bigint.

Related: `CT_GET_WELL_PLAN.md` (repo root) tracks *implementation*
fixes to the CT surface itself; its "Verification plan" section ends
with "a ctgrind / dudect run … is the right next step for actual
production attestation." This plan is that step, industrialized.

## What fixed-bigint built (inventory)

Four layers, each catching failure modes the others can't see:

- **Layer 0 — typing.** `assert_not_impl_any!` compile-fail tests.
  modmath's equivalent already exists (the `NonCt` marker + `Rem`/`Div`
  bound gating); no new work here beyond keeping it tested.
- **Layer 1 — asm-grep (`ct-verify/ct-fixtures` + `ct-verify/ct-driver`).**
  A `staticlib` fixtures crate emits one `#[no_mangle] pub extern "C"`
  symbol per (operation, carrier) pair, `black_box` at both ends. A
  host driver cross-builds it per target, disassembles with
  `llvm-objdump`, and fails on forbidden conditional-control-transfer
  mnemonics (per-arch regex tables in `ct-driver/src/mnemonics.rs`),
  with allowlists for CT-friendly conditionals (Thumb IT blocks,
  `csel`, `cmov`). Has since grown call-graph reachability (helpers
  reachable from fixtures are scanned too) and a per-helper allowlist
  for public-parameter loops. Negative-control fixtures MUST trip the
  gate — the harness self-test.
- **Layer 2 — taint (`ct-verify/ct-ctgrind`).** Host binary linking
  the same fixtures, marking secret inputs as uninitialized via
  `crabgrind` client requests, run under Valgrind memcheck. Any
  branch or memory access on tainted data is an error. Caught a real
  bug (LLVM rewriting the XOR-select idiom into a secret-flag `cmov`)
  within minutes of being wired up.
- **Layer 3 — deferred pillars.** Hardware DWT cycle counting
  (seams in place), `cargo-checkct`/binsec symbolic execution,
  `dudect` statistical timing.

CI: two dedicated workflows. `ct-verify.yml` runs the asm-grep
matrix (thumbv7em / thumbv7m / thumbv6m / riscv32imc / riscv32imac /
aarch64 / x86_64+`lzcnt,bmi1` / AVR-nightly), pinned to rustc 1.86.
`ct-ctgrind.yml` runs Valgrind taint on x86_64 **and** aarch64
runners. Both hard-fail, both separate workflow files so the check
names are distinct in the PR UI. Codegen is pinned by a workspace
`[profile.release]` (`lto = "fat"`, `codegen-units = 1`,
`opt-level = "z"`, `panic = "abort"`) so a passing PR locks the
inspection and toolchain bumps become reviewable perturbations.

## What modmath must verify (fixture inventory)

The `_ct` surface as of v0.4.1:

**Free functions**
- `montgomery::cios::cios_montgomery_mul_ct` (the wide-REDC primitive,
  via `CiosMontMulCt`)
- `basic_mont`: `compute_r_mod_n_ct`, `compute_r2_mod_n_ct`,
  `wide_redc_ct` / `strict_wide_redc_ct`,
  `wide_montgomery_mul_ct` / `strict_wide_montgomery_mul_ct`,
  `wide_montgomery_mul_acc_ct` / `strict_wide_montgomery_mul_acc_ct`,
  `basic_montgomery_mod_exp_pr_odd_ct` / `basic_montgomery_mod_exp_pr_ct`
- `inv::safegcd::safegcd_inv_ct` — plus its private leaf helpers
  (`se_select`, `se_add`, `se_neg`, `se_shr1`, `add_mod_ct`,
  `half_mod_ct`), which are exactly the mask-select idioms LLVM has
  a record of rewriting into secret-flag conditional moves

**`Field<T, Ct>` methods**
- `try_new_odd_ct` / `new_odd_ct` (precompute over a secret modulus —
  the CRT case makes the modulus itself secret)
- `reduce`, `into_raw`, `add`, `sub`, `mul`, `exp`, `mul_acc`,
  `wide_redc`, `inv_fermat`, `inv_safegcd_ct`
- `Residue::cswap`, `Residue::ct_eq`
- `exp_public_exp` is *deliberately not* CT in the exponent — exclude
  it, and note it in the fixtures crate so nobody "fixes" the omission

**Negative controls (must trip every layer)**
- `basic_mod_inv` (EEA — magnitude-dependent loop count)
- a non-`_pr` schoolbook entry (`Rem`-based reduction)
- one hand-written branchy compare, same as fixed-bigint's

**Carrier diagonals.** modmath is generic; fixtures must pick
concrete carriers. Two tiers:
- Primitive carriers `u32` / `u64` — every fixture.
- Multi-limb `FixedUInt<u32, 4, Ct>` and `FixedUInt<u32, 8, Ct>` —
  every fixture. This is the tier that actually resembles deployment,
  and it exercises the post-LTO inlined fixed-bigint primitives at
  modmath's real call sites (the stated point of the downstream
  consumer pattern).
- RSA-shaped `FixedUInt<u32, 64, Ct>` (2048-bit) — ctgrind layer
  only, and only for `wide_montgomery_mul_ct` + `safegcd_inv_ct`
  (budget check: safegcd at 2048 bits is ~10k divsteps; measure
  Valgrind wall-time before committing it to CI).

The fixtures crate is a `publish = false` workspace member, so it can
depend on `fixed-bigint` directly without violating modmath's
"fixed-bigint is dev-dep only" rule.

## The structural difference: public-bound loops

This is where a verbatim copy of fixed-bigint's asm-grep breaks, and
it needs to be decided before writing fixtures.

fixed-bigint's Ct primitives are mostly straight-line per-limb code;
its limb loops live in small helpers that the driver exempts by
symbol name (`allowed_helpers` — acceptable there because an exempted
helper like `conditional_select` contains *nothing but* the counter
loop). modmath's CT primitives **are** loops: CIOS iterates over
`word_count()`, `exp` over the exponent bit-width, `safegcd_inv_ct`
over `divsteps_total` (~2.5 × carrier bits — hundreds to thousands of
iterations, never unrolled). Under `opt-level = "z"` these compile to
`cmp; b.ne` loop-back branches *inside the fixture symbol*, which the
raw mnemonic gate flags. Exempting whole primitives by symbol name
would gut the gate — the loop body is precisely what we want scanned.

All of these loops have **public** trip counts (compile-time width,
public exponent bit-length, `divsteps_total` from carrier bits), so
they are CT-safe; the gate just can't see that statically. Handling,
in order of authority:

1. **ctgrind is the loop-tolerant authority.** Taint analysis handles
   public-bound loops correctly by construction: the loop counter is
   untainted, so its branch never errors; any *data*-dependent branch
   still does. For modmath, Layer 2 is therefore not a nice-to-have
   on top of asm-grep — it is the primary full-function gate, and it
   should land in the same phase as asm-grep, not deferred the way
   fixed-bigint sequenced it.
2. **Asm-grep stays strict where it can be.** Loop-free fixtures —
   the safegcd leaf helpers, `cswap`, `ct_eq`, single-row CIOS ops,
   and most fixtures at primitive carriers where `word_count() == 1`
   collapses the loops — get the zero-forbidden-mnemonic gate,
   unchanged from fixed-bigint.
3. **Loopy fixtures get branch-site baselining.** For fixtures whose
   bodies legitimately contain counter loops, the driver records the
   set of conditional-branch sites (symbol + offset + mnemonic) in a
   committed per-target baseline; CI fails on any *new* site. Pinned
   toolchain + pinned profile make codegen deterministic, which is
   what makes this stable — the same property fixed-bigint pins for.
   A regression that adds a data-dependent branch can't hide, and a
   toolchain bump that moves offsets becomes an explicit baseline-
   regeneration PR. This is a driver extension (~1 new module +
   a `--bless` flag), the only genuinely new code in the plan.

The calibration step in Phase 2 (build the fixtures, look at what
actually contains loop-backs at each carrier) decides which fixture
goes in which bucket — don't guess it up front.

## Stepwise plan

Each step is one PR, independently landable, ordered so the harness
is self-testing before it gates anything.

### Step 0 — decisions (no code)

Two decisions gate the fixture crate:

1. **Exposure of private CT leaf helpers.** The safegcd internals
   (`se_add`, `se_neg`, `se_shr1`, `se_select`, `add_mod_ct`,
   `half_mod_ct`) are the highest-value asm-grep targets but are
   private. Options: (a) `#[doc(hidden)] pub mod ct_internals`
   re-exports gated behind a non-default `ct-verify` cargo feature —
   recommended; costless when the feature is off, and honest that
   these are inspection seams, not API; (b) rely on their inlining
   into `safegcd_inv_ct`'s fixture and inspect only the whole — loses
   the ability to localize a violation to a helper.
2. **Loop policy** — confirm the branch-site-baseline design above
   (or pick straight allowlisting and accept the coverage hole; not
   recommended).

### Step 1 — workspace scaffolding

- Add `ct-verify/ct-fixtures`, `ct-verify/ct-driver`,
  `ct-verify/ct-ctgrind` to `[workspace] members` (all
  `publish = false`, `version = 0.0.0`).
- Add the pinned `[profile.release]` block at the workspace root,
  copied from fixed-bigint verbatim. Profiles don't publish with the
  crate, but they do change local `--release` builds and benches —
  flag that in the PR description.
- Pin rustc 1.86 in the future workflows (matches both repos' MSRV).

### Step 2 — ct-fixtures crate + local asm-grep

- Copy the `ct_fix_*` macro layer from
  `fixed-bigint-rs/ct-verify/ct-fixtures/src/lib.rs` (the `black_box`
  discipline, the `panic-handler` feature dance for staticlib vs
  rlib, `link_anchor`). The discipline is non-negotiable: every
  fixture goes through the macros; a fixture that bypasses them is a
  false-pass and gets rejected in review.
- Copy `ct-driver` wholesale — `mnemonics.rs` and `target.rs` are
  architecture facts and shared verbatim per fixed-bigint's own
  downstream guidance; `main.rs`/`parse.rs`/`report.rs` need only
  crate-name and fixture-prefix adjustments. The copy inherits the
  reachability closure and helper allowlist for free. (fixed-bigint's
  docs say to factor a shared driver crate once 2–3 consumers have
  copied it; modmath is consumer #2 — note the duplication in the
  PR, extract later, don't block on it.)
- Write the *loop-free tier* of fixtures only: safegcd leaf helpers
  (per Step 0 exposure), `cswap`, `ct_eq`, primitive-carrier
  diagonals, plus all three negative controls.
- Run locally on host + `thumbv7m-none-eabi`; iterate until positive
  fixtures pass and negative controls all trip. **Calibrate:** build
  the remaining (loopy) fixtures unregistered and record which
  bodies contain loop-backs per target — this fixes the bucket
  assignment for Step 3.

### Step 3 — loopy fixtures + branch-site baselining

- Implement the baseline mechanism in ct-driver (`--bless` to
  regenerate, committed JSON per target under
  `ct-verify/baselines/`).
- Add the multi-limb `FixedUInt` diagonals and the loopy fixtures:
  CIOS mul, wide-REDC family, `exp`, `safegcd_inv_ct`,
  `try_new_odd_ct`.
- Add one negative control *in the loopy tier*: a variable-trip-count
  loop (`while g != 0`-style) whose exit branch is data-dependent —
  the baseline mechanism must catch it as a new site, proving the
  loopy tier still self-tests.

### Step 4 — ct-verify.yml CI matrix

Copy fixed-bigint's workflow with its accumulated fixes baked in:
- Same target rows and priorities (thumbv7em, thumbv7m, thumbv6m,
  riscv32imc, riscv32imac, aarch64, x86_64, AVR-on-nightly).
- x86_64 row keeps `-C target-feature=+lzcnt,+bmi1` — baseline x86_64
  lowers `leading_zeros()` to `test/je/bsr`, a real data-dependent
  branch inside the std intrinsic that no `black_box` in our code can
  remove. Don't re-learn this.
- Per-target RUSTFLAGS via `CARGO_TARGET_<TRIPLE>_RUSTFLAGS` so flags
  don't leak into the host driver build; AVR needs bare `nightly` +
  `rust-src` + `-C target-cpu=atmega328p`.
- Hard-fail, JSON report uploaded as artifact, separate workflow file
  so the check name is distinct.

### Step 5 — ct-ctgrind (taint layer)

- Copy the ct-ctgrind harness (fixture registry via `inventory`,
  `crabgrind` behind a Linux cfg with stubs elsewhere, error-counter
  delta protocol with `--error-exitcode=0`).
- **Secret designation:** taint operand values, exponents, *and the
  modulus* (CRT gives secret moduli; `try_new_odd_ct` documents
  CT-over-secret-modulus). Loop bounds and carrier width stay
  untainted — bit-width of `T` is public by construction.
- **Asymmetric-taint fixtures from day one:** tainted secret ×
  untainted constant operand pairs. This is the configuration that
  exposed fixed-bigint's XOR-select→`cmov` rewrite; symmetric taint
  misses it. Target every mask-select site in safegcd and the CIOS
  conditional-subtract.
- **Both x86_64 and aarch64 runner rows from day one**
  (`ubuntu-latest` + `ubuntu-24.04-arm`). fixed-bigint started
  x86_64-only and got a compiler-heuristic false-pass on aarch64 for
  its trouble.
- Include the 2048-bit RSA-shaped fixtures here if the Valgrind
  wall-time budget allows (measure first).

### Step 6 — panic-free audit

Small companion member (fixed-bigint's `panic-free-audit` pattern):
`#[no_mangle]` wrappers over the CT surface at a deployment-shaped
instantiation, cross-built with the pinned release profile, then
`cargo nm` asserts no `panic_*`/`unwind` symbols trace back to them.
For CT code a reachable panic is both a DoS edge and a timing oracle.
`CT_GET_WELL_PLAN.md` already wanted exactly this for
`try_new_odd_ct`. Cheap; can fold into Step 2's PR if it stays small.

### Step 7 — documentation

`notes/CT_VERIFY.md` in this repo mirroring fixed-bigint's: the
layered model, how to run locally, the loop-policy rationale (the one
genuinely novel part — write it down before it's forgotten), known
target subtleties, and the negative-control contract. Update
`CLAUDE.md`'s command list with the local driver invocations.

## Future pillars (document now, build later)

Same list as fixed-bigint's, same rationale, one reordering note:

- **`cargo-checkct` (binsec symbolic execution)** — the only tool
  that sees the select-rewrite family on Thumb, where ctgrind can't
  run and IT blocks are allowlisted by asm-grep. Thumb is the actual
  deployment target for modmath's downstream consumers, so this is
  the highest-value deferred pillar. Start `continue-on-error: true`.
- **`dudect-bencher`** — host statistical timing; useful smoke when
  `subtle` or the toolchain bumps.
- **Hardware DWT cycle counting** — put the `ct_seam_begin/end` weak
  symbols into the fixtures crate now (they're free); the self-hosted
  Cortex-M runner waits for a concrete consumer demanding it.

## Honest limits

Inherited from fixed-bigint's write-up and still true here: even all
layers together don't cover cache-timing on shared caches, secret-
dependent prefetch, or speculative side effects. Thumb IT-block
predication passes asm-grep by design and has no taint coverage until
checkct lands. The goal is orthogonal coverage that must all be
defeated at once, plus consciously audited gaps — not a proof.

Two modmath-specific caveats to keep in view:

- The gate verifies *fixture instantiations*, not the generic code.
  A carrier whose own primitive ops branch (a third-party bigint that
  never made CT claims) invalidates everything above it; the fixtures
  pin fixed-bigint and primitives, which are the carriers with CT
  contracts. The typestate layer is what keeps other carriers out of
  the `_ct` surface at compile time.
- `divsteps_total` and `word_count()` being public is an *argument*,
  not something any layer verifies. A future change that made a loop
  bound depend on a secret value would pass ctgrind's untainted-
  counter assumption if the bound is computed outside the tainted
  region. The branch-site baseline catches the codegen shape change;
  review has to catch the semantic one.
