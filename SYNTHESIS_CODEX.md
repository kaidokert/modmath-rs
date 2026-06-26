# CIOS Trait Placement — Synthesis

This synthesizes the four `ANALYSIS_*` documents and the original
`cios_trait_placement.md`. The analyses agree on the mechanical problem but
split on where the new boundary should land.

The core issue is not one trait. `MulAccOps` currently mixes three layers:

| Layer | Content | Proper owner |
|---|---|---|
| L1: word arithmetic | primitive carrying/widening word ops | `const-num-traits` already owns this |
| L2: bigint row/access contract | word count/access plus CIOS row kernels | contested |
| L3: Montgomery algorithm/capability | CIOS double-loop and `CiosMontMul` traits | `modmath` |

The placement decision is about L2. L1 should stay in `const-num-traits`; L3
should stay in `modmath`.

## Consensus

All analyses converge on these points:

- `MulAccOps` is misnamed. The current methods are CIOS row operations, not a
  general multiply-accumulate abstraction.
- The algorithm body should not move into `fixed-bigint` unless the project is
  intentionally making `fixed-bigint` the only optimized backend.
- A slice-first `Limbs` API is attractive for non-CT heap bigints but cannot be
  the primary contract while the CT path matters. Plain indexed slices make the
  side-channel story too easy to violate.
- `GetWordOutput` is clunky. The current associated-type discriminant hides an
  important CT/NCT distinction.
- `Copy + Default` and static `word_count()` encode the fixed-size backend too
  strongly. If this surface is being changed, it is worth checking whether those
  constraints can be moved off the core trait.

## Main Disagreement

The real disagreement is whether limb access / row traits are low-level enough
for `const-num-traits`.

The “put it in `const-num-traits`” view argues:

- both crates already depend on it;
- limb access is structural, not modular arithmetic;
- no new crate is needed;
- CT/NCT personalities already show up in the ecosystem.

The opposing view argues:

- `const-num-traits` should stop at numeric primitives and word arithmetic;
- limb-array access is representation-level, not a num-traits-style operation;
- CIOS row kernels are algorithm contracts, not primitive numeric traits;
- putting this there invites later SOS/FIOS/NTT/polynomial-kernel drift.

I find the second argument stronger. `const-num-traits` is the right home for
word-level operations such as `CarryingMul`, `BorrowingSub`, `WideningMul`, and
similar primitive atoms. It is not the right home for “a bigint backend exposes
rows suitable for CIOS.” That is already one layer too high.

So: do not put `CiosRowOps` in `const-num-traits`. Be very cautious even with a
generic `WordAccess` there; it looks structural, but it is representation
exposure, not a numeric law or operation.

## Trait Shape

The minimum shape cleanup is Option G:

```rust
pub trait CiosRowOps: Sized {
    type Word;

    fn word_count(&self) -> usize;

    fn mul_acc_row(
        scalar: Self::Word,
        multiplicand: &Self,
        acc: &mut Self,
        carry_in: Self::Word,
    ) -> Self::Word;

    fn mul_acc_shift_row(
        scalar: Self::Word,
        multiplicand: &Self,
        acc: &mut Self,
        acc_hi: Self::Word,
    ) -> Self::Word;
}

pub trait WordAccess: CiosRowOps {
    fn word(&self, i: usize) -> Option<Self::Word>;
}

pub trait WordAccessCt: CiosRowOps {
    fn word_ct(&self, i: usize) -> subtle::CtOption<Self::Word>;
}
```

That is better than the current `GetWordOutput` associated type because the CT
contract appears in the bound instead of being hidden in an associated type.

However, before committing to separate `WordAccess` / `WordAccessCt`, run the
“discriminant probe” from `ANALYSIS_CLAUDE.md`: in CIOS, the word index appears
to be a public loop counter in `0..word_count`. If every lookup is guaranteed
in-bounds and the index is public, the fallible `Option`/`CtOption` may not be
load-bearing for the algorithm. In that case the best shape may be:

```rust
pub trait CiosRowOps: Sized {
    type Word;

    fn word_count(&self) -> usize;
    fn word(&self, i: usize) -> Self::Word; // precondition: i < word_count()

    fn mul_acc_row(...);
    fn mul_acc_shift_row(...);
}
```

That would collapse the two near-identical CIOS bodies into one while preserving
CT, provided the CT implementation of `word` does not branch on secret data. This
is a threat-model question, not just an API-style question. Verify it against the
actual CT backend before choosing.

## Placement Options That Remain Live

After rejecting `const-num-traits`, fixed-bigint ownership, and slice-first
limbs, there are three live placements.

### 1. Status Quo Plus Shape Cleanup

Keep the trait in `fixed-bigint`, but rename/split it. `modmath` still imports
one fixed-bigint trait family.

This is the smallest change. It does not meet the decoupling goal, but it
removes the worst naming and CT-contract opacity. It is a good first step if the
project is not ready to commit to an extraction.

### 2. Option H: Traits in `modmath`, Algorithm Stays in `modmath`

This is the missing option identified by Claude’s analysis.

`modmath` owns `CiosRowOps` plus the CIOS algorithm and blanket impls.
`fixed-bigint` optionally depends on `modmath` and implements the row traits for
its types. `modmath` can keep `fixed-bigint` as a dev-dependency for tests.

This gives:

- no `fixed_bigint::*` references in the `modmath` library;
- low backend implementor cost, because backends implement row ops rather than
  the whole algorithm;
- no new crate;
- a dependency cycle only in the optional/dev test configuration.

This dominates the original Option B, because Option B moved the algorithm into
`fixed-bigint` and made future backends reimplement too much.

### 3. Option C: New `mont-traits` / `bigint-algo-traits` Crate

A neutral sibling crate owns the L2 contracts. `modmath` owns the algorithm.
`fixed-bigint` implements the contracts.

This gives the cleanest graph and ages best if more shared algorithm contracts
are coming: SOS, FIOS, REDC variants, NTT kernels, polynomial multiplication, or
other bigint backend hooks.

The cost is real but small: one more crate, one more release cadence, and a bit
more coordination while the trait surface is still moving.

## Recommended Path

Use a two-stage path.

### Stage 1: Refactor Shape In Place

Do this now:

- Rename `MulAccOps` to `CiosRowOps`.
- Split `GetWordOutput` into explicit `WordAccess` / `WordAccessCt`, unless the
  discriminant probe proves a single infallible public-index `word()` is sound.
- Keep the CIOS algorithm body in `modmath`.
- Do not move anything into `const-num-traits`.
- While touching the trait, test whether `word_count(&self)` can replace static
  `word_count()`, and whether `Copy + Default` can move to specific impls or
  helper constructors.

This stage improves the API without forcing the crate-boundary decision.

### Stage 2: Extract When There Is a Second Consumer

If another backend or another shared algorithm family appears, extract L2 into a
new neutral crate.

Name candidates:

- `mont-traits` if the scope is strictly Montgomery contracts;
- `bigint-algo-traits` if the expected scope includes NTT, polynomial kernels,
  or other limb-row algorithm hooks.

At that point, use Option C. The extraction should be mechanical if Stage 1 kept
L2 cleanly separated from both fixed-bigint internals and modmath algorithm code.

## If Immediate Decoupling Is Required

If the immediate goal is “modmath library must not depend on fixed-bigint,” use
Option H rather than creating a tiny crate immediately:

- Move the row/access traits into `modmath`.
- Keep the algorithm and blanket impls in `modmath`.
- Add an optional `modmath` feature in `fixed-bigint` that implements those
  traits for `FixedUInt`.
- Accept the dev/feature test cycle as the temporary cost.

This is the best two-crate answer. It is not as clean as a neutral crate, but it
keeps backend implementor cost low and avoids premature crate creation.

## Final Decision Rule

- If this is a one-off CIOS cleanup: Stage 1 in place, then stop.
- If immediate decoupling is mandatory: Option H.
- If this is the beginning of a shared bigint-algorithm trait layer: Stage 1,
  then Option C with a new neutral crate.

Default recommendation: Stage 1 now. Do not use `const-num-traits`. Extract to a
new neutral crate only when the second backend or second algorithm justifies the
crate boundary.
