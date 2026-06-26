# CIOS Trait Placement — Codex Analysis

`MulAccOps` should not live in `const-num-traits`. It is not a numeric primitive; it is an algorithm-specific limb-row contract for CIOS Montgomery multiplication. Putting it in `const-num-traits` would repeat the typestate “zoo” problem in trait form: useful, but at the wrong layer.

My recommendation is:

1. Rename/split the trait shape now.
2. Keep the algorithm body reusable from `modmath`.
3. Choose either a new tiny traits crate or a temporary status quo, depending on how much future shared algorithm work is expected.

## Best Long-Term Shape

Use Option C + G: a small sibling crate, but with a clearer trait split.

A shape like this is the right direction:

```rust
pub trait CiosRowOps: Sized + Copy + Default {
    type Word: Copy + PartialOrd;

    fn word_count() -> usize;

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

Then `modmath` owns the generic CIOS algorithm over:

- `T: CiosRowOps + WordAccess` for NCT
- `T: CiosRowOps + WordAccessCt` for CT

`fixed-bigint` implements those traits. `modmath` no longer imports `fixed_bigint`, and `fixed-bigint` does not need to depend on `modmath`.

This is the cleanest architecture if more shared backend-facing contracts are expected: SOS/FIOS row ops, Montgomery reduction variants, NTT limb kernels, or other bigint arithmetic hooks.

## What Not To Do

Do not move the algorithm into `fixed-bigint` unless the explicit decision is that `fixed-bigint` is the only backend that gets the optimized CIOS path. That makes `modmath` cleaner on paper, but it forces every future backend to reimplement the whole algorithm instead of just the row primitives.

Do not use a slice-based `Limbs` trait as the primary contract. It is attractive for non-CT heap bigints, but it breaks the current CT design by exposing indexed limb access as a normal slice. A separate NCT-only fast path can exist later, but it should not replace the CT-safe contract.

Do not put `CiosMontMul` in `fixed-bigint`. That makes the modular math trait identity live in the bigint crate and makes other bigint backends second-class.

## Pragmatic Alternative

If a new crate is premature, keep the status quo but do the shape cleanup in `fixed-bigint`:

- Rename `MulAccOps` to `CiosRowOps`.
- Split `get_word` into `WordAccess` and `WordAccessCt`.
- Leave `modmath` temporarily depending on `fixed_bigint::{CiosRowOps, WordAccess, WordAccessCt}`.

That gives the main conceptual cleanup without introducing a crate. Later, moving those three traits into `mont-traits` is mostly mechanical.

This is better than moving the algorithm into `fixed-bigint`, because it preserves the reusable algorithm body in `modmath`. Moving the algorithm into `fixed-bigint` decouples crates but raises backend implementor cost too much.

## Recommended Two-Stage Path

Stage 1, low risk:

- Keep placement as-is.
- Rename `MulAccOps` to `CiosRowOps`.
- Split `GetWordOutput` into `WordAccess` / `WordAccessCt`.
- Keep the CIOS algorithm in `modmath`.

Stage 2, when a second backend or second Montgomery variant appears:

- Extract `CiosRowOps`, `WordAccess`, `WordAccessCt`, and possibly `CiosMontMul` / `CiosMontMulCt` convenience traits into `mont-traits`.
- Keep the algorithm body in `modmath`.

This avoids creating a tiny crate prematurely, but keeps the design pointed at the right final state.

## Final Call

Do Stage 1 now. Plan for Option C + G once there is either another backend or another shared algorithm family.
