# CIOS trait placement — where does `MulAccOps` belong?

A design-exploration document. The question: modmath currently depends on
`fixed_bigint::MulAccOps` in production code (exactly one site, in
`modmath/src/montgomery/cios.rs`). Everything else can be migrated to
`const_num_traits::*` directly. This last coupling is the structural
one. Where should the trait live, and what should its shape be?

This document lays out the current shape, why it looks the way it does,
and every plausible way of moving or splitting the traits between
crates. No recommendation — just the full exploration.

---

## 1. What CIOS needs from the limb representation

Coarsely Integrated Operand Scanning Montgomery multiplication
interleaves the multiply and the REDC reduction into a single
double-loop. Per outer iteration, it needs to:

1. Read the i-th word of operand `a` (scalar).
2. Compute `acc += a[i] * b` as a fused multiply-accumulate row over
   the limb array of `b`, returning a carry-out word.
3. Compute a reduction factor `m = acc[0] * n_prime_0 (mod word)`.
4. Compute `acc = (acc + m * modulus) >> word_bits` as a second fused
   row (the "shift-row").
5. After `n` iterations, conditionally subtract the modulus if the
   accumulator overflowed.

Steps 1, 2, and 4 require word-level access to the limb array of the
bigint. The algorithm is generic over the bigint type but specific to
"types whose limbs can be addressed individually."

## 2. The current trait

Lives in `fixed-bigint` at `src/mul_acc_ops.rs`:

```rust
pub trait MulAccOps: Sized + Copy + Default {
    type Word: Copy + PartialOrd;
    type GetWordOutput;

    fn word_count() -> usize;
    fn get_word(&self, i: usize) -> Self::GetWordOutput;
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
```

modmath consumes it in `cios.rs` through two algorithm bodies
(`cios_montgomery_mul` for Nct, `cios_montgomery_mul_ct` for Ct) and
two convenience traits (`CiosMontMul`, `CiosMontMulCt`) whose blanket
impls are gated on `T: MulAccOps + …`.

## 3. Three quirks of the current shape

### 3a. The `GetWordOutput` discriminant

`get_word(&self, i: usize)` returns an associated type, not a concrete
type. The two FixedUInt personalities resolve it differently:

- `FixedUInt<W, N, Nct>`: `GetWordOutput = Option<W>`. O(1) array
  index; `None` on out-of-bounds.
- `FixedUInt<W, N, Ct>`: `GetWordOutput = subtle::CtOption<W>`. O(N)
  branchless scan with constant-time selection; the discriminant is
  exposed only through `subtle`-backed combinators.

**Why it's this way:** for the CT path, exposing `Option<W>` would be a
side-channel leak — pattern matching on `Some`/`None` branches on the
position. `CtOption` exists precisely to prevent that. So the trait
*cannot* uniformly return `Option<W>`; the CT contract requires the
return type to vary.

**Why it's clunky:** modmath ends up with two algorithm bodies that
differ in one character — `?` vs `.into_option()?`. The bodies cannot
be unified because the trait method returns two different types. The
"clever" part is that the CT contract is encoded at the type level,
not in naming convention; the clunk is the cost of that correctness.

### 3b. CIOS-specific terminology in a generically-named trait

The trait is called `MulAccOps` but its actual methods —
`mul_acc_row` (Phase 1 of CIOS) and `mul_acc_shift_row` (Phase 2 of
CIOS, the right-shift-by-one-word fold) — are CIOS-specific. The
"shift" in `mul_acc_shift_row` is the CIOS post-multiply word-shift,
not a generic operation.

So the trait is functionally `CiosRowOps`, not a generic
multiply-accumulate primitive. A consumer wanting a different
Montgomery variant (SOS, FIOS) cannot use this trait directly — they
would need different row primitives.

This is a layering concern, not a correctness one. The trait works
fine for CIOS; it just doesn't generalize.

### 3c. Fixed-limb-count assumption

`word_count()` is a static method (no `&self`). This bakes in the
assumption that the bigint's word count is known at the type level.
Works for `FixedUInt<W, N, _>` (N is a const generic); does **not**
work for heap-allocated bigints (`num_bigint::BigUint`,
`crypto_bigint::BoxedUint`, etc.) where the word count is a runtime
property.

The `Self: Default + Copy` bound reinforces the same assumption —
`Default` means "zero of statically-known size," `Copy` rules out
heap-allocated storage entirely.

So this trait fundamentally excludes a large class of bigint backends
from being CIOS providers. That's a design choice (fast path for
fixed-size); it's not necessarily wrong, but it should be acknowledged.

## 4. The placement dilemma

modmath defines `CiosMontMul` and `CiosMontMulCt` and provides blanket
impls that delegate to `MulAccOps`. To remove modmath's
`fixed_bigint::*` dependency, we need to break the relationship
between modmath's `CiosMontMul` and fixed-bigint's `MulAccOps`. The
algorithm body (currently in modmath's `cios.rs`) and the trait
definitions can move independently. The placement options for **the
trait definitions** are:

- **modmath**: where `CiosMontMul`/`CiosMontMulCt` already live.
- **fixed-bigint**: where `MulAccOps` already lives.
- **const-num-traits**: the const-friendly numeric primitives crate.
- **A new dedicated crate**: e.g., `mont-traits`,
  `bigint-mont-traits`, `cios-traits`.

The placement options for **the algorithm body** are:

- **modmath**: where it is now (`cios.rs` has both algorithm bodies).
- **fixed-bigint**: alongside `MulAccOps`.
- **A new dedicated crate**: alongside the trait definitions.

These two questions are independent, which is why the option space is
larger than it looks.

## 5. Evaluation dimensions

| Dimension | What to look for |
|---|---|
| **modmath ↔ fixed-bigint coupling** | Does modmath's lib still `use fixed_bigint::*`? |
| **Cargo-level risk** | Are there dep cycles? Are they in the "dev-dep cycle is documented and supported" category, or more fragile? |
| **Number of crates to maintain** | One? Two? Three? Each is a separate version/release cadence. |
| **Identity / scope clarity** | Does each crate have a coherent stated purpose, or are they grab-bags? |
| **Future extensibility** | If you add SOS, REDC variants, NTT, polynomial mul later, does this design accommodate them? |
| **CT contract preservation** | Does the `GetWordOutput` discriminant survive? Or does the refactor force exposing limbs in a way that breaks CT? |
| **Backend implementor cost** | What does a new bigint backend have to write to plug in? One trait? Several? The whole algorithm? |
| **Static vs. dynamic limb count** | Does the design exclude heap-allocated bigints, or accommodate them? |

## 6. The options

### Option A — Status quo

Do nothing. modmath `use`s `fixed_bigint::MulAccOps` in its production
code. Everything works; the coupling is documented as the price of the
CIOS abstraction.

- **modmath ↔ fixed-bigint:** coupled (one `use` site in `cios.rs`).
- **Cargo risk:** none.
- **Crates maintained:** existing two.
- **Identity:** modmath is a modular-math library; fixed-bigint owns
  the limb-level row primitives; the coupling reflects the actual
  layering.
- **Future extensibility:** any new bigint backend that wants modmath's
  CIOS must impl fixed-bigint's `MulAccOps`. That's awkward — they'd
  be implementing a trait from fixed-bigint without using fixed-bigint.
- **CT contract:** preserved.
- **Cost:** zero.

This is the honest "do nothing" option. The only argument against it
is the original "make modmath really generic just over traits" goal.

### Option B — Trait stays in modmath, blanket impl removed

modmath drops the `impl<T: MulAccOps + ...> CiosMontMul for T` blanket.
The algorithm body in `cios.rs` is deleted from modmath and moved to
fixed-bigint. fixed-bigint gains an optional `modmath` feature; under
that feature, fixed-bigint provides `impl modmath::CiosMontMul for
FixedUInt<W, N, Nct>` and the Ct equivalent, writing the algorithm
internally against its own `MulAccOps`.

- **modmath ↔ fixed-bigint:** modmath's lib has **zero**
  `fixed_bigint::*` references. fixed-bigint's lib has `modmath` as
  an optional dep, gated behind the feature.
- **Cargo risk:** dev-dep cycle. modmath's tests need
  `fixed-bigint --features modmath`, which pulls modmath. Cargo
  documents and supports dev-dep cycles (the modmath lib is built
  first, fixed-bigint with feature is built second using that lib,
  tests link both). Not brittle, but unusual.
- **Crates maintained:** existing two.
- **Identity:** modmath is the consumer-of-abstractions; fixed-bigint
  is one such abstraction provider. Direction-inverted from today.
- **Future extensibility:** a new bigint backend wanting modmath's CIOS
  implements `CiosMontMul`/`CiosMontMulCt` directly (the whole
  algorithm). High bar; the backend can't reuse modmath's algorithm
  body because it lives in fixed-bigint.
- **CT contract:** preserved (the impl is private to each backend).
- **Cost:** code move from modmath to fixed-bigint;
  fixed-bigint adds an optional feature; modmath adds dev-dep cycle.

### Option C — New `mont-traits` (or similar) sibling crate

Extract `CiosMontMul`, `CiosMontMulCt`, and the underlying `MulAccOps`
(or whatever lower-level trait the algorithm needs) into a new tiny
crate. Both modmath and fixed-bigint depend on it. The algorithm body
can live in either modmath or fixed-bigint depending on the
backend-implementor cost preference.

- **modmath ↔ fixed-bigint:** decoupled. Both depend on `mont-traits`.
- **Cargo risk:** none. Linear dep graph.
- **Crates maintained:** **three**. The new crate is structurally
  tiny (~50–100 lines) but is a separate release/version cadence.
- **Identity:** each crate has a clear stated purpose. modmath = math
  algorithms; fixed-bigint = fixed-size bigints; mont-traits =
  Montgomery algorithm contracts. Clean.
- **Future extensibility:** ideal home for SOS, FIOS, REDC variants,
  carry-chain primitives, etc. — they all go in `mont-traits`.
  Particularly attractive if more bigint algorithm traits are
  expected.
- **CT contract:** preserved.
- **Cost:** new crate to maintain; one more line in two Cargo.toml
  files.

### Option D — Add as feature-gated module in `const-num-traits`

Put `CiosMontMul`/`CiosMontMulCt` (and possibly `MulAccOps`) in
`const-num-traits` behind a feature flag (e.g.,
`features = ["mont"]`). Both modmath and fixed-bigint already depend
on `const-num-traits`; they enable the feature.

- **modmath ↔ fixed-bigint:** decoupled. Both depend on
  `const-num-traits` (which they already do).
- **Cargo risk:** none.
- **Crates maintained:** existing two + const-num-traits (which is
  already maintained).
- **Identity:** **scope mismatch.** `const-num-traits` is a port of
  `num-traits` plus const-friendly redesigns plus Personality machinery
  — primitive-level traits (`Add`, `BorrowingSub`, `CarryingMul`,
  `Ct`/`Nct`). Putting a Montgomery algorithm trait there muddles what
  the crate is for. The README has to explain why a "numeric
  primitives" crate has a CIOS module.
- **Future extensibility:** if you add SOS, FIOS, REDC, NTT, etc.,
  they all accumulate inside `const-num-traits`, exacerbating the
  identity drift.
- **CT contract:** preserved.
- **Cost:** functionally equivalent to Option C from the downstream
  consumer's perspective. The cost is identity-over-time, not
  function.

### Option E — Move trait and algorithm into fixed-bigint, modmath imports

`CiosMontMul`/`CiosMontMulCt` move into fixed-bigint alongside
`MulAccOps`. modmath's `field.rs` writes
`use fixed_bigint::{CiosMontMul, CiosMontMulCt}` and bounds its
`Field<T, P>::mul` on those traits.

- **modmath ↔ fixed-bigint:** modmath has a direct `use fixed_bigint::*`
  for these traits. Coupling is **stronger** than the status quo, not
  weaker. Does not solve the stated goal.
- **Cargo risk:** none.
- **Crates maintained:** existing two.
- **Identity:** fixed-bigint becomes "bigints + bigint algorithm
  traits," which is reasonable for a fixed-size bigint library that
  has opinions about what algorithms its types support.
- **Future extensibility:** new bigint backends don't get modmath's
  `Field<T, P>` because they can't impl traits defined in
  fixed-bigint (they'd have to fork fixed-bigint or have the trait
  re-exported).
- **CT contract:** preserved.
- **Cost:** trivial code move.

This option is honest about the coupling but doesn't achieve the
stated decoupling goal. Included for completeness.

### Option F — Slice-based `Limbs` trait, drop the discriminant

Replace `MulAccOps` with a simpler:

```rust
pub trait Limbs {
    type Word: Copy + WordOps;
    fn as_limbs(&self) -> &[Self::Word];
    fn as_limbs_mut(&mut self) -> &mut [Self::Word];
}
```

modmath's CIOS becomes a slice algorithm: `cios<W: WordOps>(a: &[W],
b: &[W], m: &[W], n_prime_0: W) -> ...`. The `GetWordOutput` dance
goes away. The CIOS-specific terminology goes away. Heap bigints can
implement it (their `as_limbs` returns the slice they own).

- **modmath ↔ fixed-bigint:** can be decoupled by placing `Limbs` in
  any of the homes above (modmath, fixed-bigint, mont-traits,
  const-num-traits).
- **Cargo risk:** depends on home.
- **Crates maintained:** depends on home.
- **Identity:** `Limbs` is a clean, general trait — "things with
  exposed word-sized limbs."
- **Future extensibility:** any limb-based algorithm (CIOS, SOS, FIOS,
  REDC, NTT, polynomial mul) becomes a generic slice algorithm.
  Maximally reusable.
- **CT contract:** **breaks.** `as_limbs(&self) -> &[Word]` exposes the
  limb array by reference. Indexing into a `&[Word]` by position is
  variable-time; the CT bigint backends use `CtOption` and constant-
  time scans precisely to avoid this. A slice-based API is
  fundamentally incompatible with CT bigints.
- **Cost:** would require giving up CT-via-fixed-bigint, which
  contradicts the modmath Ct path's reason to exist.

So **Option F is rejected for the CT story**, but worth surfacing
because it would be the cleanest design absent the CT constraint. A
hybrid is possible (slice trait for Nct, discriminant trait for Ct)
but doubles the surface.

### Option G — Split MulAccOps into smaller traits

Decompose:

```rust
pub trait CiosRowOps: Sized {
    type Word: Copy + WordOps;
    fn word_count() -> usize;
    fn mul_acc_row(scalar: Self::Word, mult: &Self, acc: &mut Self,
                   carry: Self::Word) -> Self::Word;
    fn mul_acc_shift_row(scalar: Self::Word, mult: &Self, acc: &mut Self,
                         acc_hi: Self::Word) -> Self::Word;
}

pub trait WordAccess: CiosRowOps {
    fn word(&self, i: usize) -> Option<Self::Word>;
}

pub trait WordAccessCt: CiosRowOps {
    fn word_ct(&self, i: usize) -> subtle::CtOption<Self::Word>;
}
```

modmath's `cios_montgomery_mul` bounds on `CiosRowOps + WordAccess`;
`cios_montgomery_mul_ct` bounds on `CiosRowOps + WordAccessCt`. The
`GetWordOutput` discriminant is replaced by two separate traits — the
CT contract is now explicit at the trait level.

- **Cargo / coupling concerns:** orthogonal to splitting; combine with
  any of A–E for placement.
- **Identity:** the rename `MulAccOps` → `CiosRowOps` resolves the
  layering complaint (3b).
- **Future extensibility:** SOS / FIOS could add their own row traits
  alongside `CiosRowOps`.
- **CT contract:** preserved, and now explicit.
- **Cost:** breaks every existing impl. fixed-bigint's `FixedUInt`
  impls would split into three (`CiosRowOps + WordAccess` for Nct,
  `CiosRowOps + WordAccessCt` for Ct). modmath's two `cios_*` function
  bodies stay (still two-character difference in the unwrap), but the
  bound is now self-documenting.

This option is orthogonal to placement — it's about trait *shape*. It
can apply on top of any of A–E.

## 7. Combinations

The trait-shape question (current shape vs. split per G vs. slice-based
per F) is orthogonal to the placement question (which crate owns the
trait). Pure combinations:

| | Current shape | Split (G) | Slice (F) |
|---|---|---|---|
| **A** (status quo) | today | shape-only refactor | rejected (CT) |
| **B** (modmath owns trait) | trait moves to modmath | + split | rejected (CT) |
| **C** (new crate) | trait moves to new crate | + split | rejected (CT) |
| **D** (const-num-traits) | trait moves to cnt | + split | rejected (CT) |
| **E** (fixed-bigint owns) | today's structure | + split | rejected (CT) |

So the practical decision is two-dimensional:

1. **Where does the trait live?** A, B, C, D, or E.
2. **Does the trait get split?** No (keep current shape) or yes (G).

These can be decided independently.

## 8. Open questions

- **Is this a one-shot or the start of a pattern?** If CIOS is the only
  algorithm trait this ecosystem will ever share, the lightweight
  options (A, D, or "B without the algorithm body move") look better.
  If REDC variants, NTT, polynomial mul, carry-chain primitives are
  all coming, then C (sibling crate) ages better — it's the natural
  home for "Montgomery algorithm contracts."
- **Who pays the implementor cost?** Option B raises the bar
  (full algorithm per backend); Option A keeps it low (just impl the
  row primitives). Option C can go either way depending on whether
  the algorithm body lives in `mont-traits` or in the consumer.
- **Heap-bigint backends.** Today they're excluded by `MulAccOps`'s
  `Self: Copy + Default`. None of A–E address this directly. Option F
  would have solved it but is rejected for CT reasons. A "Limbs for
  Nct" + "discriminant trait for Ct" hybrid could re-open it, at the
  cost of two parallel surfaces.
- **CIOS-specific terminology.** Option G's rename
  (`MulAccOps` → `CiosRowOps`) is a clarity win regardless of where
  the trait lives. It can land in any of A–E as a separate change.

## 9. What "purist clean" looks like

Stripped of pragmatism, the cleanest design is:

- `mont-traits` sibling crate owns the trait definitions
  (`CiosRowOps`, `WordAccess`, `WordAccessCt`, `CiosMontMul`,
  `CiosMontMulCt`), per Option C + G.
- Algorithm body lives in the consumer (modmath) so any backend can
  plug in by implementing only the row primitives.
- modmath has zero coupling to fixed-bigint in its lib.
- fixed-bigint has zero coupling to modmath in its lib.
- Each crate has one stated purpose.

The pragmatic costs of this: one new crate to maintain, slightly more
ceremony to update three Cargo.toml files in lockstep when the trait
surface changes, and the burden of justifying that the new crate
exists at all when it's structurally tiny.

The pragmatic alternative — **Option B without a new crate** — gives
the same modmath-side decoupling at the cost of inverting the
dependency direction (fixed-bigint depends on modmath, with cargo
dev-dep cycle in tests). Functionally equivalent to the purist option
for downstream consumers; aesthetically less clean.

Both work. The decision turns on whether "one more small crate" or
"inverted dep direction with a documented cycle" feels more
maintainable in your specific situation.
