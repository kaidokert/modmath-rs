# CIOS trait placement — synthesis of four analyses

Factors `ANALYSIS_CLAUDE.md`, `ANALYSIS_CODEX.md`, `ANALYSIS_GEMINI.md`,
`ANALYSIS_GEMINI2.md` against the exploration (`cios_trait_placement.md`).
Captured 2026-06-17.

Four independent analyses converged on a large shared core and split on exactly
one axis. This document states the consensus, resolves the split with the one
new idea that reframes it, and gives a unified staged plan.

---

## 1. Where all four agree (the settled core)

These are unanimous (or 3-of-4 with the fourth silent, not opposed). Treat them
as decided:

1. **`MulAccOps` is not a numeric primitive — it is a CIOS-specific limb-row
   contract.** Putting the *algorithm* trait in `const-num-traits` repeats the
   typestate "zoo" problem in trait form (CLAUDE, CODEX explicit; GEMINI/GEMINI2
   agree the algorithm layer is "domain-high").
2. **Do the split + rename now (Option G).** `MulAccOps` → `CiosRowOps`, and
   replace the `GetWordOutput` associated-type discriminant with explicit word-
   access traits. The CT contract moves from a hidden assoc-type trick to a
   self-documenting bound. **Unanimous.**
3. **The CIOS algorithm body stays in modmath.** Do *not* move it into
   fixed-bigint (Option B / E) — that forces every future backend to re-derive
   CIOS instead of implementing only the row primitives. **Unanimous, emphatic.**
4. **Do not put `CiosMontMul` in fixed-bigint**, and do not adopt a slice-based
   `Limbs` trait as the *primary* contract (it breaks the CT story). **Unanimous.**
5. **Generalize for heap backends while you're in there:** make
   `word_count(&self)` an instance method and relax/relocate the `Copy + Default`
   bound (CLAUDE, GEMINI, GEMINI2 explicit; CODEX silent). Cheap door-opener; do
   it now, don't refactor to force it.

The exploration's Options D-as-stated (algorithm trait in cnt), E (fixed-bigint
owns the high trait), and F-as-primary (slice `Limbs`) are **rejected by
consensus**.

## 2. The one new idea that reframes the problem (from GEMINI)

My layering called the contested middle a single layer (L2). Gemini correctly
splits it in two, and that split dissolves most of the remaining argument:

| Layer | What it is | Home |
|---|---|---|
| **L1 — word arithmetic** | `CarryingMul`/`CarryingAdd`/`WideningMul` | **cnt** (already there) |
| **L2a — word/limb *access*** | "expose my i-th limb" — a *representation* primitive | **cnt** ← the new insight |
| **L2b — CIOS *row ops*** | `mul_acc_row` / `mul_acc_shift_row` — CIOS-specific | modmath **or** new crate |
| **L3 — algorithm + capability** | the double-loop, `CiosMontMul`/`CiosMontMulCt` | **modmath** |

**Why L2a belongs in cnt and L2b does not** — and why this is *not* a violation
of the scope discipline that killed Option D: limb access is a representation
capability, exactly the same category as cnt's existing `ToBytes`/`FromBytes`/
`NumBytes`. "How do you read your words" is structural; "how does CIOS fold a
row" is an algorithm. Option D was wrong for the *algorithm* trait; it is *right*
for the *access* trait. The CT variant (`WordAccessCt`, returning
`subtle::CtOption`) rides cnt's existing `ct` feature — no new dependency
surface.

This bisects the GEMINI2-vs-(CLAUDE/CODEX) dispute about Option D: **the
legitimate half of D (structural access → cnt) and the illegitimate half (CIOS
algorithm → not cnt) were being argued as one thing.** Separated, both sides are
right about their half.

It also *shrinks the contested trait*: once L2a is in cnt, the only thing whose
home is still open is L2b — two row-op methods. A whole new crate for two
methods is heavy; that tilts the short-term answer (see §5).

## 3. Pre-step that may shrink it further (from CLAUDE, unique)

Before committing to *any* `WordAccess` / `WordAccessCt` split, probe whether the
discriminant is load-bearing at all:

> In CIOS the word index is a **public loop counter** (`0..word_count`), never
> secret. A guaranteed-in-bounds `fn word(&self, i: usize) -> Self::Word` is
> data-oblivious for *both* personalities — no `Option`, no `CtOption`. The
> Some/None side-channel argument applies to a *fallible, position-revealing*
> lookup at the API boundary, not to an in-bounds read inside the algorithm.

If this holds against the real CT threat model, L2a collapses to a **single
infallible accessor** (Tier A, no `subtle` dependency, no Nct/Ct split), and the
two near-identical algorithm bodies in `cios.rs` merge into one. That is a bigger
simplification than Option G itself, and it makes the cnt placement of L2a
trivially clean. **Verify first** — the CT model is the maintainer's to confirm.
If it does *not* hold, keep the two-trait split (Gemini/Codex shape), with
`WordAccess` Tier B / `WordAccessCt` Tier A.

## 4. Resolving the secondary disagreements

**Naming — `CiosRowOps` (CLAUDE, CODEX) vs generic `LimbRowOps`/`LimbShiftRowOps`
(GEMINI2).** Side with the honest CIOS-specific name. The exploration (§3b) is
explicit that the "shift" in `mul_acc_shift_row` is the *CIOS* post-multiply
word-shift, not a generic op; SOS/FIOS have different row structures. A generic
name over-promises reuse that the trait can't deliver. If SOS/FIOS arrive, they
get their own `SosRowOps` — and *that* is what makes them first-class, not a
misleadingly-broad name. Reuse is served by the layering (L1/L2a shared), not by
renaming L2b.

**Heap accumulator as `&mut [Word]` (GEMINI2).** This is a *defensible* slice use
— unlike rejected Option F, it exposes only the mutable *accumulator* (scratch),
indexed by *public* positions, not the secret *operands*. It's CT-safe so long as
no branch is taken on the values. Treat it as the concrete mechanism for the §1.5
heap generalization if/when heap backends are actually targeted; don't let it
reshape the row-op signature pre-emptively.

**Timing — CODEX's two-stage path vs go-straight-to-final.** Codex's instinct
(don't mint a crate prematurely) is right, but Codex's Stage 1 leaves `CiosRowOps`
*in fixed-bigint*, so modmath still `use`s `fixed_bigint::*` — it does the shape
cleanup without achieving the stated decoupling. The fix is to combine Codex's
staging with my Option H (below): decouple *now* without a crate, extract to a
crate *later*.

## 5. The one residual hinge: where L2b lives

After §2 moves L2a to cnt, only L2b's home is open. Three live answers:

- **H — modmath owns L2b** (my formalization of "B without moving the
  algorithm"): modmath defines the tiny `CiosRowOps`, owns the algorithm + the
  blanket impl; fixed-bigint impls `CiosRowOps` under an optional `modmath`
  feature. modmath lib drops `fixed_bigint` **now**. Cost: dev-dep cycle
  (modmath tests ↔ fixed-bigint feature), no new crate.
- **C — new `mont-traits` crate owns L2b**: linear dep graph, +1 crate. Best if
  more shared algorithm contracts are coming.
- **CODEX Stage-1 — L2b stays in fixed-bigint**: lowest risk, but **does not
  decouple yet** (modmath still imports fixed-bigint).

The structural split changes the calculus: because L2b is now just two methods
and L2a is already low, **H's dep-inversion carries almost nothing**, and a whole
crate for two methods (C) looks premature *today*. So the synthesis ordering is:

> **H now → C later.** Decouple immediately via H (faithful to the stated goal,
> no premature crate); extract `CiosRowOps` (+ `CiosMontMul` siblings) into
> `mont-traits` when a second backend or a second Montgomery variant (SOS/FIOS)
> actually appears. By then L2a is already in cnt and L2b is already cleanly
> separated, so the extraction is mechanical — exactly Codex's Stage 2.

This honors all four: Gemini's structural split, my H + decoupling-now, Codex's
"don't mint the crate early / extract later," and everyone's G + heap fixes.

The only judgment call left for the maintainer: **is decoupling wanted *now*?**
If yes → H (accept the dev-cycle). If the team is fine deferring decoupling until
a second backend exists → Codex Stage-1 in-place is lower-risk. Since decoupling
modmath from fixed-bigint is the *explicit objective of the whole exercise*, the
synthesis recommends **H**.

## 6. Unified plan

**Phase 1 — shape + structural move (low risk, unanimous core):**
1. Run the §3 discriminant probe. Decide: single infallible accessor, or the
   `WordAccess`/`WordAccessCt` split.
2. Move L2a (word/limb access) into `const-num-traits` as a representation trait,
   feature-gated (e.g. `limb`); CT variant under the existing `ct` feature.
   Reuses L1 (`CarryingMul`/`CarryingAdd`) already present.
3. Rename `MulAccOps` → `CiosRowOps`; drop `GetWordOutput`; make
   `word_count(&self)` an instance method; relax `Copy + Default`.
4. Keep both algorithm bodies in modmath (or one body if the probe collapsed
   them).

**Phase 2 — decouple now (Option H):**
5. modmath owns `CiosRowOps` + the algorithm + the blanket
   `impl<T: CiosRowOps + WordAccess> CiosMontMul for T`.
6. fixed-bigint impls `CiosRowOps` under an optional `modmath` feature.
7. modmath lib drops all `fixed_bigint::*` references; fixed-bigint becomes a
   modmath **dev**-dependency for tests.

**Phase 3 — extract when justified (Option C, deferred):**
8. On the second backend or second Montgomery variant, lift `CiosRowOps`
   (+ `CiosMontMul`/`CiosMontMulCt`) into `mont-traits`. Mechanical.

## 7. Decision matrix

| Question | Synthesis answer | Backing |
|---|---|---|
| Algorithm trait in cnt? | **No** | CLAUDE, CODEX (GEMINI/GEMINI2 agree it's domain) |
| Word *access* in cnt? | **Yes** (representation primitive, like `ToBytes`) | GEMINI (new insight) |
| Split + rename `MulAccOps`? | **Yes → `CiosRowOps`** | unanimous (G) |
| Generic vs CIOS name for row ops? | **CIOS-specific** | CLAUDE, CODEX; exploration §3b |
| Keep `GetWordOutput` discriminant? | **Probe first; likely drop** | CLAUDE (unique) |
| Algorithm body home? | **modmath** | unanimous |
| Move algorithm to fixed-bigint? | **No** | unanimous |
| Decouple now without a crate? | **Yes — Option H** | CLAUDE |
| New `mont-traits` crate? | **Later, when a 2nd algorithm/backend lands** | CODEX staging + CLAUDE |
| Heap backends? | **Door-open now** (`word_count(&self)`, drop `Copy+Default`); `&mut [Word]` acc when actually targeted | GEMINI, GEMINI2 |

## 8. One-line summary

Split `MulAccOps` three ways: **word arithmetic and limb access go down into
`const-num-traits`** (representation primitives), **the CIOS row ops become a tiny
`CiosRowOps` owned by modmath now (Option H, decoupling immediately) and extracted
to `mont-traits` later (Option C, when a second algorithm justifies it)**, and the
**algorithm body never leaves modmath** — after first checking whether public loop
indices let you delete the CT discriminant entirely.
