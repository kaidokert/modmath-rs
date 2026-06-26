# CIOS trait placement — analysis & recommendation

Companion to `cios_trait_placement.md` (the option-space exploration). That
document deliberately stops at "no recommendation." This one reframes the
problem, disposes of the dead options, formalizes the option the exploration
buried, and lands on a recommendation. Captured 2026-06-17.

---

## 1. Reframe: three layers, not one trait

`MulAccOps` conflates three distinct layers. Naming them makes the placement
answer fall out:

| Layer | What it is | Who implements it | Where it belongs |
|---|---|---|---|
| **L1 — word arithmetic** | carrying mul/add on the limb word (`CarryingMul`, `CarryingAdd`, `WideningMul`) | the word type (primitive) | **already in `const-num-traits`** (`ops/bigint.rs`) |
| **L2 — limb/row contract** | `word_count`, word access, the two `mul_acc_*_row` ops | the bigint backend (fixed-bigint) | the contested layer |
| **L3 — algorithm + capability** | the CIOS double-loop, `CiosMontMul`/`CiosMontMulCt` | nobody re-implements it — it *is* the algorithm | **modmath** |

L1 already lives in the right place. That matters twice over:

- the L2 row ops are, in principle, *default-implementable* from L1 + word
  access — a backend overrides them only for hand-tuning / CT, and
- it sets the clean rule that kills Option D (below): cnt owns the *primitive*
  layer; the *algorithm* contract must not drift in beside it.

The whole dilemma is really "where does **L2** live, and what shape does it
have," with L1 and L3 already pinned.

## 2. The decision-crux the exploration leaves implicit

Strip it down and there is one structural law:

> To decouple two crates that **share a trait**, the trait must sit somewhere
> both can name without one depending on the other. With only two crates,
> *someone* must depend on the other — a cycle. A third (neutral) crate is
> precisely the thing that breaks the cycle.

So the entire decision reduces to:

> **two crates ⇒ a dependency cycle somewhere, vs three crates ⇒ linear graph +
> one more crate to maintain.**

Everything else (implementor cost, CT contract, extensibility) is *equal*
across the live options as long as the algorithm body stays in modmath. Stating
the crux this way is what makes the choice answerable.

## 3. Dispositions of the exploration's options

- **D (put L2/L3 in `const-num-traits`) — reject**, on the same scope
  discipline established for that crate in the typestate work. A CIOS Montgomery
  contract is a *domain* trait; the structural-low / domain-high filter sends it
  *up*, not into the primitive-mirror crate. The clean line is exactly the
  L1/L2 split: cnt already owns the word arithmetic; the row/algorithm contract
  must not accrete beside it. (Same reasoning that sent `Prime` to modmath, not
  cnt.)
- **E (fixed-bigint owns the high trait) — reject.** It *strengthens*
  modmath→fixed-bigint coupling; fails the stated goal. Dismiss.
- **F (slice `Limbs`) — reject for CT, as the exploration says** — but the
  discriminant probe in §5 recovers most of F's cleanliness *without* breaking
  CT.
- **A (status quo) — the honest baseline.** Loses only to the decoupling goal.

That leaves the real contest: **B vs C** — and the exploration under-specifies
the best member of the B family.

## 4. New option: H — "B without moving the algorithm"

Option B does the right decoupling but then makes a wrong move: it relocates the
*algorithm body* into fixed-bigint, forcing every future backend to re-derive
CIOS. Section 9 and Open-Question #2 of the exploration gesture at "B without the
algorithm body move" but never formalize it. It is the actual sweet spot, so
name it:

**Option H.** modmath owns the L2 contract (`CiosRowOps` +
`WordAccess`/`WordAccessCt`), the L3 algorithm body, *and* the blanket
`impl<T: CiosRowOps + WordAccess> CiosMontMul for T`. fixed-bigint gains an
optional `modmath` feature and under it writes only
`impl modmath::CiosRowOps for FixedUInt<…>`.

- **modmath lib:** *zero* `fixed_bigint::*` references. ✓ (the goal)
- **New-backend cost:** implement only the row trait → `CiosMontMul` for free.
  ✓ (this is B's fatal flaw, fixed)
- **Cargo:** fixed-bigint takes `modmath` as an *optional* (non-dev) dep;
  modmath takes fixed-bigint as a *dev* dep for tests. Both publish fine; the
  cycle materializes only in a feature-on test build, which Cargo supports.
  Same Cargo cost as B, none of B's extensibility loss.

**H strictly dominates B.** The only remaining axis is H (cycle, two crates) vs
C (no cycle, three crates) — i.e. exactly the §2 crux.

## 5. Trait shape: endorse G, and push it further

Independent of placement, **do the split + rename (Option G)**:
`MulAccOps` → `CiosRowOps`, and replace the `GetWordOutput` associated-type
discriminant with explicit `WordAccess` (`Option`) / `WordAccessCt` (`CtOption`)
traits. This fixes the misnomer (exploration §3b) and lifts the CT contract from
a hidden associated-type trick to a self-documenting bound. You break impls
anyway; do it once.

But **probe whether the discriminant is load-bearing at all before splitting**:

> In CIOS the word index is a **public loop counter** (`0..word_count`), never
> secret. So an in-bounds `fn word(&self, i: usize) -> Self::Word` is
> data-oblivious for *both* personalities — no `Option`, no `CtOption`. The
> Some/None side-channel argument (§3a) applies to a *fallible,
> position-revealing* lookup at the API boundary, not to a guaranteed-in-bounds
> read inside the algorithm.

If that holds against the real CT threat model, don't split `WordAccess` in two
— **delete the discriminant entirely**, collapse the two near-identical
algorithm bodies into one, and recover most of Option F's cleanliness *without*
breaking CT (CT-ness then lives only in the row-op impls, which are already
oblivious). That is a bigger simplification than G itself. Verify this first; it
may make the `WordAccess`/`WordAccessCt` split moot.

*(Flagging, not asserting — the CT threat model is yours to confirm.)*

## 6. Bonus axis no option fixes: heap backends (§3c)

While rewriting the trait, make `word_count(&self)` an instance method and drop
the `Copy + Default` bound (replace with what the accumulator genuinely needs —
likely `Clone` + an explicit zero-acc constructor). That alone admits
`BoxedUint` / `BigUint` as CIOS providers, at near-zero cost *at trait-definition
time*. If accumulator construction turns out to entangle with
`Default`-means-statically-sized, don't chase it — document the fixed-size
restriction and move on. Cheap to leave the door open; not worth a refactor to
force it.

## 7. Recommendation

Two clean exits; the hinge is the exploration's Open-Question #1 (one-shot vs
pattern):

- **If more shared algorithm contracts are coming** (SOS / FIOS / REDC / NTT …)
  — and modmath is a *modular-arithmetic* library, so those read as roadmap, not
  hypotheticals — go **C + G**: a `mont-traits` (or `bigint-algo-traits`)
  sibling crate owns L2's contract; algorithm body stays in modmath; cnt keeps
  L1. Linear graph, each crate one purpose; the new crate earns its keep the
  moment the second algorithm lands.
- **If CIOS is genuinely the last shared algorithm you foresee** — go **H + G**:
  no new crate, same modmath-side decoupling, low backend cost, dev-dep cycle as
  the only tax.

**Default lean: C + G.** The signals point at "pattern": a growing math library,
the exploration's own future-extensibility column favoring C in 3 of 4 rows, and
cycle-avoidance that compounds with each added contract. Extracting a *clean,
already-separated* module into a crate later is mechanical; starting inside a
cycle and unwinding it later is not. If CIOS is confirmed one-shot, switch to
**H + G** without hesitation.

In **every** variant:

1. Keep the algorithm body in **modmath**.
2. Apply **G** (split + rename `MulAccOps` → `CiosRowOps`).
3. Run the **discriminant probe** (§5) before committing to the
   `WordAccess`/`WordAccessCt` split — it may collapse to one trait + one
   algorithm body.
4. Make `word_count` take `&self` and drop `Copy + Default` while you're in
   there (§6).

## 8. Summary table

| Option | modmath→fixed-bigint coupling | Cargo | Crates | New-backend cost | Verdict |
|---|---|---|---|---|---|
| A status quo | coupled (1 site) | none | 2 | impl row trait | baseline |
| B | decoupled | dev cycle | 2 | **whole algorithm** | dominated by H |
| C new `mont-traits` | decoupled | linear | **3** | impl row trait | **rec. if pattern** |
| D const-num-traits | decoupled | none | 2(+cnt) | impl row trait | reject (scope) |
| E fixed-bigint owns | *stronger* | none | 2 | n/a | reject (goal) |
| F slice `Limbs` | any | any | any | impl `Limbs` | reject (CT) |
| **H** (new here) | decoupled | dev cycle | 2 | impl row trait | **rec. if one-shot** |
| G split/rename | — | — | — | — | apply on top of any |
