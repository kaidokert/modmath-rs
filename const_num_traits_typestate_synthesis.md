# const-num-traits typestate — synthesis (3 reviews + assessment)

Consolidates `const_num_traits_typestate_review_{1,2,3}.md` with the placement
logic from the design conversation (structural-low / domain-high; the by-value +
associated-`Output` convention; CT tiers). Captured 2026-06-16. **Revised after a
final-review pass on this synthesis** (`const_num_traits_synthesis_review_{1,2,3}.md`)
— all three green-lit the direction; the changes below fold in their amendments
and resolve the two splits (D3 `Into` and ship-order).

**Design driver (north star):** *"What's a perfectly-fallen-out-of-the-sky
num-traits API that allows the most beneficial, optimized implementations down
the consumer chain — for those who opt in — assuming nothing but `core::`
exists?"* Everything is derived from that: only `core` is assumed; the measure of
a feature is how much optimization it unlocks for an opt-in downstream
implementer; layering puts each capability where the widest set of consumers can
reach it. No legacy num-traits quirk and no external library (e.g.
`crypto_bigint`'s `Odd`/`NonZero`) is a design input — bridging to them is
downstream *implementation*, settled after the right API exists, never shaping
it.

Two findings independently verified against the source before writing:
- **`Parity for &T` gap is a real compile bug.** `src/ops/parity.rs` impls
  `Parity` only for the 12 primitives (`parity_impl!`), no `&T` impl, no blanket.
  So `Odd::from_ref`'s `for<'a> &'a T: Parity` bound is unsatisfiable for the
  primitive case the sketch claims it serves.
- **Signed division is not fully infallible.** `src/ops/rounding.rs` documents
  its own `MIN / -1` overflow. `DivNonZero` is clean only for unsigned.

---

## 1. Verdict

The sketch's **core instinct is right and worth shipping**: `NonZero` as a
bridge to `core::num::NonZero` (not a duplicate type), with a consuming op
(`DivNonZero`) that exploits the invariant. All three reviews agree on that.

But the sketch ships three real defects, over-includes `Odd`/`Even` (which have
no consuming op at the right layer), and under-includes the two candidates that
actually clear the bar — `PowerOfTwo` and the sign typestates. After fixes, the
crate-resident set is **shorter** than any single review proposed, because two
filters cut the list hard (below).

## 2. The bar (unifying the three reviews + placement logic)

A typestate earns a place **in const-num-traits** only if it clears all three:

1. **Consuming op** (R2/R3's meta-rule). Some op must *exploit* the invariant to
   delete a branch, an `Option`, or a division. A proof with only
   `new`/`get`/`Deref` is dead weight. `NonZero→DivNonZero` passes; `Odd`/`Even`
   as sketched fail.
2. **Right layer** (the conversation's structural-low / domain-high split). The
   consuming op must be a *numeric primitive* op, not a *domain* op. `PowerOfTwo`
   → mask/shift = numeric ✓. `Odd` → Montgomery setup = domain ✗ (belongs up in
   modmath/bignum-traits). R3 already applies this to `Normalized` ("belongs in
   the bignum crate").
3. **Respects conventions.** By-value + associated `Output`; CT tiers (branching
   constructors leak — Tier B); feature-gated; const via `c0nst!`.

This bar is the scope discipline that keeps const-num-traits a core-mirror rather
than a typestate zoo. It is *more* restrictive than the union of the three
reviews, deliberately.

## 3. Must-fix defects in the sketch (all confirmed)

| # | Defect | Source | Fix |
|---|---|---|---|
| D1 | `Odd::from_ref` bound `&T: Parity` unsatisfiable for primitives (no such impl exists) — **compile bug** | R1/R2/R3, verified | Add `impl<T: Parity + Copy> Parity for &T` to `parity.rs` **or** give `from_ref` a `T: Parity + Copy` bound and deref internally (`(*value).is_odd()`). Ship *both* constructors: Copy `new` + `&`-impl'd `from_ref`. **Honest contract (R2/R3):** the blanket covers `Copy` `T`; a non-`Copy` bignum still needs its own `impl Parity for &Big`, so document `from_ref` as "zero-copy for any `T` that provides *borrowed* parity," not a blanket guarantee. |
| D2 | `DivNonZero` returns bare `Self` — violates the by-value + `Output` convention; non-Copy backends can't `impl … for &T` | R3, verified vs CLAUDE.md | Add a fresh `type Output` (can't reuse `Div::Output` — the divisor is `Self::NonZero`, not `Self`). |
| D3 | `type NonZero` unconstrained → trait inert generically (only `from_nonzero` recovers the value) | R1/R2/R3 (split → resolved) | **Do NOT bound `type NonZero: Into<Self>`** (R1 endorsed it, but R3's catch is decisive): `From<NonZero<T>> for T` is *not* a const impl — only `NonZero::get` is `const fn` — so a `[const] Into` bound would poison every generic **const** consumer, defeating the crate's reason to exist. Instead leave `type NonZero` unconstrained and reach the value through a crate-owned **`const fn` accessor** (or a sealed `NonZeroRepr<Base = Self>` companion with a const `get`). And note the generic meaning is mostly carried by the **consuming op** anyway — `DivNonZero::div_nonzero(self, Self::NonZero)` takes the proof directly, so `Into` was solving a non-problem (R2). |
| D4 | `DivNonZero` not infallible for **signed** (`MIN / -1`) | R2, verified | Document the residual overflow, or split unsigned (total) vs signed (still `Option` / needs a `NonMin` co-proof). Unsigned-only to start is cleanest. |
| D5 | `MaybeNonZero: Zero` drags in `Add` (`Zero: Add<Self>`) — unwanted for CT/marker types | R3 | Don't require `Zero`; test zero-ness inside `into_nonzero` via the inherent path. (Or extract an `IsZero` atom — bigger change, defer.) |
| D6 | `typestate` feature absent from `Cargo.toml` | R2 | Add `typestate = []` (crate-owned stable types ⇒ not `nightly-std`). |
| D7 | Primitive `DivNonZero` payoff oversold | R3 | Reframe doc: for primitives `x / nz.get()` already elides the panic via core's niche — codegen win ≈ 0. The real wins are **(a)** the `Option`-free *API*, **(b)** bignum backends' hand-written no-branch path, **(c)** generic code. |
| D8 | CT: branching `Odd::new`/`NonZero::new`/**`PowerOfTwo::new`** leak (all Tier B — popcount/low-bit + bool) | R1/R2/R3 (split → resolved) | Plain constructors are the documented **Tier-B baseline** and **not a day-one blocker** (R2). Under the `ct` feature, provide `CtOption` constructors where matching `ct` primitives exist — a hard requirement for *that* surface (R1), reached incrementally. Reconciles R1's "hard requirement" with R2's "don't block the core API": CT is itself opt-in, so the baseline ships without it. |
| D9 | Naming: `MaybeNonZero` reads Option-shaped; `get(self)` diverges from core | R2/R3 | `HasNonZero`/`NonZeroBridge`; use `into_inner`/`as_inner`, reserve `get` for the core-style accessor. |

## 4. Candidate disposition (the two filters applied)

| Candidate | Consuming op? | Layer | Disposition |
|---|---|---|---|
| **`PowerOfTwo<T>`** (unsigned) | ✓ mask-mod, shift-div, align-up, is-multiple | numeric | **Ship (1st)** — R1+R3 rank it first (Tier-C→Tier-A for generic consumers; all predicates exist: `IsPowerOfTwo`, `ilog2`, `MultipleOf`/`NextMultipleOf`). **Pre-draft decisions (R3):** (a) **store the exponent `k`** (or `k`+mask), *not* the value — so div/mul = `>>k`/`<<k`, mask-mod = `& ((1<<k)-1)`, nothing recomputed; consequence: **no `Deref<Target=T>`** (you'd deref a reconstructed temp), `get()` reconstructs `1<<k` — deliberately *unlike* Odd/NonZero. (b) the consuming ops are **new surface** — a `PowerOfTwoOps` trait / inherent methods taking the proof (`rem_pow2(self, PowerOfTwo<Self>)`, …), **not** a speed-up of blanket `MultipleOf` (no stable specialization). |
| **`NonZero` bridge + `DivNonZero`/`RemNonZero`** | ✓ division | numeric | **Ship (2nd)** (R2 argues 1st — see §8), fix D1–D9. Unsigned-total first; signed later. |
| **`NonNegative<T>` / `Positive<T>`** (signed-only) | ✓ infallible→unsigned cast, branch-free `abs`, total `isqrt` | numeric | **Ship (3rd).** Reuses `sign.rs`, `convert.rs`, `CheckedIsqrt`. **Semantics (R2):** `Positive = >0`, `NonNegative = >=0`, **signed-only** — on unsigned `Positive ≡ NonZero` and `NonNegative ≡ all`, so neither earns a separate type there. |
| **`Odd<T>` / `Even<T>`** | ✗ in-crate (Odd's consumer is Montgomery; Even's is exact-halving) | **domain** | **Do NOT put consuming API in const-num-traits.** Odd's definition *may* sit low next to `Parity`, but its ops are modmath/bignum-traits'. `Even` → **cut** unless paired with `DivExact`-by-2 (marginal). |
| `Prime<T>` / `OddPrime<T>` | ✓ but inverse existence = domain | **domain** | **modmath**, not here (R1 proposed it for the crate — rejected on the layer filter; matches modmath's own `typestate_synthesis.md`). |
| `Finite<T>` | ✓ Ord-from-PartialOrd | numeric but **float axis** | **Defer/low.** Different axis (float guards); doesn't compose with the int typestate set. Revisit if/when a float-decomposition slice grows. |
| `BitIndex<T>` / `ShiftAmount<T>` (`n < BITS`) | ✓ removes `Option` from checked shifts, panic from funnel | numeric | **Defer.** Real, but the proof is *width-dependent* (awkward ergonomics). Good later candidate. |
| `Normalized<T>` (high bit set) | ✓ Knuth division | backend-specific | **bignum crate**, not here (R3 agrees). |
| `NonMin<T>` / `NonMax<T>` | thin (±1 can't-overflow; `MIN` abs/neg) | numeric | **Low.** Mainly useful as the signed co-proof for D4. Park until `DivNonZero`-signed needs it. |

## 5. The proof lattice (what makes it a *system*, not three newtypes)

R3's key addition, and the one most aligned with the north star. The families
aren't independent — they **refine**, via cheap `const fn` narrowings that let
stronger proofs *reuse* weaker proofs' ops instead of duplicating them:

- **`PowerOfTwo<T>` ⊂ `NonZero<T>`** — a power of two is never zero. So div/rem by
  a `PowerOfTwo` is infallible *by construction*: it narrows to `NonZero` with no
  recheck and reuses `DivNonZero`. (And its own shift/mask ops are even cheaper.)
- **`Positive<T>` ⊂ `NonNegative<T>` → unsigned cast** — `>0` implies `>=0`; a
  `NonNegative` signed value casts to unsigned with no `Option`.
- **Constructive answer to D4** — signed total division is `NonZero ∧ NonMin`,
  *composed from the lattice* rather than special-cased: a divisor carrying a
  "not `-1`"/`NonMin` fact removes `MIN / -1` by construction. This is the real
  role of the otherwise-parked `NonMin<T>`.

Ship the refinement paths as `const fn` (`PowerOfTwo → NonZero`,
`Positive → NonNegative → unsigned`). This is what turns the "proof trifecta"
into a coherent system — maximal reuse of optimization hooks downstream, which is
exactly what the north star asks for.

## 6. Cross-cutting design constraints

- **CT personality.** Every branching constructor needs a `ct`-feature
  `CtOption` sibling; plain ctors documented non-CT. This is a hard requirement
  for a crypto/modmath consumer, not a nicety.
- **Scope discipline.** The union of the three reviews would add ~9 typestates.
  The bar in §2 cuts that to **three** crate-resident families (NonZero bridge,
  PowerOfTwo, sign) + Odd's bare definition. That shrinkage *is the
  recommendation* — resist the zoo.

**Not a design input — interop is implementation.** The API is designed as if
only `core::` exists, to maximize the optimizations an opt-in downstream
implementer can provide. The existence of `crypto_bigint`'s `Odd`/`NonZero` (or
any other library's wrappers) does **not** shape the abstraction. Coexisting
with such types in a downstream that uses both is mechanical bridging
(`RefCast`/re-export/`From`) that happens *after* the right API exists; it is a
footnote, never a driver.

## 7. Where I differ from the reviewers

- **R1 over-includes** (`Prime`, `Finite` in const-num-traits). Both fail the
  layer filter — `Prime` is domain (modmath), `Finite` is a float-axis guard.
  R1's enthusiasm is right about value, wrong about *home*.
- **R1 vs R3 on the `DivNonZero` primitive payoff** — R3 is correct (codegen win
  ≈ 0 for primitives via the niche). But I'd keep R1's broader point: the
  `Option`-free *API* and the generic/backend abstraction are real wins, so the
  trait is still worth it — just market it honestly (D7).
- **External libraries are not a design driver.** An earlier pass (mine
  included) treated "align with `crypto_bigint` vs. own it" as a gating decision.
  That was a category error: it conflated *what the maximally-correct abstraction
  is* (the design question, answered from `core` + correctness + layering) with
  *how it coexists downstream* (plumbing). Design the ideal opt-in API; bridging
  is implementation, settled later and independently.
- **Odd: not a "give it a consumer" fix (R3), a placement fix.** Odd *has* a
  consumer — Montgomery — but it's domain, so the fix isn't "add an op in
  const-num-traits," it's "Odd's ops live up in modmath/bignum-traits; only its
  definition (maybe) sits low." Subtle but it changes where the work goes.

## 8. Recommended sequence

**Order split (resolved): PowerOfTwo first.** R1 + R3 rank `PowerOfTwo` first
(highest generic payoff — Tier-C→Tier-A — and it sidesteps NonZero's signed +
const-accessor complications); R2 argued `NonZero` first (simpler bridge, fewer
moving parts). Resolution: **PowerOfTwo first** on the 2-of-3 + payoff weight,
with R2's escape hatch honored — if real downstream `NonZero` div/rem call sites
are waiting and no `PowerOfTwo` ones are, swap 1↔2. It's demand-driven, not
dogma; both are early.

1. **`PowerOfTwo<T>`** — decide representation (exponent[+mask], §4) first, then
   the `PowerOfTwoOps` consuming ops (`rem_pow2`/`div_pow2`/`align_up`) in the
   `c0nst!` / by-value / associated-`Output` style. Highest payoff; all
   predicates already exist. (Reviewer 3 offered to draft this.)
2. **`NonZero` bridge + unsigned `DivNonZero`/`RemNonZero`**, fixing D1–D9.
   Unsigned-total only; signed deferred behind a `NonMin` co-proof or `Option`.
3. **`NonNegative<T>` / `Positive<T>`** + infallible cast / branch-free abs /
   total isqrt.
4. Park `BitIndex`/`ShiftAmount`, `NonMin`/`NonMax`, `Finite`. Send `Prime` to
   modmath, `Normalized` to the bignum crate.

Interop with any existing library (`crypto_bigint`, etc.) is *not* a step here —
it's downstream bridging, addressed only if and when a consumer needs both.

Across all of it: feature-gate `typestate`; `ct`-feature CtOption constructors;
fix the `Parity for &T` gap (D1) regardless, since `from_ref` is load-bearing for
every non-Copy wrapper, not just `Odd`.
</content>
