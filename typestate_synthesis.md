# modmath-rs typestate roadmap (synthesis)

Consolidates input from five external design-review documents
(`codex_typestate_1.md`, `gemini_typestate_1.md`, `gemini_typestate_2.md`,
a follow-up gemini pass adding `Borrow<T>` consolidation and
`Coprime<A, M>` framing, and a follow-up codex pass adding 20+ more
patterns including modulus-shape typestates, algorithm-family split,
and CIOS infallibility), prior in-flight decisions on this branch, and
architectural questions about where typestate wrappers should live.

## 0. Two axes of typestate

Worth distinguishing up front, because the early lists conflated them:

- **Input-shape wrappers** — `Odd<T>`, `NonZero<T>`, `Reduced<T>`, `Wide<T>`, `Canonical<T>`. Encode "this value satisfies a structural property." Easy to construct, easy to reason about, mostly improve safety and remove redundant checks.
- **Algebraic / security-contract wrappers** — `Prime<T>`, `Unit<'m, T>`, `NonZeroResidue<'m, T>`, `Public<T>`, `Secret<T>`, `SpecialModulus<T>` family. Encode "this value satisfies a property that unlocks an algorithm" or "this value has a timing/secrecy contract." Construction is harder (sometimes requires `unsafe new_unchecked` for static constants), but they unlock real new APIs (Fermat inverse, infallible division, type-checked CT discipline).

Both are valuable. The algebraic axis tends to have higher leverage per wrapper added; the input-shape axis is broader and more incremental.

---

## 1. Conceptual framing: traits vs. typestates

The two are orthogonal. They compose, but they live in different
axes of the design.

| Axis | Example | Role | Where backend optimization lives |
|---|---|---|---|
| **Trait** | `Parity`, `PartialOrd`, `Rem`, `WideMul`, `CiosMontMul` | Behavioral interface. Implementations *do* the runtime work. | Inside each backend's impl. |
| **Typestate wrapper** | `Odd<T>`, `Reduced<T>`, `MontForm<T>`, `NonZero<T>` | Static proof that a runtime check already happened. `repr(transparent)`, zero cost. | None — the wrapper itself has nothing to specialize. |

The runtime check (e.g. `is_odd()`) happens **once** at typestate
construction. Everything downstream that accepts `&Odd<T>` is a
constant-time pointer reinterpretation; backends have nothing to
contribute there.

This dissolves the user's stated concern about backends being unable
to provide optimized implementations — they specialize the
**underlying trait** (already happens via num-integer feature gate),
not the wrapper.

## 2. Where the typestate wrappers should live

Three plausible homes:

### Option A — in `modmath` (current direction)
- **Pros:** self-contained; fast iteration; ships with the consumer crate.
- **Cons:** other crypto crates (curve25519-dalek, RustCrypto stack, RSA implementations) can't share these types without depending on modmath, which couples them to a specific Montgomery implementation they may not want.
- **Implication for backend optimization:** none. Backends impl the traits, not the wrappers.

### Option B — in `num-traits` (or a fork)
- **Pros:** natural home for "generic number abstractions"; widely depended on.
- **Cons:** num-traits is conservative about adding surface; PR cycle is slow; upstream is unlikely to accept opinionated typestate wrappers like `MontForm<T, P>`.
- **A fork** creates ecosystem fragmentation — anything depending on the fork can't easily compose with code expecting upstream num-traits.

### Option C — a new dedicated crate (e.g. `num-typestate` or similar)
- **Pros:** clean separation; can ship without waiting on num-traits; reusable across modmath, RustCrypto, and other downstream crates.
- **Cons:** adoption problem — needs to reach critical mass to be worth the dep; if it's just-for-modmath, A is simpler.
- **Sweet spot if:** several downstream consumers explicitly want the wrappers (rsa-heapless, krabipqc, an ed25519 layer, etc.) and modmath isn't the right umbrella for all of them.

**Recommendation:** Start with **Option A** (modmath) for the wrappers
that are tightly bound to Montgomery semantics (`MontForm`, `MontParams`).
For the more universal wrappers (`Odd`, `NonZero`, `Reduced`, `Wide`),
keep the option open to extract them into Option C later. They're
small enough that a future extraction is mechanical (rename + republish).

## 3. Wrapper inventory (consensus across all input docs)

Marked **★** for high-consensus, immediately actionable; **○** for
plausible but lower-priority or context-dependent.

| Wrapper | Invariant | Status | Notes |
|---|---|---|---|
| ★ `Odd<T>` | `T: Parity, value.is_odd()` | **Previously abandoned** — see §5 | All three docs propose this. Existing branch `feat/odd-modulus-accessor` at `335795d` (local). |
| ★ `NonZero<T>` | `value != 0` | New | Pre-req for `Unit<T>` and `OddPrime<T>`. Mirrors stdlib's `NonZeroU32`/`NonZeroU64`/etc., so the pattern is already familiar; the generic version extends it to bigint backends. Removes the `if modulus == zero { return None }` checks at the top of every Montgomery free function. |
| ★ `Reduced<'m, T>` | `value < modulus` | **Previously deferred** — see §5 | All three propose; replaces the `_pr` suffix at the schoolbook surface. |
| ★ `MontForm<T, P>` | "in Montgomery form, branded by some modulus" | Existing `Residue` is this at the Field level — proposal is to push it down into the free-function surface | Codex §3, gemini_1 §3. Removes the silent-correctness risk in `basic::montgomery::wide::*` calls. |
| ★ `MontParams<T>` | precomputed Montgomery context: modulus + n_prime + r_mod_n + r2_mod_n | Existing `Field<T, P>` already stores these — proposal is to make `MontParams` a separable type so consumers can have the parameters without the full Field | Codex Priority 1+6. Lets standalone Montgomery functions skip per-call precompute. |
| ★ `Wide<T>` | named `lo`/`hi` instead of tuple | Easy, isolated change | Codex §4, gemini_1 §4. Removes `(hi, lo)` swap-bug class. |
| ○ `Prime<T>` / `OddPrime<T>` | n is prime (resp. odd prime) | New | gemini_2 dedicated doc; codex §13. Unlocks infallible inverse, Fermat path. Construction is `unsafe new_unchecked` for known constants — runtime primality testing is impractical at crypto sizes. |
| ○ `Unit<'m, T>` / `NonZeroResidue<'m, T>` | coprime to modulus (resp. nonzero in a prime field) | New | Codex §7. Infallible inverse where invertibility is statically known. |
| ○ `Coprime<A, M>` | `gcd(a, m) == 1` (pairwise relation, not branded to a context) | New | Gemini follow-up. Different shape from `Unit<'m, T>`: two-parameter wrapper proving a *pair* is coprime rather than "this residue is a unit in some modulus context". Right shape for RSA key-gen paths where `e` is chosen coprime to `φ(N)` — proves `e^{-1} mod φ(N)` is infallible. Could be a constructor variant of `Unit<'m, T>` rather than a separate type. |
| ○ `PublicExp<T>` / `SecretExp<T>` (or `Public<T>` / `Secret<T>` for any operand) | timing-contract on operand | New | Codex §8 frames it as exponent-specific; gemini follow-up generalizes to any operand. Compiler-checked replacement for the current `exp` vs `exp_public_exp` naming convention. The generalized form (`Secret<T>`) can also refuse to compile `_pr` / variable-time entrypoints when fed secret data — pushes side-channel safety into the type system rather than developer discipline. Composes with the `subtle` / `secrecy` ecosystem. |
| ○ `Loose<'m, T>` / `WideLoose<'m, T>` | `value < 2n`, `value < 4n` (lazy reduction) | New | Codex §10. Skips cond-subtract in addition-heavy code (curve formulas). Requires backend slack. |
| ○ `Canonical<'m, T>` vs `Raw<T>` | guards signed `%` returning negative | New | Codex §14. Mostly a correctness guard, not a perf win. Codex follow-up notes the `_pr` APIs actually require `[0, m)` (canonical) not just congruence — so `Canonical<'m, T>` may be the better name than `Reduced<'m, T>` for the fast-path wrapper, with `Congruent<'m, T>` reserved for mere modular equivalence. |
| ○ `NonTrivialModulus<T>` / `ModulusGtOne<T>` | `modulus > 1` | New | Codex follow-up §3. Many exponentiation paths special-case `modulus == 1` with a top-of-function branch (`if modulus == one { return zero }`). A context proving `n > 1` lets every operation skip that branch. Likely subsumed into `OddModulus<T>` if that constructor also rejects 1 (Montgomery modulo 1 is useless anyway). |
| ○ `SpecialModulus<T>` family: `PseudoMersenne`, `Solinas`, `Crandall`, `Mersenne`, `Goldilocks` | modulus has a special form enabling specialized reduction | New | Codex follow-up §10. Curve25519's `2^255 - 19`, Ed448's `2^448 - 2^224 - 1`, NTT-friendly primes, etc. These have reduction paths *much* faster than generic Montgomery. Belongs in specialization layers (the field.rs docs already point at this — `lazy_add`, Solinas reduction). Shared vocabulary lets multiple consumer crates compose. |
| ○ `MontParams<T, RKind>` with `RKind = FullWidthR` or `RGreaterThanN` | algorithm-family marker on the Montgomery context | New | Codex follow-up §11. The crate has two Montgomery algorithm families — `R > N` (uses `compute_params`, `to_mont`, etc.) and `R = 2^W` (uses `wide::redc`, `wide::mul`). Currently distinguished by function/module naming. A type parameter would prevent accidentally mixing parameters from one family with functions from the other. |
| ○ `OverflowSafe<T>` / `SmallModulus<T>` | modulus is small enough that intermediate products don't overflow | New | Codex follow-up §12. Some `R > N` helpers (`basic_from_montgomery`, `basic_montgomery_mul`) warn about overflow for large moduli. Typestate would gate the fast-but-overflow-prone helpers behind proof that the modulus fits the regime. Recommendation in the same doc: favor `FullWidthR` for full-width moduli by default. |
| ○ `InfallibleCios<T>` / `CiosReady<T>` | CIOS limb indices are statically valid | New | Codex follow-up §14. CIOS currently returns `Option<T>` because `get_word(i)?` can fail in principle. In practice, indices come from `0..T::word_count()` and are statically valid. A typestate-backed infallible variant would remove the `Option` *and* the `expect("CIOS mul cannot fail...")` in `Field::mul`. Composes with the `n_prime_0` caching opportunity. |
| ○ `OneResidue<'f, T, P>` / `ZeroResidue<'f, T, P>` | identity values | New | Codex follow-up §15. Less central but useful: `mul_by_one` is a no-op, `mul_by_zero` short-circuits, sparse polynomial / NTT code can avoid Montgomery operations on identity terms. Probably wait until a consumer has a clear use case. |
| ○ Generative `SameModulus` brand | residues from one specific `Field` instance, not just any Field with the same `P` | New | Codex follow-up §17. The current lifetime brand allows mixing residues from two `Field` instances built in the same scope with the same personality. A generative brand (`PhantomData<fn(&'brand ()) -> &'brand ()>` + closure scope) closes that gap. The module docs already flag this as a known limitation. |

## 4. New ideas that aren't just typestate wrappers

Concrete optimizations and API shapes that don't require a typestate decision to land:

- **Cache `n_prime_0` inside `Field`** (codex §11, Priority 3). `Field::mul`'s CIOS path calls `n_prime.get_word(0)?` on every call. Storing the word at `Field::new` time removes the per-call extract and lets the trait expose an infallible variant — meaningful for the exponentiation hot loop. Composes with the future `InfallibleCios` typestate but doesn't require it.
- **`MontParams` as a separable, precomputed context** (codex §6, reinforced by the codex follow-up and gemini follow-up §1). The standalone `basic_montgomery_mod_mul`, `basic_montgomery_mod_exp`, and CT variants currently recompute `n_prime`, `r_mod_n`, `r2_mod_n` *every call*. Even without the wrapper question, splitting `Field`'s precompute out into a reusable struct closes that performance trap.
- **Dedicated `Field::square` / `Field::square_n<const N>`** (codex follow-up §7). Multiplication and squaring currently share the same CIOS path. A dedicated squaring entrypoint enables backend-specific optimization (fewer cross products, symmetric limb mults) and is especially relevant for Fermat-inverse and sqrt-by-exponentiation addition chains, where repeated squaring dominates.
- **Windowed exponentiation tables** (codex follow-up §6). The current `exp` is bit-by-bit square-and-multiply. Fixed-window / sliding-window with a precomputed `OddPowerTable<'f, T, P, const W>` is meaningfully faster for the same base. The const-generic table size is the typestate part, but the addition of the windowed exp method is the substantive change.
- **API consolidation via `Borrow<T>`** — see §5 below (separate section because it's structurally different from typestate wrappers).

## 5. API consolidation via `Borrow<T>` (gemini follow-up §2)

**Not a typestate.** A structural-refactoring proposal: collapse the three flavor modules (`basic_`, `constrained_`, `strict_`) into single generic functions using `core::borrow::Borrow<T>` to abstract over owned vs. borrowed operands.

```rust
pub fn mod_add<T, A, B, M>(a: A, b: B, m: M) -> T
where
    A: Borrow<T>, B: Borrow<T>, M: Borrow<T>,
    for<'a> &'a T: core::ops::Add<&'a T, Output = T>
                 + core::ops::Rem<&'a T, Output = T>,
{
    let sum = a.borrow() + b.borrow();
    sum % m.borrow()
}
```

The idea is that the compiler monomorphizes one definition into the three current forms — call with `u32, u32, u32` and it's the `basic` path, call with `&BigUint, &BigUint, &BigUint` and it's the `strict` path.

**Assessment:**

- **Pro:** real code reduction — three flavor modules collapse to one per operation. Maintenance burden drops significantly.
- **Pro:** consumers stop having to choose a flavor module up front. The choice is implicit in the argument types.
- **Con:** the bound `for<'a> &'a T: Add<&'a T, Output = T> + Rem<...>` is harder to read and document than the current three-way split. Backend authors may struggle to satisfy it cleanly for their types.
- **Con:** doesn't address the **orthogonal axis** of `_pr` vs. non-`_pr` (precondition-reduced vs. reducing-internally). That split remains regardless of how `Borrow<T>` collapses the flavor axis.
- **Con:** the existing split is documented as deliberate — `basic` (Copy) targets primitives, `strict` (ref) targets heap-allocated bigints, `constrained` (Clone, mixed) is the middle ground. Some backends impl `&T: Add<&T>` cleanly; some only impl `T: Add<T>` for `T: Copy`. The unified-bound approach may force backends to expand their trait implementations.
- **Composability with typestate:** if `Reduced<T>` lands and replaces `_pr`, the consolidated function takes `Reduced<T>` (or `&Reduced<T>`, or any `Borrow<Reduced<T>>`) operands. The two axes are independent and can land in either order.

**Recommendation:** evaluate after `Reduced<T>` (or `Canonical<T>`) is settled. Doing both at once compounds the redesign cost. A reasonable interim step is to prototype `mod_add` with `Borrow` bounds on a single operation, measure the bound complexity, and see if backends ship cleanly.

## 6. Prior modmath decisions on this branch (preserve the context)

**Abandoned: `Odd<T>` accessor on `Field`** (June 2026)

- Branch `feat/odd-modulus-accessor` at `335795d` (local-only, never pushed).
- Implemented a local `modmath::Odd<T>` newtype + `Field::modulus_odd() -> &Odd<T>` accessor.
- Rejected after rsa-heapless side argued it didn't actually close the downstream stack regression: the RSA crate's trait wants `&crypto_bigint::Odd<T>` specifically, and a modmath-local `Odd<T>` doesn't bridge to it without unsafe transmute in the consumer.
- **Correct fix was identified as upstream in crypto_bigint** — `Odd::from_ref(&T) -> Option<&Odd<T>>` or `derive(RefCast)`. Until that lands, modmath shipping its own `Odd<T>` doesn't help the downstream consumer.
- Branch is preserved as a starting point for revival if upstream cooperates.

**Deferred: `Reduced<T>` for the schoolbook `_pr` surface** (June 2026)

- Discussion concluded that the schoolbook `_pr` suffix could be replaced by `Reduced<T>` typestate, but:
  - krabipqc voted **unbranded** `Reduced<T>` (zero-sized type-level brand on their side, lifetime branding from modmath would be friction).
  - User chose to defer with "no Reduced<T> for now" — krabipqc immediate need was `mul_acc`, not Reduced.
- The design conclusion: when revived, `Reduced<T>` should be **unbranded** at the modmath level. Consumers add their own brand (curve, scheme, modulus) as a thin wrapper.

**Shipped on this branch:**

- v0.3.0: optional `zeroize` feature; `Residue` drops `Copy`, picks up `Drop + ZeroizeOnDrop`. `MontStorage` trait added as feature-aware bound on `Residue`'s mont field.
- v0.3.1: `wide_montgomery_mul_acc` family (free functions) + Field method sweep (`Field::mul_acc`, `Field::wide_redc`, `Field::inv_fermat`, `Field::inv_eea`, `ResidueCt::ct_eq`).

## 7. Recommended next-action prioritization

This consolidates codex's prioritization with our prior decisions and the
new perf-only opportunities.

**Tier 1 — additive, no breaking change, ship as 0.3.x**

1. **Cache `n_prime_0` in `Field`.** Hot-loop win, no API change. Internal struct field + infallible internal CIOS entrypoint. Sets up the future `InfallibleCios` typestate without requiring it. *(codex §11, follow-up §13)*
2. **`Field::square` / `Field::square_n<const N>` entrypoints.** Currently squaring goes through CIOS multiply with `a == b`. Dedicated entrypoint enables backend-specific specialization (symmetric limb products) and is on the critical path for Fermat-inverse and sqrt addition chains. Purely additive. *(codex follow-up §7)*
3. **`Wide<T>` newtype.** Cosmetic + safety; removes `(hi, lo)` swap-bug class. Either change `WideMul::wide_mul` to return `Wide<T>` (breaking) or add a `Wide<T>` view type that wraps `(T, T)` losslessly (additive). *(codex §4, gemini_1 §4)*
4. **Wait on `Odd<T>` until crypto_bigint upstream cooperates.** No work here.

**Tier 2 — additive but larger surface, 0.4.0 candidate**

5. **`MontParams<T>` as separable precomputed context.** Refactor `Field` to wrap a `MontParams` internally; expose `MontParams::new()` and standalone-function entrypoints that take `&MontParams`. Eliminates the per-call recompute in `basic_montgomery_mod_*` functions — the most consensus-flagged perf trap in the codebase. *(codex §6, follow-up §1; gemini follow-up §1)*
6. **`Canonical<'m, T>` (unbranded) + retiring `_pr` suffix.** Replace `basic_mod_add_pr` family with `basic_mod_add` taking `Canonical<T>` operands and returning `Canonical<T>`. Breaking on the schoolbook surface. Decision deferred from June. *(Naming note: codex follow-up §8 argues `Canonical` is more accurate than `Reduced` because the `_pr` APIs require `[0, m)`, not just congruence.)*
7. **`NonZero<T>` for moduli.** Mirrors stdlib `NonZeroUxx`. Removes the `if modulus == zero { return None }` top-of-function check across Montgomery entrypoints. Composes with `MontParams<T>` (which wraps an `Odd<NonZero<T>>` internally). *(gemini follow-up §3, codex follow-up §16)*

**Tier 3 — design-significant, 0.4.0 or later**

8. **`MontForm<T, P>` at the free-function level.** Pushes the brand currently on `Residue` down to `basic::montgomery::wide::*`. Removes the raw-vs-Mont silent-correctness risk. *(codex §3, gemini_1 §3)*
9. **`Public<T>` / `Secret<T>` (or `PublicExp` / `SecretExp`).** Make `Field::exp` vs `Field::exp_public_exp` compiler-checked rather than naming-convention. The generalized `Secret<T>` form can also gate variable-time `_pr` entrypoints. *(codex §8, follow-up §4; gemini follow-up §5)*
10. **`OddPrime<T>` + `PrimeField<T, P>`.** Adds infallible inverse / Fermat path for prime moduli. Construction is `unsafe new_unchecked` for static constants; runtime primality testing is impractical at crypto sizes. *(gemini_2; codex §13, follow-up §1)*
11. **`MontParams<T, RKind>` algorithm-family split.** Const-or-marker parameter distinguishing `R > N` from `R = 2^W` Montgomery. Prevents accidentally mixing helpers from one family with parameters from the other. *(codex follow-up §11)*
12. **API consolidation via `Borrow<T>`.** Collapse `basic`/`constrained`/`strict` flavor modules. Evaluate after `Canonical<T>` (item 6) is settled — the two changes compound. *(gemini follow-up §2; see §5 for assessment)*

**Tier 4 — speculative / context-dependent / specialization layers**

13. **`Loose<'m, T, K>` / `WideLoose<'m, T, K>`** for lazy reduction in addition-heavy curve code. Requires the backing type to have headroom; only useful in specialization layers (e.g. Curve25519's `Field<Loose>`). *(codex §10, follow-up §9)*
14. **`Unit<'m, T>` / `NonZeroResidue<'m, T>` / `Coprime<A, M>`** for infallible inverse on known-coprime values. The three are related; pick one consistent shape. *(codex §7, follow-up §2; gemini follow-up §4)*
15. **`SpecialModulus<T>` family** — `PseudoMersenne<const C>`, `Solinas`, `Crandall`, `Mersenne`, `Goldilocks`. Belongs in downstream specialization crates (Curve25519's `2^255-19`, Ed448, NTT-friendly primes); modmath's role is to provide the shared vocabulary so multiple consumers compose. *(codex follow-up §10)*
16. **`OverflowSafe<T>` / `SmallModulus<T>`** to gate `R > N` helpers behind their non-overflow regime. *(codex follow-up §12)*
17. **Generative `SameModulus` brand** via `PhantomData<fn(&'brand ()) -> &'brand ()>` + closure scope. Closes the documented covariance gap that lets residues from two same-`P` Fields mix. *(codex follow-up §17)*
18. **`OneResidue` / `ZeroResidue` identity typestates** — wait until a consumer needs them. *(codex follow-up §15)*
19. **`NonTrivialModulus<T>` / `ModulusGtOne<T>`** — likely absorbed into `OddModulus<T>` if that constructor rejects 1 (Montgomery modulo 1 is useless). Not a separate wrapper. *(codex follow-up §3)*

## 8. Open architectural questions

- **Backend specialization:** confirmed irrelevant to typestate wrappers themselves. Lives on the underlying trait (`Parity`, `PartialOrd`, `Rem`). The wrapper is a static guarantee, not a behavior.
- **Wrapper hosting:** start in modmath (Option A). Extract to a dedicated crate (Option C) only if a second consumer crate explicitly wants the surface and modmath isn't a natural umbrella.
- **Branding strategy for `Reduced<T>` / `Canonical<T>` / `MontForm<T>`:** unbranded at the modmath level; consumers add their own brand. Settled by the June krabipqc discussion.
- **Construction trust model for `Prime<T>` / `OddPrime<T>`:** `unsafe new_unchecked` for cryptographic constants; optionally a `new_checked` deterministic constructor for small primitive integers and a Miller-Rabin variant for dynamic sizes.
- **`Canonical<T>` vs `Reduced<T>` naming:** codex follow-up argues `Canonical` is more accurate because the `_pr` APIs require `[0, m)` (canonical representative), not just modular equivalence. Probably resolve in favor of `Canonical` and reserve `Congruent` for mere equivalence if it's ever needed.
- **Algorithm-family split (`R > N` vs `R = 2^W`):** currently distinguished only by function naming. Adding a type parameter on `MontParams` is the cleanest fix, but it adds a generic everywhere `Field` flows. Defer until `MontParams<T>` itself lands so the type parameter has a natural home.
- **`Borrow<T>` consolidation vs. `Canonical<T>` ordering:** both change the schoolbook surface significantly. Doing both at once compounds the redesign cost. Pick one to land first.
