# Roundup — typestate, CIOS placement, and the panic-free north star

A consolidating pass over the whole arc: the const-num-traits typestate sketch,
the CIOS/`MulAccOps` trait-placement exploration, and the two
`PANIC_FREE_REQUESTS.md` shopping lists (fixed-bigint + modmath). Supersedes the
standing conclusions in `SYNTHESIS_CLAUDE.md` and `ANALYSIS_CLAUDE.md` where they
conflict. Captured 2026-06-17.

This document exists because the panic-free requirements retroactively sharpen
the design driver behind *all* of the prior work, and in doing so reverse two
conclusions I argued earlier. It records the corrected state.

---

## 0. The real north star (now measurable)

Earlier docs measured a typestate by "how much downstream optimization it
unlocks." The `PANIC_FREE_REQUESTS.md` docs replace that vague target with a
concrete, falsifiable gate:

> **Zero `panic_fmt` / `panic_bounds_check` / `rust_begin_unwind` symbols
> surviving linker DCE on embedded targets (AVR / Cortex-M / RISC-V), verified
> by `cargo nm <binary> | egrep "panic|unwind"` → empty.**

The downstream consumer is `ed25519_heapless`. Three in-crate strategies are in
play (typestate, const generics, explicit `Result`/`?`), with an explicit ban on
trading panics for UB (`unreachable_unchecked`). A handful of sites need upstream
API-shape changes — those are the two shopping lists.

Everything below is re-derived from this gate. The measure of a feature is no
longer "is it faster" but "does it delete a panic symbol (or move a check to
compile time) without introducing `unsafe`."

## 1. Two corrections to my earlier reviews

The panic-free framing forces me to walk back two positions, both in the same
direction — **fallibility-removal is itself the payoff, because the metric is
panic symbols, not cycles.**

### 1a. `Odd<T>` is the flagship typestate, not a dead newtype

My first review said Odd/Even "have no consuming op; the consumer is Montgomery,
which is domain" and recommended placing Odd's definition low but its ops high,
treating it as weak. **Wrong.** `modmath/PANIC_FREE_REQUESTS.md` Ask A shows the
consuming op precisely: `Field::from_odd_modulus(Odd<T>)` **deletes the
`Field::new(p).unwrap()` → `panic_fmt` symbol** from the AVR binary. That is the
literal, named, verifiable goal — not a codegen nicety. `Odd<T>` is the typestate
with the most concrete consumer in the entire set.

**Simplification:** the doc asks for `OddNonzero<T>`, but `Field::new` rejects
"even *or* zero," and **zero is even**, so the precondition collapses to "odd."
Odd ⇒ nonzero. Montgomery needs the modulus coprime to the radix `2^k`, i.e.
odd; that single bit is the whole precondition. **Use `Odd<T>`, not
`OddNonzero<T>`** — the "Nonzero" half is implied and the wrapper is redundant.

### 1b. The primitive `DivNonZero`/infallible-ctor payoff is real, not "≈0"

I earlier flagged (synthesis D7) that the primitive `DivNonZero` payoff is
illusory because the `NonZero` niche already elides the divide-by-zero panic.
Under the panic gate that is backwards: **even when codegen is byte-identical,
removing the `Option` return / the fallible branch removes a panic symbol.**
API-level fallibility-removal *is* the win. So `NonZero`/`Odd`/infallible
constructors rank at the top, not as "API sugar."

## 2. Typestate set — re-ranked under the panic gate

| Proof | Consuming op | Payoff under the gate | Rank |
|---|---|---|---|
| **`Odd<T>`** | `Field::from_odd_modulus` (modmath) | kills `Field::new().unwrap()` `panic_fmt` | **top** |
| **`NonZero` bridge + `DivNonZero`** | infallible division | removes `Option` + div-by-zero panic site | **top** |
| **`PowerOfTwo<T>`** | mask-mod / shift-div / align-up | branchless, no div-by-zero, no overflow panic | strong |
| `NonNegative`/`Positive` | infallible→unsigned cast, total `isqrt` | removes `try_from`/`Option` panic paths | second |
| `Even<T>` | exact-halving only | no named consumer in the panic work | **cut** |

The constructor-variant priority flips from the synthesis:

- **`const` constructor first.** `Odd::new_const(t)` that `panic!`s in *const
  context* = a compile error at the call site, no runtime panic path
  (`const { assert!(…) }` shape). Curve25519's modulus is a public compile-time
  constant, so this is the exact fit.
- **`CtOption` constructor second.** Only needed for *secret* moduli; rarer.
  Keep it `ct`-feature-gated, but it is no longer the headline.

This is the cleanest convergence of the const-trait and typestate threads:
`Odd::new_const` requires **const `Parity` + const zero-check on the backend
type**, which is exactly what const-num-traits' const-trait machinery exists to
provide. The typestate proofs are the consumer that *justifies* const `Parity`.

**Placement (unchanged, now receipted):** `Odd<T>` lives low — in
const-num-traits' `typestate` feature, built on cnt's own const `Parity` — and
passes the charter test because it is a proof over a *value-level predicate*
needing nothing about composite storage. `from_odd_modulus` lives high in
modmath. Structural-low / domain-high, with a panic-symbol receipt.

## 3. CIOS / `MulAccOps` placement — unchanged

The panic-free gate is a property of trait **bodies**, not trait **homes**
(fixed-bigint's doc is explicit: even with the const-assert, the *body* must
avoid indexed `out[i]` writes or it synthesizes `panic_bounds_check`). So none of
the placement conclusions move:

- **Reject limb access in cnt.** The charter test (a cnt trait must mirror a
  `core` inherent method on a primitive scalar) is untouched by the panic gate —
  limb access still has no `core` referent and rides with the row ops, not into
  cnt. The "is it a representation primitive like `ToBytes`?" bifurcation
  resolves *against* cnt: `ToBytes` mirrors `core` and is value-canonical;
  limb access does neither.
- **Split + rename** `MulAccOps` → `CiosRowOps` (Option G): unanimous, stands.
- **Algorithm body stays in modmath**; reject moving it to fixed-bigint.
- **Decouple now via Option H** (modmath owns the tiny `CiosRowOps` + algorithm +
  blanket impl; fixed-bigint impls it under an optional `modmath` feature),
  **extract to `bigint-traits` later (Option C)** when a second backend or
  Montgomery variant (SOS/FIOS/NTT) justifies the crate.

**One ripple:** the **discriminant probe** gets a second, independent vote.
An infallible in-bounds `word(&self, i) -> Self::Word` is now both *CT-cleaner*
and *panic-free-cleaner* — a fallible `Option`/`CtOption` accessor invites
caller-side `unwrap`/match that synthesizes panic sites. So: verify the CT
threat model, and if public loop indices make the infallible accessor safe,
collapse the `WordAccess`/`WordAccessCt` split and merge the two `cios_*`
algorithm bodies into one. Strengthens an existing recommendation; changes no
direction.

## 4. Cross-cutting convention (new, applies to both threads)

Add to the design conventions for *every* new trait method body, independent of
placement:

> **Bodies must be DCE-friendly panic-free.** Iterate via
> `iter_mut().enumerate().take(n)`, not indexed `out[i]`. Use `const { assert!(…) }`
> for size preconditions (compile-time failure, no runtime path). No internal
> `unwrap`/`expect`. No `unsafe { unreachable_unchecked() }` — the gate forbids
> trading panics for UB.

This is a *body-shape* rule, which is exactly why it perturbs no home decision.
It applies to cnt's typestate constructors and to the CIOS row-op bodies
(`mul_acc_row` over a limb array is a prime `panic_bounds_check` source) alike.

The fixed-bigint asks are the concrete instances:
- `to_le_bytes_fixed<const M>(&self, &mut [u8; M])` with
  `const { assert!(M >= N * size_of::<T>()) }` — replaces the `Result` +
  `.unwrap()` chain with a compile-time size check. (Ask A; stable today.)
- `from_le_bytes_fixed<const M>` — symmetric. (Ask C.)
- Owned `to_le_byte_array()` is the long game (needs `generic_const_exprs`); the
  per-canonical-size inherent methods (`u8×32`, `u32×16`, `u64×8`, `u64×4`) are
  the stable-today fallback.

## 5. Where each piece lives — the consolidated picture

| Layer / item | Home | Why |
|---|---|---|
| Word arithmetic (`CarryingMul`/`CarryingAdd`/`WideningMul`) | **const-num-traits** | mirrors `core` inherent methods on scalars |
| `Parity` (+ const) , `IsPowerOfTwo`, zero-check | **const-num-traits** | value-level predicates with `core` referents |
| `Odd<T>`, `NonZero` bridge, `PowerOfTwo<T>` proofs | **cnt `typestate` feature** | proofs over cnt's own predicates; `const` constructors |
| `DivNonZero`, `div_pow2`, infallible casts | **cnt `typestate` feature** | consuming ops that delete `Option`/panic sites |
| Limb/word access (`WordAccess`[`Ct`] or infallible `word`) | **with the row ops** (modmath now → `bigint-traits` later) | composite-storage hook; no `core` referent |
| `CiosRowOps` (`mul_acc_row`/`mul_acc_shift_row`) | **modmath now → `bigint-traits` later** | CIOS-specific; Option H then C |
| CIOS algorithm body, `CiosMontMul`/`CiosMontMulCt` | **modmath** | the algorithm; never leaves |
| `Field::from_odd_modulus(Odd<T>)` | **modmath** | domain consuming op; the panic-symbol kill |
| `to_le_bytes_fixed` / `from_le_bytes_fixed` (panic-free serde) | **fixed-bigint** | backend-owned representation |

## 6. The maintainer's remaining calls

Everything above is settled by the analysis except three genuine judgment calls:

1. **Decouple modmath↔fixed-bigint now, or defer?** Now → Option H (accept the
   dev-dep cycle). Defer-until-second-backend → Codex's lower-risk in-place
   Stage 1. The panic-free work is orthogonal to this; pick on the decoupling
   timeline alone. (Synthesis recommends H, since decoupling is the stated goal.)
2. **CT threat model for `word(&self, i)`.** Confirm whether public loop indices
   make the infallible accessor CT-safe. If yes, collapse the access split and
   the two algorithm bodies. This is a security-model call only the maintainer
   can sign off.
3. **`generic_const_exprs` appetite.** Determines whether `to_le_byte_array()` /
   `OddNonzero`-style owned-sized returns are on the table, or whether the
   stable-today const-generic-`M` + `const assert` shape is the ceiling for now.

## 7. One-line summary

The panic-free gate (`cargo nm` empty of panic/unwind symbols) is the real driver
behind all of it: it promotes **`Odd<T>` with a `const` constructor** to the
flagship typestate (it deletes `Field::new().unwrap()`), confirms
**fallibility-removal as the payoff** (so primitive `DivNonZero` is *not*
illusory), adds a **DCE-friendly-panic-free body** convention that is orthogonal
to placement — and therefore leaves the **CIOS placement verdict
(limb-ops out of cnt; `CiosRowOps` modmath-now / `bigint-traits`-later; algorithm
stays in modmath) unchanged**, while giving the infallible-`word` discriminant
probe a second reason to win.
