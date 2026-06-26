# const-num-traits typestate sketch — review 3 of 3

Reviews the `MaybeNonZero` / `Odd<T>` / `Even<T>` sketch proposed for the
const-num-traits fork. Saved verbatim for later synthesis. Captured 2026-06-16.

---

I have enough to review against the crate's actual conventions. Here's my assessment.

## Part 1 — Review of the proposed typestate sketch

### MaybeNonZero / DivNonZero (the bridge)

The bridge-not-duplicate instinct is right, and routing primitives to `core::num::NonZero` is the correct call — duplicating NonZero would fragment the ecosystem and lose the niche. But several things don't hold up:

**`type NonZero` is unconstrained, so the trait is nearly inert generically.** There's no bound tying `Self::NonZero` back to `Self`. A generic caller with `T: MaybeNonZero` can call `into_nonzero`/`from_nonzero` but can do nothing else with the result — no `.get()`, no arithmetic, no `Into<Self>`. For this to pull weight you need something like `type NonZero: Into<Self>` or a small sealed `NonZeroRepr<Base = Self>` companion trait. As written, only `from_nonzero` recovers the value, which makes the associated type opaque dead weight.

**`MaybeNonZero: Zero` drags in Add.** In `identities.rs:16`, `Zero: Add<Self, Output = Self>`. So requiring Zero couples a nonzero typestate bridge to an arithmetic supertrait — unwanted for CT-personality types or pure markers that want the predicate without committing to Add. You only need "what counts as zero," and there's no standalone `IsZero` atom (it lives on Zero; CT has `CtIsZero`). Worth either extracting an `IsZero` atom or just not requiring Zero and testing zero-ness inside `into_nonzero` via the inherent path.

**`DivNonZero` violates the crate's associated-Output convention.** Per CLAUDE.md's by-value sweep: every owned-Self-returning method needs an Output so non-Copy types can `impl … for &T`. `div_nonzero(self, …) -> Self` and `rem_nonzero(self, …) -> Self` return bare `Self` — a non-Copy bignum can't implement them via the `&T` rescue. And you can't reuse `Div<Self>::Output` here because the divisor is `Self::NonZero`, not `Self`. So this needs a fresh `type Output`. Flag it; otherwise it's inconsistent with every other arith trait in the crate.

**The primitive payoff is mostly illusory; the real payoff is backends.** `x / nz.get()` on a primitive already elides the divide-by-zero panic — NonZero carries an assume-style niche and LLVM knows the divisor is nonzero. So `DivNonZero` for the 12 primitives is pure API sugar, buying no codegen. The genuine win is (a) the trait abstraction for bignum backends that hand-write the no-branch long-division path, and (b) generic code. The doc comment ("infallible division — no Option, no divide-by-zero branch") oversells the primitive case; reframe it as "the abstraction backends implement; primitives inherit core's niche."

Minor: `MaybeNonZero` reads like an Option-shaped thing; it's a bridge — `HasNonZero`/`NonZeroBridge` is clearer. Const-wise this is fine (`NonZero::new`/`get` are const).

### Odd<T> / Even<T>

The `crypto_bigint::Odd` analogy is exactly the right prior art, and `from_ref` zero-cost reinterpret is the form downstream wants. But:

**`from_ref`'s bound is currently unsatisfiable for primitives — this is a real bug, not a nit.** `from_ref` requires `for<'a> &'a T: Parity`. I checked `parity.rs`: it impls Parity only for the 12 primitives (`parity_impl!(usize u8 …)`), with no `impl Parity for &T`. So `&u32: Parity` does not hold, and `Odd::<u32>::from_ref` would not compile — yet the comment calls `from_ref` "the form downstream actually wants." CLAUDE.md describes the `impl Trait for &T` rescue as the general non-Copy mechanism, but Parity never actually got that impl. You need either (a) add a blanket `impl<T: Parity + Copy> Parity for &T` to `parity.rs`, or (b) give `from_ref` the `T: Parity + Copy` bound and deref internally (`if (*value).is_odd()`). The two constructors then serve Copy (primitives) vs. genuinely-`&`-impl'd (non-Copy) types — both legitimate, but the sketch only ships the one that excludes primitives.

**The proof has no consumer — so it's just a newtype.** This is the structural weakness. NonZero earns its keep because DivNonZero consumes the invariant to drop a branch. `Odd<T>`/`Even<T>` as sketched expose only `new`/`get`/`Deref` — nothing requires oddness or exploits it. A typestate without a consuming op is dead weight. Odd should gate something real: Montgomery setup, odd-modulus inverse existence, or at minimum the parity-algebra facts. `Even<T>` is even weaker — its natural payoff is exact halving (`>> 1` with no remainder); pair it with a `HalveExact`/`DivExact`-by-2 or drop it. Ship each typestate with ≥1 op that consumes the proof, or don't ship it.

**No niche, unlike NonZero.** `Option<Odd<u32>>` is 8 bytes, not 4 — the plain newtype has no forbidden bit pattern. (crypto_bigint's Odd is the same, so this is acceptable, but the sketch shouldn't be mentally filed next to NonZero's layout win.)

**CT interaction is unaddressed.** `Odd::new` branches on `is_odd`, which CLAUDE.md marks Tier B (the parity bool is the operand's low bit). For a CT-personality bigint, a branching `Odd::new` leaks. The honest move is a CtOption-returning constructor (mirroring the `ct.rs` pattern) under the `ct` feature, and a doc note that the plain constructors are not constant-time.

**Const consistency:** `new`/`from_ref`/`get`/`Deref` should go through `c0nst!`/`pub c0nst fn` like the rest of the crate (the sketch acknowledges bodies are plain-Rust for readability — just make sure the real impl does this; `new_unchecked` already is).

## Part 2 — Other typestate opportunities in this crate

The crate already has the predicates for two strong typestates; the meta-rule is the one the sketch half-misses: **a typestate is only worth a feature when a consuming op exploits the invariant to delete a branch, an Option, or a division.** By that test:

### Strongest: PowerOfTwo<T> (builds on the existing IsPowerOfTwo)

This is the best additional candidate by a wide margin, because — unlike Even — it has a rich set of ops that genuinely get cheaper, and every piece already exists in the crate:

- `x % p` → `x & (p − 1)` (mask, no division — `ops/rounding.rs`'s MultipleOf becomes branchless)
- `x / p` → `x >> p.ilog2()` (ilog2 is already in `ops/log.rs`)
- `next_multiple_of(p)` → `(x + p − 1) & !(p − 1)` (alignment, branch-free — feeds NextMultipleOf)
- `is_multiple_of(p)` → `x & (p − 1) == 0`

It pairs directly with MultipleOf/NextMultipleOf/DivCeil, has massive real-world use (allocators, alignment, ring buffers, hashing), and the proof actually consumes into faster code. Tier A/B (mask + shift). Recommend this first.

### Strong second: sign typestates NonNegative<T> / Positive<T> (build on Signed)

`sign.rs` already gives `is_positive`/`is_negative`, and `convert.rs` has `UnsignedAbs`/`CastUnsigned`. The consuming payoffs are concrete:

- infallible signed→unsigned cast: `NonNegative<i32>` → `u32` with no try_from/Option (today's CastUnsigned is lossy-or-checked).
- branch-free abs: identity on NonNegative (drops the `if self.is_negative()` in `sign.rs:78`).
- always-defined isqrt: CheckedIsqrt (signed-only, returns Option) collapses to total on NonNegative.

Together with NonZero (division) and Odd (Montgomery), these three are exactly the proof-trifecta a fixed-bigint / modmath backend wants.

### Weaker / mention-only

- **NonMax<T>/NonMin<T>** — niche-layout siblings of NonZero (core has NonMax unstable). Same Option size win, but the numeric op story is thin (mainly +1/-1 can't-overflow). Low priority.
- **Normalized<T>** (high bit set) — real in bignum long division (Knuth normalization) and float mantissas, and the crate already does wide-multiply work in `ops/bigint.rs`. But it's backend-specific enough that it probably belongs in the bignum crate, not here.

## Bottom line

- **Fix before merge:** the `from_ref`/`Parity for &T` gap (compile bug), `DivNonZero`'s missing Output, and the unconstrained `type NonZero`.
- **Reconsider:** Even and Odd without a consuming op — give them one or cut them. `MaybeNonZero: Zero` coupling.
- **Add instead:** `PowerOfTwo<T>` (top), then NonNegative/Positive. Both reuse predicates already in the crate and have real consuming ops, which is the bar the sketch sets with NonZero→DivNonZero but doesn't meet for Odd/Even.

Want me to sketch a corrected `PowerOfTwo<T>` + its consuming ops (mask-mod, shift-div, align-up) in the crate's `c0nst!` / by-value / associated-Output style so you can compare it against the proposal?
</content>
