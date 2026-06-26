# Review of the synthesis doc — reviewer 3

Reviews `const_num_traits_typestate_synthesis.md` (not the original sketch).
Saved verbatim for the final-review pass. Captured 2026-06-16.

---

Read it. This is a genuinely strong synthesis — sharper than my original review in two places, and the scope discipline is the real win. My final commentary, focused on what's still worth deciding rather than what's already settled.

## What the synthesis got right that I didn't

**D4 (signed `MIN / -1`) is a catch I missed.** I framed `DivNonZero`'s payoff as "infallible division"; that's only true for unsigned. A nonzero signed divisor can still be `-1`, and `MIN / -1` overflows. The unsigned-total / signed-deferred split is the correct resolution, and "unsigned-only to start" is the right call. Good correction.

**The bar in §2 is the load-bearing contribution.** Consuming-op + right-layer + conventions, applied as a hard filter that shrinks the list to three families. That triad is what turns a pile of "wouldn't this be neat" wrappers into a defensible API surface, and the structural-low / domain-high cut (Prime→modmath, Normalized→bignum, Odd's ops up but its definition maybe low) is exactly the placement logic this crate needs. The "shrinkage is the recommendation" line should survive into whatever ships.

## Three things still undecided that I'd settle before drafting

1. **D3's `type NonZero: Into<Self>` collides with the crate's whole reason to exist — const.** "`From<NonZeroU32> for u32` exists" is true, but existence ≠ const-callable. If `DivNonZero` is a c0nst trait and you bound `type NonZero: [const] Into<Self>`, you need core to provide `impl const From<NonZero<T>> for T` — which I don't believe is a const trait impl yet (only `NonZero::get` is const fn). So the recommended bound risks making any generic const consumer non-const, defeating the point. Safer: leave `type NonZero` unconstrained but reach the value through a crate-owned const fn accessor (or a sealed companion with a const fn get), not a core From/Into bound. This is the one D-item I'd revise, not just implement.

2. **`PowerOfTwo<T>` representation is an unmade decision that determines the entire payoff.** The synthesis says "ship" and lists mask-mod / shift-div / align-up, but never says what the wrapper stores. Three options, materially different:
   - store the value `p` → every rem/div/align recomputes `p-1` and `p.ilog2()` (cheap bsr/dec, but not free, and ilog2 per-div);
   - store the mask `p-1` → mask-mod and is-multiple and align-up are immediate; div still needs the shift count;
   - store the exponent `k` → div/mul are `>> k`/`<< k` directly, mask is `(1<<k)-1`.

   Given the north star ("most beneficial optimized implementations downstream"), I'd lean toward storing the exponent (or exponent+mask) so the consuming ops are pure shifts/masks with nothing recomputed. This should be decided before Reviewer-3-drafts-it, because it changes the constructor, the `get()` contract, and whether `Deref<Target=T>` even makes sense (deref to the value, but store the exponent — mismatch).

3. **The consuming ops are new methods that take the proof as an argument — say so explicitly.** You can't specialize the existing blanket `MultipleOf`/`NextMultipleOf` when the divisor happens to be a power of two (no specialization on stable). So these are new surface — `fn rem_pow2(self, p: PowerOfTwo<Self>) -> …`, etc., likely a small `PowerOfTwoOps` trait or inherent methods on the wrapper. The synthesis lists the ops but is silent on their home; pin it so it doesn't get reinvented as "make MultipleOf faster," which won't work.

## One addition the synthesis is missing: a proof lattice

The three families aren't independent — they refine. A `PowerOfTwo<T>` is always nonzero; a `NonNegative<i32>` that's also nonzero is `Positive`. That's free invariant-narrowing, and it directly strengthens the consuming-op story: dividing by a `PowerOfTwo` is infallible because it narrows to `NonZero` without a recheck. So a cheap const fn refinement path (`PowerOfTwo<T>` → `NonZero<T>`, `Positive` → `NonNegative` → unsigned cast) lets the stronger proofs reuse the weaker ones' ops instead of duplicating them. This is also a partial answer to D4: if your signed divisor carries a "not -1" fact, `MIN / -1` is gone by construction. Worth one section — it's the thing that makes "proof trifecta" a system rather than three unrelated newtypes.

## Minor

- "Ship (1st)" appears on both NonZero (§4) and PowerOfTwo (§4 "concurrent / first"), and §7 sequences PowerOfTwo first. Harmless, but pick one label — I agree with §7: PowerOfTwo first (lowest risk, all predicates exist, no signed/const-Into complications).
- §5 covers it via the general rule, but worth stating outright: `PowerOfTwo::new` is Tier B (popcount + bool, per pow2.rs's own header), so it needs the same ct-feature CtOption constructor as NonZero/Odd — it's not a "branchless predicate so it's fine" case.

## Bottom line

The synthesis is sound and I'd green-light the §7 sequence essentially as written, with one substantive revision (D3 — don't bound on core's non-const Into; use a crate-owned const accessor) and two decisions made before drafting (PowerOfTwo representation; consuming-op home). Add the refinement-lattice section and it goes from "three good typestates" to a coherent proof system, which is what the north star is actually asking for.
</content>
