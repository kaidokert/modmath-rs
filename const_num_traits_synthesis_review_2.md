# Review of the synthesis doc — reviewer 2

Reviews `const_num_traits_typestate_synthesis.md` (not the original sketch).
Saved verbatim for the final-review pass. Captured 2026-06-16.

---

The synthesis is materially sharper than the original sketch. I'm going to call out only the parts I'd keep, amend, or push back on, because the document already captures the broad review well.

I'd treat the synthesis as a solid direction, with a few amendments.

**The strongest conclusions are right:**

- NonZero belongs low, but `DivNonZero` should start unsigned-only if it promises totality. Signed `MIN / -1` makes "infallible signed division" false.
- Odd is useful, but its main payoff is modular arithmetic, so the consuming ops belong in modmath/bignum layers. A bare proof wrapper can live low only if you really want shared vocabulary.
- PowerOfTwo and sign-state wrappers are better crate-level typestate candidates than Even.
- Fixing `Parity for &T` is necessary if `from_ref` is part of the API.

**The parts I'd adjust:**

1. **I would not make PowerOfTwo the first implementation by default.**
   It has great payoff, but its consuming ops need more API design than NonZero. NonZero has a clearer bridge story to core, obvious primitive aliases, and obvious div/rem consumers. I'd implement NonZero first unless you already have downstream PowerOfTwo call sites waiting.

2. **Be careful with `type NonZero: Into<Self>`.**
   This is probably fine for primitive `NonZeroU*`/`NonZeroI*`, but it bakes in ownership conversion. For bignums, `from_nonzero(nz)` may be the more flexible canonical path. I'd prefer the trait's own `from_nonzero`/`as_inner` contract over requiring `Into<Self>` unless a generic caller actually needs `nz.into()`.

3. **"Add `impl<T: Parity + Copy> Parity for &T`" is only a partial rescue.**
   It fixes primitive/reference ergonomics, but it does not make `from_ref` work for non-Copy bignums unless those backends implement `Parity` for `&Big`. That is acceptable, but the doc should say so explicitly: `from_ref` is zero-copy for types that provide borrowed parity.

4. **CT constructors should be optional, not a hard blocker for the first stable typestate slice.**
   The note is correct for crypto consumers, but requiring every typestate constructor to ship with CtOption on day one may slow down the core API. I'd document plain constructors as Tier B and add CT constructors where the `ct` feature already has matching primitives.

5. **Positive/NonNegative need signed-only semantics nailed down before naming.**
   Positive usually means `> 0`; NonNegative means `>= 0`. For unsigned types these are either redundant or subtly different for zero. I'd keep them signed-only unless there is a strong generic story.

**My final recommended order:**

1. `typestate` feature.
2. `NonZeroBridge`/`HasNonZero` without a `Zero` supertrait, using stable concrete primitive `NonZero*` aliases.
3. Unsigned `DivNonZero` / `RemNonZero` with `type Output`.
4. Add `Parity for &T` helper impl and optionally bare `Odd<T>` only if downstream needs shared proof vocabulary.
5. Design `PowerOfTwo<T>` separately with concrete consuming ops before committing names.
6. Add signed `NonNegative`/`Positive` after the first slice lands.

Net: the synthesis is directionally correct and appropriately cuts the typestate zoo. I'd just avoid letting PowerOfTwo jump ahead of the simpler NonZero bridge, and I'd keep CT and borrowed-reference support honest rather than making them implicit promises.
</content>
