# Review of the synthesis doc — reviewer 1

Reviews `const_num_traits_typestate_synthesis.md` (not the original sketch).
Saved verbatim for the final-review pass. Captured 2026-06-16.

---

The synthesis document is an exceptionally rigorous and mature distillation of the typestate design space. By establishing a strict "bar" for inclusion—specifically requiring a consuming op at the correct layer—it prevents const-num-traits from bloating into a typestate zoo and keeps it strictly focused on core-mirroring structural algebra.

Here is the final assessment of the synthesis and the proposed sequence:

## 1. The "Bar" and Scope Discipline
The two filters (Consuming Op + Layer) are the most valuable architectural contributions in the document. A typestate only justifies its existence in this foundational crate if it generically removes work (e.g., eliding branches, options, or panics in mathematical operations). Pushing domain-specific proofs like Prime, Odd (whose consumer is Montgomery reduction), and Normalized up the stack to modmath or a bignum crate is the correct structural decision.

## 2. Defect Resolution Assessment
The identified defects (D1-D9) are accurate and their proposed fixes are solid:
- **D1 (Parity for &T gap):** Spot on. Adding the blanket `impl<'a, T: Parity + Copy> Parity for &'a T` is a necessary prerequisite for any `from_ref` zero-cost proof constructor.
- **D2 & D3 (Output types and constraints):** Introducing `type Output` for the consuming operations (DivNonZero) is critical to support non-Copy bignum backends without violating the crate's operational conventions. Constraining the internal wrapper (e.g., `type NonZero: Into<Self>`) ensures generic code can actually unwrap the proof.
- **D4 (MIN / -1 signed overflow):** Restricting DivNonZero to unsigned types initially is the safest and cleanest path. Signed total division requires a composite proof (NonZero + NonMin), which is too complex for an initial typestate rollout.
- **D8 (Constant-Time Constructors):** Providing CtOption returning constructors gated behind the ct feature is a hard requirement given the crate's commitment to Tier A/B constant-time operations.

## 3. The Recommended Sequence
The sequence is prioritized correctly based on the generic payoff:
1. **PowerOfTwo<T>:** This is undeniably the strongest candidate. It turns Tier-C operations (division/modulo) into Tier-A operations (shift/mask) for generic downstream consumers. All the underlying predicates already exist, making it a low-risk, high-reward first target.
2. **NonZero Bridge:** Renaming to HasNonZero (D9) and focusing on unsigned DivNonZero/RemNonZero provides the next best algebraic payoff by removing Option from the generic division API.
3. **Positive<T> / NonNegative<T>:** Useful for branch-free abs and total isqrt, successfully reusing the existing sign.rs and convert.rs infrastructure.

## Final Conclusion
The document provides a clear, actionable roadmap. By defining interop as downstream implementation rather than an upstream design driver, the synthesis ensures the API remains perfectly suited for generic core emulation. You are cleared to proceed with implementing PowerOfTwo<T> as the trailblazer, followed by the HasNonZero bridge.
</content>
