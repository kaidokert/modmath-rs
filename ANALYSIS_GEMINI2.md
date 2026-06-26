# Analysis: CIOS Trait Placement & Refactoring

This document provides a technical recommendation for resolving the coupling between \`modmath-rs\` and \`fixed-bigint-rs\` regarding the \`MulAccOps\` trait, as explored in \`cios_trait_placement.md\`.

## 1. Recommended Strategy: The "Clean Primitives" Approach

To achieve the goal of decoupling while improving the architectural integrity of the ecosystem, I recommend a combination of **Option G (Trait Splitting)** and **Option D (Relocation to \`const-num-traits\`)**, with additional refinements for future-proofing.

### A. Refactor the Trait Shape (Evolution of Option G)
The current \`MulAccOps\` is overloaded and specifically tied to CIOS. Splitting it makes the Constant-Time (CT) contract explicit and removes the "discriminant dance" (\`GetWordOutput\`).

- **Rename to Primitive Operations:** Instead of \`CiosRowOps\`, name the traits after the underlying limb-array operations. This makes them reusable for other Montgomery variants (SOS, FIOS).
    - \`LimbRowOps\`: For \`mul_acc_row\`.
    - \`LimbShiftRowOps\`: For \`mul_acc_shift_row\` (a fused multiply-accumulate with a word-shift).
- **Explicit Word Access:** Use \`WordAccess\` and \`WordAccessCt\` as separate traits. The CT/Nct contract is now enforced at the trait bound level rather than via associated type magic.

### B. Placement in \`const-num-traits\` (Refined Option D)
\`const-num-traits\` is the most pragmatic home for these traits.
- **No New Crates:** Avoids the maintenance tax of a \`mont-traits\` crate.
- **Existing Dependency:** Both \`modmath\` and \`fixed-bigint\` already depend on \`const-num-traits\`.
- **Personality Alignment:** \`const-num-traits\` already manages the \`Ct\`/\`Nct\` personalities. It is the natural home for traits defining how limbs interact with those personalities.
- **Implementation:** Gate these traits behind a feature flag (e.g., \`features = ["limb-ops"]\`) to maintain a lean core for the crate.

## 2. Generalization for Heap-Bigints

To address the "fixed-limb-count" assumption (Point 3c in the exploration), the traits should be refactored to be less about "the bigint type" and more about "an object that provides limb-level row operations."

### Proposed Refinement:
Instead of:
\`\`\`rust
pub trait CiosRowOps: Sized + Copy + Default { ... }
\`\`\`
Define:
\`\`\`rust
pub trait LimbArrayOps {
    type Word: WordOps;
    fn len(&self) -> usize;
    fn mul_acc_row(&self, scalar: Self::Word, acc: &mut [Self::Word], carry: Self::Word) -> Self::Word;
    // ...
}
\`\`\`
By passing the accumulator as a slice (\`&mut [Self::Word]\`), the trait becomes compatible with both fixed-size and heap-allocated bigints, provided they can expose a view of their limbs.

## 3. Implementation Roadmap

| Phase | Action | Outcome |
| :--- | :--- | :--- |
| **1. Relocate & Refactor** | Define split traits in \`const-num-traits\` under a feature flag. | Traits are decoupled from \`fixed-bigint\`. |
| **2. Implement** | Update \`fixed-bigint\` to implement the new split traits. | Backends conform to the new primitive-focused API. |
| **3. Consume** | Update \`modmath/src/montgomery/cios.rs\` to bound on the new traits. | \`modmath\` drops the \`fixed-bigint\` dependency. |

## 4. Conclusion

This approach avoids the "inverted dependency" of Option B and the "crate explosion" of Option C. It treats limb-level row operations as the primitives they are, placing them in the crate already dedicated to numeric primitives (\`const-num-traits\`), and sets the stage for supporting heap-allocated bigints in the future without breaking the Constant-Time requirements.
