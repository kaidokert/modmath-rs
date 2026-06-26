# Roundup: CIOS Trait Placement & Panic-Free Architecture

This document summarizes the final architectural decision regarding the \`MulAccOps\` refactor, incorporating insights from the CIOS exploration and the panic-free operation mandates.

## 1. The Critical Pivot: Panic-Free Mandates

The requirement for "Panic-Free" operation (as detailed in \`PANIC_FREE_REQUESTS.md\`) has shifted the architectural calculus. To achieve zero \`panic_fmt\` symbols in embedded binaries, we must prioritize **infallible interfaces** and **structural proofs**.

### The Impact on Trait Design:
1.  **Infallibility over Discriminants:** The previous debate about \`Option<W>\` vs \`CtOption<W>\` is resolved by the panic-free requirement. Any \`unwrap()\` or branching on a \`None\` variant pulls in panic infrastructure.
2.  **The Public Index Probe:** Because the CIOS loop index is a public constant, word access can and should be infallible: \`fn word(&self, i: usize) -> Word\`. This is security-equivalent to the fallible versions while being DCE-friendly.
3.  **The Structural Bridge:** Panic-free coordination requires a standardized "Representation Bridge." Fragmenting these structural traits across multiple crates makes it impossible to write generic, panic-free math code.

## 2. Final Recommendation: The "Refined Option D" Path

We will adopt a hybrid of **Option D (Relocation to \`const-num-traits\`)** and **Option G (Infallible Shape Refactor)**.

### A. Trait Shape (The "Infallible G")
- **Rename:** \`MulAccOps\` -> \`LimbRowOps\` (and \`LimbShiftRowOps\`).
- **Accessor:** \`fn word(&self, i: usize) -> Word\`.
    - Bounds checks are handled via \`const\` assertions in the implementation or by the caller's loop logic, allowing the linker to DCE the panic sites.
- **Generalization:** \`word_count(&self)\` takes \`&self\` to allow for potential heap backends, but implementations for \`fixed-bigint\` will use const-generics to ensure zero-cost operation.

### B. Trait Placement (The "Structural Bridge")
- **Location:** \`const-num-traits\` behind a feature gate (e.g., \`features = ["limb-ops"]\`).
- **Logic:** This places limb access and row primitives alongside other "Representation Primitives" like \`ToBytes\` and \`OddNonzero\`. It provides the necessary "neutral ground" for \`modmath\` and \`fixed-bigint\` to communicate without dependency cycles.
- **Identity:** \`const-num-traits\` evolves from a "primitive mirror" to a "high-performance numeric bridge," which is consistent with its role in the panic-free ecosystem.

## 3. Revised Implementation Roadmap

| Phase | Action | Goal |
| :--- | :--- | :--- |
| **1. Bridge** | Define infallible \`WordAccess\` and \`LimbRowOps\` in \`const-num-traits\`. | Establish the panic-free contract. |
| **2. Align** | Update \`fixed-bigint\` to implement these traits using DCE-friendly patterns. | Eliminate panic sites in the bigint backend. |
| **3. Decouple** | Update \`modmath\` to use the new traits and remove the \`fixed-bigint\` dependency. | Achieve a truly generic, panic-free math library. |

## 4. Conclusion

The "genuine bifurcation" between Purists and Pragmatists is resolved by the engineering reality of embedded targets. A centralized, infallible "Structural Bridge" in \`const-num-traits\` is the only path that simultaneously satisfies the decoupling goal, the security requirements, and the zero-panic mandate.
