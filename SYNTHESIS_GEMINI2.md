# Synthesis: CIOS Trait Placement & Montgomery Architecture

This document synthesizes the analysis from four separate investigations (Claude, Codex, Gemini, and Gemini 2) into a unified architectural recommendation for the \`MulAccOps\` trait and its role in the \`modmath\` / \`fixed-bigint\` ecosystem.

---

## 1. Unified Architectural Model: The Three Layers

All analyses agree that the current \`MulAccOps\` conflates concerns. We must separate them into three distinct layers:

| Layer | Responsibility | Ownership (Proposed) |
| :--- | :--- | :--- |
| **L1: Word Primitives** | Fundamental math on single limbs (\`CarryingMul\`, etc.) | \`const-num-traits\` (Core) |
| **L2: Limb-Array Contract** | Reading limbs and row-level primitives (\`WordAccess\`, \`LimbRowOps\`) | \`const-num-traits\` (Feature-gated) |
| **L3: Algorithm Logic** | The CIOS double-loop and modular orchestration | \`modmath\` |

### The "Structural Bridge" Philosophy
The consensus is that **L2** is a structural bridge. It describes *how* a bigint exposes its limbs and *how* it performs raw row-level math. To decouple \`modmath\` (the consumer) from \`fixed-bigint\` (the provider), this bridge must sit in a neutral location.

---

## 2. The Shared Recommendation: "Refined Option D/G"

The collective recommendation is a hybrid that avoids the "crate explosion" of a new repository while achieving total decoupling.

### A. Refactor the Trait Shape (The "Option G" consensus)
- **Rename & Specialize**: \`MulAccOps\` is deprecated. In its place:
    - \`WordAccess\` / \`WordAccessCt\`: Explicit limb reading for NCT/CT personalities.
    - \`LimbRowOps\` / \`LimbShiftRowOps\`: Generic row primitives (reusable for SOS, FIOS, etc.).
- **Drop the Discriminant**: Replace the associated type \`GetWordOutput\` with explicit trait bounds (\`T: WordAccess\` or \`T: WordAccessCt\`). This lifts the security contract to the type system.
- **Future-Proofing**: Remove \`Copy + Default\` and \`word_count() -> usize\` (static). Use \`&self\` for word count and slice-based or factory-based accumulators to support **heap-allocated bigints**.

### B. Relocation (The "Option D" consensus)
- **Home**: \`const-num-traits\` under a feature flag (e.g., \`features = ["limb-ops"]\`).
- **Rationale**: Both main crates already depend on \`const-num-traits\`. It is the natural home for representation primitives that allow diverse backends to "speak" to generic algorithms.

---

## 3. The "Discriminant Probe" (Strategic Validation)

A critical point raised (Claude/Gemini) is whether the \`Option\`/\`CtOption\` discriminant is even necessary for word access in CIOS.
- **The Observation**: The word index in CIOS is a public loop counter (\`0..N\`), not a secret.
- **The Potential Gain**: If an in-bounds \`fn word(&self, i: usize) -> Word\` is confirmed safe for both NCT and CT, we can collapse the two algorithm bodies into one, significantly simplifying the \`modmath\` implementation.

---

## 4. Implementation Roadmap

1.  **Stage 1 (Neutral Bridge)**: Define \`WordAccess[Ct]\` and \`LimbRowOps\` in \`const-num-traits\`. Implement them in \`fixed-bigint\`.
2.  **Stage 2 (Algorithm Migration)**: Update \`modmath/src/montgomery/cios.rs\` to bound on these neutral traits. Delete the \`fixed-bigint\` dependency.
3.  **Stage 3 (Generalization)**: Relax bounds to support \`BoxedUint\` or other heap backends by allowing runtime word counts.

## 5. Summary Table of Benefits

| Aspect | Status Quo | Synthesis Recommendation |
| :--- | :--- | :--- |
| **Coupling** | Direct \`modmath\` -> \`fixed-bigint\` | **Zero** (Neutral Bridge) |
| **CT Security** | Hidden in associated types | **Explicit** in trait bounds |
| **Extensibility** | Fixed-size only | Supports **Heap Bigints** |
| **Complexity** | 1 Overloaded Trait | 3 Atomic, Reusable Traits |
| **Maintenance** | 2 Crates | 2 Crates (+ Feature in \`cnt\`) |

**Conclusion:** By treating limb-level operations as the structural primitives they are, we can achieve a generic, high-performance ecosystem that is both cleaner to maintain and ready for future bigint backends.
