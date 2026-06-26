# Synthesis: CIOS Trait Placement & Bigint Ecosystem Architecture

This document synthesizes the analysis from four independent reviews regarding the placement and refactoring of the `MulAccOps` trait. It establishes a unified architectural roadmap for decoupling `modmath`, `fixed-bigint`, and `const-num-traits`.

---

## 1. The Unified Vision: Layered Abstraction

All reviews converge on the fact that `MulAccOps` currently conflates three distinct architectural layers. The synthesis recommendation is to formalize these layers:

| Layer | Responsibility | Current Home | Targeted Home |
| :--- | :--- | :--- | :--- |
| **L1: Primitives** | Word-level carrying arithmetic (`CarryingMul`, etc.) | `const-num-traits` | `const-num-traits` (as-is) |
| **L2: Structural** | Limb/Word access and CT-contract (`WordAccess`) | `fixed-bigint` | **`const-num-traits`** (new bridge) |
| **L3: Contract** | Row-level algorithm primitives (`CiosRowOps`) | `fixed-bigint` | **`bigint-traits`** (Option C) or **`modmath`** (Option H) |
| **L4: Algorithm** | The actual CIOS double-loop implementation | `modmath` | `modmath` (reusable body) |

---

## 2. Strategic Decisions

### 2.1 Trait Shape: The "Option G" Mandate
We must split `MulAccOps` to make the Constant-Time (CT) contract explicit and the row-ops reusable.
*   **Decomposition:** Split into `WordAccess` (structural) and `CiosRowOps` (algorithmic).
*   **CT Safety:** Replace the associated-type discriminant with explicit `WordAccess` (returning `Option`) and `WordAccessCt` (returning `CtOption`).
*   **The Discriminant Probe:** Before finalizing the split, verify if the word index in CIOS (a public loop counter) genuinely requires a `CtOption`. If it doesn't, we can collapse to a single `WordAccess` trait + one algorithm body, maximizing cleanliness.

### 2.2 Placement: Option C + H Hybrid
The consensus favors decoupling via a neutral ground, but with a nuance on "one-shot vs. pattern."

*   **Recommendation:** Move L2 (Structural) to `const-num-traits` immediately. It is a foundational primitive.
*   **Recommendation:** Place L3 (Contract) in a new **`bigint-traits`** crate (Option C) if NTT, SOS, or other algorithms are planned. If CIOS is a one-shot, keep the contract in `modmath` (Option H).
*   **Constraint:** The **algorithm body** (L4) must stay in `modmath`. Moving it to `fixed-bigint` (Option B) is rejected as it forces re-implementation for every new backend.

### 2.3 Future-Proofing (The "Heap" Bridge)
To support heap-allocated bigints (e.g., `BoxedUint`), the traits must be refactored:
*   Change `word_count()` to an instance method: `fn word_count(&self) -> usize`.
*   Relax `Copy + Default` bounds on the bigint type; instead, bound the accumulator or provide explicit constructors.

---

## 3. Implementation Roadmap

1.  **Stage 1 (Foundation):** Implement `WordAccess` and `WordAccessCt` in `const-num-traits/src/ops/bigint.rs`.
2.  **Stage 2 (Refactor):** Rename `MulAccOps` to `CiosRowOps` in `fixed-bigint`, updating it to inherit from the new `WordAccess` traits.
3.  **Stage 3 (Decouple):**
    *   If "Pattern": Create `bigint-traits` for `CiosRowOps`.
    *   If "One-shot": Move `CiosRowOps` definition to `modmath`.
4.  **Stage 4 (Cleanup):** Update `modmath/src/montgomery/cios.rs` to bound on the new traits and drop the direct `fixed-bigint` dependency.

---

## 4. Final Verdict

The "Purist Clean" path is **Option C + G**: a linear dependency graph where `const-num-traits` provides the structural glue, a new `bigint-traits` provides the algorithm contracts, and `modmath` provides the high-level algebra. This path maximizes reuse, preserves CT safety, and admits future heap-based backends with zero architectural friction.
