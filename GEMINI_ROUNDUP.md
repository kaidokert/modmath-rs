# GEMINI Roundup: The Zero-Panic Architectural Mandate

This document provides the final assessment for the `MulAccOps` relocation and typestate synthesis, factoring in the critical requirements for panic-free operation on embedded targets (AVR/Cortex-M) as requested in `PANIC_FREE_REQUESTS.md`.

## 1. The Decision: Linker Necessity Trumps "Pure" Layering

The discovery of the `PANIC_FREE_REQUESTS` from downstream consumers (`ed25519_heapless`) resolves the bifurcation between the "Clean Architecture" and "Representation Primitive" views.

To achieve zero `panic_fmt` and `panic_bounds_check` symbols in a linked binary, the ecosystem requires **infallible structural bridges**. This requirement dictates the placement of limb-access and numeric proofs.

### 1.1 The Role of `const-num-traits` (CNT)
CNT is no longer just a "numeric primitives" crate; it is the **structural foundation** that provides the typestates and representation traits needed to eliminate panics across the entire `modmath`/`fixed-bigint` stack.

*   **Limb Access (`WordAccess`)**: Belong in CNT. It is a "Representation Primitive" parallel to `ToBytes`. Standardizing this in CNT allows any bigint algorithm to be written without `Option/Result` unwraps, which are the primary source of panic symbols.
*   **Numeric Proofs (`Odd`, `NonZero`)**: Belong in CNT. These are the "Panic-Free Proofs" that enable infallible constructors in `modmath` (e.g., `Field::from_odd_modulus`).

## 2. Refining the CIOS Contract

The "discriminant probe" (the observation that CIOS uses a public loop counter) leads to a major simplification:

*   **Infallible Indexing**: Since the word index is a public counter (`0..word_count`), the `WordAccess` trait should provide an **infallible** lookup: `fn get_word(&self, i: usize) -> Self::Word`.
*   **Abolishing the Discriminant**: We can drop the `GetWordOutput` associated type entirely. Constant-time (Ct) and Non-constant-time (Nct) personalities both provide the same infallible signature; the Ct implementation simply uses a branchless scan internally.
*   **Outcome**: This collapses the two near-identical algorithm bodies in `modmath/src/montgomery/cios.rs` into a single, clean, and **panic-free** implementation.

## 3. The Final Architectural Stack

| Layer | Component | Home | Rationale |
| :--- | :--- | :--- | :--- |
| **L1: Primitives** | Word-arithmetic (`CarryingMul`) | `const-num-traits` | Fundamental atoms. |
| **L2: Structural** | `WordAccess`, `Odd`, `NonZero` | **`const-num-traits`** | The "Zero-Panic" Bridge. |
| **L3: Contract** | `CiosRowOps` (Renamed `MulAccOps`) | **`modmath`** | Domain-high algorithmic contract. |
| **L4: Algorithm** | CIOS Double-Loop | `modmath` | Reusable, panic-free body. |

---

## 4. Actionable Roadmap

1.  **Stage 1: The Foundation (CNT)**
    *   Implement `WordAccess` (infallible word lookup).
    *   Implement `Odd`, `Even`, `NonZero` typestates.
    *   All gated behind a `typestate` / `bigint` feature.

2.  **Stage 2: The Contract (modmath)**
    *   Rename `MulAccOps` to `CiosRowOps`.
    *   Update signature: `word_count(&self)` and relaxed `Copy + Default` to admit heap-bigints.
    *   Refactor the CIOS body to use infallible `get_word`, ensuring zero `panic_fmt` symbols survive DCE.

3.  **Stage 3: The Proof (fixed-bigint)**
    *   Implement the new CNT traits for `FixedUInt`.
    *   Add the `to_le_bytes_fixed<const M: usize>` and `from_le_bytes_fixed` methods to satisfy the direct panic-free requests for serialization.

## 5. Final Conclusion

By treating limb access and numeric proofs as foundational structural primitives in `const-num-traits`, we provide the high-level crates (`modmath`, `fixed-bigint`) with the alphabet they need to speak **infallibly**. This resolves the architectural tension while directly satisfying the linker-level constraints of embedded consumers.
