# CIOS Trait Placement & Typestate Synthesis Analysis

This document provides a strategic analysis of the CIOS (Coarsely Integrated Operand Scanning) trait placement and its integration with the broader `num-traits` (fork) and `modmath` ecosystem.

## 1. The Core Architectural Split

The challenge with `MulAccOps` is that it conflates **structural primitives** (how to read a limb) with **algorithmic primitives** (CIOS-specific row operations). To resolve the coupling between `modmath` and `fixed-bigint`, we must separate these concerns along the established "structural-low / domain-high" boundary.

### 1.1 Structural Layer (The "Bridge")
Accessing a limb from a multi-word integer is a structural primitive, not a domain-specific algorithm. It belongs at the lowest layer to allow any bigint backend to "speak" to any math algorithm.

**Recommendation:** Move limb access to `const-num-traits`.
*   **Trait:** `WordAccess` (Tier B) and `WordAccessCt` (Tier A).
*   **Payoff:** Standardizes how `modmath`, `fixed-bigint`, and third-party crates (like `num-bigint`) interact with limbs without requiring knowledge of the specific algorithm (CIOS, SOS, etc.).

### 1.2 Algorithmic Layer (The "Contract")
The row-multiply-accumulate operations are specific to the CIOS variant of Montgomery multiplication. These are "domain-high" operations.

**Recommendation:** Rename and relocate `MulAccOps`.
*   **New Name:** `CiosRowOps`.
*   **Placement:** Either `modmath` (as the algorithm owner) or a dedicated `bigint-traits` crate if further algorithms (NTT, SOS, FIOS) are planned.

---

## 2. Resolving the CIOS Quirks

### 2.1 The `GetWordOutput` Discriminant (Quirk 3a)
By splitting limb access into `WordAccess` (returning `Option<W>`) and `WordAccessCt` (returning `subtle::CtOption<W>`), we eliminate the need for a complex associated type discriminant.
*   `cios_montgomery_mul` bounds on `WordAccess`.
*   `cios_montgomery_mul_ct` bounds on `WordAccessCt`.
*   The type system now explicitly enforces the constant-time contract.

### 2.2 Static vs. Dynamic Size (Quirk 3c)
The requirement for `Self: Copy + Default` and a static `word_count()` method unnecessarily excludes heap-allocated bigints (e.g., `BoxedUint`).
*   **Fix:** Change to `fn word_count(&self) -> usize`.
*   **Fix:** Relax `Copy + Default` bounds where possible, or move them to the specific `FixedUInt` implementation rather than the trait definition.

---

## 3. Final Strategic Recommendation: Option H (The Structural Bridge)

This path combines the best of Option C (New Crate/Relocation) and Option G (Trait Splitting):

1.  **Decompose `MulAccOps`**: Split into `WordAccess` (Structural) and `CiosRowOps` (Algorithmic).
2.  **Limb Access in `const-num-traits`**: Implement `WordAccess` and `WordAccessCt` in `const-num-traits/src/ops/bigint.rs`. This provides the "Tier-A/B" foundation for all multi-limb types.
3.  **CIOS Traits in `bigint-traits`**: Create a lightweight `bigint-traits` crate to house `CiosRowOps`. This crate becomes the "Montgomery Contract" hub.
4.  **Zero-Coupling Goal**:
    *   `const-num-traits` is the universal dependency.
    *   `modmath` depends on `const-num-traits` and `bigint-traits`.
    *   `fixed-bigint` implements the traits from both, but `modmath` no longer needs to `use fixed_bigint::*`.

## 4. Integration with Typestate Synthesis

This relocation mirrors the philosophy applied to the `typestate` feature in `const-num-traits`:
*   **`PowerOfTwo<T>`** and **`HasNonZero`** are structural proofs (numeric primitives).
*   **`WordAccess`** is a structural proof (representation primitive).
*   All three enable **consuming operations** (`div_pow2`, `div_nonzero`, `mul_acc_row`) that provide real codegen or safety payoffs for downstream consumers.

By following this "Option H" path, `modmath` becomes truly generic, `fixed-bigint` remains a specialized provider, and `const-num-traits` serves as the high-performance, const-friendly glue that binds them.
