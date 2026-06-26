Based on a review of the `modmath-rs` codebase, there is a pervasive pattern of "primitive obsession" in the low-level free functions (using `T` for everything). While `field.rs` introduces a nice `Residue` brand at the highest level, the underlying arithmetic modules (`add.rs`, `mul.rs`, `montgomery/basic_mont.rs`, etc.) duplicate a massive amount of code to handle different states of reduction and validity.

Introducing typestate wrappers could drastically compact the API surface, eliminate runtime branches, and remove `Option` return types. Here are the highest-value opportunities:

### 1. `Odd<T>` or `OddNonZero<T>` for Montgomery Moduli

**The Problem:**
Montgomery arithmetic strictly requires an odd, non-zero modulus. Currently, functions like `basic_montgomery_mod_mul`, `basic_montgomery_mod_exp`, and `compute_montgomery_params` start with runtime checks:

```rust
if modulus == T::zero() || modulus.is_even() {
    return None;
}

```

Because of this, all complete pipeline functions return `Option<T>`, forcing the caller to `unwrap()` or handle errors on every single multiplication or exponentiation.

**The Typestate Solution:**
Introduce an `Odd<T>` wrapper.

```rust
pub struct Odd<T>(T);

impl<T: Parity + Zero> Odd<T> {
    pub fn new(val: T) -> Option<Self> {
        if val.is_even() || val == T::zero() { None } else { Some(Self(val)) }
    }
    // Safety escape hatch for constants
    pub const fn new_unchecked(val: T) -> Self { Self(val) }
}

```

**Benefits:**

* **Eliminates Options**: Functions like `basic_montgomery_mod_mul(a: T, b: T, m: &Odd<T>) -> T` no longer need to return `Option`.
* **Removes Branches**: The `is_even()` and `is_zero()` checks are removed from the hot path (like inside inner loop instantiations or pipeline calls). The validation happens exactly once when the modulus is wrapped.

### 2. `Reduced<T>` to Eliminate the `_pr` (Pre-Reduced) API Duplication

**The Problem:**
Across `add.rs`, `sub.rs`, `mul.rs`, and `exp.rs`, the API is duplicated to handle values that are already modulo $M$. For example, there is `basic_mod_add` (which performs `% m`) and `basic_mod_add_pr` (which assumes $a < m$ and $b < m$). This doubles the entire API surface (`basic_mod_exp` vs `basic_mod_exp_pr`, `strict_mod_mul` vs `strict_mod_mul_pr`, etc.).

**The Typestate Solution:**
Introduce a `Reduced<T>` wrapper to indicate a value is strictly in the range $[0, m)$.

```rust
pub struct Reduced<T> {
    val: T,
}

impl<T> Reduced<T> {
    /// Performs the modulo reduction once
    pub fn new(val: T, modulus: &T) -> Self { ... }
    pub unsafe fn new_unchecked(val: T) -> Self { Self { val } }
}

```

**Benefits:**

* **API Compaction**: You can delete all `_pr` variants. Instead of `basic_mod_add_pr`, you just implement arithmetic traits directly on `Reduced<T>`.
* **Zero-Cost Chaining**: If a user chains `(a + b) * c`, `add(Reduced, Reduced) -> Reduced` guarantees the output is already reduced, so `mul(Reduced, Reduced)` statically knows it doesn't need to apply `% m` to its inputs.

### 3. `MontForm<T>` for the Montgomery Free Functions

**The Problem:**
In `montgomery/basic_mont.rs` and others, there are functions like `basic_montgomery_mul` that expect inputs in Montgomery form (`a_mont`, `b_mont`). However, they just accept raw `T`. If a user accidentally passes a standard integer instead of a Montgomery-transformed integer, the compiler won't catch it, leading to silently incorrect math.

**The Typestate Solution:**

```rust
pub struct MontForm<T>(T);

```

Functions like `basic_to_montgomery` should return `MontForm<T>`. Functions like `basic_montgomery_mul` should accept `MontForm<T>` and return `MontForm<T>`. `basic_from_montgomery` should consume `MontForm<T>` and return `Reduced<T>`.

**Benefits:**

* **Type Safety**: Prevents domain mixing (multiplying a raw integer by a Montgomery integer).
* **Cleaner abstraction**: While `Residue<'f, T, P>` in `field.rs` does this, pushing `MontForm<T>` down into the free functions makes the underlying math primitives foolproof as well.

### 4. `Wide<T>` for Double-Width Products

**The Problem:**
In `wide_mul.rs` and throughout the Montgomery REDC code, widening multiplication returns a tuple `(T, T)` representing `(lo, hi)`. `wide_redc` takes `t_lo: T, t_hi: T`. This forces the reader to track argument order, and the compiler cannot prevent accidentally swapping `lo` and `hi`.

**The Typestate Solution:**

```rust
pub struct Wide<T> {
    pub lo: T,
    pub hi: T,
}

```

Change `WideMul::wide_mul` to return `Wide<T>`, and `wide_redc` to accept `Wide<T>`.

**Benefits:**

* **Readability & Safety**: `wide_redc(product, ...)` is much clearer than `wide_redc(lo, hi, ...)`. It completely eliminates the class of bugs where `hi` and `lo` are accidentally swapped in function parameters.

### 5. `SignedMagnitude<T>` replacing the ad-hoc `Signed<T>`

**The Problem:**
In `inv/signed.rs`, there is an internal `Signed<T>` struct used to support the Extended Euclidean Algorithm. It requires manual, repetitive canonicalization checks (e.g., ensuring `-0` becomes `+0`) throughout the `Add` and `AddAssign` trait implementations, littering the code with `if result.value == T::zero() { result.negative = false; }`.

**The Typestate Solution:**
If `NonZero<T>` was implemented, `Signed<T>` could be refactored into an enum that leverages the type system to make `-0` impossible to represent:

```rust
pub enum Signed<T> {
    Zero,
    Positive(NonZero<T>),
    Negative(NonZero<T>),
}

```

**Benefits:**

* **Eliminates Canonicalization Branches**: You completely remove the need to check and fix `-0` states after every mathematical operation, speeding up the hot loop of the modular inverse algorithm.
