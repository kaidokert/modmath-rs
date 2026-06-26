# CIOS trait redesign — draft

Concrete shape for Option H + post-probe trait collapse. Lives in modmath;
fixed-bigint impls under an optional `modmath` feature.

## 1. The new trait — `CiosRowOps`

```rust
/// Row-level operations for CIOS Montgomery multiplication.
///
/// # Preconditions
///
/// - `word(i)` is called only with `i < self.word_count()`.
/// - `i` is a public value (loop counter, constant) — never secret.
pub trait CiosRowOps: Default + Sized {
    type Word: Copy + PartialOrd;

    fn word_count(&self) -> usize;

    fn word(&self, i: usize) -> Self::Word;

    /// Phase 1: `acc += scalar * multiplicand`. Returns carry-out.
    fn mul_acc_row(
        scalar: Self::Word,
        multiplicand: &Self,
        acc: &mut Self,
        carry_in: Self::Word,
    ) -> Self::Word;

    /// Phase 2: `[acc, acc_hi] = ([acc, acc_hi] + scalar * multiplicand) >> word_bits`.
    /// Returns the carry word (0 or 1) from the fold.
    fn mul_acc_shift_row(
        scalar: Self::Word,
        multiplicand: &Self,
        acc: &mut Self,
        acc_hi: Self::Word,
    ) -> Self::Word;
}
```

**What's gone from `MulAccOps`:**
- `type GetWordOutput` discriminant — replaced by infallible `fn word(&self, i)`.
- `Self: Copy` bound — not actually needed by the algorithm; mutation goes through `&mut Self`.
- Static `word_count()` — now `&self`, admits heap-allocated bigints.
- The wholesale rename: `MulAccOps` → `CiosRowOps` (truth-in-naming; this is CIOS-specific).

**What's kept:**
- `Default` for constructing the zero accumulator. Natural; works for both fixed and heap.
- The two row ops, unchanged in shape.

## 2. The collapsed algorithm body

The two `cios_montgomery_mul` / `cios_montgomery_mul_ct` functions converge on
the inner loop (no more `?` vs `.into_option()?` discriminant). They still
diverge on the **final reduction** (branchful subtract vs. `conditional_select`),
which is a real algorithmic split, not a discriminant artifact. Extract the
shared loop into a helper.

```rust
/// Shared inner loop: produces `(acc, acc_hi)` ready for final reduction.
/// Pure function, no CT-vs-NCT split — both paths share this body.
fn cios_inner_loop<T>(
    a: &T,
    b: &T,
    modulus: &T,
    n_prime_0: T::Word,
) -> (T, T::Word)
where
    T: CiosRowOps,
    T::Word: num_traits::Zero
        + num_traits::One
        + num_traits::WrappingMul
        + num_traits::ops::overflowing::OverflowingAdd
        + core::ops::Add<Output = T::Word>
        + subtle::ConditionallySelectable,  // for the branchless carry-up step
{
    let n = a.word_count();
    let zero = <T::Word as num_traits::Zero>::zero();
    let one = <T::Word as num_traits::One>::one();
    let mut acc = T::default();
    let mut acc_hi = zero;
    let mut acc_hi2 = zero;

    for i in 0..n {
        let ai = a.word(i);  // infallible

        // Phase 1
        let carry = T::mul_acc_row(ai, b, &mut acc, zero);
        let (sum, overflow) = acc_hi.overflowing_add(&carry);
        acc_hi = sum;
        // Branchless overflow accumulation (was already CT in cios.rs).
        let overflow_word = T::Word::conditional_select(
            &zero, &one, subtle::Choice::from(overflow as u8),
        );
        acc_hi2 = acc_hi2 + overflow_word;

        // Reduction factor
        let m = acc.word(0).wrapping_mul(&n_prime_0);  // infallible

        // Phase 2
        let new_overflow = T::mul_acc_shift_row(m, modulus, &mut acc, acc_hi);
        acc_hi = acc_hi2 + new_overflow;
        acc_hi2 = zero;
    }

    (acc, acc_hi)
}

/// Variable-time CIOS — branchful final reduction.
pub fn cios_montgomery_mul<T>(
    a: &T,
    b: &T,
    modulus: &T,
    n_prime_0: T::Word,
) -> T
where
    T: CiosRowOps + PartialOrd + BorrowingSub + core::ops::Sub<Output = T>,
    T::Word: num_traits::Zero + num_traits::One + num_traits::WrappingMul
        + num_traits::ops::overflowing::OverflowingAdd
        + core::ops::Add<Output = T::Word>
        + subtle::ConditionallySelectable,
{
    debug_assert!(a < modulus, "CIOS input a must be in [0, modulus)");
    debug_assert!(b < modulus, "CIOS input b must be in [0, modulus)");

    let (mut acc, acc_hi) = cios_inner_loop(a, b, modulus, n_prime_0);
    let zero = <T::Word as num_traits::Zero>::zero();

    if acc_hi > zero || &acc >= modulus {
        let (result, _) = BorrowingSub::borrowing_sub(acc, *modulus, false);
        acc = result;
    }
    acc  // ← was Option<T>, now T directly
}

/// Constant-time CIOS — branchless final reduction.
pub fn cios_montgomery_mul_ct<T>(
    a: &T,
    b: &T,
    modulus: &T,
    n_prime_0: T::Word,
) -> T
where
    T: CiosRowOps + BorrowingSub + core::ops::Sub<Output = T>
        + subtle::ConditionallySelectable,
    T::Word: num_traits::Zero + num_traits::One + num_traits::WrappingMul
        + num_traits::ops::overflowing::OverflowingAdd
        + core::ops::Add<Output = T::Word>
        + subtle::ConstantTimeEq
        + subtle::ConditionallySelectable,
{
    use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

    let (acc, acc_hi) = cios_inner_loop(a, b, modulus, n_prime_0);
    let zero = <T::Word as num_traits::Zero>::zero();

    let (sub_result, borrow) = BorrowingSub::borrowing_sub(acc, *modulus, false);
    let acc_hi_nonzero = !acc_hi.ct_eq(&zero);
    let needs_sub = acc_hi_nonzero | !Choice::from(borrow as u8);
    T::conditional_select(&acc, &sub_result, needs_sub)  // ← was Option<T>
}
```

**Net effect:**
- One shared inner loop instead of two near-duplicates.
- Return type drops `Option<T>` → `T` (no failure mode left in either path).
- The `cios.rs`-internal `?` propagation goes away entirely.

## 3. `CiosMontMul` / `CiosMontMulCt` convenience traits

These also lose their `Option` return and the blanket impl loses the `MulAccOps`
dependency.

```rust
pub trait CiosMontMul: CiosRowOps + PartialOrd + BorrowingSub
    + core::ops::Sub<Output = Self>
{
    fn cios_mont_mul(a: &Self, b: &Self, modulus: &Self, n_prime: &Self) -> Self;
}

impl<T> CiosMontMul for T
where T: CiosRowOps + PartialOrd + BorrowingSub + core::ops::Sub<Output = T>,
      T::Word: /* same bounds as cios_montgomery_mul above */,
{
    fn cios_mont_mul(a: &Self, b: &Self, modulus: &Self, n_prime: &Self) -> Self {
        cios_montgomery_mul(a, b, modulus, n_prime.word(0))  // infallible
    }
}
```

Same shape for `CiosMontMulCt`. Both blanket impls now depend only on
`CiosRowOps` (modmath-defined) plus const-num-traits — **no `fixed_bigint::*` reference**.

## 4. `Field::mul` cleanup

```rust
// Before:
let mont = CiosMontMul::cios_mont_mul(&a.mont, &b.mont, &self.modulus, &self.n_prime)
    .expect("CIOS mul cannot fail with valid Montgomery parameters");

// After:
let mont = CiosMontMul::cios_mont_mul(&a.mont, &b.mont, &self.modulus, &self.n_prime);
```

The `expect` goes away → one less `panic_fmt` symbol in the linked binary,
per PANIC_FREE_REQUESTS.md gate.

## 5. fixed-bigint side

Under a new optional `modmath` feature in fixed-bigint:

```rust
// fixed-bigint/src/fixeduint/cios_row_ops_impl.rs (new file)
#[cfg(feature = "modmath")]
impl<T, const N: usize, P: Personality> modmath::CiosRowOps for FixedUInt<T, N, P>
where
    T: MachineWord + /* same bounds as the current MulAccOps impl */,
{
    type Word = T;

    fn word_count(&self) -> usize { N }

    fn word(&self, i: usize) -> T {
        self.array.get(i).copied().unwrap_or(T::zero())
        // Panic-free body: no `panic_bounds_check` synthesized.
        // CT-safe: i is a public value per the trait precondition.
    }

    fn mul_acc_row(...) -> T { /* reuse existing impl from MulAccOps */ }
    fn mul_acc_shift_row(...) -> T { /* reuse existing impl */ }
}
```

Both personalities (`Nct` and `Ct`) share the SAME impl now — no
discriminant split. The Ct version no longer needs the O(N) scan in `word`;
direct indexing is CT-safe under the public-index precondition.

**fixed-bigint Cargo.toml addition:**

```toml
[dependencies]
modmath = { version = "0.4", optional = true, default-features = false }

[features]
modmath = ["dep:modmath"]
```

`MulAccOps` itself stays in fixed-bigint for now (deprecated), to avoid a hard
break of any existing direct consumers. Modmath stops using it.

## 6. modmath Cargo.toml changes

```toml
[dependencies]
# Drop:
# fixed-bigint = { version = "0.4", optional = true, default-features = false }

[dev-dependencies]
fixed-bigint = { version = "0.4", features = ["modmath"], default-features = false }
```

modmath lib has zero `fixed_bigint::*` references. fixed-bigint is dev-only.
Cargo dev-dep cycle is documented and supported.

## 7. What goes away

- `use fixed_bigint::MulAccOps` in `modmath/src/montgomery/cios.rs` — gone.
- `use fixed_bigint::const_numtraits::BorrowingSub` — gone (use `const_num_traits::BorrowingSub`).
- `MulAccOps` itself — gone entirely from fixed-bigint (no deprecation window; coordinated cut).
- `GetWordOutput` associated type — gone.
- Two near-duplicate algorithm bodies — collapsed to one inner loop + two finalizes.
- `Option<T>` return from CIOS and `Field::mul`'s `expect` — gone (panic site killed).
- The `Self: Copy + Default` constraint on `MulAccOps` — relaxed to `Default + Sized`. Heap backends become possible.

## 8. Migration plan

Coordinated single-step cut across const-num-traits, fixed-bigint, and modmath.
No deprecation period; `MulAccOps` is deleted in the same release that introduces
`CiosRowOps`. Validation is end-to-end through downstream experimental branches
(ed25519, rsa, pqc) before publication.

1. **const-num-traits**: ship the typestate work already in flight (already
   independent of CIOS).
2. **modmath**: define `CiosRowOps` + the new infallible CIOS body. Drop
   `fixed-bigint` from `[dependencies]`. Add to `[dev-dependencies]` with
   `features = ["modmath"]`.
3. **fixed-bigint**: add the `modmath` feature + `CiosRowOps` impl;
   delete `MulAccOps` entirely (file + module + re-export). See
   `fixed-bigint-rs/CIOS_MIGRATION.md` for the receiving side.
4. **Cut release candidates** of all three crates together.
5. **Experimental branches** on ed25519, rsa, and pqc test against the
   release candidates end-to-end.
6. Once green on all three downstreams, **publish** the coordinated release.
