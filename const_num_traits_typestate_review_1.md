# const-num-traits typestate sketch — review 1 of 3

Reviews the `MaybeNonZero` / `Odd<T>` / `Even<T>` sketch proposed for the
const-num-traits fork (NOT the modmath typestate roadmap, though related).
Saved verbatim for later synthesis. Captured 2026-06-16.

---

The proposed sketch is an excellent addition to the num-traits (fork) ecosystem. It aligns perfectly with the crate's "Tier-based" philosophy and its use of c0nst for nightly compatibility.

## Review of the Proposed Addition

1. **MaybeNonZero as a Bridge:**
   - **Architectural Fit:** This is the correct design for a trait-centric library. By using an associated type `type NonZero`, you allow backends like crypto-bigint or num-bigint to provide their own "proof" types (e.g., `Odd<Uint<L>>`) while primitives stay mapped to `core::num::NonZero`.
   - **Payload:** The `DivNonZero` trait is the primary payoff. In generic numeric code (like Montgomery reduction), the ability to perform infallible division (Tier A/B instead of Tier C) is a significant performance and safety win.

2. **Odd<T> / Even<T>:**
   - **Transparency:** `repr(transparent)` is essential for zero-cost proofs.
   - **Borrowing:** The `from_ref` pattern using unsafe is the standard way to handle these proofs for non-Copy types (like large bigints).
   - **Consistency Note:** To make the generic `from_ref` work, you should ensure `Parity` is implemented for `&T` (e.g., `impl<'a, T: Parity> Parity for &'a T`).

---

## Other "Obvious" Typestate Opportunities

Beyond NonZero, Odd, and Even, there are several other numeric proofs that would provide direct optimization payoffs in this crate:

### 1. PowerOfTwo<T> (The Bit-Twiddler's Choice)
- **Property:** `x > 0 && (x & (x - 1)) == 0`.
- **Optimization:** Division and remainder by a PowerOfTwo are simple bit-shifts and masks (Tier A).
- **Payoff:** `pub fn div_pow2(self, p: PowerOfTwo<T>) -> T`. This turns a Tier-C division into a Tier-A bitwise op in generic code.
- **Existing Base:** Already has a predicate in `src/ops/pow2.rs`.

### 2. Positive<T> / Negative<T> (Algebraic Proofs)
- **Property:** `x > 0` or `x < 0`.
- **Optimization:**
  - `abs(Positive<T>)` is a no-op.
  - `signum(Positive<T>)` is constant 1.
- **Domain Safety:** Algorithms for log, sqrt, and pow can skip domain checks if they accept `Positive<T>`.

### 3. Unit<T> (The Multiplicative Identity)
- **Property:** `x == 1`.
- **Optimization:** Multiplication and division by Unit are no-ops. In generic unrolled loops (like matrix ops or dot products), knowing a value is 1 allows the compiler to prune entire branches of arithmetic.

### 4. Finite<T> (Floating-Point Guard)
- **Property:** `!is_nan() && !is_infinite()`.
- **Optimization:** Allows PartialOrd to act as Ord (no NaN branches).
- **Payoff:** Many algorithms (sorting, stats, interpolation) can avoid the "Tier B" branching logic required to handle exceptional float values.

### 5. Prime<T> (Marker Typestate)
- **Property:** `x` is prime.
- **Payoff:** In modular arithmetic, it guarantees the existence of a multiplicative inverse. While verifying primality is expensive, the type allows "inheriting" the proof from a generator (like a CSPRNG) to downstream generic crypto-ops.

## Implementation Recommendation
All of these should be gated behind your proposed typestate feature. Wrap the constructors and traits in `c0nst::c0nst!` so that `Odd::new(3).unwrap()` can be evaluated at compile-time on nightly, which is the "killer feature" of this specific fork.
</content>
