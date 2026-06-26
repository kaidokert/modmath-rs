# const-num-traits typestate sketch — review 2 of 3

Reviews the `MaybeNonZero` / `Odd<T>` / `Even<T>` sketch proposed for the
const-num-traits fork. Saved verbatim for later synthesis. Captured 2026-06-16.

---

## Sketch Review

**Findings:**

1. **`core::num::NonZero<Self>` is not a stable public primitive bridge.**
   The crate already imports and implements over stable concrete aliases like `NonZeroU32`, `NonZeroI64`, etc. in `src/bounds.rs:2`. A `MaybeNonZero` bridge should map primitives to those concrete aliases, probably via macro arms, unless you gate the whole thing behind a nightly std feature.

2. **NonZero does not make signed division fully infallible.**
   It removes divide-by-zero, but signed `MIN / -1` and `MIN % -1` still overflow. That affects `CheckedDiv`, `CheckedRem`, `DivCeil`, `DivFloor`, `Euclid`, `DivExact`, and `NextMultipleOf`; see the existing panics/checks in `src/ops/rounding.rs:23` and `src/ops/rounding.rs:152`. For unsigned, `DivNonZero` is genuinely infallible. For signed, the trait either needs to document remaining overflow, return Option, or use a stronger divisor proof.

3. **The proposed typestate feature is not present in Cargo.toml.**
   This is easy to add, but it matters because this crate is careful about feature split semantics. If these are crate-owned stable types, `typestate = []` is fine; if they expose unstable std names, they belong under the existing `nightly-std` style instead.

4. **Odd<T> / Even<T> are coherent, but their payoff is narrower than NonZero.**
   The `Parity` trait is already a cheap per-type atom in `src/ops/parity.rs:19`. `Odd<T>` becomes valuable where "odd modulus" or "invertible modulo 2^n" is a real precondition. `Even<T>` is less obviously useful unless paired with exact-halving/trailing-zero APIs.

5. **`from_ref` is a reasonable pattern, but it only helps non-Copy types that implement `Parity` for `&T`.**
   The sketch's bound handles that correctly. I'd add `as_inner(&self) -> &T`, maybe `into_inner(self) -> T`, and avoid naming `get(self)` unless you intentionally want to diverge from `core::num::NonZero*` expectations.

## Other Obvious Typestate Opportunities

The strongest candidates I see:

- **NonZero<T> / MaybeNonZero**
  Highest value. It can remove zero checks from division, remainder, Euclid, MultipleOf, logs, HighestOne, and LowestOne. Existing Option/panic surfaces are visible in `src/ops/log.rs:20` and `src/ops/bits.rs:270`.

- **Positive<T> and NonNegative<T> for signed integers**
  Positive makes integer logs infallible for signed values; NonNegative makes signed isqrt infallible. See `src/ops/log.rs:32` and `src/ops/sqrt.rs:33`. For unsigned, NonZero already covers the log input case.

- **PowerOfTwo<T> for unsigned integers**
  This is a natural wrapper over `src/ops/pow2.rs:15`. It enables infallible/exact ilog2, cheap masks, alignment, divide/rem-by-power-of-two specialization, and can replace repeated is_power_of_two checks.

- **BitIndex<T> or ShiftAmount<T> proving n < BITS**
  This would remove Option from checked shifts and remove panics from funnel shifts. The relevant checks are in `src/ops/bits.rs:156`, `src/ops/bits.rs:197`, and `src/ops/bits.rs:214`. This is useful but more ergonomically awkward because the proof depends on the target width.

- **NonMin<T> for signed integers**
  This would make abs/neg infallible where `CheckedAbs`/`CheckedNeg` only fail on MIN. Useful, but less broadly compelling than NonZero, Positive, or PowerOfTwo.

I would start with `MaybeNonZero` plus unsigned `DivNonZero`/`RemNonZero`, then Positive/NonNegative, then PowerOfTwo. Odd is worth adding if downstream crypto/modular arithmetic is a target; Even is probably secondary unless you also add exact-halving or valuation-style APIs.
</content>
