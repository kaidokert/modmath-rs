# Leveraging `NonZero` (and where it differs from `Odd`)

How modmath should spend const-num-traits' non-zero capabilities to delete
divide-by-zero panic paths — and why `NonZero` is **not** redundant with `Odd`
even though "odd ⇒ non-zero."

The measure is the panic-free gate, not cycles:

> `cargo nm <binary> | egrep "panic|unwind"` → empty on the embedded targets.

A `% m` by a runtime modulus emits a `panic_fmt` (divide-by-zero) symbol that
survives linker DCE. Proving the divisor non-zero at the type level deletes it.

## The non-obvious part: two reduction paths, two different proofs

modmath has two families of modular ops, and they spend *different* proofs.
This is the crux, because at the value level `Odd` and `NonZero` look
interchangeable (zero is even, so every odd value is non-zero) — but they are
cashed in at different sites:

| Path | Divides by `m`? | What it actually needs | Proof to spend |
|---|---|---|---|
| **Montgomery** (`Field<T, P>`) | **No** — multiplies by `R⁻¹` | modulus **odd** (coprime to `2^k`) | **`Odd<T>`** |
| **Basic reduction** (`basic_mod_*`, `constrained_mod_*`) | **Yes** — `a % m` | modulus **non-zero** | **`NonZero` form of `T`** |

So:

- **`Field` does not benefit from `NonZero`.** It never divides by the modulus,
  so there is no divide-by-zero to delete. Its proof is `Odd<T>`, already wired:
  `Field::new_odd(Odd<T>)` is infallible and `Field::new(m)` is just
  `Odd::new(m).map(Self::new_odd)`. `Odd` covers *both* preconditions a
  Montgomery modulus needs — odd **and** non-zero — so do **not** reach for a
  separate `NonZero` here. (field.rs already does this.)

- **The `basic_mod_*` / `constrained_mod_*` functions are where `NonZero`
  pays off.** Every `a % m` / `b % m` in `add.rs`, `sub.rs`, `exp.rs` is a
  division by the runtime modulus `T`, each carrying a divide-by-zero branch.
  A modulus proven non-zero turns those into an infallible remainder.

Mixing them up is the trap: "the modulus is `Odd`, and odd implies non-zero, so
`Field` already removed the divide-by-zero" is **false** — `Field` never had a
divide-by-zero to begin with, and the `basic_mod_*` functions that *do* are a
separate code path that `Odd` doesn't touch.

## The capability path (how to actually wire it)

const-num-traits exposes the non-zero capability as two **carrier-implementable**
traits — plain traits with an associated type, no primitive-only blanket:

```rust
pub trait HasNonZero: Sized {
    type NonZero: Copy;                                  // for primitives: core::num::NonZero<Self>
    fn into_nonzero(self) -> Option<Self::NonZero>;      // the one check, at the boundary
    fn nonzero_get(nz: Self::NonZero) -> Self;
}
pub trait DivNonZero: HasNonZero {
    type Output;
    fn div_nonzero(self, d: Self::NonZero) -> Self::Output;   // no zero-check branch
    fn rem_nonzero(self, d: Self::NonZero) -> Self::Output;
}
```

cnt implements both for the 12 primitives. For a bignum carrier the impls live
**downstream** — this is the structural-low / domain-high split: the *capability*
is in cnt, the *carrier impl* is in fixed-bigint, the *consumer* is in modmath.

1. **fixed-bigint** implements `HasNonZero` + `DivNonZero` for `FixedUInt`,
   choosing its own `type NonZero` (a `FixedUInt` newtype carrying the
   `!= 0` invariant — `core::num::NonZero` is sealed to primitives and cannot
   wrap `FixedUInt`).

2. **modmath** takes the modulus once as a non-zero proof and threads it through
   the reduction, replacing `%` with `rem_nonzero`:

   ```rust
   // before: each `% m` has a divide-by-zero panic path
   pub fn basic_mod_add_pr<T>(a: T, b: T, m: T) -> T { /* … a % m … */ }

   // after: prove non-zero once at the boundary, reductions are infallible
   pub fn basic_mod_add_nz<T>(a: T, b: T, m: T::NonZero) -> T
   where
       T: DivNonZero<Output = T> + /* … */,
   {
       let a = a.rem_nonzero(m);        // no zero-check, no panic symbol
       let b = b.rem_nonzero(m);
       /* … */
   }
   ```

   The caller builds `m: T::NonZero` once (`m.into_nonzero()?` at the parse /
   key-load boundary) and reuses it; every reduction inside is straight-line.

### The `Copy` constraint and the cfg split

`HasNonZero::NonZero: Copy` means only a `Copy` carrier can use this bridge —
`FixedUInt` qualifies, `BoxedUint` does not. That is exactly the existing
fixed-vs-alloc split: the heapless/`FixedUInt` path uses cnt's `DivNonZero`; the
alloc/`BoxedUint` path stays on `crypto_bigint`'s own `NonZero`. Don't try to
unify the two non-zero designs — cfg-gate the adapter and let them coexist.

## Bridging `Odd` → `NonZero` (only if a Montgomery modulus feeds a `%` path)

If a `Field` ever has to hand its odd modulus to a basic-reduction subroutine
(e.g. a final canonical reduction, or an EEA/gcd inverse that divides), narrow
the proof rather than re-checking:

```rust
let nz: T::NonZero = odd_modulus.into_nonzero();   // infallible: odd ⇒ non-zero
```

This `Odd::into_nonzero` does **not exist in cnt yet** — it is the in-character
addition to make when such a caller appears, mirroring `Positive::into_nonzero`
(owned, routed through `HasNonZero`, vetted-unsafe construction justified by the
`Odd` proof). Until there is a concrete caller, do not add it: a proof-narrowing
with no consumer is speculative surface.

## Checklist

- [ ] `Field` / Montgomery: keep spending `Odd<T>`. No `NonZero` here. (done)
- [ ] fixed-bigint: `impl HasNonZero + DivNonZero for FixedUInt` with a
      `FixedUInt` `NonZero` newtype.
- [ ] modmath: add `*_nz` variants of `basic_mod_add` / `sub` / `exp` (and the
      `constrained_*` forms) taking `T::NonZero`, calling `rem_nonzero`.
- [ ] Verify with `cargo nm` that the divide-by-zero `panic_fmt` is gone from the
      `*_nz` path on the embedded fixture.
- [ ] Add `Odd::into_nonzero` to cnt **only** when a Montgomery modulus must feed
      a `%`-reduction subroutine.
