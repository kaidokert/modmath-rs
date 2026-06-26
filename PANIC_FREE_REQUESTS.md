# Panic-free operation: design requests for modmath

A downstream consumer (`ed25519_heapless`) is making a pass at being
panic-free in the linked-binary sense: zero `core::panicking::panic_fmt`,
`core::panicking::panic_bounds_check`, or `rust_begin_unwind` symbols
surviving linker DCE on embedded targets (AVR / Cortex-M / RISC-V).
The technique: cross-build, then `cargo nm <binary> | egrep "panic|unwind"`
must produce empty output.

We're applying three strategies in-crate (typestate, const generics,
explicit `Result/?`) and explicitly **not** trading panics for UB via
`unsafe { unreachable_unchecked() }`. A few remaining sites need API
shape changes upstream — this document is the shopping list.

## Blocker: `Field::new` / `FieldCt::new` return `Option`

`Field<T, P>::new(modulus: T) -> Option<Self>` returns `None` when the
modulus is even or zero. For Curve25519 (`p = 2^255 − 19`), the modulus
is statically known to be odd and non-zero, but consumers still write:

```rust
let field = Curve25519FieldCt::new(p).unwrap();
```

That `unwrap` pulls `panic_fmt` into the linked binary even though the
`None` arm is unreachable in practice. There's no way to express
"`Self::new` never returns `None` for this specific input" within the
type system as-is.

### Ask A — `OddNonzero<T>` typestate + infallible `from_odd_modulus`

Add a wrapper that encodes the validity check at the type level:

```rust
/// `T` known to be odd and non-zero.
pub struct OddNonzero<T>(T);

impl<T: /* whatever Field::new needs to verify */> OddNonzero<T> {
    /// Runtime check; returns None for even or zero.
    pub fn new(t: T) -> Option<Self>;

    /// `const fn` variant. Panics in const context → compile error
    /// at the call site if the input is invalid. No runtime panic.
    pub const fn new_const(t: T) -> Self;
}

impl<T, P> Field<T, P> {
    /// Infallible constructor — the typestate proves the precondition.
    pub fn from_odd_modulus(modulus: OddNonzero<T>) -> Self;
}
```

Downstream usage becomes:

```rust
const P: OddNonzero<T> = OddNonzero::new_const(/* compile-time prime */);
let field = FieldCt::from_odd_modulus(P);   // infallible, no panic
```

For this to work, `OddNonzero::new_const` needs `T` to support the
oddness/nonzero check in a const context. For `fixed_bigint::FixedUInt`
that means the `Parity` impl and equality-with-zero need to be
const-evaluable for our backend types (separate ask to fixed-bigint).

A weaker variant that still helps: a non-const `OddNonzero::new` plus
`Field::from_odd_modulus`. Callers that build the field once at
program startup pay the panic infra for that one site instead of for
every `verify` call — much smaller blast radius.

### Ask B — `from_modulus_unchecked` (less preferred)

```rust
impl<T, P> Field<T, P> {
    /// SAFETY: caller asserts modulus is odd and non-zero.
    pub unsafe fn from_modulus_unchecked(modulus: T) -> Self;
}
```

This pushes the cost to the caller and makes the obligation explicit.
But it forces consumers to write `unsafe`, which we'd prefer to avoid
in panic-elimination work. Ask A is strictly better.

## Nice-to-have: `MontStorage` is currently a marker

`MontStorage` requires `Zeroize` when the feature is on. That's load-
bearing for hardening. No request here — just noting it's already
panic-free as far as we can tell.

## Out of scope for this doc

- `FixedUInt::to_le_bytes(out)` returning `Result` (panics on
  buffer-too-small in the unwrap chain) — that's a `fixed-bigint` ask.
- `FixedUInt::from_le_bytes` panic-on-too-many-bytes behavior — same.
- Consumer-side trait method redesign in our `UnsignedModularInt` —
  we own that.

## Verification path

Once Ask A lands in modmath, we can update `ed25519_heapless`'s
`Curve25519FieldCt::curve25519()` factory to route through
`from_odd_modulus`, drop the `Curve25519FieldCt::new(p).unwrap()` /
`.expect()` chain, and re-run the AVR `cargo nm` gate. If the only
remaining panic symbols come from `fixed-bigint::to_le_bytes`, the
modmath side is done.
