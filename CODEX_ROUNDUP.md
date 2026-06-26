# Codex Roundup — CIOS Placement, Typestate, Panic-Free APIs

This is a compact synthesis after reading the CIOS placement analyses and the
panic-free request documents in `fixed-bigint-rs` and `modmath-rs`.

## Bottom Line

The panic-free requests strengthen the case for typestate and infallible,
precondition-proven APIs. They do not change the crate-placement conclusion for
CIOS row traits.

`CiosRowOps` is still an operational bigint-algorithm backend contract, not a
`const-num-traits` primitive. Keep `const-num-traits` at the word-arithmetic and
low numeric-proof layer; keep modular-field construction and CIOS algorithm
contracts in `modmath` or a neutral bigint-algorithm traits crate.

## What Is Settled

Universal agreement across the analyses:

- `MulAccOps` is misnamed. Rename it to `CiosRowOps` or equivalent.
- The current `GetWordOutput` associated-type discriminant should go.
- The CIOS algorithm body should stay in `modmath`.
- Do not move the algorithm into `fixed-bigint`.
- Do not use slice-based `Limbs` as the primary contract while CT matters.
- Open the door to heap backends where possible: prefer `word_count(&self)` and
  avoid trait-level `Copy + Default` unless the algorithm truly requires them.
- Verify the discriminant probe before committing to `WordAccess` /
  `WordAccessCt`: if the index is public and guaranteed in-bounds, a single
  infallible `word(&self, i) -> Word` may be CT-safe and more panic-free.

## Crate Boundaries

Recommended boundary:

- `const-num-traits`: scalar/word numeric atoms and low proof vocabulary.
  Examples: `Parity`, `NonZero` bridge, `CarryingMul`, `BorrowingSub`,
  `WideningMul`.
- `modmath`: modular arithmetic APIs, field construction, CIOS algorithm body,
  and possibly the first home for `CiosRowOps`.
- `fixed-bigint`: concrete backend implementation, panic-free byte APIs, and
  implementations of the row/access traits.
- future `bigint-algo-traits` / `mont-traits`: neutral home if a second backend
  or second algorithm family appears.

The important distinction is serialization representation vs operational limb
access. `ToBytes`/`FromBytes` expose canonical value encoding. `word(i)` and
`mul_acc_shift_row` expose algorithm layout and backend strategy. The former can
belong low; the latter should not drift into `const-num-traits`.

## Panic-Free Impact

The panic-free documents add a hard API constraint: downstream code must be able
to express impossible failure arms without `.unwrap()`, `.expect()`, or caller
`unsafe`.

For `modmath`, this means:

- Add a proof type for valid Montgomery modulus, e.g. `OddNonzero<T>`.
- Add an infallible constructor such as `Field::from_odd_modulus` /
  `FieldCt::from_odd_modulus`.
- Prefer a const constructor for statically known moduli where possible:
  invalid modulus fails at compile time, not through a runtime panic path.

For `fixed-bigint`, this means:

- Add const-generic fixed-buffer byte APIs such as `to_le_bytes_fixed` and
  `from_le_bytes_fixed`.
- Ensure the method bodies are themselves panic-symbol-free. Avoid runtime range
  indexing patterns that synthesize `panic_bounds_check`; prefer iterator or
  chunk-based bodies that the symbol gate verifies.

For CIOS traits, this means:

- Prefer infallible APIs when the algorithm proves the precondition.
- If `word(i)` is infallible, document the caller precondition
  `i < word_count()` and verify implementations do not introduce panic symbols
  for valid algorithm usage.
- Do not replace panics with `unreachable_unchecked`; the design goal is safe
  panic-free API shape, not UB.

## Recommended Path

Stage 1, in place:

- Rename `MulAccOps` to `CiosRowOps`.
- Probe whether word access can be a single infallible `word(&self, i) -> Word`.
- If not, split explicitly into `WordAccess` and `WordAccessCt` instead of using
  `GetWordOutput`.
- Keep the algorithm body in `modmath`.
- Keep placement as-is unless immediate decoupling is mandatory.
- Add panic-free field-construction API using `OddNonzero<T>`.
- Add panic-free fixed-size byte APIs in `fixed-bigint`.

Stage 2, only when justified:

- If immediate decoupling is mandatory, move row/access traits into `modmath`
  and have `fixed-bigint` implement them behind an optional `modmath` feature
  (Option H).
- If a second backend or second shared algorithm appears, extract the row/access
  traits into a neutral `bigint-algo-traits` or `mont-traits` crate.

## Final Position

The new panic-free context changes the urgency and the API shape, but not the
layering call.

Do not put CIOS operational limb traits in `const-num-traits`. Do make the APIs
infallible where typestate or const generics prove the preconditions. Then verify
with the actual `cargo nm | egrep "panic|unwind"` gate, because panic-freedom is
an emitted-binary property, not just a source-level property.
