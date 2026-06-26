# Typestate synthesis — follow-up: `Canonical<T>` revival and `Prime<T>` placement

Extends `typestate_synthesis.md` (the June 2026 synthesis). Captures two
decisions raised during the v0.4 const-num-traits migration work:

1. Reviving the `_pr` → `Canonical<T>` cut that was deferred in June.
2. Where `Prime<T>` (and its sibling number-theoretic typestates) should
   live when revived: **a new `const-num-integer` leaf crate**, not
   modmath.

Both are extensions of the same architectural line the original
synthesis drew — structural-low typestates in `const-num-traits` (or its
sibling leaves), domain-high typestates in modmath. Nothing in here
overrides the original; it sharpens two open items.

## Status anchor (what's shipped since the original synthesis)

Settled and landed during the v0.4 cycle, anchoring the design context
the proposals below build on:

- **`Odd<T>` cascade through the Montgomery surface.** `Field::new_odd`,
  `basic_montgomery_mod_mul_pr_odd`, the four exp variants (`_pr_odd`,
  `_odd`, `_pr_odd_ct`, `_odd_ct`). Old `Option`-returning shapes
  preserved as thin wrappers. Re-exported at `modmath::Odd`.
- **`Field::try_new_odd_ct`** — `CtOption<Field<T, Ct>>` for secret
  modulus (RSA-CRT private-prime construction). Uses `Odd::new_ct` from
  cnt's `ct` feature; modmath enabled the feature.
- **`HasNonZero` + `DivNonZero` `*_nz` surface.** Twelve siblings across
  add/sub/mul/exp × basic/constrained/strict, taking `T::NonZero`,
  routing through `rem_nonzero` instead of `% m`. fixed-bigint shipped
  the carrier impl (`NonZeroFixedUInt` newtype) with audit confirming
  panic-deletion. Re-exported at `modmath::{HasNonZero, DivNonZero}` and
  available as `modmath::{basic,constrained,strict}::nonzero::*`.
- **`modmath-cios` leaf crate.** `CiosRowOps` trait + primitive impls
  extracted to a stable-identity leaf so the dev-dep cycle resolves at
  the rustc-metadata level, not just the cargo level. **This is the
  template `const-num-integer` should follow.**

The path-by-path panic-deletion / typestate-discharge story now reads:

| Path | Proof spent | Where the proof construction lives |
|---|---|---|
| Montgomery (`Field`, `basic_montgomery_mod_*_odd`) | `Odd<T>` | cnt |
| Montgomery, secret modulus (RSA-CRT) | `CtOption<Odd<T>>` via `Odd::new_ct` | cnt (`ct` feature) |
| Schoolbook with reduction (`basic_mod_*_nz`) | `T::NonZero` | cnt (capability) + FB (carrier) |
| Schoolbook pre-reduced (`*_pr`) | **nothing yet — doc'd precondition only** | TBD (this proposal) |
| Algebraic operations (Fermat inv, etc.) | `Prime<T>` (or `OddPrime<T>`) | TBD (this proposal) |

The last two rows are this doc's scope.

---

## Part 1: Revive `Canonical<T>` for the `_pr` surface

### Prior status

The June 2026 synthesis recorded this as **deferred**:

> Deferred: `Reduced<T>` for the schoolbook `_pr` surface (June 2026)
> - Discussion concluded that the schoolbook `_pr` suffix could be
>   replaced by `Reduced<T>` typestate, but:
>   - krabipqc voted **unbranded** `Reduced<T>` (their own brand on top,
>     lifetime branding from modmath would be friction).
>   - User chose to defer with "no Reduced<T> for now" — krabipqc
>     immediate need was `mul_acc`, not Reduced.
> - The design conclusion: when revived, `Reduced<T>` should be
>   **unbranded** at the modmath level. Consumers add their own brand
>   (curve, scheme, modulus) as a thin wrapper.

Also from the same doc, on naming:

> `Canonical` vs `Reduced` naming: codex follow-up argues `Canonical` is
> more accurate because the `_pr` APIs require `[0, m)` (canonical
> representative), not just modular equivalence. Probably resolve in
> favor of `Canonical` and reserve `Congruent` for mere equivalence if
> it's ever needed.

### Why revive now

The deferral reason ("krabipqc immediate need was `mul_acc`") is no
longer load-bearing. `wide_montgomery_mul_acc` shipped in v0.3.1. krabipqc
is on `_pr` throughout (confirmed by direct audit during the `*_nz`
discussion) and would benefit from typestate-enforcement of the
"operands in `[0, m)`" precondition — currently doc-only, silently wrong
on un-reduced input.

### Confirmed consumers

Two independent consumer signals after the `*_nz` interlude — strengthens
the proposal beyond "krabipqc would benefit":

**1. krabipqc (immediate, schoolbook hot loop).** Confirmed during the
`*_nz` audit: 100% on the `_pr` discipline, `_pr` bodies are already
`%`-free by design, `Canonical<T>` would catch the only failure mode
left (silent garbage on un-reduced input).

**2. rsa-heapless private-key path (when scope opens).** Confirmed by
the RSA agent during the v0.4 migration. The public-key surface
(verify, encrypt) is 100% Montgomery and doesn't need `Canonical<T>` —
`Residue<'f, T, P>` lifetime branding suffices. But the private-key
path is *half* Montgomery, *half* schoolbook:

- **CRT decrypt:** `c^dp mod p` and `c^dq mod q` are Montgomery, but
  the reassembly `(m1 - m2) * qinv mod p` and the full-width
  `m2 + h*q mod n` step **cross fields and drop back to schoolbook mod
  ops**. That's where `*_nz` (boundary-validate `p`, `q`, `n` once) +
  `Canonical<T>` (track `< p` / `< n` through the reassembly chain)
  land cleanly. Without them, the schoolbook reassembly is the
  panic-fmt + boundary-friction soup the public path got to avoid.
- **Keygen:** Miller-Rabin loops, `gcd` / `e⁻¹ mod φ(n)`, computing
  `dp` / `dq` / `qinv` — heavily schoolbook with strong `raw < something`
  invariants propagating through long chains. Carrying those as types
  instead of re-checking is a real win.

The RSA-private-key consumer is on the horizon (next experimental scope
on `rsa-heapless`), not in flight today. But the design pressure it
applies is concrete: `Canonical<T>` lands now for krabipqc; when the
RSA-private branch opens, it's load-bearing infrastructure rather than
"nice to have."

The framing matters: the **public-key verify/encrypt path is
Montgomery-dominant and these typestates are noise**; the
**private-key CRT/keygen path is half Montgomery, half schoolbook**,
and the schoolbook half is where `_nz` + `Canonical<T>` infrastructure
becomes load-bearing. Same answer-shape as the `*_nz` interlude.

The structural argument is also cleaner now:

- `_pr` bodies are designed to be `%`-free (cond-sub via `wrapping_sub`
  + magnitude comparison). That's the structural guarantee.
- The `Canonical<T>` typestate captures exactly the precondition those
  bodies require: `value < modulus` at entry.
- Output is *also* canonical because the cond-sub finalize produces
  values in `[0, m)` by construction. The wrapper round-trips
  naturally — `Canonical → Canonical`.

### Proposed shape (concrete)

Following the June consensus: **unbranded** at the modmath level.

```rust
// modmath::canonical (new submodule)

/// Proof that a value is in `[0, m)` for some modulus `m`.
///
/// **Unbranded** — the type does not carry modulus identity. Consumers
/// that need modulus identity wrap `Canonical<T>` in their own newtype
/// (e.g. `KyberReduced<u32>(Canonical<u32>)`). modmath does not
/// presume the brand strategy; per the June 2026 synthesis, that
/// decision belongs to the consumer.
///
/// Discharges the documented precondition on the `_pr` schoolbook
/// surface (`basic_mod_*_pr`, `constrained_mod_*_pr`, `strict_mod_*_pr`):
/// today the API requires `a < m && b < m` and silently produces
/// nonsense on out-of-range input. With `Canonical<T>` operands, that
/// precondition becomes a compile-time guarantee.
#[repr(transparent)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Canonical<T>(T);

impl<T> Canonical<T> {
    /// # Safety
    /// `value < modulus` for the modulus this will be used with.
    pub const unsafe fn new_unchecked(value: T) -> Self { Self(value) }

    /// The proven-canonical value.
    pub fn get(self) -> T { self.0 }
}

impl<T: PartialOrd> Canonical<T> {
    /// `Some` iff `value < modulus`. One runtime check at the trust
    /// boundary; subsequent `_canonical` calls take and return
    /// `Canonical<T>` straight-line, no re-check.
    pub fn new(value: T, modulus: &T) -> Option<Self> {
        if &value < modulus { Some(Self(value)) } else { None }
    }
}
```

`_pr` siblings get `_canonical` (or `_c`) variants — additive, mirroring
the `Odd` / `NonZero` cuts:

```rust
pub fn basic_mod_mul_pr_canonical<T>(
    a: Canonical<T>,
    b: Canonical<T>,
    m: T,                  // raw modulus — its own invariants (odd,
                           // non-zero) are separately discharged by
                           // `Odd<T>` / `T::NonZero` when needed.
) -> Canonical<T>          // output preserves the invariant
where T: /* existing _pr bounds */
{
    let raw = basic_mod_mul_pr(a.get(), b.get(), m);
    // SAFETY: `_pr` body's cond-sub finalize produces values in [0, m).
    // The body's contract IS the proof of canonicity.
    unsafe { Canonical::new_unchecked(raw) }
}
```

Chains without re-checking, which is the actual ergonomic win:

```rust
let a = Canonical::new(raw_a, &q).expect("PQC sampler invariant");
let b = Canonical::new(raw_b, &q).expect("PQC sampler invariant");
let r = basic_mod_mul_pr_canonical(a, b, q);   // Canonical<T>
let s = basic_mod_add_pr_canonical(r, c, q);   // still Canonical<T>
```

### Open knobs (need decision before landing)

1. **Naming on the function side: `_pr_canonical` or just `_canonical`?**
   The `_pr` suffix becomes redundant once the typestate encodes the
   precondition. Two options:
   - Land as `*_pr_canonical` siblings, keep `*_pr` alongside. Maximum
     back-compat, name length feels excessive.
   - Land as `*_canonical` siblings (drop `_pr` from the new path),
     keep `*_pr` as the unchecked-precondition entry. Cleaner names,
     same back-compat. **Lean: this.**

2. **Construction side: `Canonical::new -> Option`, plus `new_ct`?**
   - Phase 1: just `Option`. No `Canonical::new_ct` until a secret-input
     caller materializes (per NONZERO_LEVERAGE's "no speculative
     surface" rule).
   - When a CT consumer appears (RSA-CRT decoded ciphertext as a
     candidate), add `Canonical::new_ct(value, &modulus) -> CtOption<Self>`
     using `ConstantTimeLess` on `&value < modulus`. **Lean: phase 1
     only.**

3. **Scope: 12 schoolbook entry points only, or also Montgomery `_pr_odd`?**
   - The Montgomery `_pr_odd` family (`basic_montgomery_mod_mul_pr_odd`,
     etc.) has the same `a < modulus && b < modulus` precondition.
     Hoisting it into `Canonical<T>` is natural.
   - But the Montgomery consumer typically goes through `Field` →
     `Residue<'f, T, P>`, which already carries the proof via lifetime
     branding. Hoisting into the low-level `_pr_odd` direct entries
     duplicates the Residue mechanism.
   - **Lean: schoolbook only, defer Montgomery `_pr_odd_canonical` until
     a low-level Montgomery consumer asks.**

4. **`Reduced<T>` deprecated alias?**
   - The June synthesis used `Reduced<T>` as the working name. Anyone
     who wrote against that name (probably nobody yet, since the cut
     was deferred) would want a `pub type Reduced<T> = Canonical<T>;`
     for one cycle.
   - **Lean: no alias. The deferral means no one shipped against
     `Reduced<T>`; ship `Canonical<T>` cleanly.**

### Why this isn't `_nz`'d redundantly

`Canonical<T>` and `T::NonZero` are orthogonal — they discharge different
preconditions:

- `T::NonZero`: modulus is non-zero (operand-shape on the modulus side,
  enables panic-free `% m`).
- `Canonical<T>`: operand is `< modulus` (operand-shape on the value
  side, enables `_pr` skipping `% m` entirely).

A fully type-state-discharged `_pr` Montgomery hot loop wants both, plus
`Odd<T>` for the modulus:

```rust
fn pqc_montgomery_step<T>(
    a: Canonical<T>,
    b: Canonical<T>,
    m: Odd<T>,
) -> Canonical<T>
{ /* infallible, no `Option`, no `%`, no panic_fmt */ }
```

Three orthogonal proofs, three independent boundary checks, one
straight-line hot loop. That's the design center.

---

## Part 2: `Prime<T>` (and friends) → `const-num-integer` leaf crate

### Why not in cnt

`Prime<T>` doesn't fit cnt's design center:

| Wrapper | Predicate cost | Domain |
|---|---|---|
| `Odd<T>`, `Even<T>` | O(1) bit check | structural-low |
| `NonZero<T>`, `PowerOfTwo<T>` | O(1) — zero check / popcount | structural-low |
| **`Prime<T>`** | **O(poly)** — Miller-Rabin, or `unsafe new_unchecked` for cryptographic constants | **algebraic / number-theoretic** |

cnt's typestates are uniform: O(1) structural predicates with cheap
total constructors. `Prime<T>` is polynomial-time or unsafe-constructor —
distinct design axis.

### Why not in modmath

`Canonical<T>` lives in modmath because it's **modulus-relative** — the
predicate (`value < modulus`) requires the modulus. `Prime<T>` is not
modulus-relative; it's a property of the value itself, like `Odd<T>`.
The "this has to live in modmath" argument doesn't apply.

Reasons against modmath specifically:

- **Wider consumer ecosystem.** Crypto libs that don't use modmath
  (rsa with its own field math, ed25519/x25519 implementations bringing
  their own scalar field, ZK-SNARK field-math layers, signature schemes,
  ECC code) would still want `Prime<T>` / `Coprime<T, U>`. Stuffing it
  in modmath ties it to modmath's release cadence and design.
- **Domain mismatch.** Modmath is the *modular arithmetic* layer.
  `Prime<T>` is a *number-theoretic* property. They cluster differently.
- **Sibling typestates.** Once `Prime<T>` lands, natural neighbors
  follow: `Coprime<T, U>` (RSA `gcd(e, λ(n)) = 1`; EEA inversion is
  infallible when modulus is coprime to value), `Composite<T>` (RSA
  modulus `n = p·q`), `OddPrime<T>` (the cryptographic-prime case).
  These don't belong in modmath either.
- **Companion operations.** Const-callable `gcd`, `lcm`, `extended_gcd`,
  Miller-Rabin (feature-gated), Jacobi / Legendre symbol — all natural
  fits for a number-theory crate. Modmath shouldn't accrete these.

### Why a separate leaf crate (`const-num-integer`)

Matches the existing stable-ecosystem layering:

```
num-traits      → num-integer      → various-math
const-num-traits → const-num-integer → modmath, crypto crates, ...
```

Anyone reaching for the const analog of `num_integer::Integer::gcd` would
look there first. The leaf-crate pattern also worked for `modmath-cios`:
single-purpose, stable-identity, zero-dep, primary tests on primitives,
carrier impls downstream.

Reversibility cost favors getting this right early:

- Put `Prime<T>` in modmath → later want `const-num-integer` → breaking
  migration for any consumer that pinned to `modmath::Prime`.
- Put it in `const-num-integer` → if uptake is narrow, just an extra
  crate to publish. Much smaller cost.

### Why not spin it up *right now*

The "no speculative surface" rule applies. Currently:

- `Field::inv_fermat` returns `Option<Residue>` because residue might be
  zero. With `Prime<T>` modulus + `NonZero<Residue>`, it's infallible.
  But: **no current downstream needs `inv_fermat` to be infallible**.
  ed25519 verify doesn't use it. krabipqc PQC arithmetic doesn't use
  it.
- ed25519's scalar-inversion uses Curve25519's scalar field — separate
  consumer from the field-prime inversion `Prime<T>` would unlock.

### The concrete consumer that *would* trigger it: RSA-heapless keygen

Refining "wait for a concrete consumer" — the next planned scope
opening on `rsa-heapless` is **keygen**, which is exactly the
`Prime<T>` / `OddPrime<T>` / Miller-Rabin trigger:

- Generate candidate `p`, `q` via deterministic CSPRNG.
- Probabilistic primality test (Miller-Rabin) on each candidate.
- A passing candidate becomes an `OddPrime<T>` proof (the witness is
  the test result; the typestate captures the conclusion).
- Compute `n = p * q` (Composite<T>), `φ(n) = (p-1)(q-1)`, validate
  `gcd(e, φ(n)) = 1` (Coprime<T, U>), compute `d = e⁻¹ mod φ(n)`.
- Per-prime CRT precompute: `dp = d mod (p-1)`, `dq = d mod (q-1)`,
  `qinv = q⁻¹ mod p` (the inversion is infallible because `q` is
  coprime to `p` — both prime).

This is the consumer that exercises the full `const-num-integer`
surface and gives the design real pressure rather than speculation.
Miller-Rabin behind a feature, `OddPrime<T>` as primary, `Coprime<T, U>`
for the EEA inversion path, const-callable `gcd` / `extended_gcd`.

**Trigger condition:** when `rsa-heapless` opens its keygen scope (not
in flight today; on the horizon per the RSA agent's report).
**Not:** "now, on speculation." The build-it-when-summoned discipline
that worked for `modmath-cios` works here too.

### Sequencing

1. **Now: land `Canonical<T>` in modmath** (Part 1). Real consumer
   (krabipqc + anyone using `_pr` directly), clear value, mature
   design.
2. **Hold `const-num-integer`** until a `Prime<T>` consumer raises a
   flag.
3. **At that point: spin up `const-num-integer`** as a leaf crate
   following the `modmath-cios` template. Initial surface:
   - `OddPrime<T>` (primary — the cryptographic-prime case is the
     common consumer, and carrying odd + prime separately gets ugly
     fast)
   - `Prime<T>` (the `Prime == 2` edge case, defer indefinitely if no
     consumer asks)
   - `Coprime<T, U>`
   - Const-callable `gcd`, `extended_gcd`, `lcm` (where feasible on
     stable; richer surface on nightly via `c0nst!`)
   - Miller-Rabin behind a feature (`primality-test` or similar)
4. Modmath consumes it: `Field::inv_fermat_prime(odd_prime: OddPrime<T>,
   nonzero_residue: NonZero<Residue>) -> Residue` becomes infallible.

### Pre-emptive design lock-in

Even before `const-num-integer` exists, one decision worth recording so
nobody re-litigates it later: **`OddPrime<T>` is the primary surface,
not `Odd<T>` + `Prime<T>` separately.**

Reasons:

- Every cryptographic prime > 2 is odd.
- The common Montgomery+Fermat consumer wants both proofs in one
  parameter. `Field::inv_fermat_prime(&self, modulus: OddPrime<T>, ...)`
  is the natural signature; carrying them separately fragments the API.
- `OddPrime: Deref<Target = Odd<T>>` (or an explicit
  `.into_odd() -> Odd<T>`) gives the structural sub-proof without
  forcing redundant typing.

`Prime<T>` alone exists for the `Prime == 2` case. Defer until a
consumer for `Prime<u8>(2)` shows up — likely never.

---

## Consumer matchup matrix (three-way confirmation)

Three downstream agents independently audited their consumption against
the new typestate surface (`Odd<T>` boundary, `try_new_odd_ct`, `*_nz`,
proposed `Canonical<T>`, future `OddPrime<T>` via `const-num-integer`).
The dividing line is sharp and consistent: **Montgomery-dominant
consumers don't benefit from the new schoolbook typestates;
schoolbook-touching consumers do.** This vindicates the scope decisions
on every recent cut.

| Consumer | Path shape | `Odd<T>` boundary | `try_new_odd_ct` | `*_nz` | `Canonical<T>` | `OddPrime<T>` (future) |
|---|---|---|---|---|---|---|
| **krabipqc** (PQC) | 100% Montgomery + `_pr` schoolbook | yes — `Field::new_odd(KYBER_Q)` | no — modulus public const | no — `_pr` already %-free | **yes — caught the silent-on-bad-input case** | no — q is fixed const |
| **ed25519 verify** | 100% Montgomery | yes — `Field::new_odd(P25519)` | no — modulus public const | no | no | no |
| **ed25519 sign (planned)** | 100% Montgomery + new scalar-q field | **yes — `Field::new_odd(SCALAR_Q)` for the new ScalarField** | no — both p and q public consts | no | no | no |
| **rsa-heapless verify/encrypt** | 100% Montgomery | yes | no — public RSA modulus | no | no | no |
| **rsa-heapless CRT decrypt (future)** | half Montgomery, half schoolbook | yes — for p, q Fields | **yes — secret p, q (CRT)** | **yes — schoolbook reassembly** | **yes — `(m1-m2)*qinv mod p` chain** | no |
| **rsa-heapless keygen (future)** | mostly schoolbook | yes | yes (Miller-Rabin tested odd primes) | yes — `gcd(e, φ(n))` etc. | yes — long invariant chains | **yes — primality testing trigger** |

The pattern that falls out:

- **Public-key, public-modulus, Montgomery-dominant** (ed25519
  verify+sign, krabipqc, RSA verify/encrypt): the `Odd<T>` cut is
  load-bearing for the construction site, everything else is noise.
  The verify-side typestate vocabulary (`Odd<T>`, Ct personality,
  `Residue<'f, T, P>` branding) is sufficient.
- **Secret-modulus / schoolbook-touching** (RSA-CRT decrypt, RSA
  keygen, krabipqc's `_pr` reductions): the `*_nz` + `Canonical<T>`
  infrastructure becomes load-bearing. `try_new_odd_ct` matters for
  secret-prime construction.
- **Keygen / number-theoretic** (RSA keygen): triggers
  `const-num-integer`. No other consumer in the current roadmap does.

This is the dividing line the original synthesis drew implicitly between
"structural-low" (cnt + sibling leaves) and "domain-high" (modmath +
consumers). The three-agent confirmation pins it concretely: every
landing decision was correctly scoped, every deferral is on the right
side of the line.

**Implication for `Canonical<T>` knob #3 (scope: schoolbook only vs.
also Montgomery `_pr_odd`):** triple-confirmed lean toward
**schoolbook only**. None of the Montgomery-dominant consumers
(krabipqc PQC, ed25519 verify+sign, RSA public-key) want
Montgomery `_pr_odd_canonical` — `Residue<'f, T, P>` lifetime branding
covers them. Adding it would be type-system theater for contracts the
Residue mechanism already proves dynamically.

## What this doesn't change

- **The `_nz` cut as shipped.** `Canonical<T>` doesn't subsume it;
  they're orthogonal. The `_nz` surface stays as-is; a `_canonical_nz`
  combined variant (operand canonical + modulus non-zero) is a
  separate decision when both proofs flow at the same call site, which
  is the schoolbook hot loop case.
- **`Odd<T>` re-export at `modmath::Odd`.** Stays. Consumers shouldn't
  have to track which leaf crate hosts which typestate.
- **`Residue<'f, T, P>` lifetime branding for Field.** Stays. It's a
  *different* mechanism — modulus-identity branding — that
  `Canonical<T>` (unbranded) does not replace. They coexist.
- **The deferral list in the original synthesis.** Items not addressed
  here (e.g. `Borrow<T>` consolidation of the schoolbook surface,
  `MontParams<T>` algorithm-family split) remain in their original
  status.

---

## Checklist

- [ ] Land `Canonical<T>` in modmath (Part 1). Decide the four open
      knobs (naming, `new_ct`, scope, `Reduced<T>` alias).
- [ ] Hold `const-num-integer` (Part 2). When a `Prime<T>` consumer
      materializes, spin it up following the `modmath-cios` template
      with `OddPrime<T>` as primary surface.
- [ ] Record the `OddPrime<T> > Prime<T>` decision in the eventual
      `const-num-integer` design doc so reviewers don't re-litigate.
