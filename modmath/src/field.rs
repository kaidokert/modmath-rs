//! `Field<T, P>` — Montgomery-form modular arithmetic context, parameterized
//! over the limb personality (`P: Personality`).
//!
//! This module owns the *generic* Montgomery surface (`mul`, `exp`, `add`,
//! `sub`, `reduce`, `into_raw`, `zero`, `one`); curve- or RSA-specific
//! specializations (`lazy_add`, Fermat `inv`, Solinas reduction, etc.) live
//! in consumer crates that wrap these types.
//!
//! ## Personality parameter
//!
//! `Field<T, P>` selects the *algorithm* at the modmath level (variable-
//! time CIOS with branching finalize vs. CT CIOS with conditional-select
//! finalize); the personality of `T` itself (e.g. `FixedUInt<W, N, Nct>`
//! vs `FixedUInt<W, N, Ct>`) selects the *limb-primitive* bodies. In
//! practice the two personalities co-vary at the construction site —
//! you build `Field<FixedUInt<W, N, Nct>, Nct>` for verify paths and
//! `Field<FixedUInt<W, N, Ct>, Ct>` for signing paths. Type aliases
//! [`FieldCt`] / [`FieldNct`] (and the matching [`ResidueCt`] /
//! [`ResidueNct`]) are the canonical per-personality spellings — they
//! read naturally at the call site and sidestep the type-inference
//! ambiguity that bare `Field::new(modulus)` hits when two `impl`
//! blocks (Nct/Ct) are both method-resolution candidates.
//!
//! Because the Nct and Ct algorithms have **different trait bounds on `T`**
//! (`T: CiosMontMul` for Nct vs `T: CiosMontMulCt + ConditionallySelectable`
//! for Ct), the per-personality method bodies live in separate `impl`
//! blocks rather than dispatching via `match P::TAG`. Only the bounds
//! shared by both algorithms (precompute, `new`, `zero`, `one`, the
//! residue brand) live in the common `impl<T, P: Personality>` block.
//!
//! ## Branding
//!
//! Each `Field<T, P>` instance is implicitly tagged by its borrow lifetime
//! `'f`, and the residues it produces carry that same brand plus the
//! personality parameter:
//!
//! ```ignore
//! let field = Field::new(modulus).unwrap();
//! let r: Residue<'_, U256, Nct> = field.reduce(&seven);
//! ```
//!
//! The borrow checker prevents a `Residue` from outliving its parent
//! `Field`, and the `P` parameter prevents an Nct residue from being fed
//! to a Ct method (and vice versa) at compile time. Covariance does not
//! prevent mixing residues from two `Field` instances built in the same
//! scope with the same `P` — a known limitation matched by ed25519's
//! `field.rs`. A generative brand
//! (`PhantomData<fn(&'f ()) -> &'f ()>` + closure) would close that gap
//! when a future consumer needs it.

use core::marker::PhantomData;

use const_num_traits::{Ct, Nct, Odd, Personality};

use crate::montgomery::basic_mont::{
    wide_montgomery_mul, wide_montgomery_mul_acc, wide_montgomery_mul_acc_ct,
    wide_montgomery_mul_ct, wide_redc, wide_redc_ct,
};
use crate::montgomery::{
    CiosMontMul, CiosMontMulCt, compute_n_prime_newton, compute_r_mod_n, compute_r2_mod_n,
    type_bit_width,
};
use crate::parity::Parity;
use crate::wide_mul::WideMul;

/// Bound on the value stored in a [`Residue`]. With the `zeroize`
/// feature this requires `T: Zeroize`; otherwise it's vacuous.
#[cfg(feature = "zeroize")]
pub trait MontStorage: zeroize::Zeroize {}
#[cfg(feature = "zeroize")]
impl<T: zeroize::Zeroize> MontStorage for T {}

/// Bound on the value stored in a [`Residue`]. With the `zeroize`
/// feature this requires `T: Zeroize`; otherwise it's vacuous.
#[cfg(not(feature = "zeroize"))]
pub trait MontStorage {}
#[cfg(not(feature = "zeroize"))]
impl<T> MontStorage for T {}

// ---------------------------------------------------------------------------
// Field<T, P>
// ---------------------------------------------------------------------------

/// Montgomery context over modulus `T`, with algorithm choice driven by
/// the personality marker `P` (defaults to [`Nct`] = variable-time fast
/// path).
///
/// Use `Field<T, Nct>` for operations on **public** data only — signature
/// verification, RSA public-key encryption, anything whose inputs are not
/// secret. Use `Field<T, Ct>` (also reachable via the [`FieldCt`] type
/// alias) for secret-handling paths.
///
/// `Clone` is a trivial 4×`T` memcpy — it does NOT re-run the
/// `compute_r_mod_n` / `compute_r2_mod_n` precompute that `new()` does.
/// Callers building a `Field` once per key and then reusing it (e.g.
/// RSA-CRT) should clone the prebuilt instance rather than calling
/// `new()` again with the same modulus.
#[derive(Clone, Debug)]
pub struct Field<T, P: Personality = Nct> {
    modulus: T,
    n_prime: T,
    r_mod_n: T,
    r2_mod_n: T,
    _p: PhantomData<fn() -> P>,
}

/// Alias for the Nct variant of [`Field`]. Equivalent to `Field<T, Nct>`
/// (matches the default personality). Provided for symmetry with
/// [`FieldCt`] and to side-step the construction-site type-ambiguity
/// pitfall — `FieldNct::new(modulus)` resolves unambiguously without
/// the type-annotation/turbofish friction of `Field::new(modulus)`.
pub type FieldNct<T> = Field<T, Nct>;

/// Alias for the Ct variant of [`Field`]. Equivalent to `Field<T, Ct>`.
/// Reads naturally at construction sites and sidesteps the type-inference
/// ambiguity that bare `Field::new(modulus)` hits — `FieldCt::new(modulus)`
/// resolves unambiguously because the alias fixes `P = Ct` at the type
/// level. Symmetric with [`FieldNct`].
pub type FieldCt<T> = Field<T, Ct>;

/// A value in `Field<T, P>`, stored implicitly in Montgomery form.
///
/// The `'f` lifetime brand ties this residue to its parent `Field`; the
/// `P` parameter ties it to the parent's algorithm personality. The
/// borrow checker rejects code that uses a residue after its `Field` is
/// dropped, or that mixes residues across personalities. See module docs
/// for the covariance caveat.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Residue<'f, T: MontStorage, P: Personality = Nct> {
    mont: T,
    _brand: PhantomData<&'f ()>,
    _p: PhantomData<fn() -> P>,
}

#[cfg(feature = "zeroize")]
impl<'f, T: MontStorage, P: Personality> zeroize::Zeroize for Residue<'f, T, P> {
    fn zeroize(&mut self) {
        self.mont.zeroize();
    }
}

#[cfg(feature = "zeroize")]
impl<'f, T: MontStorage, P: Personality> Drop for Residue<'f, T, P> {
    fn drop(&mut self) {
        self.mont.zeroize();
    }
}

#[cfg(feature = "zeroize")]
impl<'f, T: MontStorage, P: Personality> zeroize::ZeroizeOnDrop for Residue<'f, T, P> {}

/// Alias for the Nct variant of [`Residue`]. Equivalent to
/// `Residue<'f, T, Nct>`. Symmetric with [`ResidueCt`].
pub type ResidueNct<'f, T> = Residue<'f, T, Nct>;

/// Alias for the Ct variant of [`Residue`]. Equivalent to
/// `Residue<'f, T, Ct>`. Symmetric with [`ResidueNct`].
pub type ResidueCt<'f, T> = Residue<'f, T, Ct>;

// ---------------------------------------------------------------------------
// Shared impls (any P)
// ---------------------------------------------------------------------------

impl<T: MontStorage, P: Personality> Residue<'_, T, P> {
    /// Returns a reference to the underlying Montgomery-form value.
    ///
    /// **Escape hatch.** Intended for downstream specialization layers
    /// (e.g. `Curve25519Field`) that implement fast paths reading the raw
    /// limbs. General consumers should not call this — use the methods on
    /// [`Field`] instead.
    pub fn mont_value(&self) -> &T {
        &self.mont
    }
}

impl<T, P: Personality> Field<T, P> {
    /// Construct a `Field` directly from already-computed Montgomery
    /// parameters. **`const fn` — usable in const initializers.**
    ///
    /// Intended for callers whose modulus is statically known at compile
    /// time (curve constants, PQC parameters, RSA group constants for a
    /// fixed key, etc.) and who want to expose the Field as a `const`
    /// associated item or static, rather than paying the runtime
    /// [`new`](Self::new) precompute on each instantiation.
    ///
    /// **The caller is responsible for the correctness of `n_prime`,
    /// `r_mod_n`, and `r2_mod_n`** — see [`compute_n_prime_newton`],
    /// [`compute_r_mod_n`], and [`compute_r2_mod_n`] for the algorithms.
    /// Those helpers are not `const fn` today (their bodies use trait
    /// method calls on `T`), so const-initializer callers must either
    /// hand-compute the values, use a build script, or — for primitive
    /// integers — compute them in a non-const context once at startup
    /// and cache.
    ///
    /// No invariant checking is performed here. Passing inconsistent
    /// parameters produces a `Field` whose arithmetic methods return
    /// silently incorrect results.
    ///
    /// [`compute_n_prime_newton`]: crate::compute_n_prime_newton
    /// [`compute_r_mod_n`]: crate::compute_r_mod_n
    /// [`compute_r2_mod_n`]: crate::compute_r2_mod_n
    pub const fn from_precomputed(modulus: T, n_prime: T, r_mod_n: T, r2_mod_n: T) -> Self {
        Self {
            modulus,
            n_prime,
            r_mod_n,
            r2_mod_n,
            _p: PhantomData,
        }
    }
}

impl<T, P: Personality> Field<T, P>
where
    T: Copy
        + PartialEq
        + PartialOrd
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::WrappingMul
        + const_num_traits::WrappingAdd
        + const_num_traits::WrappingSub
        + const_num_traits::ops::overflowing::OverflowingAdd
        + core::ops::Add<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Mul<Output = T>
        + Parity
        + MontStorage,
{
    /// Construct a new `Field` from an already-proven-odd modulus.
    ///
    /// **Infallible.** The `Odd<T>` typestate hoists the "modulus is odd and
    /// nonzero" precondition to the caller's trust boundary — typically a
    /// single `Odd::new(p)?` (or `Odd::new(p).unwrap()` for a const modulus)
    /// at config load. No runtime check inside this constructor, and the
    /// `panic_fmt` symbol that an `unwrap()` on the old `Option` API would
    /// have synthesized stays out of the linked binary on embedded targets
    /// when the boundary check is const-evaluated.
    ///
    /// `Odd<T>` covers both the "non-zero" and "odd" halves (zero is even),
    /// so this also discharges the modulus-nonzero check that [`new`] does.
    ///
    /// [`new`]: Self::new
    pub fn new_odd(modulus: Odd<T>) -> Self {
        let modulus = modulus.get();
        let w = type_bit_width::<T>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = compute_r_mod_n(modulus, w);
        let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);
        Self {
            modulus,
            n_prime,
            r_mod_n,
            r2_mod_n,
            _p: PhantomData,
        }
    }

    /// Construct a new `Field` from an already-proven-odd modulus,
    /// using the **constant-time** precompute path.
    ///
    /// Same precompute as [`new_odd`] (Newton's iteration for `N'`,
    /// repeated modular doublings for `R mod N` and `R² mod N`), but
    /// the doubling-and-reduction loop in
    /// [`compute_r_mod_n_ct`](crate::montgomery::compute_r_mod_n_ct)
    /// avoids value-dependent branches on the modulus.
    ///
    /// Use this when `modulus` is secret (e.g. RSA-CRT private primes
    /// `p`, `q`). For public moduli (ed25519 / Curve25519 / krabipqc),
    /// [`new_odd`] is faster and equivalent.
    ///
    /// Cost vs [`new_odd`]: one extra `wrapping_sub` and one
    /// `conditional_select` per modular doubling step (`w` per
    /// precompute call). Negligible against the subsequent field
    /// operations the precompute amortizes.
    ///
    /// [`new_odd`]: Self::new_odd
    pub fn new_odd_ct(modulus: Odd<T>) -> Self
    where
        T: subtle::ConditionallySelectable + subtle::ConstantTimeLess + const_num_traits::CtIsZero,
    {
        let modulus = modulus.get();
        let w = type_bit_width::<T>();
        let n_prime = compute_n_prime_newton(modulus, w);
        let r_mod_n = crate::montgomery::compute_r_mod_n_ct(modulus, w);
        let r2_mod_n = crate::montgomery::compute_r2_mod_n_ct(r_mod_n, modulus, w);
        Self {
            modulus,
            n_prime,
            r_mod_n,
            r2_mod_n,
            _p: PhantomData,
        }
    }

    /// Construct a new `Field` over the given (odd, nonzero) `modulus`.
    ///
    /// Returns `None` if `modulus` is zero or even (Montgomery requires odd N).
    /// Thin wrapper around [`new_odd`] that performs the parity proof at
    /// runtime. Prefer [`new_odd`] in panic-sensitive paths so the modulus
    /// proof becomes a one-shot boundary check rather than a returned
    /// `Option<Self>` the caller must `.unwrap()`.
    ///
    /// [`new_odd`]: Self::new_odd
    pub fn new(modulus: T) -> Option<Self> {
        Odd::new(modulus).map(Self::new_odd)
    }

    /// Returns the modulus by reference.
    ///
    /// Returning `&T` rather than `T` avoids a memcpy of the full modulus
    /// (~256 bytes for a 2048-bit carrier) at the call site. Consumers that
    /// need a `T` by value can copy at the use point.
    pub fn modulus(&self) -> &T {
        &self.modulus
    }

    /// The additive identity (0 in Montgomery form is 0).
    pub fn zero(&self) -> Residue<'_, T, P> {
        Residue {
            mont: T::zero(),
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// The multiplicative identity (1 in Montgomery form is `R mod N`).
    pub fn one(&self) -> Residue<'_, T, P> {
        Residue {
            mont: self.r_mod_n,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Reconstruct a [`Residue`] from a raw value already in Montgomery form.
    ///
    /// **Escape hatch.** Intended for downstream specialization layers that
    /// persist or compute Montgomery-form values outside this module and
    /// need to re-attach the brand. The caller must guarantee `mont` is in
    /// `[0, modulus)` and represents some value `x` such that
    /// `mont == x * R mod modulus`.
    pub fn residue_from_mont(&self, mont: T) -> Residue<'_, T, P> {
        Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }
}

// ---------------------------------------------------------------------------
// Nct-only impls — variable-time CIOS with branching finalize
// ---------------------------------------------------------------------------

impl<T> Field<T, Nct>
where
    T: Copy
        + PartialEq
        + PartialOrd
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::WrappingMul
        + const_num_traits::WrappingAdd
        + const_num_traits::WrappingSub
        + const_num_traits::ops::overflowing::OverflowingAdd
        + core::ops::Add<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Mul<Output = T>
        + Parity
        + crate::NonCt
        + MontStorage,
{
    /// Convert a raw value `< modulus` (or arbitrary value, which is then
    /// reduced) to Montgomery form. Returns a brand-tagged [`Residue`].
    pub fn reduce(&self, raw: &T) -> Residue<'_, T, Nct>
    where
        T: WideMul,
    {
        let mont = wide_montgomery_mul(*raw, self.r2_mod_n, self.modulus, self.n_prime);
        Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Convert a [`Residue`] back to its raw `T` representative in `[0, modulus)`.
    #[allow(clippy::wrong_self_convention)]
    pub fn into_raw(&self, r: &Residue<'_, T, Nct>) -> T
    where
        T: WideMul,
    {
        wide_redc(r.mont, T::zero(), self.modulus, self.n_prime)
    }

    /// Modular addition: `(a + b) mod modulus`.
    ///
    /// **The conditional-subtract here is non-negotiable.** Future consumers
    /// with a modulus narrower than `T::BITS` may be tempted to skip the
    /// `wrapping_add` + cond-sub path (since `2 * modulus < 2^T::BITS` for
    /// them, the wraparound never fires). That's a correct optimization, but
    /// it belongs in a specialization layer (e.g. `Curve25519Field` in the
    /// ed25519 crate), not in `modmath::Field`. Patches removing the
    /// wrapping path will break RSA-CRT (full-width modulus, no slack) and
    /// any other consumer at full type width; ed25519 has slack and uses its
    /// own lazy variant in `Curve25519Field`.
    pub fn add(&self, a: &Residue<'_, T, Nct>, b: &Residue<'_, T, Nct>) -> Residue<'_, T, Nct> {
        let mont = crate::add::basic_mod_add_pr(a.mont, b.mont, self.modulus);
        Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Modular subtraction: `(a - b) mod modulus`.
    ///
    /// Same load-bearing contract as [`add`](Self::add) — the borrow-detect
    /// branch is required at full type width.
    pub fn sub(&self, a: &Residue<'_, T, Nct>, b: &Residue<'_, T, Nct>) -> Residue<'_, T, Nct> {
        let mont = crate::sub::basic_mod_sub_pr(a.mont, b.mont, self.modulus);
        Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Modular multiplication via CIOS Montgomery multiplication.
    ///
    /// CIOS interleaves multiplication and reduction in one pass (~2N² + N
    /// limb mults vs ~3N² for separate wide-mul + REDC), which dominates the
    /// inner-loop cost on constrained cores. The functional output is
    /// identical to `wide_montgomery_mul`.
    ///
    /// Marked `#[inline]` deliberately: this is the documented inner-loop
    /// wrapper for Montgomery exponentiation, the body is a single trait
    /// method call, and skipping it across crate boundaries costs ~250
    /// cycles per call under `opt-level="z"` on Cortex-M. Not blanket cargo
    /// culting — surgical on the actual hot path.
    #[inline]
    pub fn mul(&self, a: &Residue<'_, T, Nct>, b: &Residue<'_, T, Nct>) -> Residue<'_, T, Nct>
    where
        T: CiosMontMul,
    {
        let mont = CiosMontMul::cios_mont_mul(&a.mont, &b.mont, &self.modulus, &self.n_prime);
        Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Modular exponentiation via square-and-multiply.
    ///
    /// `base` is taken as a [`Residue`] (already in Montgomery form); `exp`
    /// is a raw `T`. The result is a [`Residue`] in Montgomery form.
    ///
    /// **Variable-time in `exp`.** The loop iterates `bit_length(exp)` times
    /// and branches on each bit. Do not call with a secret `exp` — use the
    /// Ct-variant `Field<T, Ct>::exp` instead.
    pub fn exp(&self, base: &Residue<'_, T, Nct>, exp: &T) -> Residue<'_, T, Nct>
    where
        T: CiosMontMul + core::ops::ShrAssign<usize>,
    {
        let mut result = self.r_mod_n;
        let mut base_var = base.mont;
        let mut exp_val = *exp;
        while exp_val > T::zero() {
            if exp_val.is_odd() {
                result =
                    CiosMontMul::cios_mont_mul(&result, &base_var, &self.modulus, &self.n_prime);
            }
            exp_val >>= 1;
            if exp_val > T::zero() {
                base_var =
                    CiosMontMul::cios_mont_mul(&base_var, &base_var, &self.modulus, &self.n_prime);
            }
        }
        Residue {
            mont: result,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Wide multiply-accumulate: `(acc.0, acc.1) += a.mont * b.mont`.
    ///
    /// Brand-tagged wrapper around [`wide_montgomery_mul_acc`]. Pair
    /// with [`Field::wide_redc`] to close the accumulator after a fused
    /// inner-product loop. See the free-function for the `N ≤ R/q`
    /// bound contract.
    pub fn mul_acc(&self, acc: (T, T), a: &Residue<'_, T, Nct>, b: &Residue<'_, T, Nct>) -> (T, T)
    where
        T: WideMul,
    {
        wide_montgomery_mul_acc(acc.0, acc.1, a.mont, b.mont)
    }

    /// Close a wide accumulator into a brand-tagged [`Residue`].
    pub fn wide_redc(&self, acc: (T, T)) -> Residue<'_, T, Nct>
    where
        T: WideMul,
    {
        let mont = wide_redc(acc.0, acc.1, self.modulus, self.n_prime);
        Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Modular inverse via Fermat's little theorem: `a^(modulus − 2)`.
    ///
    /// **Requires `modulus` to be prime.** Variable-time over the bits
    /// of `modulus − 2`. Returns `None` when `a` is the zero residue.
    pub fn inv_fermat(&self, a: &Residue<'_, T, Nct>) -> Option<Residue<'_, T, Nct>>
    where
        T: CiosMontMul + core::ops::ShrAssign<usize>,
    {
        if a.mont == T::zero() {
            return None;
        }
        let two = T::one().wrapping_add(T::one());
        let exp_val = self.modulus.wrapping_sub(two);
        Some(self.exp(a, &exp_val))
    }

    /// Modular inverse via extended Euclidean GCD on the raw Mont
    /// value, then rebrand to Mont form via two Mont multiplies by
    /// `R^2 mod N`.
    ///
    /// Works for any odd modulus (composite is fine). Variable-time —
    /// do not call with secret inputs; use [`Self::inv_fermat`] for CT
    /// paths. Returns `None` when `a` is not coprime to modulus.
    pub fn inv_eea(&self, a: &Residue<'_, T, Nct>) -> Option<Residue<'_, T, Nct>>
    where
        T: WideMul + core::ops::Div<Output = T> + core::ops::Sub<Output = T>,
    {
        if a.mont == T::zero() {
            return None;
        }
        let raw_inv = crate::inv::basic_mod_inv(a.mont, self.modulus)?;
        // raw_inv = (a*R)^{-1} = a^{-1} * R^{-1} (residue form).
        // Two Mont mults by R^2 lift it back to a^{-1} * R = Mont(a^{-1}).
        let step1 = wide_montgomery_mul(raw_inv, self.r2_mod_n, self.modulus, self.n_prime);
        let mont = wide_montgomery_mul(step1, self.r2_mod_n, self.modulus, self.n_prime);
        Some(Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        })
    }
}

// ---------------------------------------------------------------------------
// Ct-only impls — CT CIOS with conditional-select finalize
// ---------------------------------------------------------------------------

impl<'f, T> Residue<'f, T, Ct>
where
    T: subtle::ConditionallySelectable + MontStorage,
{
    /// Conditionally swap two residues in constant time.
    ///
    /// If `choice` is set, `a` and `b` exchange Montgomery-form values;
    /// otherwise both are left unchanged. The operation is branchless.
    ///
    /// This is the primitive used by Montgomery ladders (x25519 scalar
    /// multiplication, RSA blinded exponentiation). It is the **only**
    /// residue swap that should appear in such a ladder; `std::mem::swap`
    /// is not guaranteed to be branchless.
    pub fn cswap(choice: subtle::Choice, a: &mut Self, b: &mut Self) {
        T::conditional_swap(&mut a.mont, &mut b.mont, choice);
    }
}

impl<'f, T> Residue<'f, T, Ct>
where
    T: subtle::ConstantTimeEq + MontStorage,
{
    /// Constant-time equality on the underlying Montgomery values.
    ///
    /// Use in place of derived `PartialEq` on Ct paths where the
    /// equality outcome must not leak through timing (ML-KEM
    /// decapsulation tag check, ed25519 signature verification).
    pub fn ct_eq(&self, other: &Self) -> subtle::Choice {
        self.mont.ct_eq(&other.mont)
    }
}

impl<T> Field<T, Ct>
where
    T: Copy
        + PartialEq
        + PartialOrd
        + const_num_traits::Zero
        + const_num_traits::One
        + const_num_traits::WrappingMul
        + const_num_traits::WrappingAdd
        + const_num_traits::WrappingSub
        + const_num_traits::ops::overflowing::OverflowingAdd
        + core::ops::Add<Output = T>
        + core::ops::Sub<Output = T>
        + core::ops::Mul<Output = T>
        + Parity
        + MontStorage,
{
    /// Construct a `Field<T, Ct>` from a **secret** modulus without a
    /// value-dependent branch on the parity check.
    ///
    /// `Odd::new_ct` performs the parity check via [`CtParity`], producing a
    /// masked [`subtle::CtOption`] rather than a control-flow branch. The
    /// precompute (`compute_n_prime_newton`, `compute_r_mod_n`,
    /// `compute_r2_mod_n`) runs unconditionally — its inputs are the secret
    /// modulus's value, but the operations are constant-time word arithmetic
    /// over the existing CT trait surface, and the `CtOption` wrapper
    /// branchlessly masks the result if the modulus turned out to be even.
    ///
    /// Intended for the **RSA-CRT private-key path** where `p` and `q` are
    /// secret primes. Public-modulus / verify-side callers should use
    /// [`Field::new_odd`] instead — the secret-aware code path is strictly
    /// more expensive on platforms with branch prediction.
    ///
    /// Collapses the boundary check at the consumer:
    ///
    /// ```ignore
    /// // Old shape, panics on a secret-derived branch:
    /// let field = Field::<_, Ct>::new(secret_p).expect("p is odd prime");
    ///
    /// // New shape, masked:
    /// let field = Field::<_, Ct>::try_new_odd_ct(secret_p);
    /// let result = field.map(|f| /* CT-sensitive ops */ );
    /// ```
    ///
    /// [`CtParity`]: const_num_traits::CtParity
    pub fn try_new_odd_ct(modulus: T) -> subtle::CtOption<Self>
    where
        T: const_num_traits::CtParity
            + const_num_traits::CtIsZero
            + subtle::ConditionallySelectable
            + subtle::ConstantTimeLess,
    {
        // Mask the parity check (no branch on the secret modulus). The
        // precompute below uses the CT path ([`Self::new_odd_ct`]) so
        // no value-dependent branches on the modulus value either —
        // every step is `subtle::Choice`-masked. `CtOption::new(_,
        // choice)` discards the result via the standard masked-`Some`
        // pattern if the modulus turned out to be even.
        let is_odd = modulus.ct_is_odd();
        // SAFETY: when `is_odd` is unset the wrapped `Odd` proof carries a
        // false predicate, but the resulting `Field` is unreachable through
        // the `CtOption` mask. No body downstream consumes the proof except
        // via the masked output.
        let proof = unsafe { Odd::new_unchecked(modulus) };
        let field = Self::new_odd_ct(proof);
        subtle::CtOption::new(field, is_odd)
    }

    /// Convert a raw value to Montgomery form. Constant-time finalize.
    pub fn reduce(&self, raw: &T) -> Residue<'_, T, Ct>
    where
        T: WideMul
            + subtle::ConditionallySelectable
            + subtle::ConstantTimeLess
            + const_num_traits::CtIsZero,
    {
        let mont = wide_montgomery_mul_ct(*raw, self.r2_mod_n, self.modulus, self.n_prime);
        Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Convert a [`Residue`] back to raw form. Constant-time finalize.
    #[allow(clippy::wrong_self_convention)]
    pub fn into_raw(&self, r: &Residue<'_, T, Ct>) -> T
    where
        T: WideMul
            + subtle::ConditionallySelectable
            + subtle::ConstantTimeLess
            + const_num_traits::CtIsZero,
    {
        wide_redc_ct(r.mont, T::zero(), self.modulus, self.n_prime)
    }

    /// Modular addition — constant-time finalize.
    ///
    /// See `Field<T, Nct>::add` for the load-bearing comment about why the
    /// wrapping cond-sub path is non-negotiable.
    pub fn add(&self, a: &Residue<'_, T, Ct>, b: &Residue<'_, T, Ct>) -> Residue<'_, T, Ct>
    where
        T: subtle::ConditionallySelectable + subtle::ConstantTimeLess,
    {
        let sum = a.mont.wrapping_add(b.mont);
        let sub = sum.wrapping_sub(self.modulus);
        // Carry from wrapping: sum < a means wraparound occurred.
        let carry = sum.ct_lt(&a.mont);
        // Result >= modulus when !(sum < modulus).
        let ge_m = !sum.ct_lt(&self.modulus);
        let needs_sub = carry | ge_m;
        let mont = T::conditional_select(&sum, &sub, needs_sub);
        Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Modular subtraction — constant-time finalize.
    ///
    /// Same contract as `Field<T, Nct>::sub`.
    pub fn sub(&self, a: &Residue<'_, T, Ct>, b: &Residue<'_, T, Ct>) -> Residue<'_, T, Ct>
    where
        T: subtle::ConditionallySelectable + subtle::ConstantTimeLess,
    {
        let diff = a.mont.wrapping_sub(b.mont);
        let corrected = diff.wrapping_add(self.modulus);
        // borrow == (a < b)
        let borrow = a.mont.ct_lt(&b.mont);
        let mont = T::conditional_select(&diff, &corrected, borrow);
        Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Modular multiplication via CIOS — constant-time finalize.
    ///
    /// See `Field<T, Nct>::mul` for the rationale on CIOS vs. wide-REDC and
    /// the `#[inline]` justification.
    #[inline]
    pub fn mul(&self, a: &Residue<'_, T, Ct>, b: &Residue<'_, T, Ct>) -> Residue<'_, T, Ct>
    where
        T: CiosMontMulCt,
    {
        let mont = CiosMontMulCt::cios_mont_mul_ct(&a.mont, &b.mont, &self.modulus, &self.n_prime);
        Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Modular exponentiation — constant-time over `exp`.
    ///
    /// Implements a fixed-iteration Montgomery ladder over all
    /// `bit_length(T)` bits of the exponent. Both square and multiply are
    /// performed every iteration; the result is selected branchlessly. Loop
    /// count does not depend on `exp`; per-iteration timing does not depend
    /// on the bit pattern.
    pub fn exp(&self, base: &Residue<'_, T, Ct>, exp: &T) -> Residue<'_, T, Ct>
    where
        T: CiosMontMulCt
            + const_num_traits::CtIsZero
            + subtle::ConditionallySelectable
            + subtle::ConstantTimeEq
            + core::ops::Shr<usize, Output = T>
            + core::ops::BitAnd<Output = T>,
    {
        let w = type_bit_width::<T>();
        let one = T::one();
        let mut result = self.r_mod_n;

        for i in (0..w).rev() {
            // Always square.
            result =
                CiosMontMulCt::cios_mont_mul_ct(&result, &result, &self.modulus, &self.n_prime);
            // Always compute the conditional product.
            let multiplied =
                CiosMontMulCt::cios_mont_mul_ct(&result, &base.mont, &self.modulus, &self.n_prime);
            // Select based on bit i of exp.
            let bit_t = (*exp >> i) & one;
            let choice = !bit_t.ct_is_zero();
            result = T::conditional_select(&result, &multiplied, choice);
        }
        Residue {
            mont: result,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Modular exponentiation — constant-time over the base, **variable-time
    /// over the exponent**. Use when the exponent is public.
    ///
    /// This is the right primitive for several common cryptographic shapes:
    ///
    /// - **RSA encrypt / verify** — `m^e mod n` with the secret message `m`
    ///   and the public exponent `e` (typically 65537). Saves `bit_length(T)
    ///   - bit_length(e)` squarings vs. the fixed-iteration ladder, which is
    ///   ~2031 squarings at 2048-bit modulus when `e = 65537`.
    /// - **Curve25519 Fermat inverse** — `a^(p-2) mod p` where `p - 2` is the
    ///   curve constant `2^255 - 21`. The exponent is public; the base is
    ///   the secret intermediate `Z`. Skip-on-zero square-and-multiply
    ///   matches the ~252-of-255 bits set without spending the per-bit
    ///   `conditional_select` cost of the fixed-iteration ladder.
    /// - **Curve25519 square root** — `a^((p+3)/8) mod p`, same shape.
    ///
    /// The squarings and multiplications themselves go through CT primitives
    /// ([`cios_montgomery_mul_ct`](crate::montgomery::cios::cios_montgomery_mul_ct)),
    /// so the
    /// base and intermediate Montgomery values do not leak through timing.
    /// What DOES leak is the bit pattern of `exp` — which is fine by
    /// construction: the caller asserts the exponent is public.
    ///
    /// **Do not call with a secret exponent.** Use [`exp`](Self::exp)
    /// instead, which is a fixed-iteration Montgomery ladder.
    pub fn exp_public_exp(&self, base: &Residue<'_, T, Ct>, exp: &T) -> Residue<'_, T, Ct>
    where
        T: CiosMontMulCt + core::ops::Shr<usize, Output = T> + core::ops::BitAnd<Output = T>,
    {
        let w = type_bit_width::<T>();
        let one = T::one();
        let zero = T::zero();

        // Find the position of the highest set bit (1-indexed: hi == top + 1).
        // This loop and the rest of the function leak `bit_length(exp)`,
        // which is the documented contract — `exp` is public.
        let mut hi = w;
        while hi > 0 {
            if (*exp >> (hi - 1)) & one != zero {
                break;
            }
            hi -= 1;
        }

        if hi == 0 {
            // exp == 0: return 1 in Montgomery form.
            return Residue {
                mont: self.r_mod_n,
                _brand: PhantomData,
                _p: PhantomData,
            };
        }

        // The top bit is set, so result starts at `base` (base^1 contribution
        // for the 2^(hi-1) term). Then iterate over the remaining bits.
        let mut result = base.mont;
        for i in (0..hi - 1).rev() {
            // Square.
            result =
                CiosMontMulCt::cios_mont_mul_ct(&result, &result, &self.modulus, &self.n_prime);
            // Multiply only when the bit is set — branch on a public value.
            if (*exp >> i) & one != zero {
                result = CiosMontMulCt::cios_mont_mul_ct(
                    &result,
                    &base.mont,
                    &self.modulus,
                    &self.n_prime,
                );
            }
        }

        Residue {
            mont: result,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Wide multiply-accumulate (CT carry).
    ///
    /// Brand-tagged wrapper around [`wide_montgomery_mul_acc_ct`]. Pair
    /// with [`Field::wide_redc`] (CT variant) to close the accumulator.
    /// See the free-function for the `N ≤ R/q` bound contract.
    pub fn mul_acc(&self, acc: (T, T), a: &Residue<'_, T, Ct>, b: &Residue<'_, T, Ct>) -> (T, T)
    where
        T: WideMul + subtle::ConditionallySelectable + subtle::ConstantTimeLess,
    {
        wide_montgomery_mul_acc_ct(acc.0, acc.1, a.mont, b.mont)
    }

    /// Close a wide accumulator (CT finalize) into a brand-tagged
    /// [`Residue`].
    pub fn wide_redc(&self, acc: (T, T)) -> Residue<'_, T, Ct>
    where
        T: WideMul
            + subtle::ConditionallySelectable
            + subtle::ConstantTimeLess
            + const_num_traits::CtIsZero,
    {
        let mont = wide_redc_ct(acc.0, acc.1, self.modulus, self.n_prime);
        Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        }
    }

    /// Modular inverse via Fermat: `a^(modulus − 2)` through the fixed-
    /// iteration CT Montgomery ladder.
    ///
    /// **Requires `modulus` to be prime.** Constant-time over `a`'s
    /// bits and zero-ness via the fixed `T::BITS`-iteration Montgomery
    /// ladder in [`Self::exp`]. The loop count depends only on the
    /// carrier type's bit width, not on `modulus - 2`'s significant
    /// bit count or `a`'s value. Returns `CtOption::None`-masked for
    /// the zero residue.
    ///
    /// Cost: one full ladder over every bit of `T` (e.g. 256
    /// square-and-multiply iterations for a 256-bit carrier over a
    /// Curve25519 scalar field), regardless of whether `modulus - 2`
    /// occupies the full carrier width. For composite moduli (RSA
    /// `n = p·q`) where Fermat doesn't apply, use
    /// [`Self::inv_safegcd_ct`] instead.
    pub fn inv_fermat(&self, a: &Residue<'_, T, Ct>) -> subtle::CtOption<Residue<'_, T, Ct>>
    where
        T: CiosMontMulCt
            + const_num_traits::CtIsZero
            + subtle::ConditionallySelectable
            + subtle::ConstantTimeEq
            + core::ops::Shr<usize, Output = T>
            + core::ops::BitAnd<Output = T>,
    {
        let a_is_nonzero = !a.mont.ct_is_zero();
        let two = T::one().wrapping_add(T::one());
        let exp_val = self.modulus.wrapping_sub(two);
        let result = self.exp(a, &exp_val);
        subtle::CtOption::new(result, a_is_nonzero)
    }

    /// Constant-time modular inverse via Bernstein-Yang divsteps.
    /// **Works for any modulus** — composite (RSA `n = p·q`) or prime —
    /// unlike [`inv_fermat`] which assumes a prime modulus.
    ///
    /// Returns `CtOption::None` masked when `gcd(value, modulus) != 1`
    /// (no inverse exists). Failure timing is independent of input
    /// magnitudes.
    ///
    /// The modulus may occupy the full carrier width (MSB set in
    /// `T`) — a 2048-bit RSA modulus works in an exact 2048-bit
    /// carrier. The algorithm's signed intermediates are carried in
    /// an internally widened representation, so no headroom bit is
    /// required of `T`.
    ///
    /// Used by RSA private-key blinding, where the modulus is the
    /// composite `n = p·q` and Fermat's little theorem doesn't apply.
    /// See the `inv::safegcd` module source for the algorithm and
    /// full precondition list.
    ///
    /// [`inv_fermat`]: Self::inv_fermat
    pub fn inv_safegcd_ct(&self, a: &Residue<'_, T, Ct>) -> subtle::CtOption<Residue<'_, T, Ct>>
    where
        T: CiosMontMulCt
            + WideMul
            + subtle::ConditionallySelectable
            + subtle::ConstantTimeLess
            + const_num_traits::CtIsZero
            + modmath_cios::CiosRowOps
            + core::ops::Shr<usize, Output = T>
            + core::ops::Shl<usize, Output = T>
            + core::ops::BitOr<Output = T>,
        <T as modmath_cios::CiosRowOps>::Word: const_num_traits::CtParity,
    {
        // The value in the Residue is in Montgomery form. To get the
        // Montgomery form of the inverse:
        //   a.mont           = value · R mod n
        //   raw_inv          = safegcd(a.mont, n) = a.mont⁻¹ mod n
        //                    = (value · R)⁻¹ mod n
        //                    = value⁻¹ · R⁻¹ mod n
        //   wanted: inv.mont = value⁻¹ · R mod n
        //                    = raw_inv · R² mod n
        // Computing raw_inv · R² mod n via Mont multiplications requires
        // **two** multiplications by R², not one:
        //   m1 = REDC(raw_inv · R²) = raw_inv · R mod n  (= value⁻¹ raw)
        //   m2 = REDC(m1 · R²)      = m1 · R mod n       (= value⁻¹ · R = inv.mont)
        // The first multiplication "converts raw_inv into something that
        // multiplied by R again gives the desired Mont form". The
        // second multiplication does that final · R step. Equivalent
        // to one multiplication by R³, but we only have R² cached.
        let inv_raw = crate::inv::safegcd::safegcd_inv_ct(&a.mont, &self.modulus);
        // Extract the raw inverse, defaulting to zero when safegcd
        // reports `None`. The two REDCs run unconditionally on the
        // extracted value — under the CtOption mask any garbage they
        // produce on the failure path is discarded.
        let inv_exists = inv_raw.is_some();
        let raw_inv = inv_raw.unwrap_or(T::zero());
        let m1 = wide_montgomery_mul_ct(raw_inv, self.r2_mod_n, self.modulus, self.n_prime);
        let mont = wide_montgomery_mul_ct(m1, self.r2_mod_n, self.modulus, self.n_prime);
        let residue = Residue {
            mont,
            _brand: PhantomData,
            _p: PhantomData,
        };
        subtle::CtOption::new(residue, inv_exists)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use fixed_bigint::FixedUInt;
    use subtle::Choice;
    #[cfg(feature = "zeroize")]
    use zeroize::Zeroize;

    // Field<T, P> requires the right combination of T-bounds for the chosen
    // P (CiosMontMul for Nct, CiosMontMulCt for Ct), which in practice means
    // T must be a FixedUInt of the matching personality. Small tests use
    // FixedUInt<u8, 2> aliases for tight ranges; larger tests use the U128
    // family. Cross-personality tests bridge values via `.into()` (Nct → Ct)
    // and `.forget_ct()` (explicit Ct → Nct).
    type U16 = FixedUInt<u8, 2>;
    type U16Ct = FixedUInt<u8, 2, Ct>;
    type U128Ct = FixedUInt<u32, 4, Ct>;

    fn u16(n: u16) -> U16 {
        U16::from(n)
    }

    fn u16ct(n: u16) -> U16Ct {
        U16Ct::from(n)
    }

    #[test]
    fn round_trip_small() {
        let f: Field<U16> = Field::new(u16(13)).unwrap();
        for raw in 0u16..13 {
            let r = f.reduce(&u16(raw));
            assert_eq!(f.into_raw(&r), u16(raw), "round trip failed for {raw}");
        }
    }

    #[test]
    fn new_odd_matches_new() {
        // The infallible Odd-typestate constructor and the runtime-checked
        // `Option`-returning one must agree on the precompute (modulus,
        // n_prime, r_mod_n, r2_mod_n) for the same modulus value.
        let m = u16(13);
        let modulus_odd = Odd::new(m).expect("13 is odd");
        let from_odd: Field<U16> = Field::new_odd(modulus_odd);
        let from_opt: Field<U16> = Field::new(m).unwrap();
        assert_eq!(from_odd.modulus(), from_opt.modulus());
        // Round-trip through Field::mul under each to confirm the precompute
        // tables match observably.
        let a = from_odd.reduce(&u16(7));
        let b = from_odd.reduce(&u16(5));
        let via_odd = from_odd.into_raw(&from_odd.mul(&a, &b));
        let a2 = from_opt.reduce(&u16(7));
        let b2 = from_opt.reduce(&u16(5));
        let via_opt = from_opt.into_raw(&from_opt.mul(&a2, &b2));
        assert_eq!(via_odd, via_opt);
        assert_eq!(via_odd, u16(35 % 13));
    }

    #[test]
    fn new_rejects_even_and_zero() {
        // Wrapper preserves the rejection semantics of the old API.
        assert!(Field::<U16>::new(u16(0)).is_none());
        assert!(Field::<U16>::new(u16(12)).is_none()); // even
        assert!(Field::<U16>::new(u16(13)).is_some()); // odd
    }

    /// `try_new_odd_ct` produces a `CtOption<Field<T, Ct>>` whose
    /// `Some`-ness tracks `T::ct_is_odd`. The precompute runs
    /// unconditionally; the parity check is masked, not branched. Test
    /// on `u32` (which impls `CtParity` directly) since that's the
    /// straightforward case — the RSA-CRT consumer pattern will be on a
    /// bigint type, but the contract we're pinning here is the
    /// modmath-side adapter.
    #[test]
    fn try_new_odd_ct_masks_parity() {
        // Even modulus → `None`-masked.
        let even = Field::<u32, Ct>::try_new_odd_ct(12);
        assert_eq!(even.is_some().unwrap_u8(), 0);

        // Zero is even → `None`-masked.
        let zero = Field::<u32, Ct>::try_new_odd_ct(0);
        assert_eq!(zero.is_some().unwrap_u8(), 0);

        // Odd modulus → `Some` with a usable Field.
        let odd = Field::<u32, Ct>::try_new_odd_ct(13);
        assert_eq!(odd.is_some().unwrap_u8(), 1);
        let field: Field<u32, Ct> = odd.unwrap();
        // Same precompute as the infallible boundary constructor:
        let baseline = Field::<u32, Ct>::new_odd(Odd::new(13u32).unwrap());
        assert_eq!(field.modulus(), baseline.modulus());
    }

    /// `Field::new_odd_ct` (the CT precompute path) must produce
    /// identical precompute values to `Field::new_odd` (the
    /// variable-time path) for every modulus. Pins the contract that
    /// `mod_double_ct` / `mod_exp2_ct` are CT-equivalent, not just
    /// "CT but different output."
    #[test]
    fn new_odd_ct_precompute_matches_new_odd() {
        for m in [3u32, 5, 7, 11, 13, 97, 65521, 0x7FFF_FFE7] {
            let modulus = Odd::new(m).unwrap();
            let f_nct = Field::<u32, Ct>::new_odd(modulus);
            let f_ct = Field::<u32, Ct>::new_odd_ct(modulus);
            assert_eq!(f_nct.modulus(), f_ct.modulus(), "modulus mismatch at m={m}");
            assert_eq!(f_nct.n_prime, f_ct.n_prime, "n_prime mismatch at m={m}");
            assert_eq!(f_nct.r_mod_n, f_ct.r_mod_n, "r_mod_n mismatch at m={m}");
            assert_eq!(f_nct.r2_mod_n, f_ct.r2_mod_n, "r2_mod_n mismatch at m={m}");
        }
    }

    /// CT precompute on multi-limb FixedUInt produces identical
    /// output to the variable-time precompute. The actual RSA-CRT
    /// shape — the precompute is what would silently produce
    /// wrong results if `mod_double_ct` had a bug.
    #[test]
    fn new_odd_ct_precompute_matches_new_odd_fixed_bigint() {
        // 128-bit odd modulus (composite, RSA-CRT-shape)
        let m = U128Ct::from(0xFFFF_FFFF_FFFF_FFE7u64);
        let modulus = Odd::new(m).unwrap();
        let f_nct = Field::<U128Ct, Ct>::new_odd(modulus);
        let f_ct = Field::<U128Ct, Ct>::new_odd_ct(modulus);
        assert_eq!(f_nct.modulus(), f_ct.modulus());
        assert_eq!(f_nct.n_prime, f_ct.n_prime);
        assert_eq!(f_nct.r_mod_n, f_ct.r_mod_n);
        assert_eq!(f_nct.r2_mod_n, f_ct.r2_mod_n);
    }

    #[test]
    fn add_sub_mul_small() {
        let f: Field<U16> = Field::new(u16(13)).unwrap();
        for a_raw in 0u16..13 {
            for b_raw in 0u16..13 {
                let a = f.reduce(&u16(a_raw));
                let b = f.reduce(&u16(b_raw));
                assert_eq!(f.into_raw(&f.add(&a, &b)), u16((a_raw + b_raw) % 13));
                assert_eq!(
                    f.into_raw(&f.sub(&a, &b)),
                    u16((a_raw + 13 - b_raw) % 13),
                    "sub failed for {a_raw}, {b_raw}"
                );
                assert_eq!(f.into_raw(&f.mul(&a, &b)), u16((a_raw * b_raw) % 13));
            }
        }
    }

    #[test]
    fn zero_one_identity_small() {
        let f: Field<U16> = Field::new(u16(13)).unwrap();
        let z = f.zero();
        let o = f.one();
        assert_eq!(f.into_raw(&z), u16(0));
        assert_eq!(f.into_raw(&o), u16(1));
        // a + 0 = a, a * 1 = a
        for raw in 0u16..13 {
            let a = f.reduce(&u16(raw));
            assert_eq!(f.into_raw(&f.add(&a, &z)), u16(raw));
            assert_eq!(f.into_raw(&f.mul(&a, &o)), u16(raw));
        }
    }

    #[test]
    fn exp_small() {
        let f: Field<U16> = Field::new(u16(13)).unwrap();
        // 7^5 mod 13 = 11
        let base = f.reduce(&u16(7));
        let result = f.exp(&base, &u16(5));
        assert_eq!(f.into_raw(&result), u16(11));
        // x^0 = 1
        let r0 = f.exp(&base, &u16(0));
        assert_eq!(f.into_raw(&r0), u16(1));
    }

    #[test]
    fn ct_round_trip_small() {
        let f = FieldCt::new(u16ct(13)).unwrap();
        for raw in 0u16..13 {
            let r = f.reduce(&u16ct(raw));
            assert_eq!(f.into_raw(&r), u16ct(raw));
        }
    }

    #[test]
    fn ct_matches_nct_small() {
        let f: Field<U16> = Field::new(u16(13)).unwrap();
        let fc = FieldCt::new(u16ct(13)).unwrap();
        for a_raw in 0u16..13 {
            for b_raw in 0u16..13 {
                let a = f.reduce(&u16(a_raw));
                let b = f.reduce(&u16(b_raw));
                let ac = fc.reduce(&u16ct(a_raw));
                let bc = fc.reduce(&u16ct(b_raw));

                assert_eq!(
                    f.into_raw(&f.add(&a, &b)),
                    fc.into_raw(&fc.add(&ac, &bc)).forget_ct()
                );
                assert_eq!(
                    f.into_raw(&f.sub(&a, &b)),
                    fc.into_raw(&fc.sub(&ac, &bc)).forget_ct()
                );
                assert_eq!(
                    f.into_raw(&f.mul(&a, &b)),
                    fc.into_raw(&fc.mul(&ac, &bc)).forget_ct()
                );
            }
        }
        // exp cross-check
        let base = f.reduce(&u16(7));
        let base_ct = fc.reduce(&u16ct(7));
        for e in 0u16..20 {
            assert_eq!(
                f.into_raw(&f.exp(&base, &u16(e))),
                fc.into_raw(&fc.exp(&base_ct, &u16ct(e))).forget_ct()
            );
        }
    }

    #[test]
    fn ct_cswap_small() {
        let f = FieldCt::new(u16ct(13)).unwrap();
        let mut a = f.reduce(&u16ct(3));
        let mut b = f.reduce(&u16ct(7));
        // choice = 0: no swap
        ResidueCt::cswap(Choice::from(0), &mut a, &mut b);
        assert_eq!(f.into_raw(&a), u16ct(3));
        assert_eq!(f.into_raw(&b), u16ct(7));
        // choice = 1: swap
        ResidueCt::cswap(Choice::from(1), &mut a, &mut b);
        assert_eq!(f.into_raw(&a), u16ct(7));
        assert_eq!(f.into_raw(&b), u16ct(3));
    }

    /// Under personality, the safe NCT → CT bridge requires converting the
    /// underlying T's personality (free `.into()` from fixed-bigint), then
    /// constructing a fresh `Field<_, Ct>` on the Ct-typed modulus. Same-T
    /// `Field<T, Nct> -> Field<T, Ct>` conversion is degenerate (the per-P
    /// impl blocks have disjoint trait bounds, so the result is a methodless
    /// variant); this test documents the actual bridge pattern.
    #[test]
    fn nct_to_ct_upgrade_small() {
        let f: Field<U16> = Field::new(u16(13)).unwrap();
        let modulus_ct: U16Ct = (*f.modulus()).into();
        let fc = FieldCt::new(modulus_ct).unwrap();
        let a = fc.reduce(&u16ct(7));
        let b = fc.reduce(&u16ct(5));
        assert_eq!(fc.into_raw(&fc.mul(&a, &b)), u16ct(9)); // 35 mod 13 = 9
    }

    #[test]
    fn exp_public_exp_matches_ct_exp_small() {
        // For every (base, exp) pair, exp_public_exp must produce the same
        // result as the fixed-iteration ladder exp.
        let f = FieldCt::new(u16ct(13)).unwrap();
        let base = f.reduce(&u16ct(7));
        for e in 0u16..32 {
            let via_ladder = f.exp(&base, &u16ct(e));
            let via_pub = f.exp_public_exp(&base, &u16ct(e));
            assert_eq!(
                f.into_raw(&via_ladder),
                f.into_raw(&via_pub),
                "exp_public_exp mismatch at e={e}"
            );
        }
    }

    #[test]
    fn exp_public_exp_matches_ct_exp_u128() {
        // Same cross-check at FixedUInt<u32, 4> sizes against a few
        // characteristic exponents: 0, 1, small, a value with both low and
        // high set bits.
        let modulus = !U128Ct::from(0u64) - U128Ct::from(58u64);
        let f = FieldCt::new(modulus).unwrap();
        let base = f.reduce(&U128Ct::from(0xDEAD_BEEF_u64));
        let exps = [
            U128Ct::from(0u64),
            U128Ct::from(1u64),
            U128Ct::from(7u64),
            U128Ct::from(65537u64), // RSA-style public exponent
            U128Ct::from(0xCAFE_BABEu64),
        ];
        for e in &exps {
            let via_ladder = f.exp(&base, e);
            let via_pub = f.exp_public_exp(&base, e);
            assert_eq!(
                f.into_raw(&via_ladder),
                f.into_raw(&via_pub),
                "exp_public_exp mismatch at e={e:?}"
            );
        }
    }

    #[test]
    fn brand_round_trip_fixed_bigint_u128() {
        // A larger odd modulus.
        let modulus = !U128Ct::from(0u64) - U128Ct::from(58u64);
        let f = FieldCt::new(modulus).unwrap();
        let raw = U128Ct::from(0xDEAD_BEEF_u64);
        let r = f.reduce(&raw);
        assert_eq!(f.into_raw(&r), raw);
    }

    /// `inv_safegcd_ct` round-trip on a prime modulus. The CT
    /// composite-modulus inverse is the load-bearing primitive for RSA
    /// blinding; here we test on a prime (smaller test surface) and
    /// verify `inv * value ≡ 1 mod modulus`.
    #[test]
    fn inv_safegcd_ct_round_trip_prime_modulus() {
        let f = FieldCt::new(u16ct(13)).unwrap();
        for raw_val in 1u16..13 {
            let r = f.reduce(&u16ct(raw_val));
            let inv = f.inv_safegcd_ct(&r);
            assert_eq!(
                inv.is_some().unwrap_u8(),
                1,
                "expected inverse for {raw_val} mod 13"
            );
            let inv_residue = inv.unwrap();
            let product = f.mul(&r, &inv_residue);
            assert_eq!(
                f.into_raw(&product),
                u16ct(1),
                "{raw_val} * inv != 1 mod 13"
            );
        }
    }

    /// `inv_safegcd_ct` on a composite modulus — the RSA blinding case.
    /// Confirms the algorithm works when the modulus is `p·q`, not
    /// prime, where Fermat inversion would fail.
    #[test]
    fn inv_safegcd_ct_composite_modulus() {
        // n = 3 * 5 = 15. Coprime values: 1, 2, 4, 7, 8, 11, 13, 14.
        let f = FieldCt::new(u16ct(15)).unwrap();
        for &raw_val in &[1u16, 2, 4, 7, 8, 11, 13, 14] {
            let r = f.reduce(&u16ct(raw_val));
            let inv = f.inv_safegcd_ct(&r);
            assert_eq!(
                inv.is_some().unwrap_u8(),
                1,
                "expected inverse for {raw_val} mod 15"
            );
            let product = f.mul(&r, &inv.unwrap());
            assert_eq!(
                f.into_raw(&product),
                u16ct(1),
                "{raw_val} * inv != 1 mod 15"
            );
        }
        // Non-coprime values: safegcd returns None.
        for &raw_val in &[3u16, 5, 6, 9, 10, 12] {
            let r = f.reduce(&u16ct(raw_val));
            let inv = f.inv_safegcd_ct(&r);
            assert_eq!(
                inv.is_some().unwrap_u8(),
                0,
                "expected None for non-coprime {raw_val} mod 15"
            );
        }
    }

    /// `inv_safegcd_ct` with a modulus that occupies the full carrier
    /// width (MSB set) — the exact-width-carrier case (a 2048-bit RSA
    /// modulus in a 2048-bit `T`), scaled down to a 16-bit carrier.
    #[test]
    fn inv_safegcd_ct_full_width_modulus() {
        // U16Ct = FixedUInt<u8, 2, Ct> — 16-bit carrier holding the
        // odd 16-bit modulus 0xFFFD = 13 · 71² (MSB set, composite).
        let modulus = u16ct(0xFFFD);
        let f = FieldCt::new(modulus).unwrap();
        for raw_val in [1u16, 2, 7, 0xBEEF, 0xFFFC] {
            let r = f.reduce(&u16ct(raw_val));
            let inv = f.inv_safegcd_ct(&r);
            assert_eq!(
                inv.is_some().unwrap_u8(),
                1,
                "expected inverse for {raw_val:#x} mod 0xFFFD"
            );
            let product = f.mul(&r, &inv.unwrap());
            assert_eq!(
                f.into_raw(&product),
                u16ct(1),
                "{raw_val:#x} * inv != 1 mod 0xFFFD"
            );
        }
        // Non-coprime (shares factor 13) still masks to None.
        let r = f.reduce(&u16ct(13));
        assert_eq!(f.inv_safegcd_ct(&r).is_some().unwrap_u8(), 0);
    }

    /// `inv_safegcd_ct` on a larger RSA-CRT-shaped composite modulus.
    /// n = p · q with small primes p, q. Confirms the algorithm runs
    /// correctly on a multi-limb FixedUInt and at sizes more
    /// representative of the RSA blinding workload than the toy
    /// `mod 15` case (still small enough that we can exhaustively
    /// check inv * value ≡ 1).
    #[test]
    fn inv_safegcd_ct_composite_modulus_u128() {
        // n = (2^32 + 7) · (2^24 + 7) — RSA-CRT-shape two-prime
        // composite, ~52 bits. safegcd handles composites; the result
        // works for any coprime value.
        let n_raw: u64 = 0x1_0000_0007 * 0x100_0007u64; // = 4503599644606465
        let modulus = U128Ct::from(n_raw);
        let f = FieldCt::new(modulus).unwrap();

        // A handful of values coprime to n. (0xDEAD_BEEF deliberately
        // omitted — it shares factor 11 with this n.)
        let test_vals = [
            U128Ct::from(1u64),
            U128Ct::from(2u64),
            U128Ct::from(3u64),
            U128Ct::from(0xCAFE_BABEu64),
            U128Ct::from(0xFEED_FACEu64),
        ];
        for v in test_vals {
            let r = f.reduce(&v);
            let inv = f.inv_safegcd_ct(&r);
            assert_eq!(
                inv.is_some().unwrap_u8(),
                1,
                "expected inverse to exist for v={:?}",
                v
            );
            let product = f.mul(&r, &inv.unwrap());
            assert_eq!(f.into_raw(&product), U128Ct::from(1u64));
        }
    }

    #[cfg(feature = "zeroize")]
    #[test]
    fn residue_zeroize_wipes_mont_small() {
        fn assert_zeroize_on_drop<T: zeroize::ZeroizeOnDrop>(_: &T) {}
        let f = FieldCt::new(u16ct(13)).unwrap();
        let mut r = f.reduce(&u16ct(7));
        assert_zeroize_on_drop(&r);
        assert_ne!(*r.mont_value(), u16ct(0));
        r.zeroize();
        assert_eq!(*r.mont_value(), u16ct(0));
    }

    #[test]
    fn residue_from_mont_escape_hatch_small() {
        // Round-trip via the escape hatch: reduce -> mont_value -> residue_from_mont.
        let f: Field<U16> = Field::new(u16(13)).unwrap();
        for raw in 0u16..13 {
            let r = f.reduce(&u16(raw));
            let mont = *r.mont_value();
            let r2 = f.residue_from_mont(mont);
            assert_eq!(f.into_raw(&r2), u16(raw));
        }
    }

    /// Documented limitation: covariance allows mixing residues across two
    /// distinct Field instances built in the same scope. Asserts current
    /// behavior (the compiler does NOT reject this) so a future generative
    /// brand can be observed as a hardening change.
    #[test]
    fn covariance_mixes_residues_documented_limitation() {
        let f1: Field<U16> = Field::new(u16(13)).unwrap();
        let f2: Field<U16> = Field::new(u16(13)).unwrap();
        let r1 = f1.reduce(&u16(5));
        // f2 accepting r1 compiles today. This is a documented limitation; a
        // generative brand would make this a type error.
        let _ = f2.into_raw(&r1);
    }

    /// Personality demonstration: the same `Field` type signature
    /// parameterized differently (`<_, Nct>` vs `<_, Ct>`) computes the
    /// same modular arithmetic, with the personality choice driving which
    /// algorithm (variable-time branch vs CT conditional-select finalize)
    /// the compiler routes to via the per-P impl blocks.
    ///
    /// Also exercises the residue type discipline: a `Residue<_, _, Nct>`
    /// passed to a `Field<_, Ct>` method would be a compile error
    /// (different `P` parameter), and vice versa. Cross-personality
    /// comparison goes through `.forget_ct()` rather than a same-type
    /// `assert_eq!`.
    #[test]
    fn field_p_personality_cross_check_small() {
        // Same modulus value, two personalities.
        let m_nct = u16(13);
        let m_ct: U16Ct = m_nct.into();

        let f_nct: Field<U16, Nct> = Field::new(m_nct).unwrap();
        let f_ct: Field<U16Ct, Ct> = Field::new(m_ct).unwrap();

        // Pick a non-trivial product and exponentiation.
        let a_nct = f_nct.reduce(&u16(7));
        let b_nct = f_nct.reduce(&u16(5));
        let a_ct = f_ct.reduce(&u16ct(7));
        let b_ct = f_ct.reduce(&u16ct(5));

        // Multiplication agrees across personalities.
        let mul_nct = f_nct.into_raw(&f_nct.mul(&a_nct, &b_nct));
        let mul_ct = f_ct.into_raw(&f_ct.mul(&a_ct, &b_ct));
        assert_eq!(mul_nct, mul_ct.forget_ct());

        // Exponentiation agrees (with different algorithms underneath:
        // f_nct.exp is variable-time square-and-multiply, f_ct.exp is
        // fixed-iteration ladder).
        let exp_nct = f_nct.into_raw(&f_nct.exp(&a_nct, &u16(11)));
        let exp_ct = f_ct.into_raw(&f_ct.exp(&a_ct, &u16ct(11)));
        assert_eq!(exp_nct, exp_ct.forget_ct());
    }

    /// The `FieldNct<T>` alias side-steps the construction-site type-
    /// ambiguity that bare `Field::new(modulus)` hits — no type annotation
    /// or turbofish required, because the alias fixes `P = Nct` at the
    /// type level (mirroring how `FieldCt::new` fixes `P = Ct`).
    ///
    /// Symmetric `ResidueNct` alias also exists for downstream consumers
    /// who want symmetric naming. Used here just for the type spelling.
    #[test]
    fn field_nct_alias_resolves_without_annotation() {
        let f = FieldNct::new(u16(13)).unwrap();
        let r: ResidueNct<'_, U16> = f.reduce(&u16(7));
        assert_eq!(f.into_raw(&r), u16(7));
        let two = f.reduce(&u16(2));
        assert_eq!(f.into_raw(&f.mul(&r, &two)), u16(14 % 13));
    }

    /// `from_precomputed` is `const fn` and usable in a const initializer.
    /// This is the constructor static-modulus consumers (PQC, embedded RSA
    /// with a baked key, etc.) reach for when they want to expose a `Field`
    /// as a `const` associated item rather than paying the runtime
    /// `Field::new` precompute.
    ///
    /// Demonstrated here over `u32` (primitive) — `Field<u32, Nct>` is
    /// methodless because `u32` doesn't impl `CiosMontMul` (MulAccOps is
    /// FixedUInt-only), but `from_precomputed` itself works for any
    /// `T: Copy`. The intended consumer path is a downstream Mont-newtype
    /// wrapper that calls modmath's standalone `wide_montgomery_mul[_ct]`
    /// free functions, using `f.modulus()` to read the static modulus.
    #[test]
    fn from_precomputed_const_construction_u32() {
        // Hand-computed Montgomery params for modulus 13 at word width 32:
        //   n_prime  = -13^-1 mod 2^32 = 0x4EC4EC4F
        //   r_mod_n  = 2^32 mod 13     = 9
        //   r2_mod_n = (2^32)^2 mod 13 = 3
        const F: Field<u32, Nct> = Field::from_precomputed(13u32, 0x4EC4EC4F, 9, 3);
        assert_eq!(*F.modulus(), 13u32);
        // The struct fields are accessible to anyone in the same crate
        // through Field::modulus(); downstream consumers driving the Mont
        // newtype pattern will pull modulus + n_prime + r/r2 via a Modulus
        // trait extension on their own side and call modmath's standalone
        // primitives. This test just proves the const-context construction
        // path is real.
    }

    /// `Field::mul_acc` + `Field::wide_redc` from a zero accumulator
    /// must equal `Field::mul` on the same operands.
    #[test]
    fn field_mul_acc_round_trip_small() {
        let f: Field<U16> = Field::new(u16(13)).unwrap();
        for a_raw in 0u16..13 {
            for b_raw in 0u16..13 {
                let a = f.reduce(&u16(a_raw));
                let b = f.reduce(&u16(b_raw));
                let direct = f.mul(&a, &b);
                let via_acc = f.wide_redc(f.mul_acc((u16(0), u16(0)), &a, &b));
                assert_eq!(f.into_raw(&direct), f.into_raw(&via_acc));
            }
        }
    }

    /// Dot product through `Field::mul_acc` + single `Field::wide_redc`
    /// must equal the direct residue-domain sum of products.
    #[test]
    fn field_mul_acc_dot_product_small() {
        let f: Field<U16> = Field::new(u16(13)).unwrap();
        let pairs: &[(u16, u16)] = &[(2, 3), (5, 7), (11, 4), (1, 12)];
        let mut acc = (u16(0), u16(0));
        for &(a_raw, b_raw) in pairs {
            let a = f.reduce(&u16(a_raw));
            let b = f.reduce(&u16(b_raw));
            acc = f.mul_acc(acc, &a, &b);
        }
        let result = f.wide_redc(acc);
        let expected: u16 = pairs
            .iter()
            .fold(0u16, |s, &(a, b)| (s + (a * b) % 13) % 13);
        assert_eq!(f.into_raw(&result), u16(expected));
    }

    /// `a * inv_fermat(a) == 1` for every nonzero residue at prime
    /// modulus 13; zero returns `None`.
    #[test]
    fn field_inv_fermat_small() {
        let f: Field<U16> = Field::new(u16(13)).unwrap();
        for raw in 1u16..13 {
            let a = f.reduce(&u16(raw));
            let inv = f.inv_fermat(&a).unwrap();
            assert_eq!(
                f.into_raw(&f.mul(&a, &inv)),
                u16(1),
                "fermat fails at {raw}"
            );
        }
        assert!(f.inv_fermat(&f.zero()).is_none());
    }

    /// Same contract as inv_fermat but via EEA path; cross-checks the
    /// two methods agree at every nonzero residue.
    #[test]
    fn field_inv_eea_small() {
        let f: Field<U16> = Field::new(u16(13)).unwrap();
        for raw in 1u16..13 {
            let a = f.reduce(&u16(raw));
            let inv_e = f.inv_eea(&a).unwrap();
            let inv_f = f.inv_fermat(&a).unwrap();
            assert_eq!(f.into_raw(&f.mul(&a, &inv_e)), u16(1), "eea fails at {raw}");
            assert_eq!(
                f.into_raw(&inv_e),
                f.into_raw(&inv_f),
                "fermat/eea disagree at {raw}"
            );
        }
        assert!(f.inv_eea(&f.zero()).is_none());
    }

    /// Ct variant of `Field::mul_acc` + `Field::wide_redc` must agree
    /// with `Field::mul`.
    #[test]
    fn field_mul_acc_ct_round_trip_small() {
        let f = FieldCt::new(u16ct(13)).unwrap();
        for a_raw in 0u16..13 {
            for b_raw in 0u16..13 {
                let a = f.reduce(&u16ct(a_raw));
                let b = f.reduce(&u16ct(b_raw));
                let direct = f.mul(&a, &b);
                let via_acc = f.wide_redc(f.mul_acc((u16ct(0), u16ct(0)), &a, &b));
                assert_eq!(f.into_raw(&direct), f.into_raw(&via_acc));
            }
        }
    }

    /// Ct `inv_fermat` must satisfy `a * inv(a) == 1` at prime modulus.
    #[test]
    fn field_inv_fermat_ct_small() {
        let f = FieldCt::new(u16ct(13)).unwrap();
        for raw in 1u16..13 {
            let a = f.reduce(&u16ct(raw));
            let inv = f.inv_fermat(&a).into_option().unwrap();
            assert_eq!(
                f.into_raw(&f.mul(&a, &inv)),
                u16ct(1),
                "ct fermat fails at {raw}"
            );
        }
        assert!(f.inv_fermat(&f.zero()).into_option().is_none());
    }

    /// `ResidueCt::ct_eq` matches `PartialEq` outcomes on representative
    /// inputs (true and false cases).
    #[test]
    fn residue_ct_eq_small() {
        let f = FieldCt::new(u16ct(13)).unwrap();
        let a = f.reduce(&u16ct(7));
        let b = f.reduce(&u16ct(7));
        let c = f.reduce(&u16ct(8));
        let eq_ab: bool = a.ct_eq(&b).into();
        let eq_ac: bool = a.ct_eq(&c).into();
        assert!(eq_ab);
        assert!(!eq_ac);
    }
}
