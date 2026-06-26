Review: typestate-wrapper opportunities in modmath-rs

I did not make code changes, commit, or open a PR, because this was a review/analysis request rather than an implementation request.
Commands used

    ✅ pwd && rg --files -g 'AGENTS.md' -g 'Cargo.toml' -g '*.rs'

    ✅ find .. -name AGENTS.md -print

    ✅ sed -n '1,220p' modmath/src/lib.rs && sed -n '1,240p' modmath/src/parity.rs && sed -n '1,240p' modmath/src/inv.rs

    ✅ sed -n '1,260p' modmath/src/mul.rs && sed -n '1,260p' modmath/src/wide_mul.rs && sed -n '1,220p' modmath/src/add.rs && sed -n '1,220p' modmath/src/sub.rs

    ✅ sed -n '1,260p' modmath/src/exp.rs && sed -n '1,260p' modmath/src/montgomery/basic_mont.rs && sed -n '1,260p' modmath/src/montgomery/strict_mont.rs

    ✅ rg -n "pub fn .*(_pr|mont|redc|params|inv|exp|mul|add|sub)|if .*None|return None|is_even|is_odd|modulus ==|while .*<=|Option|Result" modmath/src

    ✅ sed -n '220,660p' modmath/src/field.rs && sed -n '760,1010p' modmath/src/montgomery/basic_mont.rs && sed -n '520,750p' modmath/src/montgomery/basic_mont.rs

    ✅ nl -ba modmath/src/{add.rs,sub.rs,mul.rs,exp.rs,inv.rs,parity.rs,field.rs} | sed -n '1,260p'

    ✅ nl -ba modmath/src/montgomery/basic_mont.rs | sed -n '169,990p'

    ✅ nl -ba modmath/src/montgomery/cios.rs | sed -n '1,210p'

High-level finding

This crate already has an informal typestate split in a few places:

    _pr functions encode “pre-reduced” as a naming convention and documentation precondition, but not as a type. For example, basic_mod_add_pr, constrained_mod_add_pr, and strict_mod_add_pr require operands already below the modulus and therefore avoid %/Rem bounds.

    Field<T, P> and Residue<'f, T, P> already encode Montgomery form, field association, and constant-time personality at the type level.

    Field::new still validates “odd and nonzero modulus” dynamically, returning Option<Self>.

    CIOS Montgomery multiplication documents that its inputs must already be Montgomery-form residues in [0, modulus), but the lower-level function still accepts raw &T and returns Option<T> because limb access can theoretically fail.

The largest opportunity is to turn these documented preconditions into lightweight wrappers:

NonZero<T>
Odd<T>
Reduced<'m, T>
MontForm<'m, T, P>
Wide<T>

That would let the crate expose faster, narrower APIs without repeatedly proving the same facts at runtime.
1. NonZero<T> and Odd<T> for Montgomery moduli
Current cost / issue

Montgomery APIs repeatedly check that the modulus is nonzero and odd:

if modulus == T::zero() || modulus.is_even() {
    return None;
}

This pattern appears in the wide-REDC Montgomery mod_mul, mod_exp, CT exponentiation, and field construction paths.

This also forces Option<T> returns in places where, after construction, failure is impossible.
Typestate opportunity

Introduce something like:

pub struct NonZero<T>(T);
pub struct Odd<T>(NonZero<T>);
pub struct MontModulus<T> {
    n: Odd<T>,
    n_prime: T,
    r_mod_n: T,
    r2_mod_n: T,
}

Construction remains fallible once:

impl<T> Odd<T> {
    pub fn new(value: T) -> Option<Self> { ... }
}

But all hot-path Montgomery operations receive an Odd<T> or MontModulus<T> and become infallible:

pub fn montgomery_mod_mul_checked<T>(a: T, b: T, n: T) -> Option<T>;
pub fn montgomery_mod_mul_odd<T>(a: Reduced<T>, b: Reduced<T>, n: &MontModulus<T>) -> T;

Expected benefits

    Remove repeated modulus == 0 || modulus.is_even() checks from every operation after context construction.

    Remove Option from optimized Montgomery entrypoints.

    Remove repeated Parity bounds from inner APIs that no longer inspect modulus parity.

    Make APIs clearer: a Montgomery context is valid by construction.

Best target files

    modmath/src/field.rs

    modmath/src/montgomery/basic_mont.rs

    modmath/src/montgomery/constrained_mont.rs

    modmath/src/montgomery/strict_mont.rs

The biggest immediate win is Field::new: once a Field exists, self.modulus is known odd and nonzero, so all downstream operations should be able to use infallible internal entrypoints.
2. Reduced<'m, T> for pre-reduced residues
Current cost / issue

The crate already has _pr functions for pre-reduced operands:

    basic_mod_add_pr

    constrained_mod_add_pr

    strict_mod_add_pr

    basic_mod_sub_pr

    basic_mod_mul_pr

    basic_mod_exp_pr

    Montgomery _pr variants

These functions remove % and Rem bounds, but only by convention/documentation. For example, basic_mod_add reduces both operands, then delegates to basic_mod_add_pr.
The _pr variant assumes a < m and b < m.

Similarly, basic_mod_mul reduces both operands before calling basic_mod_mul_pr.
The _pr variant documents a < m and b < m and avoids the Rem bound.
Typestate opportunity

Add:

pub struct Reduced<'m, T> {
    value: T,
    modulus: PhantomData<&'m T>,
}

or, if the modulus itself needs to be carried:

pub struct Reduced<'m, T> {
    value: T,
    modulus: &'m NonZero<T>,
}

Then provide:

impl<'m, T> Reduced<'m, T> {
    pub fn new(value: T, modulus: &'m NonZero<T>) -> Self {
        Self { value: value % modulus.get(), ... }
    }

    /// # Safety
    /// Caller promises value < modulus.
    pub unsafe fn new_unchecked(value: T, modulus: &'m NonZero<T>) -> Self;
}

Operations then become:

pub fn add_reduced<'m, T>(
    a: Reduced<'m, T>,
    b: Reduced<'m, T>,
) -> Reduced<'m, T>;

Expected benefits

    Collapse basic_mod_add / basic_mod_add_pr style pairs into safer layered APIs.

    Keep the fast path but make invalid usage unrepresentable.

    Eliminate repeated % reductions across chained arithmetic.

    Avoid accidental mixing of residues from different moduli if the lifetime/brand is chosen carefully.

    Potentially remove Rem trait bounds from many downstream APIs.

Best target files

    modmath/src/add.rs

    modmath/src/sub.rs

    modmath/src/mul.rs

    modmath/src/exp.rs

    modmath/src/montgomery/basic_mont.rs

This is likely the broadest, safest typestate improvement because the crate already exposes _pr variants whose entire contract is “pre-reduced.”
3. MontForm<'m, T, P> as a lower-level public wrapper
Current state

Residue<'f, T, P> already represents a value in Montgomery form, branded by a Field<T, P>.
It is used by Field::add, Field::sub, Field::mul, and exponentiation.

However, lower-level Montgomery functions still expose raw T parameters:

pub fn wide_montgomery_mul<T>(a_mont: T, b_mont: T, modulus: T, n_prime: T) -> T

This assumes a_mont and b_mont are Montgomery residues, but the type does not encode that.

CIOS has the same shape: docs say both inputs must be Montgomery-form and in range, but the signature accepts raw &T.
Typestate opportunity

Add a reusable lower-level wrapper distinct from Field::Residue:

pub struct MontForm<'m, T, P = Nct> {
    value: Reduced<'m, T>,
    _p: PhantomData<P>,
}

or make Residue the canonical MontForm and expose it through lower-level APIs:

pub type MontForm<'m, T, P> = Residue<'m, T, P>;

Then specialized operations can require MontForm instead of raw T:

pub fn cios_mul<'m, T, P>(
    a: &MontForm<'m, T, P>,
    b: &MontForm<'m, T, P>,
    ctx: &MontCtx<'m, T, P>,
) -> MontForm<'m, T, P>;

Expected benefits

    Remove debug-only assertions from CIOS as the primary correctness mechanism. Current CIOS only debug_assert!s a < modulus and b < modulus.

    Avoid accidental raw-vs-Montgomery confusion.

    Allow infallible CIOS wrappers for valid MontForm values.

    Make “already converted to Montgomery form” composable outside Field.

Important note

Field already does much of this. The real opportunity is to push that typestate down into the lower-level Montgomery APIs so users do not have to choose between:

    safe high-level Field, and

    fast but convention-heavy raw Montgomery helpers.

4. Wide<T> / DoubleWidth<T> for REDC inputs
Current state

Wide REDC currently takes (t_lo, t_hi) as two raw T values:

pub fn wide_redc<T>(t_lo: T, t_hi: T, modulus: T, n_prime: T) -> T

The docs state that (t_lo, t_hi) represent a double-width value and that the input must be less than N * R.

wide_montgomery_mul computes this pair from a_mont.wide_mul(&b_mont) and immediately calls wide_redc.
Typestate opportunity

Add:

pub struct Wide<T> {
    lo: T,
    hi: T,
}

Then:

pub trait WideMul {
    fn wide_mul(&self, rhs: &Self) -> Wide<Self>
    where
        Self: Sized;
}

pub fn wide_redc<T>(t: Wide<T>, modulus: Odd<T>, n_prime: T) -> Reduced<T>;

Potential stronger variants:

pub struct ProductWide<T>(Wide<T>); // known product of two reduced words
pub struct RedcInput<T>(Wide<T>);   // known < N * R

Expected benefits

    Prevent swapped (hi, lo) argument bugs.

    Make the REDC precondition explicit.

    Allow optimized REDC paths that assume “product of two reduced Montgomery residues” rather than arbitrary double-width input.

    Make wide multiplication and reduction compose more naturally.

Best target files

    modmath/src/wide_mul.rs

    modmath/src/montgomery/basic_mont.rs

This is mostly a correctness/clarity improvement, but it can also help specialization: a ProductWide<Reduced<T>> can justify skipping certain defensive checks or debug assertions.
5. Odd<T> for Newton/Hensel n_prime computation
Current state

compute_n_prime_newton assumes the modulus is odd because Newton inversion modulo 2^W only works for odd moduli. It does not check this itself; callers generally check before using it.

Field::new validates the modulus, then computes n_prime using Newton.

basic_montgomery_mod_mul_pr, basic_montgomery_mod_exp_pr, and CT variants also validate then recompute n_prime.
Typestate opportunity

Change the internal API from:

pub fn compute_n_prime_newton<T>(modulus: T, w: usize) -> T

to an additional overload:

pub fn compute_n_prime_newton_odd<T>(modulus: Odd<T>, w: usize) -> T;

or, better, have MontModulus<T> construction compute it once:

pub struct MontParams<T> {
    modulus: Odd<T>,
    n_prime: T,
    r_mod_n: T,
    r2_mod_n: T,
}

Expected benefits

    Remove repeated oddness checks before each Newton computation.

    Document the mathematical precondition in the type system.

    Encourage precomputation reuse instead of mod_mul / mod_exp recomputing n_prime, r_mod_n, and r2_mod_n every call.

This is especially relevant because standalone basic_montgomery_mod_mul_pr recomputes all Montgomery parameters even when called repeatedly with the same modulus.
Field already avoids that by storing modulus, n_prime, r_mod_n, and r2_mod_n.
6. MontParams<T> / Field as the primary optimized entrypoint
Current cost / issue

The free Montgomery functions are convenient, but several optimized-looking functions still do full per-call precomputation:

let w = type_bit_width::<T>();
let n_prime = compute_n_prime_newton(modulus, w);
let r_mod_n = compute_r_mod_n(modulus, w);
let r2_mod_n = compute_r2_mod_n(r_mod_n, modulus, w);

This occurs in basic_montgomery_mod_mul_pr, basic_montgomery_mod_exp_pr, and CT variants.

Field<T, P> already stores these four values.
Typestate opportunity

Make a precomputed context the explicit optimized API:

pub struct MontCtx<T, P = Nct> {
    modulus: Odd<T>,
    n_prime: T,
    r_mod_n: T,
    r2_mod_n: T,
    _p: PhantomData<P>,
}

Then standalone functions become thin convenience wrappers:

pub fn mod_mul(a: T, b: T, modulus: T) -> Option<T> {
    let ctx = MontCtx::new(modulus)?;
    ctx.mul_raw(a, b)
}

And optimized callers use:

ctx.mul_reduced(a, b)
ctx.exp_public_exp(base, exp)
ctx.exp_secret_exp(base, exp)

Expected benefits

    Avoid repeated precomputation.

    Remove repeated validity checks after context creation.

    Centralize invariants.

    Make fast entrypoints discoverable.

Field<T, P> already nearly fits this role, so the API decision is whether to:

    extend Field as the canonical context, or

    introduce a smaller lower-level MontCtx.

Given the existing implementation, I would prefer extending Field unless there is a desire to keep field.rs behind the wide-mul feature.
7. NonZero<T> / Unit<T> / Invertible<T> for modular inverse
Current state

*_mod_inv returns Option<T> because not all values are invertible. The extended Euclidean loop returns None when gcd(a, modulus) > 1.

This is mathematically necessary for arbitrary rings. But in common fields — prime modulus, odd Montgomery modulus with known coprime inputs, etc. — callers often know invertibility.
Typestate opportunity

Introduce one or more wrappers:

pub struct Unit<'m, T>(Reduced<'m, T>);      // gcd(value, modulus) == 1
pub struct NonZeroResidue<'m, T>(Reduced<'m, T>); // for prime fields, enough
pub struct PrimeModulus<T>(Odd<NonZero<T>>);

Then:

pub fn inv_unit<'m, T>(x: Unit<'m, T>) -> Reduced<'m, T>;
pub fn inv_prime_nonzero<'m, T>(x: NonZeroResidue<'m, T>, p: PrimeModulus<T>) -> Reduced<'m, T>;

Expected benefits

    Expose infallible inverse APIs where invertibility is already known.

    In prime fields, NonZeroResidue can eliminate the None case.

    For Montgomery fields, this could enable optimized Fermat inversion entrypoints using p - 2 when modulus is statically prime.

Caveat

A generic NonZero<T> is not enough for invertibility under composite modulus. For example, 6 mod 8 is nonzero but not invertible. So the right typestate is either:

    Unit<'m, T>: proven coprime to modulus, or

    NonZeroResidue<'m, T> plus PrimeModulus<T>.

8. Exponent<T> typestates: public, secret, fixed-width, nonzero
Current state

There are several exponentiation variants with different leakage/performance contracts:

    Schoolbook basic_mod_exp_pr branches on exponent bits and loop count.

    Montgomery basic_montgomery_mod_exp_pr also branches on exp.is_odd() and stops when exp == 0.

    CT Montgomery exponentiation iterates over every bit and uses conditional selection.

    Field<T, Ct>::exp is fixed-iteration and constant-time over the exponent.

    Field<T, Ct>::exp_public_exp deliberately leaks public exponent bit length and bit pattern for speed.

Typestate opportunity

Add exponent wrappers:

pub struct PublicExp<T>(T);
pub struct SecretExp<T>(T);
pub struct NonZeroExp<T>(T);
pub struct FixedWidthExp<T>(T);

Then expose:

pub fn exp_public(&self, base: &MontForm<T>, exp: PublicExp<T>) -> MontForm<T>;
pub fn exp_secret(&self, base: &MontForm<T>, exp: SecretExp<T>) -> MontForm<T>;
pub fn exp_nonzero_public(&self, base: &MontForm<T>, exp: NonZeroExp<PublicExp<T>>) -> MontForm<T>;

Expected benefits

    Make the security/performance choice explicit at the type level.

    Avoid accidental use of variable-time exponentiation with secret exponents.

    A NonZeroExp<T> can remove the hi == 0 branch in public exponent code. Current exp_public_exp searches for the highest bit and has a special exp == 0 return path.

    Common public exponent 65537 could receive a compact addition-chain or fixed-shape specialization.

Suggested concrete optimized entrypoints

pub fn exp_u65537(&self, base: &Residue<'_, T, P>) -> Residue<'_, T, P>;
pub fn square(&self, x: &Residue<'_, T, P>) -> Residue<'_, T, P>;
pub fn square_n<const N: usize>(&self, x: &Residue<'_, T, P>) -> Residue<'_, T, P>;

This would be useful for RSA verify/encrypt and prime-field inversion chains.
9. Even<T> / Odd<T> for parity checks and bit loops
Current state

The crate has a generic Parity trait. Without num-integer, the default implementation clones and computes self & one == one, which can be expensive for big integers.

Many loops call is_odd() repeatedly:

    double-and-add multiplication.

    constrained/strict multiplication.

    exponentiation.

    wide multiplication fallback.

Typestate opportunity

Odd<T> and Even<T> are mostly useful for inputs known once, not for loop variables whose parity changes each iteration. Good targets:

    Montgomery modulus must be Odd<T>.

    Prime field modulus is odd except p = 2.

    Algorithms specialized for odd exponent or even exponent.

Less useful:

    The shifting loop in multiplication/exponentiation, because b or exp changes every iteration.

Expected benefits

    Remove repeated modulus parity checks.

    Allow specialized algorithms for odd modulus.

    Avoid requiring Parity on internal functions that only need oddness of a pre-validated modulus.

Not recommended

Do not wrap each loop iteration’s b/exp into Odd/Even; that would add overhead and not compact the math. Use bit iteration APIs instead.

A better direction for loops is a trait like:

trait Bits {
    fn bit(&self, i: usize) -> bool;
    fn bit_len(&self) -> usize;
}

This could replace repeated is_odd + shift cloning for big integers.
10. Reduced<T> with lazy bounds: Reduced2<T>, Reduced4<T>, etc.
Current state

Field::add always canonicalizes with a conditional subtract. The code explicitly warns that skipping the cond-sub is only valid for moduli with slack and belongs in specialization layers.
Typestate opportunity

Introduce lazy-reduction wrappers:

pub struct Reduced<'m, T>;   // 0 <= x < n
pub struct Loose<'m, T>;     // 0 <= x < 2n
pub struct WideLoose<'m, T>; // 0 <= x < 4n

Then:

pub fn add_lazy(a: Reduced<T>, b: Reduced<T>) -> Loose<T>;
pub fn normalize(x: Loose<T>) -> Reduced<T>;
pub fn mul_requires_reduced(a: Reduced<T>, b: Reduced<T>) -> Reduced<T>;

Expected benefits

    Fewer conditional subtracts in addition-heavy code.

    Better compactness for curve-field formulas, where multiple additions/subtractions occur before a multiplication.

    Keeps generic Field::add correct while allowing downstream specialized layers to encode their slack assumptions.

Caveat

This is only safe when the backing type has enough headroom for the chosen loose bound. That requires either:

pub struct HasSlack<const LIMBS: usize, const BITS: usize>;

or backend-specific typestates. The existing comment in Field::add is correct: the generic full-width RSA use case cannot skip the conditional subtract.
11. Replace internal expect("CIOS mul cannot fail...") with infallible typestate-backed CIOS
Current state

Field::mul calls CiosMontMul::cios_mont_mul(...).expect("CIOS mul cannot fail with valid Montgomery parameters").

This is semantically true if:

    Field parameters are valid,

    n_prime has at least one word,

    CIOS loop indices are valid,

    residues are in range and Montgomery form.

But the trait still returns Option<Self>.
Typestate opportunity

Provide an infallible internal trait for fixed-size integer backends whose limb indexing is statically valid:

pub trait CiosMontMulInfallible {
    fn cios_mont_mul_infallible(
        a: &Self,
        b: &Self,
        modulus: &Self,
        n_prime_0: Self::Word,
    ) -> Self;
}

Or encode the first word:

pub struct NPrime0<W>(W);

and precompute it in MontCtx.
Expected benefits

    Remove Option and expect from the hot path.

    Avoid repeated n_prime.get_word(0)? in every multiplication. The trait currently extracts word 0 on each call.

    Store n_prime_0 directly in Field/MontCtx when CIOS is available.

This is a concrete speed opportunity: CIOS multiplication is the exponentiation inner loop, and Field::exp calls CIOS for every multiply/square.
12. PrecomputedR<T> / RModN<T> / R2ModN<T>
Current state

Field stores r_mod_n and r2_mod_n.
from_precomputed accepts raw T values and documents that callers are responsible for correctness.
Typestate opportunity

Add wrappers:

pub struct NPrime<T>(T);
pub struct RModN<T>(T);
pub struct R2ModN<T>(T);

pub struct MontParams<T> {
    modulus: Odd<T>,
    n_prime: NPrime<T>,
    r_mod_n: RModN<T>,
    r2_mod_n: R2ModN<T>,
}

Then:

pub const fn from_precomputed(params: MontParams<T>) -> Self;

Expected benefits

    Reduce accidental parameter-order bugs in from_precomputed.

    Make build-script generated constants safer.

    Allow storing n_prime_0 as part of params for CIOS.

This is more of a safety/maintainability win than a runtime speed win, except for the possible n_prime_0 caching.
13. PrimeModulus<T> for field-specific optimizations
Current state

Field<T, P> only requires an odd nonzero modulus, not primality.

That is correct for Montgomery arithmetic over odd composite moduli, including RSA. But prime fields have extra opportunities:

    Every nonzero residue is invertible.

    Fermat inversion x^(p-2) is valid.

    Square roots may use exponent formulas depending on p mod 4 or p mod 8.

    Addition chains can be modulus-specific.

Typestate opportunity

Add:

pub struct PrimeField<T, P> {
    field: Field<T, P>,
}

or:

pub struct PrimeModulus<T>(Odd<T>);

Construction can be:

    unchecked for compile-time known primes,

    probabilistic/deterministic checked for primitive sizes,

    delegated to downstream crates for curve constants.

Expected benefits

    Infallible inverse for NonZeroResidue.

    Specialized inverse/sqrt exponentiation chains.

    Smaller APIs for prime-field consumers without compromising RSA/composite support.

Caveat

Do not make Field itself assume primality. RSA needs odd composite moduli.
14. Canonical<T> vs Raw<T> for signed / unsigned reduction
Current state

reduce_mod assumes unsigned types because % on signed types may produce negative values, violating [0, modulus).
Typestate opportunity

Add:

pub struct Canonical<'m, T>(Reduced<'m, T>);
pub struct Raw<T>(T);

Or split constructors:

Reduced::from_unsigned_rem(raw, modulus)
Reduced::from_signed_rem(raw, modulus)

Expected benefits

    Make the signedness/canonical-representative issue explicit.

    Avoid accidentally feeding negative signed remainders into Montgomery conversion.

    Allow Reduced<T> to guarantee [0, modulus) independent of backend % behavior.

This is more correctness than speed, but it helps ensure optimized entrypoints are sound.
Prioritized implementation plan
Priority 1: Odd<NonZero<T>> / MontParams<T>

This gives the biggest immediate win for Montgomery APIs:

    one-time validation,

    no repeated Option,

    no repeated parity checks,

    better precomputation reuse.

Start with Field::new and basic_montgomery_mod_* paths.

Recommended shape:

pub struct OddModulus<T> {
    value: T,
}

pub struct MontParams<T> {
    modulus: OddModulus<T>,
    n_prime: T,
    n_prime_0: Option<<T as MulAccOps>::Word>, // when available
    r_mod_n: T,
    r2_mod_n: T,
}

Priority 2: Reduced<'m, T>

This formalizes the existing _pr APIs and removes repeated % in chained operations.

Start by wrapping the _pr family:

    add_pr

    sub_pr

    mul_pr

    exp_pr

    Montgomery mod_mul_pr

    Montgomery mod_exp_pr

The crate already has the exact fast implementations; the task is mostly API layering.
Priority 3: CIOS infallible/cached-n_prime_0 path

CIOS currently returns Option<T> and extracts n_prime.get_word(0)? per call.
In Field::mul, this is immediately unwrapped with expect.

For exponentiation, this is hot. Cache n_prime_0 in the context if possible.
Priority 4: MontForm<'m, T, P> lower-level API

Residue already covers this at the high level.
The opportunity is to make lower-level wide-REDC/CIOS APIs safer and infallible without forcing users into all of Field.
Priority 5: exponent typestates

Add PublicExp<T> and SecretExp<T> to prevent variable-time exponentiation misuse and enable optimized public-exponent entrypoints. This maps cleanly onto existing Field<T, Ct>::exp vs exp_public_exp.
Concrete API sketch

pub struct NonZero<T>(T);
pub struct Odd<T>(NonZero<T>);

pub struct Modulus<T> {
    n: Odd<T>,
}

pub struct MontParams<T> {
    modulus: Modulus<T>,
    n_prime: T,
    r_mod_n: T,
    r2_mod_n: T,
}

pub struct Reduced<'m, T> {
    value: T,
    _modulus: PhantomData<&'m Modulus<T>>,
}

pub struct MontForm<'m, T, P = Nct> {
    value: Reduced<'m, T>,
    _p: PhantomData<P>,
}

pub struct Wide<T> {
    lo: T,
    hi: T,
}

pub struct PublicExp<T>(T);
pub struct SecretExp<T>(T);
pub struct NonZeroExp<T>(T);

Fast-path examples:

impl<T> MontParams<T> {
    pub fn new(modulus: T) -> Option<Self>;
}

impl<'m, T> Reduced<'m, T> {
    pub fn reduce(value: T, params: &'m MontParams<T>) -> Self;
    pub unsafe fn new_unchecked(value: T, params: &'m MontParams<T>) -> Self;
}

impl<'m, T, P> MontForm<'m, T, P> {
    pub fn from_reduced(x: Reduced<'m, T>, params: &'m MontParams<T>) -> Self;
    pub fn into_reduced(self, params: &'m MontParams<T>) -> Reduced<'m, T>;
}

impl<T, P> Field<T, P> {
    pub fn from_params(params: MontParams<T>) -> Self;
    pub fn params(&self) -> &MontParams<T>;
}

Summary of best opportunities
Typestate	Best use	Primary benefit
NonZero<T>	generic modular arithmetic modulus	remove division-by-zero/error branches
Odd<T>	Montgomery modulus	eliminate repeated even-modulus checks
Reduced<'m, T>	_pr APIs	remove %, Rem bounds, and invalid pre-reduced calls
MontForm<'m, T, P>	lower-level Montgomery APIs	avoid raw/Montgomery mixups and make CIOS safer
Wide<T>	REDC inputs	prevent (lo, hi) misuse and enable REDC-specific specialization
MontParams<T>	repeated Montgomery operations	cache n_prime, r_mod_n, r2_mod_n, maybe n_prime_0
PublicExp<T> / SecretExp<T>	exponentiation	make timing contract explicit and select optimized algorithm
Unit<'m, T>	modular inverse	infallible inverse where gcd is known 1
PrimeModulus<T>	prime fields	infallible nonzero inverse and addition-chain optimizations
Loose<'m, T>	curve-style formulas	fewer conditional reductions in addition-heavy code

The most compelling path is: OddModulus<T> → MontParams<T> → Reduced<'m, T> → MontForm<'m, T, P>. That sequence matches the crate’s existing architecture, formalizes current _pr contracts, and should remove the most repeated checks and conversions without disrupting the generic schoolbook APIs.
What about Prime<T> ?
Yes — Prime<T> is a very good typestate, but it should be separate from Odd<T>

A Prime<T> wrapper would be valuable, but it encodes a stronger and more expensive invariant than the ones discussed earlier:

Odd<T>       // n is odd
NonZero<T>   // n != 0
Prime<T>     // n is prime

For this crate, Prime<T> would be most useful not for basic Montgomery multiplication itself, but for field-level operations that become faster, infallible, or more specialized when the modulus is prime.

The current Field<T, P> only requires an odd nonzero modulus; Field::new rejects zero and even moduli, computes Montgomery parameters, and then accepts the modulus as valid for Montgomery arithmetic.
That is correct because Montgomery arithmetic works over odd composite moduli too, including RSA-style moduli. So Prime<T> should be an additional layer, not a replacement for Field<T, P>.
Where Prime<T> helps most
1. Infallible inverse for nonzero residues

The biggest semantic win is modular inverse.

Currently, inverse returns Option<T> because not every element is invertible modulo an arbitrary modulus. The strict, constrained, and basic inverse functions all return Option<T>.

That Option is necessary for composite moduli:

6 mod 8 is nonzero, but not invertible.

But in a prime field:

Every nonzero residue mod p is invertible.

So with typestates:

Prime<T>
Reduced<'p, T>
NonZeroResidue<'p, T>

you can expose:

impl<T, P> PrimeField<T, P> {
    pub fn inv_nonzero(&self, x: NonZeroResidue<'_, T, P>) -> Residue<'_, T, P>;
}

No Option.

This is strictly stronger than NonZero<T> alone. NonZero<T> tells you x != 0; Prime<T> tells you that x != 0 implies gcd(x, p) == 1.
2. Fermat inverse: x^(p - 2) mod p

A Prime<T> modulus enables Fermat’s little theorem:

x⁻¹ ≡ x^(p-2) mod p, for x != 0 and p prime

The crate already has optimized Montgomery exponentiation machinery and even calls out Curve25519 Fermat inversion as a valid use case for public-exponent CT-base exponentiation.

That suggests a natural API:

impl<T, P> PrimeField<T, P> {
    pub fn inv_fermat(&self, x: NonZeroResidue<'_, T, P>) -> Residue<'_, T, P>;
}

For generic prime fields, this would compute p - 2 and use existing exponentiation.

For well-known static primes, downstream crates could provide hand-optimized addition chains:

impl Curve25519Field {
    pub fn inv(&self, x: NonZeroResidue<'_, U25519>) -> Residue<'_, U25519> {
        // fixed addition chain for p - 2
    }
}

This is especially attractive because the existing Field<T, Ct>::exp_public_exp already has the right contract: constant-time over the base but variable-time over a public exponent.
3. Infallible division

Once you have prime modulus + nonzero denominator, division can be infallible:

pub fn div(
    &self,
    numerator: &Residue<'_, T, P>,
    denominator: NonZeroResidue<'_, T, P>,
) -> Residue<'_, T, P>;

Without Prime<T>, this is unsound for composite moduli because a nonzero denominator may still not be a unit.

So the proper generic split is:

// Works for any modulus, but fallible.
fn checked_inv(x: Reduced<'m, T>, modulus: &Modulus<T>) -> Option<Reduced<'m, T>>;

// Works for any modulus, infallible only after proving gcd(x, n) == 1.
fn inv_unit(x: Unit<'m, T>) -> Reduced<'m, T>;

// Works for prime modulus, because nonzero implies unit.
fn inv_prime_nonzero(x: NonZeroResidue<'m, T>, p: Prime<T>) -> Reduced<'m, T>;

4. Square root and quadratic-residue APIs

Prime fields also enable specialized square-root operations:

pub fn sqrt(&self, x: Residue<'_, T, P>) -> Option<Residue<'_, T, P>>;
pub fn sqrt_checked_qr(&self, x: QuadraticResidue<'_, T, P>) -> Residue<'_, T, P>;

Possible typestates:

Prime<T>
PrimeMod4Eq3<T>
PrimeMod8Eq5<T>
QuadraticResidue<'p, T>

For primes with special congruence classes:

    If p % 4 == 3, then sqrt(x) = x^((p + 1) / 4).

    If p % 8 == 5, there are also specialized formulas.

    For general odd prime fields, Tonelli–Shanks applies.

This would be outside the current crate’s core, but the Field abstraction is a natural base.
5. Stronger “field” semantics distinct from Montgomery “ring” semantics

Right now the type is named Field<T, P>, but the constructor accepts any odd nonzero modulus.
Mathematically, that is a ring for composite moduli, not always a field.

That naming is understandable from a crypto ergonomics perspective, but a Prime<T> typestate would let the crate distinguish:

MontRing<T, P>      // odd modulus, Montgomery arithmetic valid
PrimeField<T, P>   // prime modulus, field laws valid

Possible compatibility-friendly approach:

pub type Field<T, P = Nct> = MontRing<T, P>;

pub struct PrimeField<T, P = Nct> {
    inner: Field<T, P>,
    _prime: Prime<T>,
}

Or, less disruptive:

pub struct Field<T, P = Nct> {
    // current type
}

pub struct PrimeField<T, P = Nct> {
    field: Field<T, P>,
}

Then PrimeField adds APIs that require primality:

impl<T, P> PrimeField<T, P> {
    pub fn inv_nonzero(...);
    pub fn div_nonzero(...);
    pub fn fermat_inv(...);
    pub fn legendre(...);
    pub fn sqrt(...);
}

Construction: checked, probable, and unchecked

The hard part is not the wrapper — it is proving primality.

I would recommend three constructors:

impl<T> Prime<T> {
    /// Deterministic for small primitive integer types.
    pub fn new_checked(n: T) -> Option<Self>;

    /// Probabilistic / Miller-Rabin style, for large dynamic integers.
    pub fn new_probable_prime(n: T, rounds: usize) -> Option<Self>;

    /// For static known primes, curve constants, test vectors, generated params.
    ///
    /// # Safety
    /// Caller promises n is prime.
    pub const unsafe fn new_unchecked(n: T) -> Self;
}

For cryptographic code, unsafe fn new_unchecked is often the most practical API because primes are usually constants from standards. But it should be explicit and well documented.

A less scary alternative is:

pub const fn new_trusted(n: T) -> Self;

with docs saying “this does not test primality,” but unsafe better communicates that invalid construction can make later infallible APIs mathematically wrong.
Relationship to Odd<T>

Prime<T> should imply Odd<T> for all primes except 2.

This raises an API design question: do you want to support the field modulo 2?

For this crate’s Montgomery paths, modulus 2 is not usable because Montgomery requires odd modulus and current Field::new rejects even modulus.

So I would probably define:

pub struct OddPrime<T>(T);

rather than plain Prime<T> for Montgomery contexts.

Possible hierarchy:

pub struct Prime<T>(T);      // includes 2
pub struct OddPrime<T>(T);   // excludes 2, valid for Montgomery

Then:

impl<T> TryFrom<Prime<T>> for OddPrime<T>;
impl<T> From<OddPrime<T>> for Odd<T>;
impl<T> From<OddPrime<T>> for NonZero<T>;

For this crate, OddPrime<T> is more directly useful than Prime<T>.
Best concrete design for this crate

I would introduce:

pub struct OddModulus<T> {
    value: T,
}

pub struct OddPrime<T> {
    value: T,
}

pub struct MontParams<T> {
    modulus: OddModulus<T>,
    n_prime: T,
    r_mod_n: T,
    r2_mod_n: T,
}

pub struct PrimeField<T, P = Nct> {
    inner: Field<T, P>,
    prime: OddPrime<T>,
}

Then:

impl<T, P> Field<T, P> {
    pub fn new_odd(modulus: OddModulus<T>) -> Self;
}

impl<T, P> PrimeField<T, P> {
    pub fn new_prime(prime: OddPrime<T>) -> Self;

    pub fn inv_nonzero(
        &self,
        x: NonZeroResidue<'_, T, P>,
    ) -> Residue<'_, T, P>;

    pub fn div_nonzero(
        &self,
        a: &Residue<'_, T, P>,
        b: NonZeroResidue<'_, T, P>,
    ) -> Residue<'_, T, P>;
}

This keeps the generic odd-modulus Montgomery use case intact while adding true field-only APIs.
What Prime<T> should not be used for
Do not require Prime<T> for Montgomery multiplication

Montgomery multiplication only requires an odd modulus, not a prime modulus. The crate’s current Montgomery machinery correctly accepts odd composite moduli; tests and docs even mention RSA-style reuse and full-width moduli. Field stores Montgomery parameters and works over any odd nonzero modulus.

Requiring Prime<T> for Montgomery multiplication would unnecessarily exclude RSA and other ring uses.
Do not treat Prime<T> as proof that a residue is nonzero

Primality of the modulus does not mean a residue is nonzero. You still need:

NonZeroResidue<'p, T>

for infallible inverse/division.
Do not make dynamic primality testing mandatory

For fixed cryptographic primes, runtime primality testing is wasted code size and startup time. Support unchecked/static construction.
Speed/compactness wins from Prime<T>
Opportunity	Why Prime<T> helps	Speed/compactness impact
Inverse	nonzero residue is always invertible	remove Option, avoid gcd failure branch
Division	denominator nonzero implies denominator is a unit	infallible div_nonzero
Fermat inverse	x^(p-2) valid	reuse Montgomery exponentiation; specialize addition chains
Square root	prime-field algorithms available	specialized formulas for p mod 4/8
Legendre symbol	only defined cleanly for odd primes	field-specific APIs
Addition chains	known prime constants imply known exponents	fewer branches/loops for fixed primes
API correctness	distinguish ring from field	prevents using field-only operations under composite modulus
My recommendation

Add Prime<T> only if it is paired with more specific wrappers:

OddPrime<T>
PrimeField<T, P>
NonZeroResidue<'p, T, P>
Unit<'m, T, P>

The most useful stack would be:

OddPrime<T>
    -> PrimeField<T, P>
        -> Residue<'p, T, P>
        -> NonZeroResidue<'p, T, P>
            -> inv_nonzero() // infallible

Prime<T> alone is not enough to speed up ordinary add/mul; those only need Reduced<T> and Odd<T>. But Prime<T> becomes very powerful for inversion, division, square roots, known-exponent chains, and true field-only APIs.
