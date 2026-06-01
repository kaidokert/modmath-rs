//! R > N Montgomery arithmetic for `Clone`-bound types with
//! mixed by-value / by-reference parameters.
//!
//! Textbook Montgomery: pick `R` as the smallest power of 2 greater
//! than the modulus, precompute `R⁻¹ mod N` and `N' = -N⁻¹ mod R`,
//! transform values into Montgomery form, multiply in the Mont domain,
//! transform back. Use [`compute_params`] to derive the parameters,
//! [`to_mont`] / [`from_mont`] for representation transforms, and
//! [`mul`] for in-Mont multiplication. The full-pipeline operations
//! [`mod_mul`] and [`mod_exp`] wrap precompute + transform + arith.
//!
//! ## Distinction between `mul` and `mod_mul`
//!
//! - [`mul`] — takes Mont-form values and Montgomery params; returns
//!   a Mont-form value. The in-domain multiplication step.
//! - [`mod_mul`] — takes raw values and the modulus only; internally
//!   computes Mont params, transforms inputs, multiplies, transforms
//!   back. The complete `a * b mod m` pipeline.
//!
//! ## N' computation method
//!
//! [`compute_params`] uses a default algorithm. To pick explicitly
//! between trial-search / extended-Euclidean / Hensel's lifting, use
//! [`compute_params_with_method`], [`mod_mul_with_method`], or
//! [`mod_exp_with_method`] with a [`crate::NPrimeMethod`] choice.

#[doc(inline)]
pub use crate::montgomery::constrained_mont::constrained_compute_montgomery_params as compute_params;
#[doc(inline)]
pub use crate::montgomery::constrained_mont::constrained_compute_montgomery_params_with_method as compute_params_with_method;
#[doc(inline)]
pub use crate::montgomery::constrained_mont::constrained_from_montgomery as from_mont;
#[doc(inline)]
pub use crate::montgomery::constrained_mont::constrained_montgomery_mod_exp as mod_exp;
#[doc(inline)]
pub use crate::montgomery::constrained_mont::constrained_montgomery_mod_exp_with_method as mod_exp_with_method;
#[doc(inline)]
pub use crate::montgomery::constrained_mont::constrained_montgomery_mod_mul as mod_mul;
#[doc(inline)]
pub use crate::montgomery::constrained_mont::constrained_montgomery_mod_mul_with_method as mod_mul_with_method;
#[doc(inline)]
pub use crate::montgomery::constrained_mont::constrained_montgomery_mul as mul;
#[doc(inline)]
pub use crate::montgomery::constrained_mont::constrained_to_montgomery as to_mont;
