//! Montgomery arithmetic for `Copy`-bound types.
//!
//! Two algorithm families live under this module:
//!
//! - **R > N path**: textbook Montgomery construction. Pick `R` as the
//!   smallest power of 2 greater than the modulus, precompute
//!   `R⁻¹ mod N` and `N' = -N⁻¹ mod R`, transform values into
//!   Montgomery form, multiply in the Mont domain, transform back.
//!   Building blocks: [`compute_params`], [`to_mont`], [`from_mont`],
//!   [`mul`].
//!
//! - **Wide-REDC path** (`R = 2^W` exactly the word width): the
//!   complete-pipeline operations [`mod_mul`] and [`mod_exp`] use this
//!   internally. Overflow-free for full-width moduli; no caller-visible
//!   Montgomery params required.
//!
//! For the high-level type-safe surface, see
//! [`modmath::Field`](crate::Field) — Montgomery context with
//! lifetime-branded residues.
//!
//! ## N' computation method
//!
//! [`compute_params`] uses a default algorithm internally. To pick
//! between trial-search / extended-Euclidean / Hensel's lifting, use
//! [`compute_params_with_method`] with an explicit
//! [`crate::NPrimeMethod`] choice. (The complete pipelines [`mod_mul`]
//! and [`mod_exp`] use wide-REDC which mandates Newton; method
//! selection is meaningful only on the R>N building-block path.)

// R > N path building blocks
#[doc(inline)]
pub use crate::montgomery::basic_mont::basic_compute_montgomery_params as compute_params;
#[doc(inline)]
pub use crate::montgomery::basic_mont::basic_compute_montgomery_params_with_method as compute_params_with_method;
#[doc(inline)]
pub use crate::montgomery::basic_mont::basic_from_montgomery as from_mont;
#[doc(inline)]
pub use crate::montgomery::basic_mont::basic_montgomery_mul as mul;
#[doc(inline)]
pub use crate::montgomery::basic_mont::basic_to_montgomery as to_mont;

// Wide-REDC complete-pipeline operations
#[doc(inline)]
pub use crate::montgomery::basic_mont::basic_montgomery_mod_exp as mod_exp;
#[doc(inline)]
pub use crate::montgomery::basic_mont::basic_montgomery_mod_exp_odd as mod_exp_odd;
#[doc(inline)]
pub use crate::montgomery::basic_mont::basic_montgomery_mod_mul as mod_mul;
#[doc(inline)]
pub use crate::montgomery::basic_mont::basic_montgomery_mod_mul_odd as mod_mul_odd;

/// Wide-REDC primitives at the `R = 2^W` working width — by-value
/// operands, `Copy`-bound `T`. These are the building blocks the
/// [`mod_mul`] / [`mod_exp`] pipelines call internally; consumers
/// reach for them directly when implementing their own Montgomery
/// pipelines (PQC's `Mont` newtype pattern, RSA's `ModMathParams`).
///
/// For the reference-based variant suitable for non-`Copy` bigint
/// backends, see [`modmath::strict::montgomery::wide`](crate::strict::montgomery::wide).
pub mod wide {
    #[doc(inline)]
    pub use crate::montgomery::basic_mont::wide_montgomery_mul as mul;
    #[doc(inline)]
    pub use crate::montgomery::basic_mont::wide_montgomery_mul_acc as mul_acc;
    #[doc(inline)]
    pub use crate::montgomery::basic_mont::wide_redc as redc;

    /// Constant-time finalize. Branchless conditional-subtract via
    /// `subtle::ConditionallySelectable`.
    pub mod ct {
        #[doc(inline)]
        pub use crate::montgomery::basic_mont::wide_montgomery_mul_acc_ct as mul_acc;
        #[doc(inline)]
        pub use crate::montgomery::basic_mont::wide_montgomery_mul_ct as mul;
        #[doc(inline)]
        pub use crate::montgomery::basic_mont::wide_redc_ct as redc;
    }
}
