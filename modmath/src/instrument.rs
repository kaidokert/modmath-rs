//! Division/remainder instrumentation counters.
//!
//! When the `instrument` feature is enabled, every `%`, `%=`, and `/` operation
//! in the modular arithmetic functions increments a static atomic counter.
//! This lets you find out exactly where division budget is being spent.
//!
//! Usage from an application:
//! ```ignore
//! modmath::instrument::reset_all();
//! // ... do work ...
//! modmath::instrument::dump_nonzero();
//! ```

use core::sync::atomic::{AtomicU64, Ordering};

// ---- add.rs ----
/// strict_mod_add: a.rem_assign(m)
pub static STRICT_ADD_A: AtomicU64 = AtomicU64::new(0);
/// strict_mod_add: b % m
pub static STRICT_ADD_B: AtomicU64 = AtomicU64::new(0);
/// constrained_mod_add: &a % m
pub static CONSTRAINED_ADD_A: AtomicU64 = AtomicU64::new(0);
/// constrained_mod_add: b % m
pub static CONSTRAINED_ADD_B: AtomicU64 = AtomicU64::new(0);
/// basic_mod_add: a % m
pub static BASIC_ADD_A: AtomicU64 = AtomicU64::new(0);
/// basic_mod_add: b % m
pub static BASIC_ADD_B: AtomicU64 = AtomicU64::new(0);

// ---- sub.rs ----
/// strict_mod_sub: a.rem_assign(m)
pub static STRICT_SUB_A: AtomicU64 = AtomicU64::new(0);
/// strict_mod_sub: b % m
pub static STRICT_SUB_B: AtomicU64 = AtomicU64::new(0);
/// constrained_mod_sub: &a % m
pub static CONSTRAINED_SUB_A: AtomicU64 = AtomicU64::new(0);
/// constrained_mod_sub: b % m
pub static CONSTRAINED_SUB_B: AtomicU64 = AtomicU64::new(0);
/// basic_mod_sub: a % m
pub static BASIC_SUB_A: AtomicU64 = AtomicU64::new(0);
/// basic_mod_sub: b % m
pub static BASIC_SUB_B: AtomicU64 = AtomicU64::new(0);

// ---- mul.rs ----
/// strict_mod_mul: a.rem_assign(m)
pub static STRICT_MUL_A: AtomicU64 = AtomicU64::new(0);
/// strict_mod_mul: b % m
pub static STRICT_MUL_B: AtomicU64 = AtomicU64::new(0);
/// constrained_mod_mul: a.rem_assign(m)
pub static CONSTRAINED_MUL_A: AtomicU64 = AtomicU64::new(0);
/// constrained_mod_mul: b % m
pub static CONSTRAINED_MUL_B: AtomicU64 = AtomicU64::new(0);
/// basic_mod_mul: a % m
pub static BASIC_MUL_A: AtomicU64 = AtomicU64::new(0);
/// basic_mod_mul: b % m
pub static BASIC_MUL_B: AtomicU64 = AtomicU64::new(0);

// ---- exp.rs ----
/// strict_mod_exp: T::one() % modulus
pub static STRICT_EXP_ONE: AtomicU64 = AtomicU64::new(0);
/// strict_mod_exp: base.rem_assign(modulus)
pub static STRICT_EXP_BASE: AtomicU64 = AtomicU64::new(0);
/// constrained_mod_exp: base.rem_assign(modulus)
pub static CONSTRAINED_EXP_BASE: AtomicU64 = AtomicU64::new(0);
/// constrained_mod_exp: T::one() % modulus
pub static CONSTRAINED_EXP_ONE: AtomicU64 = AtomicU64::new(0);
/// basic_mod_exp: T::one() % modulus
pub static BASIC_EXP_ONE: AtomicU64 = AtomicU64::new(0);
/// basic_mod_exp: base = base % modulus
pub static BASIC_EXP_BASE: AtomicU64 = AtomicU64::new(0);

// ---- inv.rs ----
/// strict_mod_inv: &r / &new_r (loop body)
pub static STRICT_INV_DIV: AtomicU64 = AtomicU64::new(0);
/// constrained_mod_inv: &r / &new_r (loop body)
pub static CONSTRAINED_INV_DIV: AtomicU64 = AtomicU64::new(0);
/// basic_mod_inv: r / new_r (loop body, via Signed)
pub static BASIC_INV_DIV: AtomicU64 = AtomicU64::new(0);

// ---- strictest.rs ----
/// strictest_mod_add: a % m
pub static STRICTEST_ADD_A: AtomicU64 = AtomicU64::new(0);
/// strictest_mod_add: b % m
pub static STRICTEST_ADD_B: AtomicU64 = AtomicU64::new(0);
/// strictest_mod_sub: a % m
pub static STRICTEST_SUB_A: AtomicU64 = AtomicU64::new(0);
/// strictest_mod_sub: b % m
pub static STRICTEST_SUB_B: AtomicU64 = AtomicU64::new(0);
/// strictest_mod_mul: a % m
pub static STRICTEST_MUL_A: AtomicU64 = AtomicU64::new(0);
/// strictest_mod_mul: b % m
pub static STRICTEST_MUL_B: AtomicU64 = AtomicU64::new(0);

// ---- external (ed25519 lazy_field.rs) ----
/// lazy_mod_mul: product % p (direct Rem op in ed25519 hot path)
pub static LAZY_MUL_REM: AtomicU64 = AtomicU64::new(0);
/// force_reduce: a % p (direct Rem op in ed25519)
pub static LAZY_FORCE_REDUCE: AtomicU64 = AtomicU64::new(0);

// ---- montgomery/basic_mont.rs ----
/// reduce_mod: val % modulus
pub static MONT_REDUCE_MOD: AtomicU64 = AtomicU64::new(0);
/// compute_n_prime_trial_search: (modulus * n_prime) % r (loop)
pub static MONT_NPRIME_TRIAL: AtomicU64 = AtomicU64::new(0);
/// compute_n_prime_hensels_lifting: (modulus * n_prime + 1) % target_mod (loop)
pub static MONT_NPRIME_HENSEL_LOOP: AtomicU64 = AtomicU64::new(0);
/// compute_n_prime_hensels_lifting: final check (modulus * n_prime) % r
pub static MONT_NPRIME_HENSEL_FINAL: AtomicU64 = AtomicU64::new(0);
/// compute_n_prime_extended_euclidean: uses basic_mod_inv internally
pub static MONT_NPRIME_EUCLID: AtomicU64 = AtomicU64::new(0);

/// Helper macro to collect all counters for iteration.
macro_rules! all_counters {
    ($callback:ident) => {
        $callback!(
            STRICT_ADD_A,
            STRICT_ADD_B,
            CONSTRAINED_ADD_A,
            CONSTRAINED_ADD_B,
            BASIC_ADD_A,
            BASIC_ADD_B,
            STRICT_SUB_A,
            STRICT_SUB_B,
            CONSTRAINED_SUB_A,
            CONSTRAINED_SUB_B,
            BASIC_SUB_A,
            BASIC_SUB_B,
            STRICT_MUL_A,
            STRICT_MUL_B,
            CONSTRAINED_MUL_A,
            CONSTRAINED_MUL_B,
            BASIC_MUL_A,
            BASIC_MUL_B,
            STRICT_EXP_ONE,
            STRICT_EXP_BASE,
            CONSTRAINED_EXP_BASE,
            CONSTRAINED_EXP_ONE,
            BASIC_EXP_ONE,
            BASIC_EXP_BASE,
            STRICT_INV_DIV,
            CONSTRAINED_INV_DIV,
            BASIC_INV_DIV,
            STRICTEST_ADD_A,
            STRICTEST_ADD_B,
            STRICTEST_SUB_A,
            STRICTEST_SUB_B,
            STRICTEST_MUL_A,
            STRICTEST_MUL_B,
            LAZY_MUL_REM,
            LAZY_FORCE_REDUCE,
            MONT_REDUCE_MOD,
            MONT_NPRIME_TRIAL,
            MONT_NPRIME_HENSEL_LOOP,
            MONT_NPRIME_HENSEL_FINAL,
            MONT_NPRIME_EUCLID
        );
    };
}

/// Reset all counters to zero.
pub fn reset_all() {
    macro_rules! do_reset {
        ($($name:ident),+) => {
            { $( $name.store(0, Ordering::Relaxed); )+ }
        };
    }
    all_counters!(do_reset);
}

/// Get total division count across all counters.
pub fn total_count() -> u64 {
    let mut total = 0u64;
    macro_rules! do_sum {
        ($($name:ident),+) => {
            { $( total += $name.load(Ordering::Relaxed); )+ }
        };
    }
    all_counters!(do_sum);
    total
}

/// Return a Vec of (name, count) for all non-zero counters.
/// Requires std feature.
#[cfg(feature = "std")]
pub fn nonzero_counts() -> Vec<(&'static str, u64)> {
    let mut result = Vec::new();
    macro_rules! do_collect {
        ($($name:ident),+) => {
            { $(
                let v = $name.load(Ordering::Relaxed);
                if v > 0 {
                    result.push((stringify!($name), v));
                }
            )+ }
        };
    }
    all_counters!(do_collect);
    result
}

/// Print all non-zero counters to stdout.
/// Requires std feature.
#[cfg(feature = "std")]
pub fn dump_nonzero() {
    let counts = nonzero_counts();
    if counts.is_empty() {
        println!("[instrument] No division/remainder operations recorded.");
        return;
    }
    println!("[instrument] Division/remainder operation counts:");
    let mut total = 0u64;
    for (name, count) in &counts {
        println!("  {:<30} {:>10}", name, count);
        total += count;
    }
    println!("  {:<30} {:>10}", "TOTAL", total);
}
