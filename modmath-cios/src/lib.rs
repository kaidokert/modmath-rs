//! Row-op trait surface for CIOS Montgomery multiplication.
//!
//! [`CiosRowOps`] is the minimal interface a multi-limb integer type
//! exposes so a CIOS (Coarsely Integrated Operand Scanning) Montgomery
//! multiplication loop can drive it: limb access plus the two fused
//! row kernels. The algorithm itself lives with the consumer; bigint
//! backends implement the trait on their own type.
//!
//! Single-word degenerate impls for `u8`/`u16`/`u32`/`u64` are
//! provided here, treating each primitive as a 1-limb value.

#![no_std]

/// Row-level operations for CIOS Montgomery multiplication.
///
/// Implementors provide infallible limb access and the two CIOS row
/// kernels; the multiplication loop driving them is generic over this
/// trait.
///
/// # Preconditions
///
/// - `word(i)` is called only with `i < self.word_count()`.
/// - `i` is a public value (loop counter, constant) — never secret.
///
/// Under these preconditions, an infallible accessor is CT-safe for
/// both variable-time and constant-time backends. A well-formed CIOS
/// driver upholds both by construction — `i` is only ever a public
/// loop counter (`0..word_count()`) or the constant `0`.
pub trait CiosRowOps: Default + Sized {
    type Word: Copy + PartialOrd;

    fn word_count(&self) -> usize;

    /// Infallible. Caller guarantees `i < self.word_count()` and `i`
    /// is public.
    fn word(&self, i: usize) -> Self::Word;

    /// Phase 1: `acc += scalar * multiplicand`. Returns carry-out.
    fn mul_acc_row(
        scalar: Self::Word,
        multiplicand: &Self,
        acc: &mut Self,
        carry_in: Self::Word,
    ) -> Self::Word;

    /// Phase 2: `[acc, acc_hi] = ([acc, acc_hi] + scalar * multiplicand) >> word_bits`.
    /// Returns the carry word (0 or 1) from the fold.
    fn mul_acc_shift_row(
        scalar: Self::Word,
        multiplicand: &Self,
        acc: &mut Self,
        acc_hi: Self::Word,
    ) -> Self::Word;
}

/// Single-word degenerate impl for primitive integers — treats `$narrow`
/// as a 1-limb "bigint" with `Word = $narrow`. Uses the natural `$wide`
/// companion for the row-op widening multiply.
macro_rules! impl_cios_row_ops_primitive {
    ($narrow:ty, $wide:ty, $bits:expr) => {
        impl CiosRowOps for $narrow {
            type Word = $narrow;

            fn word_count(&self) -> usize {
                1
            }

            fn word(&self, _i: usize) -> $narrow {
                *self
            }

            fn mul_acc_row(
                scalar: $narrow,
                multiplicand: &$narrow,
                acc: &mut $narrow,
                carry_in: $narrow,
            ) -> $narrow {
                let p = (scalar as $wide) * (*multiplicand as $wide)
                    + (*acc as $wide)
                    + (carry_in as $wide);
                *acc = p as $narrow;
                (p >> $bits) as $narrow
            }

            fn mul_acc_shift_row(
                scalar: $narrow,
                multiplicand: &$narrow,
                acc: &mut $narrow,
                acc_hi: $narrow,
            ) -> $narrow {
                let product = (scalar as $wide) * (*multiplicand as $wide);
                let total = (*acc as $wide) + product;
                let total_hi = total >> $bits;
                let (sum, carry) = total_hi.overflowing_add(acc_hi as $wide);
                *acc = sum as $narrow;
                ((sum >> $bits) as $narrow) + if carry { 1 } else { 0 }
            }
        }
    };
}

impl_cios_row_ops_primitive!(u8, u16, 8);
impl_cios_row_ops_primitive!(u16, u32, 16);
impl_cios_row_ops_primitive!(u32, u64, 32);
impl_cios_row_ops_primitive!(u64, u128, 64);

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! primitive_smoke {
        ($name_count:ident, $name_word:ident, $name_phase1:ident,
         $name_phase1_carry:ident, $name_phase2:ident, $t:ty) => {
            #[test]
            fn $name_count() {
                let x: $t = 0xAB;
                assert_eq!(CiosRowOps::word_count(&x), 1);
            }

            #[test]
            fn $name_word() {
                let x: $t = 0xAB;
                assert_eq!(CiosRowOps::word(&x, 0), 0xAB);
            }

            #[test]
            fn $name_phase1() {
                // acc = 0 + 3*4 = 12, no carry.
                let mult: $t = 3;
                let mut acc: $t = 0;
                let carry = <$t as CiosRowOps>::mul_acc_row(4, &mult, &mut acc, 0);
                assert_eq!(acc, 12);
                assert_eq!(carry, 0);
            }

            #[test]
            fn $name_phase1_carry() {
                // 2 * MAX + MAX = 3 * MAX = 3 * (2^bits - 1).
                // Low word = (3*2^bits - 3) mod 2^bits = MAX - 2; carry = 2.
                let mult: $t = <$t>::MAX;
                let mut acc: $t = <$t>::MAX;
                let carry = <$t as CiosRowOps>::mul_acc_row(2, &mult, &mut acc, 0);
                assert_eq!(acc, <$t>::MAX - 2);
                assert_eq!(carry, 2);
            }

            #[test]
            fn $name_phase2() {
                // 1-limb shift: acc' = (acc + s*m) >> bits, returns the fold.
                // acc = 0, s = 1, m = 0 → result word 0, carry 0.
                let mult: $t = 0;
                let mut acc: $t = 0;
                let carry = <$t as CiosRowOps>::mul_acc_shift_row(1, &mult, &mut acc, 0);
                assert_eq!(acc, 0);
                assert_eq!(carry, 0);
            }
        };
    }

    primitive_smoke!(
        u8_word_count,
        u8_word,
        u8_phase1,
        u8_phase1_carry,
        u8_phase2,
        u8
    );
    primitive_smoke!(
        u16_word_count,
        u16_word,
        u16_phase1,
        u16_phase1_carry,
        u16_phase2,
        u16
    );
    primitive_smoke!(
        u32_word_count,
        u32_word,
        u32_phase1,
        u32_phase1_carry,
        u32_phase2,
        u32
    );
    primitive_smoke!(
        u64_word_count,
        u64_word,
        u64_phase1,
        u64_phase1_carry,
        u64_phase2,
        u64
    );
}
