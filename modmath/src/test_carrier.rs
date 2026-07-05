//! Non-`Copy` test carrier for the constrained/strict flavors.
//!
//! The heap-backend dev-deps that would naturally pin "constrained and
//! strict work without `Copy`" are disabled during the const-num-traits
//! migration, so nothing else instantiates those flavors with a
//! non-`Copy` type. `NcU64` wraps `u64` with `Clone` but not `Copy`;
//! if a `Copy` bound sneaks back into either flavor, these tests stop
//! compiling.

use const_num_traits::{HasPersonality, Nct};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct NcU64(pub u64);

impl HasPersonality for NcU64 {
    type P = Nct;
}

impl const_num_traits::Zero for NcU64 {
    fn zero() -> Self {
        NcU64(0)
    }
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
    fn set_zero(&mut self) {
        self.0 = 0;
    }
}

impl const_num_traits::One for NcU64 {
    fn one() -> Self {
        NcU64(1)
    }
    fn set_one(&mut self) {
        self.0 = 1;
    }
    fn is_one(&self) -> bool {
        self.0 == 1
    }
}

impl const_num_traits::Parity for &NcU64 {
    fn is_odd(self) -> bool {
        self.0 & 1 == 1
    }
    fn is_even(self) -> bool {
        self.0 & 1 == 0
    }
}

impl const_num_traits::ops::wrapping::WrappingAdd for NcU64 {
    fn wrapping_add(self, v: Self) -> Self {
        NcU64(self.0.wrapping_add(v.0))
    }
}

impl const_num_traits::ops::wrapping::WrappingSub for NcU64 {
    fn wrapping_sub(self, v: Self) -> Self {
        NcU64(self.0.wrapping_sub(v.0))
    }
}

impl const_num_traits::ops::overflowing::OverflowingAdd for NcU64 {
    fn overflowing_add(self, v: Self) -> (Self, bool) {
        let (r, o) = self.0.overflowing_add(v.0);
        (NcU64(r), o)
    }
}

impl const_num_traits::ops::overflowing::OverflowingSub for NcU64 {
    fn overflowing_sub(self, v: Self) -> (Self, bool) {
        let (r, o) = self.0.overflowing_sub(v.0);
        (NcU64(r), o)
    }
}

impl core::ops::Add for NcU64 {
    type Output = NcU64;
    fn add(self, rhs: Self) -> Self {
        NcU64(self.0 + rhs.0)
    }
}

impl core::ops::Sub for NcU64 {
    type Output = NcU64;
    fn sub(self, rhs: Self) -> Self {
        NcU64(self.0 - rhs.0)
    }
}

impl core::ops::Mul for NcU64 {
    type Output = NcU64;
    fn mul(self, rhs: Self) -> Self {
        NcU64(self.0 * rhs.0)
    }
}

impl core::ops::Rem<&NcU64> for &NcU64 {
    type Output = NcU64;
    fn rem(self, rhs: &NcU64) -> NcU64 {
        NcU64(self.0 % rhs.0)
    }
}

impl core::ops::Rem<&NcU64> for NcU64 {
    type Output = NcU64;
    fn rem(self, rhs: &NcU64) -> NcU64 {
        NcU64(self.0 % rhs.0)
    }
}

impl core::ops::RemAssign<&NcU64> for NcU64 {
    fn rem_assign(&mut self, rhs: &NcU64) {
        self.0 %= rhs.0;
    }
}

impl core::ops::Shl<usize> for NcU64 {
    type Output = NcU64;
    fn shl(self, rhs: usize) -> Self {
        NcU64(self.0 << rhs)
    }
}

impl core::ops::Shr<usize> for NcU64 {
    type Output = NcU64;
    fn shr(self, rhs: usize) -> Self {
        NcU64(self.0 >> rhs)
    }
}

impl core::ops::ShrAssign<usize> for NcU64 {
    fn shr_assign(&mut self, rhs: usize) {
        self.0 >>= rhs;
    }
}

impl core::ops::BitAnd<&NcU64> for &NcU64 {
    type Output = NcU64;
    fn bitand(self, rhs: &NcU64) -> NcU64 {
        NcU64(self.0 & rhs.0)
    }
}

macro_rules! ncu64_binop {
    ($trait:ident, $method:ident, $op:tt) => {
        impl core::ops::$trait<NcU64> for &NcU64 {
            type Output = NcU64;
            fn $method(self, rhs: NcU64) -> NcU64 {
                NcU64(self.0 $op rhs.0)
            }
        }
        impl core::ops::$trait<&NcU64> for NcU64 {
            type Output = NcU64;
            fn $method(self, rhs: &NcU64) -> NcU64 {
                NcU64(self.0 $op rhs.0)
            }
        }
        impl core::ops::$trait<&NcU64> for &NcU64 {
            type Output = NcU64;
            fn $method(self, rhs: &NcU64) -> NcU64 {
                NcU64(self.0 $op rhs.0)
            }
        }
    };
}

ncu64_binop!(Add, add, +);
ncu64_binop!(Sub, sub, -);
ncu64_binop!(Mul, mul, *);
ncu64_binop!(Div, div, /);

impl core::ops::Div for NcU64 {
    type Output = NcU64;
    fn div(self, rhs: Self) -> Self {
        NcU64(self.0 / rhs.0)
    }
}

impl core::ops::AddAssign<&NcU64> for NcU64 {
    fn add_assign(&mut self, rhs: &NcU64) {
        self.0 += rhs.0;
    }
}

impl core::ops::DivAssign<&NcU64> for NcU64 {
    fn div_assign(&mut self, rhs: &NcU64) {
        self.0 /= rhs.0;
    }
}

#[cfg(test)]
mod tests {
    use super::NcU64;

    const M: u64 = 998_244_353; // prime
    const A: u64 = 123_456_789_012;
    const B: u64 = 987_654_321_098;

    #[test]
    fn schoolbook_constrained_matches_u64() {
        let m = NcU64(M);
        assert_eq!(
            crate::add::constrained_mod_add(NcU64(A), &NcU64(B), &m).0,
            crate::add::constrained_mod_add(A, &B, &M)
        );
        assert_eq!(
            crate::sub::constrained_mod_sub(NcU64(A), &NcU64(B), &m).0,
            crate::sub::constrained_mod_sub(A, &B, &M)
        );
        assert_eq!(
            crate::mul::constrained_mod_mul(NcU64(A), &NcU64(B), &m).0,
            crate::mul::constrained_mod_mul(A, &B, &M)
        );
        assert_eq!(
            crate::exp::constrained_mod_exp(NcU64(A), &NcU64(1000), &m).0,
            crate::exp::constrained_mod_exp(A, &1000u64, &M)
        );
    }

    #[test]
    fn schoolbook_strict_matches_u64() {
        let m = NcU64(M);
        assert_eq!(
            crate::add::strict_mod_add(NcU64(A), &NcU64(B), &m).0,
            crate::add::strict_mod_add(A, &B, &M)
        );
        assert_eq!(
            crate::sub::strict_mod_sub(NcU64(A), &NcU64(B), &m).0,
            crate::sub::strict_mod_sub(A, &B, &M)
        );
        assert_eq!(
            crate::mul::strict_mod_mul(NcU64(A), &NcU64(B), &m).0,
            crate::mul::strict_mod_mul(A, &B, &M)
        );
        assert_eq!(
            crate::exp::strict_mod_exp(NcU64(A), &NcU64(1000), &m).0,
            crate::exp::strict_mod_exp(A, &1000u64, &M)
        );
    }

    #[test]
    fn montgomery_constrained_matches_u64() {
        use crate::montgomery::constrained_mont::{
            constrained_montgomery_mod_exp, constrained_montgomery_mod_mul,
        };
        assert_eq!(
            constrained_montgomery_mod_mul(NcU64(A), &NcU64(B), &NcU64(M)).map(|v| v.0),
            constrained_montgomery_mod_mul(A, &B, &M)
        );
        assert_eq!(
            constrained_montgomery_mod_exp(NcU64(A), &NcU64(1000), &NcU64(M)).map(|v| v.0),
            constrained_montgomery_mod_exp(A, &1000u64, &M)
        );
    }

    #[test]
    fn montgomery_strict_matches_u64() {
        use crate::montgomery::strict_mont::{
            strict_montgomery_mod_exp, strict_montgomery_mod_mul,
        };
        assert_eq!(
            strict_montgomery_mod_mul(NcU64(A), &NcU64(B), &NcU64(M)).map(|v| v.0),
            strict_montgomery_mod_mul(A, &B, &M)
        );
        assert_eq!(
            strict_montgomery_mod_exp(NcU64(A), &NcU64(1000), &NcU64(M)).map(|v| v.0),
            strict_montgomery_mod_exp(A, &1000u64, &M)
        );
    }
}
