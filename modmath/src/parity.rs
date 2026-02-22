/// Trait for checking parity (odd/even) of integer types.
///
/// When the `num-integer` feature is enabled, this delegates to
/// `num_integer::Integer::is_odd` which can be optimized per-type
/// (e.g. fixed-bigint checks a single byte instead of all limbs).
///
/// Without the feature, falls back to `self & one == one` which
/// touches all limbs via BitAnd.
pub trait Parity {
    fn is_odd(&self) -> bool;
    fn is_even(&self) -> bool {
        !self.is_odd()
    }
}

#[cfg(not(feature = "num-integer"))]
impl<T> Parity for T
where
    T: core::ops::BitAnd<Output = T> + num_traits::One + PartialEq + Clone,
{
    fn is_odd(&self) -> bool {
        self.clone() & T::one() == T::one()
    }
}

#[cfg(feature = "num-integer")]
impl<T> Parity for T
where
    T: num_integer::Integer,
{
    fn is_odd(&self) -> bool {
        num_integer::Integer::is_odd(self)
    }
}
