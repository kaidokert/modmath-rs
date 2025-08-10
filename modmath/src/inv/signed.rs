use core::cmp::Ordering;

/// Partial signed implementation for modular inverse calculations
///
/// This type represents a signed number using a magnitude-sign representation.
///
/// # Type Parameter Requirements
///
/// **IMPORTANT**: The type parameter `T` **MUST** be an unsigned/non-negative type
/// (such as `u32`, `u64`, `BigUint`, etc.). Using signed types (like `i32`, `i64`)
/// will break the internal invariant that `value` represents a non-negative magnitude.
///
/// # Representation Invariant
///
/// - `value`: Always represents the **absolute magnitude** (non-negative)
/// - `negative`: The sign bit (`false` = positive, `true` = negative)
/// - Zero is canonicalized to have `negative = false`
///
/// # Examples
///
/// ```rust,ignore
/// // Note: This example uses private API for demonstration
/// use crate::inv::signed::Signed;
///
/// // Correct usage with unsigned types
/// let positive_five = Signed::new(5u32, false);  // +5
/// let negative_three = Signed::new(3u32, true);  // -3
/// let zero = Signed::new(0u32, false);           // +0 (canonical)
///
/// // DO NOT use with signed types - breaks invariants!
/// // let bad_example = Signed::new(-5i32, false); // ❌ WRONG!
/// ```
///
/// Only here to support mod_inv minimal operations.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Signed<T> {
    /// The absolute magnitude of the number (must be non-negative)
    value: T,
    /// The sign: false for positive/zero, true for negative
    negative: bool,
}

impl<T> From<T> for Signed<T> {
    /// Creates a positive Signed<T> from an unsigned value.
    ///
    /// # Requirements
    /// `T` must be an unsigned type (u32, u64, BigUint, etc.)
    fn from(value: T) -> Self {
        Self {
            value,
            negative: false,
        }
    }
}

impl<T> Signed<T> {
    /// Creates a new Signed<T> with the specified magnitude and sign.
    ///
    /// # Requirements
    /// - `value` must represent a non-negative magnitude
    /// - `T` must be an unsigned type (u32, u64, BigUint, etc.)
    ///
    /// # Arguments
    /// - `value`: The absolute magnitude (must be non-negative)
    /// - `negative`: The sign (false = positive, true = negative)
    ///
    /// # Zero Canonicalization
    /// If `value` is zero, `negative` will automatically be set to `false`
    /// regardless of the input, maintaining the invariant that zero is always positive.
    /// This requires `T: PartialEq + num_traits::Zero`.
    pub fn new(value: T, negative: bool) -> Self
    where
        T: PartialEq + num_traits::Zero,
    {
        let is_zero = value == T::zero();
        Self {
            value,
            negative: if is_zero { false } else { negative },
        }
    }

    /// Internal constructor that bypasses zero canonicalization for performance.
    /// Used internally where we know the canonicalization is already handled.
    ///
    /// # Safety
    /// This should only be used when the caller ensures zero canonicalization
    /// is handled elsewhere or when creating non-zero values.
    pub(crate) fn new_unchecked(value: T, negative: bool) -> Self {
        Self { value, negative }
    }

    /// Returns the inner magnitude value.
    ///
    /// Note: This returns only the magnitude; use the `negative` field
    /// to determine the actual sign of the number.
    pub fn into_inner(self) -> T {
        self.value
    }
}

impl<T> PartialOrd for Signed<T>
where
    T: PartialOrd,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self.negative, other.negative) {
            (false, true) => Some(Ordering::Greater),
            (true, false) => Some(Ordering::Less),
            (true, true) => other.value.partial_cmp(&self.value), // Flip comparison for negative numbers
            (false, false) => self.value.partial_cmp(&other.value),
        }
    }
}

impl<T> core::ops::Add<&Signed<T>> for Signed<T>
where
    for<'a> T: core::ops::Add<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::cmp::PartialOrd
        + PartialEq
        + num_traits::Zero,
    for<'a> &'a T: core::ops::Sub<T, Output = T>,
{
    type Output = Self;

    fn add(self, other: &Signed<T>) -> Self::Output {
        // can i delegate this implementation ??
        let mut result = match (self.negative, other.negative) {
            (false, false) => Self::new_unchecked(self.value + &other.value, false),
            (true, true) => Self::new_unchecked(self.value + &other.value, true),
            (false, true) | (true, false) => {
                if self.value >= other.value {
                    Self::new_unchecked(self.value - &other.value, self.negative)
                } else {
                    Self::new_unchecked(&other.value - self.value, other.negative)
                }
            }
        };

        // Canonicalize zero: ensure zero is always represented as positive
        if result.value == T::zero() {
            result.negative = false;
        }

        result
    }
}

impl<'a, T> core::ops::AddAssign<&'a T> for Signed<T>
where
    T: core::ops::AddAssign<&'a T>
        + core::ops::SubAssign<&'a T>
        + PartialEq
        + PartialOrd
        + num_traits::Zero
        + Clone,
    for<'b> &'b T: core::ops::Sub<T, Output = T>,
{
    fn add_assign(&mut self, other: &'a T) {
        match self.negative {
            false => self.value += other,
            true => {
                // Handle negative case: -self.value + other
                if &self.value >= other {
                    // |self| >= other, result stays negative: -(|self| - other)
                    self.value -= other;
                } else {
                    // |self| < other, result becomes positive: other - |self|
                    self.value = other - self.value.clone();
                    self.negative = false; // Flip sign to positive
                }
            }
        }

        // Canonicalize zero: ensure zero is always represented as positive
        if self.value == T::zero() {
            self.negative = false;
        }
    }
}

impl<T> core::ops::Add<&T> for Signed<T>
where
    for<'a> T: core::ops::Add<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::cmp::PartialOrd
        + num_traits::Zero
        + PartialEq,
    for<'a> &'a T: core::ops::Sub<T, Output = T>,
{
    type Output = Self;

    fn add(self, other: &T) -> Self::Output {
        let mut result = match self.negative {
            false => Self::new_unchecked(self.value + other, self.negative),
            true => {
                if self.value >= *other {
                    Self::new_unchecked(self.value - other, self.negative)
                } else {
                    Self::new_unchecked(other - self.value, !self.negative)
                }
            }
        };

        // Canonicalize zero: ensure zero is always represented as positive
        if result.value == T::zero() {
            result.negative = false;
        }

        result
    }
}

impl<T> core::ops::Add for Signed<T>
where
    T: core::ops::Add<Output = T> + core::ops::Sub<Output = T> + core::cmp::PartialOrd,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        match (self.negative, rhs.negative) {
            (false, false) => Self::new_unchecked(self.value + rhs.value, false),
            (true, true) => Self::new_unchecked(self.value + rhs.value, true),
            (false, true) | (true, false) => {
                if self.value >= rhs.value {
                    Self::new_unchecked(self.value - rhs.value, self.negative)
                } else {
                    Self::new_unchecked(rhs.value - self.value, rhs.negative)
                }
            }
        }
    }
}

impl<T> core::ops::Sub for Signed<T>
where
    T: core::ops::Add<Output = T> + core::ops::Sub<Output = T> + core::cmp::PartialOrd,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        // a - b = a + (-b)
        // To negate a number, just flip its sign
        self + Self::new_unchecked(rhs.value, !rhs.negative)
    }
}

/// Multiplication of two Signed<T> values.
///
/// # Requirements
/// `T` must be an unsigned type to maintain the invariant that `value` represents
/// a non-negative magnitude. Using signed types will break this assumption.
impl<T> core::ops::Mul for Signed<T>
where
    T: core::ops::Mul<Output = T>,
{
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        Self {
            value: self.value * other.value,
            negative: self.negative != other.negative,
        }
    }
}

/// Multiplication of Signed<T> with an unsigned value T.
///
/// # Requirements
/// Both `T` and `other` must be unsigned types. The `other` parameter should
/// represent a non-negative value to maintain correctness.
impl<T> core::ops::Mul<T> for Signed<T>
where
    T: core::ops::Mul<Output = T>,
{
    type Output = Self;

    fn mul(self, other: T) -> Self::Output {
        Self {
            value: self.value * other,
            negative: self.negative,
        }
    }
}

/// Multiplication of Signed<T> with a reference to an unsigned value &T.
///
/// # Requirements
/// Both `T` and `other` must be unsigned types. The `other` parameter should
/// represent a non-negative value to maintain correctness.
impl<'a, T> core::ops::Mul<&'a T> for Signed<T>
where
    T: core::ops::Mul<&'a T, Output = T> + 'a,
{
    type Output = Self;

    fn mul(self, other: &'a T) -> Self::Output {
        Self {
            value: self.value * other,
            negative: self.negative,
        }
    }
}

impl<T> core::ops::Div for Signed<T>
where
    T: core::ops::Div<Output = T>,
{
    type Output = Self;

    fn div(self, other: Self) -> Self::Output {
        Self {
            value: self.value / other.value,
            negative: self.negative != other.negative,
        }
    }
}

#[cfg(test)]
mod signed_tests {
    use super::Signed;

    #[test]
    fn test_signed_test_greater_than() {
        fn check(a: u32, a_neg: bool, b: u32, b_neg: bool, expected: bool) {
            let a = Signed::new(a, a_neg);
            let b = Signed::new(b, b_neg);
            let a_greater_than_b = a > b;
            assert_eq!(a_greater_than_b, expected);
        }

        // 1 , 1 , a_greater_than_b = false
        check(1, false, 1, false, false);

        // 1 , 2 , a_greater_than_b = false
        check(1, false, 2, false, false);
        // -1 , 2 , a_greater_than_b = false
        check(1, true, 2, false, false);
        // 1, -2 , a_greater_than_b = false
        check(1, false, 2, true, true);
        // -1, -2 , a_greater_than_b = false
        check(1, true, 2, true, true);

        // 2 , 1 , a_greater_than_b = true
        check(2, false, 1, false, true);
        // -2 , 1 , a_greater_than_b = false
        check(2, true, 1, false, false);
        // 2, -1 , a_greater_than_b = true
        check(2, false, 1, true, true);
        // -2, -1 , a_greater_than_b = false
        check(2, true, 1, true, false);
    }

    #[test]
    fn test_signed_add() {
        fn check(a: u32, a_neg: bool, b: u32, b_neg: bool, expected_val: u32, expected_sign: bool) {
            let a = Signed::new(a, a_neg);
            let b = Signed::new(b, b_neg);
            let result = a + b;

            assert_eq!(result.value, expected_val);
            assert_eq!(result.negative, expected_sign);
        }
        // 1 + 2 = 3
        check(1, false, 2, false, 3, false);
        // 1 + (-2) = -1
        check(1, false, 2, true, 1, true);
        // -1 + 2 = 1
        check(1, true, 2, false, 1, false);
        // -1 + (-2) = -3
        check(1, true, 2, true, 3, true);
        // -2 + 1 = -1
        check(2, true, 1, false, 1, true);
        // -2 + (-1) = -3
        check(2, true, 1, true, 3, true);
        // 2 + 1 = 3
        check(2, false, 1, false, 3, false);
        // 2 + (-1) = 1
        check(2, false, 1, true, 1, false);
    }

    #[test]
    fn test_signed_sub() {
        fn check(a: u32, a_neg: bool, b: u32, b_neg: bool, expected_val: u32, expected_sign: bool) {
            let a = Signed::new(a, a_neg);
            let b = Signed::new(b, b_neg);
            let result = a - b;

            assert_eq!(result.value, expected_val);
            assert_eq!(result.negative, expected_sign);
        }
        // 1 - 2 = -1
        check(1, false, 2, false, 1, true);
        // 1 - (-2) = 3
        check(1, false, 2, true, 3, false);
        // -1 - 2 = -3
        check(1, true, 2, false, 3, true);
        // -1 - (-2) = 1
        check(1, true, 2, true, 1, false);
        // -2 - 1 = -3
        check(2, true, 1, false, 3, true);
        // -2 - (-1) = -1
        check(2, true, 1, true, 1, true);
        // 2 - 1 = 1
        check(2, false, 1, false, 1, false);
        // 2 - (-1) = 3
        check(2, false, 1, true, 3, false);
    }

    #[test]
    fn test_signed_div() {
        fn check(a: u32, a_neg: bool, b: u32, b_neg: bool, expected_val: u32, expected_sign: bool) {
            let a = Signed::new(a, a_neg);
            let b = Signed::new(b, b_neg);
            let result = a / b;

            assert_eq!(result.value, expected_val);
            assert_eq!(result.negative, expected_sign);
        }

        check(1, false, 2, false, 0, false);
        check(1, false, 2, true, 0, true);
        check(1, true, 2, false, 0, true);
    }

    #[test]
    fn test_signed_from_and_into_inner() {
        // Test From<T> implementation
        let signed_from_42 = Signed::from(42u32);
        assert_eq!(signed_from_42.value, 42u32);
        assert_eq!(signed_from_42.negative, false);

        // Test into_inner() method
        let signed_value = Signed::new(123u32, true);
        assert_eq!(signed_value.into_inner(), 123u32);
    }

    #[test]
    fn test_signed_new_zero_canonicalization() {
        // Test that new() automatically canonicalizes zero to positive
        let zero_positive = Signed::new(0u32, false);
        assert_eq!(zero_positive.value, 0u32);
        assert_eq!(zero_positive.negative, false);

        let zero_negative_input = Signed::new(0u32, true);
        assert_eq!(zero_negative_input.value, 0u32);
        assert_eq!(zero_negative_input.negative, false); // Should be canonicalized to positive

        // Test with non-zero values - should preserve the sign
        let positive_five = Signed::new(5u32, false);
        assert_eq!(positive_five.value, 5u32);
        assert_eq!(positive_five.negative, false);

        let negative_three = Signed::new(3u32, true);
        assert_eq!(negative_three.value, 3u32);
        assert_eq!(negative_three.negative, true); // Should preserve negative for non-zero
    }

    #[test]
    fn test_signed_add_with_reference() {
        // Test Add<&Signed<T>> implementation
        let a = Signed::new(5u32, false); // +5
        let b = Signed::new(3u32, true); // -3
        let result = a + &b; // +5 + (-3) = +2
        assert_eq!(result.value, 2u32);
        assert_eq!(result.negative, false);

        // Test case where other > self in mixed signs
        let a = Signed::new(2u32, false); // +2
        let b = Signed::new(5u32, true); // -5
        let result = a + &b; // +2 + (-5) = -3
        assert_eq!(result.value, 3u32);
        assert_eq!(result.negative, true);

        // Test both negative
        let a = Signed::new(4u32, true); // -4
        let b = Signed::new(6u32, true); // -6
        let result = a + &b; // -4 + (-6) = -10
        assert_eq!(result.value, 10u32);
        assert_eq!(result.negative, true);

        // Test zero canonicalization: positive + negative with equal magnitude
        let a = Signed::new(5u32, false); // +5
        let b = Signed::new(5u32, true); // -5
        let result = a + &b; // +5 + (-5) = 0
        assert_eq!(result.value, 0u32);
        assert_eq!(result.negative, false); // Fixed: zero should be canonicalized to positive

        // Test zero canonicalization: negative + positive with equal magnitude
        let a = Signed::new(7u32, true); // -7
        let b = Signed::new(7u32, false); // +7
        let result = a + &b; // -7 + 7 = 0
        assert_eq!(result.value, 0u32);
        assert_eq!(result.negative, false); // Fixed: zero should be canonicalized to positive

        // Additional zero canonicalization test cases for Add<&Signed<T>>

        // Test zero canonicalization with different values: +12 + (-12) = 0
        let a = Signed::new(12u32, false); // +12
        let b = Signed::new(12u32, true); // -12
        let result = a + &b; // +12 + (-12) = 0
        assert_eq!(result.value, 0u32);
        assert_eq!(result.negative, false); // Zero should be canonicalized to positive

        // Test zero canonicalization with different values: -25 + (+25) = 0
        let a = Signed::new(25u32, true); // -25
        let b = Signed::new(25u32, false); // +25
        let result = a + &b; // -25 + 25 = 0
        assert_eq!(result.value, 0u32);
        assert_eq!(result.negative, false); // Zero should be canonicalized to positive

        // Test zero canonicalization with larger values: +100 + (-100) = 0
        let a = Signed::new(100u32, false); // +100
        let b = Signed::new(100u32, true); // -100
        let result = a + &b; // +100 + (-100) = 0
        assert_eq!(result.value, 0u32);
        assert_eq!(result.negative, false); // Zero should be canonicalized to positive

        // Test zero canonicalization with edge case: +1 + (-1) = 0
        let a = Signed::new(1u32, false); // +1
        let b = Signed::new(1u32, true); // -1
        let result = a + &b; // +1 + (-1) = 0
        assert_eq!(result.value, 0u32);
        assert_eq!(result.negative, false); // Zero should be canonicalized to positive

        // Test zero canonicalization with maximum value: +u32::MAX + (-u32::MAX) = 0
        let a = Signed::new(u32::MAX, false); // +u32::MAX
        let b = Signed::new(u32::MAX, true); // -u32::MAX
        let result = a + &b; // +u32::MAX + (-u32::MAX) = 0
        assert_eq!(result.value, 0u32);
        assert_eq!(result.negative, false); // Zero should be canonicalized to positive
    }

    #[test]
    fn test_signed_add_with_t_reference() {
        // Test Add<&T> implementation - positive Signed + T
        let a = Signed::new(7u32, false); // +7
        let result = a + &3u32; // +7 + 3 = +10
        assert_eq!(result.value, 10u32);
        assert_eq!(result.negative, false);

        // Test negative Signed + T where |signed| >= T
        let a = Signed::new(8u32, true); // -8
        let result = a + &3u32; // -8 + 3 = -5
        assert_eq!(result.value, 5u32);
        assert_eq!(result.negative, true);

        // Test negative Signed + T where |signed| < T (should flip sign)
        let a = Signed::new(2u32, true); // -2
        let result = a + &5u32; // -2 + 5 = +3
        assert_eq!(result.value, 3u32);
        assert_eq!(result.negative, false);
    }

    #[test]
    fn test_signed_add_assign() {
        // Test AddAssign<&T> with positive Signed
        let mut a = Signed::new(10u32, false); // +10
        a += &3u32; // +10 += 3 = +13
        assert_eq!(a.value, 13u32);
        assert_eq!(a.negative, false);

        // Test AddAssign<&T> with negative Signed
        let mut a = Signed::new(15u32, true); // -15
        a += &5u32; // -15 += 5 = -10 (subtracts for negative)
        assert_eq!(a.value, 10u32);
        assert_eq!(a.negative, true);

        // Test exact cancellation case (this should work without underflow)
        let mut a = Signed::new(7u32, true); // -7
        a += &7u32; // -7 += 7 = 0 (7 - 7 = 0)
        assert_eq!(a.value, 0u32);
        assert_eq!(a.negative, false); // Fixed: zero is canonicalized to positive

        // Test another valid case where subtraction doesn't underflow
        let mut a = Signed::new(20u32, true); // -20
        a += &15u32; // -20 += 15 = -5 (20 - 15 = 5, stays negative)
        assert_eq!(a.value, 5u32);
        assert_eq!(a.negative, true);

        // Test with zero addition to positive
        let mut a = Signed::new(5u32, false); // +5
        a += &0u32; // +5 += 0 = +5
        assert_eq!(a.value, 5u32);
        assert_eq!(a.negative, false);

        // Test with zero addition to negative
        let mut a = Signed::new(8u32, true); // -8
        a += &0u32; // -8 += 0 = -8 (8 - 0 = 8, stays negative)
        assert_eq!(a.value, 8u32);
        assert_eq!(a.negative, true);

        // Extended test cases for sign flips and zero canonicalization:

        // These tests verify the fixed AddAssign implementation.
        // The fixed implementation handles sign flips properly and
        // canonicalizes zero to always be positive.

        // Test case: Zero result should be canonicalized to positive
        let mut zero_case = Signed::new(5u32, true); // -5
        zero_case += &5u32; // -5 + 5 should be +0 (canonicalized)
        assert_eq!(zero_case.value, 0u32);
        assert_eq!(zero_case.negative, false); // Fixed: zero is canonicalized to positive

        // Test case: Multiple operations that should result in zero
        let mut multiple_zero = Signed::new(3u32, true); // -3
        multiple_zero += &3u32; // -3 + 3 = 0
        assert_eq!(multiple_zero.value, 0u32);
        assert_eq!(multiple_zero.negative, false); // Fixed: zero is canonicalized to positive

        // Test case: Large exact cancellation
        let mut large_cancel = Signed::new(100u32, true); // -100
        large_cancel += &100u32; // -100 + 100 = 0
        assert_eq!(large_cancel.value, 0u32);
        assert_eq!(large_cancel.negative, false); // Fixed: zero is canonicalized to positive

        // Additional edge cases with zero
        let mut zero_to_zero = Signed::new(0u32, true); // -0 (invalid state)
        zero_to_zero += &0u32; // -0 + 0 = 0
        assert_eq!(zero_to_zero.value, 0u32);
        assert_eq!(zero_to_zero.negative, false); // Fixed: zero is canonicalized to positive

        // Test cases that should now work: sign flips from negative to positive
        let mut sign_flip1 = Signed::new(8u32, true); // -8
        sign_flip1 += &12u32; // -8 + 12 = +4
        assert_eq!(sign_flip1.value, 4u32);
        assert_eq!(sign_flip1.negative, false); // Should become positive

        let mut sign_flip2 = Signed::new(3u32, true); // -3
        sign_flip2 += &7u32; // -3 + 7 = +4
        assert_eq!(sign_flip2.value, 4u32);
        assert_eq!(sign_flip2.negative, false); // Should become positive

        let mut sign_flip3 = Signed::new(1u32, true); // -1
        sign_flip3 += &10u32; // -1 + 10 = +9
        assert_eq!(sign_flip3.value, 9u32);
        assert_eq!(sign_flip3.negative, false); // Should become positive
    }

    #[test]
    fn test_signed_add_assign_sign_flip_fixed() {
        // This test verifies the fixed AddAssign implementation
        // When adding to a negative number causes sign flip, it should work correctly
        let mut a = Signed::new(8u32, true); // -8
        a += &12u32; // -8 + 12 = +4
        assert_eq!(a.value, 4u32);
        assert_eq!(a.negative, false); // Should be positive
    }

    #[test]
    fn test_signed_add_assign_sign_flip_negative_to_positive_fixed() {
        // Test case that flips from negative to positive - now working correctly
        // -3 + 7 should equal +4
        let mut a = Signed::new(3u32, true); // -3
        a += &7u32; // Should become +4
        assert_eq!(a.value, 4u32);
        assert_eq!(a.negative, false); // Should be positive
    }

    #[test]
    fn test_signed_add_assign_small_negative_large_positive_fixed() {
        // Test case: small negative + large positive should become positive
        // -1 + 10 should equal +9
        let mut a = Signed::new(1u32, true); // -1
        a += &10u32; // Should become +9
        assert_eq!(a.value, 9u32);
        assert_eq!(a.negative, false); // Should be positive
    }

    #[test]
    fn test_signed_add_assign_edge_case_maximum_flip_fixed() {
        // Test case: maximum negative flip
        // -1 + u32::MAX should become +(u32::MAX - 1)
        let mut a = Signed::new(1u32, true); // -1
        a += &u32::MAX; // Should become +(u32::MAX - 1)
        assert_eq!(a.value, u32::MAX - 1);
        assert_eq!(a.negative, false); // Should be positive
    }

    #[test]
    fn test_signed_add_t_reference_correct_behavior() {
        // This test verifies the correct behavior of Add<&T> implementation
        // Both sign flips and zero canonicalization should work properly

        // Case 1: -8 + 12 should equal +4
        let a = Signed::new(8u32, true); // -8
        let result = a + &12u32; // This works with Add<&T>
        assert_eq!(result.value, 4u32);
        assert_eq!(result.negative, false); // Correctly becomes positive

        // Case 2: -3 + 7 should equal +4
        let a = Signed::new(3u32, true); // -3
        let result = a + &7u32;
        assert_eq!(result.value, 4u32);
        assert_eq!(result.negative, false); // Correctly becomes positive

        // Case 3: -1 + 10 should equal +9
        let a = Signed::new(1u32, true); // -1
        let result = a + &10u32;
        assert_eq!(result.value, 9u32);
        assert_eq!(result.negative, false); // Correctly becomes positive

        // Case 4: Zero canonicalization should work correctly now
        let a = Signed::new(5u32, true); // -5
        let result = a + &5u32;
        assert_eq!(result.value, 0u32);
        assert_eq!(result.negative, false); // Fixed: zero is now canonicalized to positive

        // Case 5: Test edge case where negative value equals the addend
        let a = Signed::new(7u32, true); // -7
        let result = a + &7u32; // -7 + 7 = 0
        assert_eq!(result.value, 0u32);
        assert_eq!(result.negative, false); // Fixed: zero is now canonicalized to positive

        // Case 6: Show that positive + positive works correctly
        let a = Signed::new(3u32, false); // +3
        let result = a + &5u32; // +3 + 5 = +8
        assert_eq!(result.value, 8u32);
        assert_eq!(result.negative, false); // Correctly positive
    }

    #[test]
    fn test_signed_multiplication() {
        // Test Mul<Self> - positive * positive
        let a = Signed::new(4u32, false); // +4
        let b = Signed::new(3u32, false); // +3
        let result = a * b; // +4 * +3 = +12
        assert_eq!(result.value, 12u32);
        assert_eq!(result.negative, false);

        // Test Mul<Self> - positive * negative
        let a = Signed::new(5u32, false); // +5
        let b = Signed::new(2u32, true); // -2
        let result = a * b; // +5 * -2 = -10
        assert_eq!(result.value, 10u32);
        assert_eq!(result.negative, true);

        // Test Mul<Self> - negative * negative
        let a = Signed::new(6u32, true); // -6
        let b = Signed::new(4u32, true); // -4
        let result = a * b; // -6 * -4 = +24
        assert_eq!(result.value, 24u32);
        assert_eq!(result.negative, false);

        // Test Mul<T> - positive Signed * T
        let a = Signed::new(7u32, false); // +7
        let result = a * 2u32; // +7 * 2 = +14
        assert_eq!(result.value, 14u32);
        assert_eq!(result.negative, false);

        // Test Mul<T> - negative Signed * T
        let a = Signed::new(3u32, true); // -3
        let result = a * 5u32; // -3 * 5 = -15
        assert_eq!(result.value, 15u32);
        assert_eq!(result.negative, true);

        // Test Mul<&T> - positive Signed * &T
        let a = Signed::new(8u32, false); // +8
        let result = a * &3u32; // +8 * 3 = +24
        assert_eq!(result.value, 24u32);
        assert_eq!(result.negative, false);

        // Test Mul<&T> - negative Signed * &T
        let a = Signed::new(9u32, true); // -9
        let result = a * &2u32; // -9 * 2 = -18
        assert_eq!(result.value, 18u32);
        assert_eq!(result.negative, true);
    }

    #[test]
    fn test_signed_partial_cmp_edge_cases() {
        use core::cmp::Ordering;

        // Test equal values with same sign
        let a = Signed::new(5u32, false);
        let b = Signed::new(5u32, false);
        assert_eq!(a.partial_cmp(&b), Some(Ordering::Equal));

        let a = Signed::new(5u32, true);
        let b = Signed::new(5u32, true);
        assert_eq!(a.partial_cmp(&b), Some(Ordering::Equal));

        // Test negative vs positive (different magnitudes)
        let a = Signed::new(10u32, true); // -10
        let b = Signed::new(1u32, false); // +1
        assert_eq!(a.partial_cmp(&b), Some(Ordering::Less));

        // Test positive vs negative (different magnitudes)
        let a = Signed::new(1u32, false); // +1
        let b = Signed::new(10u32, true); // -10
        assert_eq!(a.partial_cmp(&b), Some(Ordering::Greater));

        // Test both negative with different magnitudes (flipped comparison)
        let a = Signed::new(3u32, true); // -3
        let b = Signed::new(7u32, true); // -7
        assert_eq!(a.partial_cmp(&b), Some(Ordering::Greater)); // -3 > -7
    }
}
