use core::cmp::Ordering;

/// Partial signed implementation
/// Only here to support mod_inv minimal operations

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Signed<T> {
    value: T,
    negative: bool,
}

impl<T> From<T> for Signed<T> {
    fn from(value: T) -> Self {
        Self {
            value,
            negative: false,
        }
    }
}

impl<T> Signed<T> {
    pub fn new(value: T, negative: bool) -> Self {
        Self { value, negative }
    }
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
        + core::cmp::PartialOrd,
    for<'a> &'a T: core::ops::Sub<T, Output = T>,
{
    type Output = Self;

    fn add(self, other: &Signed<T>) -> Self::Output {
        // can i delegate this implementation ??
        match (self.negative, other.negative) {
            (false, false) => Self::new(self.value + &other.value, false),
            (true, true) => Self::new(self.value + &other.value, true),
            (false, true) | (true, false) => {
                if self.value >= other.value {
                    Self::new(self.value - &other.value, self.negative)
                } else {
                    Self::new(&other.value - self.value, other.negative)
                }
            }
        }
    }
}

impl<'a, T> core::ops::AddAssign<&'a T> for Signed<T>
where
    T: core::ops::AddAssign<&'a T> + core::ops::SubAssign<&'a T>,
{
    fn add_assign(&mut self, other: &'a T) {
        match self.negative {
            false => self.value += other,
            true => self.value -= other,
        }
    }
}

impl<T> core::ops::Add<&T> for Signed<T>
where
    for<'a> T: core::ops::Add<&'a T, Output = T>
        + core::ops::Sub<&'a T, Output = T>
        + core::cmp::PartialOrd,
    for<'a> &'a T: core::ops::Sub<T, Output = T>,
{
    type Output = Self;

    fn add(self, other: &T) -> Self::Output {
        match self.negative {
            false => Self::new(self.value + other, self.negative),
            true => {
                if self.value >= *other {
                    Self::new(self.value - other, self.negative)
                } else {
                    Self::new(other - self.value, !self.negative)
                }
            }
        }
    }
}

impl<T> core::ops::Add for Signed<T>
where
    T: core::ops::Add<Output = T> + core::ops::Sub<Output = T> + core::cmp::PartialOrd,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        match (self.negative, rhs.negative) {
            (false, false) => Self::new(self.value + rhs.value, false),
            (true, true) => Self::new(self.value + rhs.value, true),
            (false, true) | (true, false) => {
                if self.value >= rhs.value {
                    Self::new(self.value - rhs.value, self.negative)
                } else {
                    Self::new(rhs.value - self.value, rhs.negative)
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
        self + Self::new(rhs.value, !rhs.negative)
    }
}

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
