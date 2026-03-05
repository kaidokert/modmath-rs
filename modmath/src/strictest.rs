/// Fast modular addition assuming inputs are already reduced to [0, m)
/// SAFETY: Caller must ensure a, b ∈ [0, m)
pub fn strictest_mod_add_no_reduce<T>(a: &T, b: &T, m: &T) -> T
where
    T: core::cmp::PartialOrd
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub,
{
    debug_assert!(a < m, "Input a must be reduced: a < m");
    debug_assert!(b < m, "Input b must be reduced: b < m");
    
    let (sum, overflow) = a.overflowing_add(b);
    if &sum >= m || overflow {
        sum.overflowing_sub(m).0
    } else {
        sum
    }
}

/// Safe modular addition that reduces inputs first
pub fn strictest_mod_add<T>(a: &T, b: &T, m: &T) -> T
where
    T: core::cmp::PartialOrd
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub,
    for<'b> &'b T: core::ops::Rem<&'b T, Output = T>,
{
    // Reduce operands first to prevent overflow
    let a_reduced = a % m;
    let b_reduced = b % m;
    strictest_mod_add_no_reduce(&a_reduced, &b_reduced, m)
}

/// Fast modular subtraction assuming inputs are already reduced to [0, m)
/// SAFETY: Caller must ensure a, b ∈ [0, m)
pub fn strictest_mod_sub_no_reduce<T>(a: &T, b: &T, m: &T) -> T
where
    T: core::cmp::PartialOrd
        + num_traits::ops::overflowing::OverflowingSub,
{
    debug_assert!(a < m, "Input a must be reduced: a < m");
    debug_assert!(b < m, "Input b must be reduced: b < m");
    
    if a >= b {
        // No underflow case: simply subtract
        let (diff, _) = a.overflowing_sub(b);
        diff
    } else {
        // Underflow case: compute m - (b - a)
        // This is equivalent to (a - b) mod m in positive representation
        let (b_minus_a, _) = b.overflowing_sub(a);
        let (result, _) = m.overflowing_sub(&b_minus_a);
        result
    }
}

/// Safe modular subtraction that reduces inputs first
pub fn strictest_mod_sub<T>(a: &T, b: &T, m: &T) -> T
where
    T: core::cmp::PartialOrd
        + num_traits::ops::overflowing::OverflowingSub
        + num_traits::ops::overflowing::OverflowingAdd,
    for<'b> &'b T: core::ops::Rem<&'b T, Output = T>,
{
    // Reduce operands first
    let a_reduced = a % m;
    let b_reduced = b % m;
    strictest_mod_sub_no_reduce(&a_reduced, &b_reduced, m)
}

pub fn strictest_mod_mul<T>(a: &T, b: &T, m: &T) -> T
where
    T: Clone
        + num_traits::Zero
        + num_traits::One
        + num_traits::ops::overflowing::OverflowingAdd
        + num_traits::ops::overflowing::OverflowingSub
        + core::cmp::PartialOrd
        + core::ops::Shr<usize, Output = T>,
    for<'b> &'b T: core::ops::Rem<&'b T, Output = T> + core::ops::BitAnd<Output = T>,
{
    // Reduce operands first to prevent overflow
    let mut a_cur = a % m;
    let mut b_cur = b % m;
    let mut result = T::zero();
    let one = T::one();
    
    // Binary multiplication algorithm (peasant multiplication)
    while b_cur > T::zero() {
        // If current bit of b is set, add current a to result
        if &b_cur & &one == one {
            let (sum, overflow) = result.overflowing_add(&a_cur);
            result = if &sum >= m || overflow {
                sum.overflowing_sub(m).0
            } else {
                sum
            };
        }
        
        // Double a_cur (mod m) and halve b_cur
        let (doubled, overflow) = a_cur.overflowing_add(&a_cur);
        a_cur = if &doubled >= m || overflow {
            doubled.overflowing_sub(m).0
        } else {
            doubled
        };
        b_cur = b_cur >> 1;
    }
    result
}

#[cfg(test)]
mod strictest_mod_add_tests {
    use super::*;

    #[test]
    fn test_strictest_mod_add() {
        // 10 + 20 = 30, 30 % 30 = 0
        let a = 10;
        let b = 20;
        let m = 30;
        let result = strictest_mod_add(&a, &b, &m);
        assert_eq!(result, 0);
    }
    
    #[test]
    fn test_strictest_mod_sub() {
        // 10 - 20 = -10, -10 % 30 = 20 (in positive modular arithmetic)
        let a = 10;
        let b = 20;
        let m = 30;
        let result = strictest_mod_sub(&a, &b, &m);
        assert_eq!(result, 20);
    }

    #[test]
    fn test_strictest_mod_mul() {
        // 10 * 20 = 200, 200 % 30 = 20
        let a = 10;
        let b = 20;
        let m = 30;
        let result = strictest_mod_mul(&a, &b, &m);
        assert_eq!(result, 20);
    }

    #[test]
    fn test_strictest_mod_add_edge_cases() {
        // Test with zero
        assert_eq!(strictest_mod_add(&0u32, &5u32, &7u32), 5);
        assert_eq!(strictest_mod_add(&5u32, &0u32, &7u32), 5);
        
        // Test with values larger than modulus
        assert_eq!(strictest_mod_add(&15u32, &20u32, &7u32), 0); // (1 + 6) % 7 = 0
        
        // Test overflow prevention
        assert_eq!(strictest_mod_add(&6u32, &6u32, &7u32), 5); // (6 + 6) % 7 = 5
    }

    #[test]
    fn test_strictest_mod_sub_edge_cases() {
        // Test normal subtraction
        assert_eq!(strictest_mod_sub(&10u32, &3u32, &13u32), 7);
        
        // Test underflow case
        assert_eq!(strictest_mod_sub(&3u32, &10u32, &13u32), 6); // (3 - 10) % 13 = -7 % 13 = 6
        
        // Test with zero
        assert_eq!(strictest_mod_sub(&5u32, &0u32, &7u32), 5);
        assert_eq!(strictest_mod_sub(&0u32, &5u32, &7u32), 2); // (0 - 5) % 7 = -5 % 7 = 2
        
        // Test with equal values
        assert_eq!(strictest_mod_sub(&5u32, &5u32, &7u32), 0);
    }

    #[test]
    fn test_strictest_mod_mul_edge_cases() {
        // Test with zero
        assert_eq!(strictest_mod_mul(&0u32, &5u32, &7u32), 0);
        assert_eq!(strictest_mod_mul(&5u32, &0u32, &7u32), 0);
        
        // Test with one
        assert_eq!(strictest_mod_mul(&1u32, &5u32, &7u32), 5);
        assert_eq!(strictest_mod_mul(&5u32, &1u32, &7u32), 5);
        
        // Test large values: 15%7=1, 17%7=3, 1*3=3
        assert_eq!(strictest_mod_mul(&15u32, &17u32, &7u32), 3);
        
        // Test powers of 2
        assert_eq!(strictest_mod_mul(&4u32, &4u32, &7u32), 2); // (4 * 4) % 7 = 16 % 7 = 2
    }
}
