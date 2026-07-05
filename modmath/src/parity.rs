//! Parity (odd/even) trait — re-exported from `const-num-traits`.
//!
//! The trait lives in the const-num-traits fork so it is usable from
//! const contexts. This module keeps the historical `crate::parity::Parity`
//! path stable for the rest of modmath.
pub use const_num_traits::Parity;
