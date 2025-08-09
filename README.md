### Modular math functions

[![Build](https://github.com/kaidokert/modmath-rs/actions/workflows/rust.yaml/badge.svg)](https://github.com/kaidokert/modmath-rs/actions/workflows/rust.yaml)
[![crate](https://img.shields.io/crates/v/modmath.svg)](https://crates.io/crates/modmath)
[![documentation](https://docs.rs/modmath/badge.svg)](https://docs.rs/modmath/)
[![Coverage Status](https://coveralls.io/repos/github/kaidokert/modmath-rs/badge.svg?branch=main)](https://coveralls.io/github/kaidokert/modmath-rs?branch=main)


Yet another mod math implementation, but written for _traits_. All functions
are free functions that are constrainted by `core::ops::` and `num_traits::`
traits.

Implements:
- Unsigned modular addition and subtraction
- Unsigned modular multiplication
- Unsigned modular exponentiation
- Unsigned modular inverse
- Unsigned modular Montgomery multiply
- Unsigned modular Montgomery exponentiation

The code isn't intended to be fast or efficient, just as generic as possible
to work with multiple implementations.

Note: While `const` traits are not yet stable and commonplace, this cannot
be verify efficient. In almost all real world code you'll want to direcyl use
crates that implement big integers with `const` functions.

### Tested with

Tested with [`num-bigint`](https://crates.io/crates/num-bigint), [`crypto-bigint`](https://crates.io/crates/crypto-bigint), [`bnum`](https://crates.io/crates/bnum), [`ibig`](https://crates.io/crates/ibig)
and [`fixed-bigint`](https://crates.io/crates/fixed-bigint) crates.
