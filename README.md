### Modular math functions

[![Build](https://github.com/kaidokert/modmath-rs/actions/workflows/rust.yaml/badge.svg)](https://github.com/kaidokert/modmath-rs/actions/workflows/rust.yaml)
[![crate](https://img.shields.io/crates/v/modmath.svg)](https://crates.io/crates/modmath)
[![documentation](https://docs.rs/modmath/badge.svg)](https://docs.rs/modmath/)
[![Coverage Status](https://coveralls.io/repos/github/kaidokert/modmath-rs/badge.svg?branch=main)](https://coveralls.io/github/kaidokert/modmath-rs?branch=main)


Yet another mod math implementation, but written for _traits_. Everything
is constrained by `core::ops::` and
[`const-num-traits`](https://crates.io/crates/const-num-traits) traits.

Implements:
- Unsigned modular addition and subtraction
- Unsigned modular multiplication
- Unsigned modular exponentiation
- Unsigned modular inverse
- Unsigned modular Montgomery multiply
- Unsigned modular Montgomery exponentiation
- Constant-time variants of the Montgomery and inverse paths, behind a
  `Ct` typestate

The code isn't intended to be fast or efficient, just as generic as possible
to work with multiple implementations.

Note: While `const` traits are not yet stable and commonplace, this cannot
be very efficient. In almost all real world code you'll want to directly use
crates that implement big integers with `const` functions.
