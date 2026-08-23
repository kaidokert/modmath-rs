//! Operation-level modular-arithmetic accelerator traits.
//!
//! These traits complement the value-oriented [`Field`](crate::Field) API.
//! They fit singleton peripherals and remote secure elements that accept whole
//! operations, where making a hardware register window impersonate an ordinary
//! copyable bigint would obscure ownership, failure, and transfer costs.

/// A backend failure with portable meaning.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum AcceleratorError {
    /// The modulus is zero or otherwise rejected by the operation.
    InvalidModulus,
    /// An input encoding is not accepted by the backend.
    InvalidOperand,
    /// `out` cannot hold the modulus-width result.
    OutputSize,
    /// The selected peripheral or implementation lacks this operation.
    Unsupported,
    /// The operation failed, timed out, or reported a hardware fault.
    Backend,
}

impl core::fmt::Display for AcceleratorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::InvalidModulus => "invalid modulus",
            Self::InvalidOperand => "invalid modular operand",
            Self::OutputSize => "output must match the modulus width",
            Self::Unsupported => "operation is not supported by this accelerator",
            Self::Backend => "modular accelerator failed",
        })
    }
}

impl core::error::Error for AcceleratorError {}

/// Compute a complete modular multiplication operation.
///
/// Inputs and output are unsigned, canonical big-endian byte strings. Leading
/// zeroes are accepted. `out.len()` must equal `modulus_be.len()` and the result
/// is left-padded to that exact width. Implementations reduce operands that are
/// at least as wide as the modulus.
pub trait ModularMultiply {
    type Error;

    fn modular_multiply(
        &mut self,
        lhs_be: &[u8],
        rhs_be: &[u8],
        modulus_be: &[u8],
        out: &mut [u8],
    ) -> Result<(), Self::Error>;
}

/// Compute a complete modular exponentiation operation.
///
/// Inputs and output use the same encoding contract as [`ModularMultiply`].
/// The exponent may be zero. This base trait makes no timing claim and is
/// suitable for public exponents and public verification inputs.
pub trait ModularExponentiate {
    type Error;

    fn modular_exponentiate(
        &mut self,
        base_be: &[u8],
        exponent_be: &[u8],
        modulus_be: &[u8],
        out: &mut [u8],
    ) -> Result<(), Self::Error>;
}

/// Marker for modular exponentiation suitable for secret exponents.
///
/// # Safety
///
/// Implementing this trait asserts that secret-dependent timing, memory
/// access, peripheral command selection, and externally observable transfer
/// behavior have been reviewed for the backend's documented threat model.
/// It does not assert resistance to power or electromagnetic analysis unless
/// the implementation explicitly documents that property.
pub unsafe trait SecretModularExponentiate: ModularExponentiate {}

/// Compute a complete multiplicative inverse operation.
///
/// `value_be` and `modulus_be` are unsigned big-endian values and the result is
/// modulus-width. A non-invertible value must return an error rather than an
/// arbitrary output.
pub trait ModularInvert {
    type Error;

    fn modular_invert(
        &mut self,
        value_be: &[u8],
        modulus_be: &[u8],
        out: &mut [u8],
    ) -> Result<(), Self::Error>;
}

/// Check the common output-width invariant before touching a peripheral.
pub fn validate_output(modulus_be: &[u8], out: &[u8]) -> Result<(), AcceleratorError> {
    if modulus_be.is_empty() || modulus_be.iter().all(|byte| *byte == 0) {
        return Err(AcceleratorError::InvalidModulus);
    }
    if out.len() != modulus_be.len() {
        return Err(AcceleratorError::OutputSize);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    struct U128Accelerator;

    fn decode(bytes: &[u8]) -> Result<u128, AcceleratorError> {
        if bytes.len() > 16 {
            return Err(AcceleratorError::InvalidOperand);
        }
        let mut encoded = [0u8; 16];
        encoded[16 - bytes.len()..].copy_from_slice(bytes);
        Ok(u128::from_be_bytes(encoded))
    }

    fn encode(value: u128, out: &mut [u8]) {
        let encoded = value.to_be_bytes();
        out.copy_from_slice(&encoded[encoded.len() - out.len()..]);
    }

    impl ModularExponentiate for U128Accelerator {
        type Error = AcceleratorError;

        fn modular_exponentiate(
            &mut self,
            base_be: &[u8],
            exponent_be: &[u8],
            modulus_be: &[u8],
            out: &mut [u8],
        ) -> Result<(), Self::Error> {
            validate_output(modulus_be, out)?;
            let modulus = decode(modulus_be)?;
            let mut base = decode(base_be)? % modulus;
            let mut exponent = decode(exponent_be)?;
            let mut result = 1 % modulus;
            while exponent != 0 {
                if exponent & 1 == 1 {
                    result = result * base % modulus;
                }
                base = base * base % modulus;
                exponent >>= 1;
            }
            encode(result, out);
            Ok(())
        }
    }

    #[test]
    fn mock_backend_obeys_big_endian_width_contract() {
        let mut out = [0xaa; 2];
        U128Accelerator
            .modular_exponentiate(&[0, 4], &[0, 13], &[0x01, 0xf1], &mut out)
            .unwrap();
        assert_eq!(out, [0x01, 0xbd]); // 4^13 mod 497 = 445
    }

    #[test]
    fn output_and_modulus_are_validated_portably() {
        assert_eq!(
            validate_output(&[0, 0], &mut [0; 2]),
            Err(AcceleratorError::InvalidModulus)
        );
        assert_eq!(
            validate_output(&[3], &mut [0; 2]),
            Err(AcceleratorError::OutputSize)
        );
    }
}
