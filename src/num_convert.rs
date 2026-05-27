//! Conversions between primitives and [`Rational`](crate::rational::Rational).
//!
//! - Infallible [`From`]: `i8`–`i32`, `i64`, `u8`–`u32`
//! - [`TryFrom`]: `i128`, `isize`, `u64`, `u128`, `usize`, `f32`, `f64` (and the same for [`Expr`](crate::expr::Expr))

use crate::expr::{const_, Expr};
use crate::rational::Rational;

#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum RationalConvertError {
    #[error("value does not fit in an i64 rational coefficient")]
    Overflow,
    #[error("non-finite float")]
    NonFinite,
}

/// Denominator used when converting inexact `f32`/`f64` to [`Rational`].
pub const FLOAT_RATIONAL_DENOM: i64 = 1_000_000_000;

fn rational_from_i64(n: i64) -> Rational {
    Rational::from_integer(n)
}

fn try_i64<T>(n: T) -> Result<i64, RationalConvertError>
where
    i64: TryFrom<T>,
{
    i64::try_from(n).map_err(|_| RationalConvertError::Overflow)
}

macro_rules! impl_try_from_int_for_rational {
    ($($ty:ty),* $(,)?) => {
        $(
            impl TryFrom<$ty> for Rational {
                type Error = RationalConvertError;
                fn try_from(n: $ty) -> Result<Self, Self::Error> {
                    Ok(rational_from_i64(try_i64(n)?))
                }
            }
        )*
    };
}

/// Infallible `From` only for types whose entire range fits in `i64`.
macro_rules! impl_from_int_for_rational_and_expr {
    ($($ty:ty),* $(,)?) => {
        $(
            impl From<$ty> for Rational {
                fn from(n: $ty) -> Self {
                    rational_from_i64(i64::from(n))
                }
            }
            impl From<$ty> for Expr {
                fn from(n: $ty) -> Self {
                    const_(Rational::from(n))
                }
            }
        )*
    };
}

#[cfg(not(feature = "bigint"))]
macro_rules! impl_try_from_int_for_expr {
    ($($ty:ty),* $(,)?) => {
        $(
            impl TryFrom<$ty> for Expr {
                type Error = RationalConvertError;
                fn try_from(n: $ty) -> Result<Self, Self::Error> {
                    Ok(const_(Rational::try_from(n)?))
                }
            }
        )*
    };
}

impl_try_from_int_for_rational!(i128, isize, u64, u128, usize);

impl_from_int_for_rational_and_expr!(i8, i16, i32, i64, u8, u16, u32);

#[cfg(not(feature = "bigint"))]
impl_try_from_int_for_expr!(i128, isize, u64, u128, usize);

fn try_from_float_bits(v: f64) -> Result<Rational, RationalConvertError> {
    if !v.is_finite() {
        return Err(RationalConvertError::NonFinite);
    }
    if v.fract() == 0.0 && v >= i64::MIN as f64 && v <= i64::MAX as f64 {
        return Ok(rational_from_i64(v as i64));
    }
    let n = (v * FLOAT_RATIONAL_DENOM as f64).round() as i64;
    Ok(Rational::new(n, FLOAT_RATIONAL_DENOM))
}

impl TryFrom<f64> for Rational {
    type Error = RationalConvertError;
    fn try_from(v: f64) -> Result<Self, Self::Error> {
        try_from_float_bits(v)
    }
}

impl TryFrom<f32> for Rational {
    type Error = RationalConvertError;
    fn try_from(v: f32) -> Result<Self, Self::Error> {
        try_from_float_bits(f64::from(v))
    }
}

#[cfg(not(feature = "bigint"))]
impl TryFrom<f64> for Expr {
    type Error = RationalConvertError;
    fn try_from(v: f64) -> Result<Self, Self::Error> {
        Ok(const_(Rational::try_from(v)?))
    }
}

#[cfg(not(feature = "bigint"))]
impl TryFrom<f32> for Expr {
    type Error = RationalConvertError;
    fn try_from(v: f32) -> Result<Self, Self::Error> {
        Ok(const_(Rational::try_from(v)?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn i128_overflow() {
        assert!(Rational::try_from(i128::MAX).is_err());
    }

    #[test]
    fn f32_constant() {
        let r = Rational::try_from(0.5f32).unwrap();
        assert!((r.to_f64() - 0.5).abs() < 1e-9);
    }
}
