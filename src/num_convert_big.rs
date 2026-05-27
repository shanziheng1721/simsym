//! Wide integer → [`Expr`](crate::expr::Expr) and float → [`BigRational`](crate::rational_big::BigRational) (`bigint` feature).

use crate::constant::{big_const, Constant};
use crate::expr::Expr;
use crate::num_convert::RationalConvertError;
use crate::rational::Rational;
use crate::rational_big::BigRational;

macro_rules! impl_from_wide_int_for_expr {
    ($($ty:ty),* $(,)?) => {
        $(
            impl From<$ty> for Constant {
                fn from(n: $ty) -> Self {
                    Constant::from_big(BigRational::from(n))
                }
            }
            impl From<$ty> for Expr {
                fn from(n: $ty) -> Self {
                    big_const(BigRational::from(n))
                }
            }
        )*
    };
}

impl_from_wide_int_for_expr!(i128, isize, u64, u128, usize);

/// `TryFrom` for floats → [`BigRational`] (same semantics as [`Rational::try_from`]).
impl TryFrom<f64> for BigRational {
    type Error = RationalConvertError;
    fn try_from(v: f64) -> Result<Self, Self::Error> {
        if !v.is_finite() {
            return Err(RationalConvertError::NonFinite);
        }
        if v.fract() == 0.0 {
            if let Some(n) = BigRational::try_float_integer(v) {
                return Ok(n);
            }
        }
        let den = num_bigint::BigInt::from(crate::num_convert::FLOAT_RATIONAL_DENOM);
        let n = num_bigint::BigInt::from((v * crate::num_convert::FLOAT_RATIONAL_DENOM as f64).round() as i64);
        Ok(BigRational::new(n, den))
    }
}

impl TryFrom<f32> for BigRational {
    type Error = RationalConvertError;
    fn try_from(v: f32) -> Result<Self, Self::Error> {
        BigRational::try_from(f64::from(v))
    }
}

impl TryFrom<f64> for Expr {
    type Error = RationalConvertError;
    fn try_from(v: f64) -> Result<Self, Self::Error> {
        if let Ok(r) = Rational::try_from(v) {
            return Ok(crate::constant::const_(r));
        }
        Ok(big_const(BigRational::try_from(v)?))
    }
}

impl TryFrom<f32> for Expr {
    type Error = RationalConvertError;
    fn try_from(v: f32) -> Result<Self, Self::Error> {
        Expr::try_from(f64::from(v))
    }
}
