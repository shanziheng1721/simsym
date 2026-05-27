//! Arbitrary-precision rationals behind the `bigint` feature.

pub use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, ToPrimitive, Zero};
use std::fmt;

use crate::rational::Rational;

/// Exact rational with [`BigInt`] numerator and denominator.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BigRational(Ratio<BigInt>);

impl BigRational {
    pub fn new(num: impl Into<BigInt>, den: impl Into<BigInt>) -> Self {
        Self(Ratio::new(num.into(), den.into()))
    }

    pub fn from_integer(n: impl Into<BigInt>) -> Self {
        Self(Ratio::from_integer(n.into()))
    }

    pub fn inner(&self) -> &Ratio<BigInt> {
        &self.0
    }

    pub fn numer(&self) -> &BigInt {
        self.0.numer()
    }

    pub fn denom(&self) -> &BigInt {
        self.0.denom()
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.0.is_one()
    }

    pub fn is_integer(&self) -> bool {
        self.0.is_integer()
    }

    pub fn to_f64(&self) -> f64 {
        let n = self.0.numer().to_f64().expect("finite numerator");
        let d = self.0.denom().to_f64().expect("finite denominator");
        n / d
    }

    pub fn to_f32(&self) -> f32 {
        self.to_f64() as f32
    }
}

pub fn big_rational(num: impl Into<BigInt>, den: impl Into<BigInt>) -> BigRational {
    BigRational::new(num, den)
}

impl fmt::Display for BigRational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_integer() {
            write!(f, "{}", self.0.numer())
        } else {
            write!(f, "{}/{}", self.0.numer(), self.0.denom())
        }
    }
}

impl From<Rational> for BigRational {
    fn from(r: Rational) -> Self {
        Self::new(
            BigInt::from(r.numer()),
            BigInt::from(r.denom()),
        )
    }
}

impl BigRational {
    pub(crate) fn try_float_integer(v: f64) -> Option<Self> {
        if v.fract() != 0.0 {
            return None;
        }
        let n = BigInt::from(v as i64);
        if v >= i64::MIN as f64 && v <= i64::MAX as f64 {
            return Some(Self::from_integer(n));
        }
        None
    }
}

macro_rules! impl_from_int_for_big_rational {
    ($($ty:ty),* $(,)?) => {
        $(
            impl From<$ty> for BigRational {
                fn from(n: $ty) -> Self {
                    Self::from_integer(BigInt::from(n))
                }
            }
        )*
    };
}

impl_from_int_for_big_rational!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize,
);
