use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

#[cfg(feature = "bigint")]
use num_traits::ToPrimitive;
use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// Exact rational constant (reduced `i64` ratio).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Rational(Ratio<i64>);

impl Rational {
    pub fn new(num: i64, den: i64) -> Self {
        Self(Ratio::new(num, den))
    }

    pub fn from_integer(n: i64) -> Self {
        Self(Ratio::from_integer(n))
    }

    pub fn to_f64(self) -> f64 {
        *self.0.numer() as f64 / *self.0.denom() as f64
    }

    pub fn to_f32(self) -> f32 {
        self.to_f64() as f32
    }

    pub fn is_zero(self) -> bool {
        self.0.is_zero()
    }

    pub fn is_one(self) -> bool {
        self.0.is_one()
    }

    pub fn is_integer(self) -> bool {
        self.0.is_integer()
    }

    pub fn is_positive(self) -> bool {
        self.0.is_positive()
    }

    pub fn as_integer(self) -> Option<i64> {
        if self.is_integer() {
            Some(*self.0.numer())
        } else {
            None
        }
    }

    pub fn inner(self) -> Ratio<i64> {
        self.0
    }

    pub fn numer(self) -> i64 {
        *self.0.numer()
    }

    pub fn denom(self) -> i64 {
        *self.0.denom()
    }

    /// Convert a [`BigRational`](crate::rational_big::BigRational) when numerator and denominator fit in `i64`.
    #[cfg(feature = "bigint")]
    pub fn try_from_big(b: &crate::rational_big::BigRational) -> Option<Self> {
        let n = b.numer().to_i64()?;
        let d = b.denom().to_i64()?;
        Some(Self::new(n, d))
    }
}

/// Build a rational from numerator and denominator.
pub fn rational(num: i64, den: i64) -> Rational {
    Rational::new(num, den)
}

/// Shorthand for [`Rational::from`] on `i32`.
pub fn rational_from_i32(n: i32) -> Rational {
    Rational::from(n)
}

pub fn rat_pow_int(base: Rational, exp: i64) -> Result<Rational, RationalPowError> {
    if exp == 0 {
        return Ok(Rational::one());
    }
    if base.is_zero() {
        if exp > 0 {
            return Ok(Rational::zero());
        }
        return Err(RationalPowError::ZeroToNegativePower);
    }
    let abs_exp: i32 = exp
        .unsigned_abs()
        .try_into()
        .map_err(|_| RationalPowError::Overflow)?;
    let inner = if exp > 0 {
        base.0.pow(abs_exp)
    } else {
        base.0.pow(abs_exp).recip()
    };
    Ok(Rational(inner))
}

impl Rational {
    pub fn zero() -> Self {
        Self(Ratio::zero())
    }

    pub fn one() -> Self {
        Self(Ratio::one())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum RationalPowError {
    #[error("zero to a negative power")]
    ZeroToNegativePower,
    #[error("exponent overflow")]
    Overflow,
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_integer() {
            write!(f, "{}", self.0.numer())
        } else {
            write!(f, "{}/{}", self.0.numer(), self.0.denom())
        }
    }
}

impl Add for Rational {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl Sub for Rational {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl Mul for Rational {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self(self.0 * rhs.0)
    }
}

impl Div for Rational {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        Self(self.0 / rhs.0)
    }
}

impl AddAssign for Rational {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl SubAssign for Rational {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
    }
}

impl Neg for Rational {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

#[cfg(feature = "serde")]
impl serde::Serialize for Rational {
    fn serialize<S: serde::Serializer>(&self, ser: S) -> Result<S::Ok, S::Error> {
        (self.0.numer(), self.0.denom()).serialize(ser)
    }
}

#[cfg(feature = "serde")]
impl<'de> serde::Deserialize<'de> for Rational {
    fn deserialize<D: serde::Deserializer<'de>>(de: D) -> Result<Self, D::Error> {
        let (n, d): (i64, i64) = serde::Deserialize::deserialize(de)?;
        Ok(Rational::new(n, d))
    }
}
