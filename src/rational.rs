use num_rational::Ratio;
use num_traits::{One, Signed, Zero};
use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

/// Exact rational constant (reduced `i64` ratio).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
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
}

/// Build a rational from numerator and denominator.
pub fn rational(num: i64, den: i64) -> Rational {
    Rational::new(num, den)
}

pub fn rational_from_i32(n: i32) -> Rational {
    Rational::from_integer(i64::from(n))
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

impl From<i32> for Rational {
    fn from(n: i32) -> Self {
        rational_from_i32(n)
    }
}

impl From<i64> for Rational {
    fn from(n: i64) -> Self {
        Rational::from_integer(n)
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
