//! Expression constants: [`Rational`] or, with the `bigint` feature, arbitrary-precision [`BigRational`].

use crate::expr::Expr;
use crate::rational::Rational;
use std::fmt;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Constant {
    /// Coefficient that fits in `i64` numerator/denominator.
    Small(Rational),
    /// Wide exact rational (`bigint` feature only).
    #[cfg(feature = "bigint")]
    Wide(crate::rational_big::BigRational),
}

impl Constant {
    pub fn from_rational(r: Rational) -> Self {
        Self::Small(r)
    }

    #[cfg(feature = "bigint")]
    pub fn from_big(b: crate::rational_big::BigRational) -> Self {
        if let Some(r) = Rational::try_from_big(&b) {
            Self::Small(r)
        } else {
            Self::Wide(b)
        }
    }

    /// If this constant fits in [`Rational`], return it.
    pub fn try_as_rational(&self) -> Option<Rational> {
        match self {
            Self::Small(r) => Some(*r),
            #[cfg(feature = "bigint")]
            Self::Wide(b) => Rational::try_from_big(b),
        }
    }

    pub fn is_zero(&self) -> bool {
        match self {
            Self::Small(r) => r.is_zero(),
            #[cfg(feature = "bigint")]
            Self::Wide(b) => b.is_zero(),
        }
    }

    pub fn is_one(&self) -> bool {
        match self {
            Self::Small(r) => r.is_one(),
            #[cfg(feature = "bigint")]
            Self::Wide(b) => b.is_one(),
        }
    }

    pub fn to_f64(&self) -> f64 {
        match self {
            Self::Small(r) => r.to_f64(),
            #[cfg(feature = "bigint")]
            Self::Wide(b) => b.to_f64(),
        }
    }
}

impl fmt::Display for Constant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Small(r) => write!(f, "{r}"),
            #[cfg(feature = "bigint")]
            Self::Wide(b) => write!(f, "{b}"),
        }
    }
}

impl From<Rational> for Constant {
    fn from(r: Rational) -> Self {
        Self::from_rational(r)
    }
}

#[cfg(feature = "bigint")]
impl From<crate::rational_big::BigRational> for Constant {
    fn from(b: crate::rational_big::BigRational) -> Self {
        Self::from_big(b)
    }
}

pub fn constant(c: Constant) -> Expr {
    Expr::from_kind(crate::expr::ExprKind::Const(c))
}

pub fn const_(r: Rational) -> Expr {
    constant(Constant::from_rational(r))
}

#[cfg(feature = "bigint")]
pub fn big_const(b: crate::rational_big::BigRational) -> Expr {
    constant(Constant::from_big(b))
}
