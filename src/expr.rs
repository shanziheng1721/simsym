use crate::constant::Constant;
use crate::rational::Rational;
use crate::symbol::Symbol;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::rc::Rc;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ExprKind {
    Const(Constant),
    Var(Symbol),
    Add(Expr, Expr),
    Sub(Expr, Expr),
    Mul(Expr, Expr),
    Div(Expr, Expr),
    Neg(Expr),
    Pow(Expr, Expr),
    Sin(Expr),
    Cos(Expr),
    Tan(Expr),
    Cot(Expr),
    Sec(Expr),
    Csc(Expr),
    Asin(Expr),
    Acos(Expr),
    Atan(Expr),
    Acot(Expr),
    Asec(Expr),
    Acsc(Expr),
    Sinh(Expr),
    Cosh(Expr),
    Tanh(Expr),
    Coth(Expr),
    Sech(Expr),
    Csch(Expr),
    Asinh(Expr),
    Acosh(Expr),
    Atanh(Expr),
    Acoth(Expr),
    Asech(Expr),
    Acsch(Expr),
    Exp(Expr),
    Ln(Expr),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Expr(Rc<ExprKind>);

impl Expr {
    pub fn from_kind(kind: ExprKind) -> Self {
        Expr(Rc::new(kind))
    }

    pub fn kind(&self) -> &ExprKind {
        &self.0
    }

    /// Unwrap the AST node, moving children when this `Expr` is the sole `Rc` owner.
    pub(crate) fn into_kind(self) -> ExprKind {
        match Rc::try_unwrap(self.0) {
            Ok(k) => k,
            Err(rc) => (*rc).clone(),
        }
    }

    pub fn const_(c: Rational) -> Self {
        Self::from_kind(ExprKind::Const(Constant::from_rational(c)))
    }

    pub fn var(s: Symbol) -> Self {
        Self::from_kind(ExprKind::Var(s))
    }

    pub fn pow(self, exp: impl Into<Expr>) -> Self {
        crate::expr::pow(self, exp.into())
    }

    /// Algebraic simplification (requires the `simplify` feature; otherwise returns `self`).
    pub fn simplify(self) -> Self {
        crate::simplify::simplify(self)
    }

    /// Symbolic derivative (requires the `diff` feature).
    #[cfg(feature = "diff")]
    pub fn diff(self, var: Symbol) -> Self {
        crate::calculus::diff::diff(self, var)
    }

    /// Derivative without final simplification (requires the `diff` feature).
    #[cfg(feature = "diff")]
    pub fn diff_without_simplify(self, var: Symbol) -> Self {
        crate::calculus::diff::diff_without_simplify(self, var)
    }

    /// Symbolic integration (requires the `integrate` feature).
    #[cfg(feature = "integrate")]
    pub fn integrate(self, var: Symbol) -> Result<Self, crate::calculus::integrate::IntegrateError> {
        crate::calculus::integrate::integrate(self, var)
    }

    /// Partial derivatives (requires the `diff` feature).
    #[cfg(feature = "diff")]
    pub fn gradient(self, vars: &[Symbol]) -> Vec<Self> {
        crate::calculus::multivar::gradient(self, vars)
    }

    /// Second partial derivatives (requires the `diff` feature).
    #[cfg(feature = "diff")]
    pub fn hessian(self, vars: &[Symbol]) -> Vec<Vec<Self>> {
        crate::calculus::multivar::hessian(self, vars)
    }

    pub fn eval(
        self,
        env: &[(Symbol, Rational)],
    ) -> Result<Rational, crate::eval::EvalError> {
        crate::eval::eval(&self, env)
    }

    pub fn eval_f64(
        self,
        env: &[(Symbol, f64)],
    ) -> Result<f64, crate::eval::EvalError> {
        crate::eval::eval_f64(&self, env)
    }

    pub fn eval_f32(
        self,
        env: &[(Symbol, f32)],
    ) -> Result<f32, crate::eval::EvalError> {
        crate::eval::eval_f32(&self, env)
    }

    /// Definite integral via antiderivative or numeric fallback (requires `integrate`).
    #[cfg(feature = "integrate")]
    pub fn integrate_definite(
        self,
        var: Symbol,
        a: Rational,
        b: Rational,
        env: &[(Symbol, Rational)],
    ) -> Result<Rational, crate::calculus::numeric::DefiniteIntegralError> {
        crate::calculus::numeric::integrate_definite(self, var, a, b, env)
    }
}

// --- shared constructors (used by ops and macros) ---

pub fn const_(c: Rational) -> Expr {
    crate::constant::const_(c)
}

pub fn var(s: Symbol) -> Expr {
    Expr::var(s)
}

pub fn add(a: Expr, b: Expr) -> Expr {
    Expr::from_kind(ExprKind::Add(a, b))
}

pub fn sub(a: Expr, b: Expr) -> Expr {
    Expr::from_kind(ExprKind::Sub(a, b))
}

pub fn mul(a: Expr, b: Expr) -> Expr {
    Expr::from_kind(ExprKind::Mul(a, b))
}

pub fn div(a: Expr, b: Expr) -> Expr {
    Expr::from_kind(ExprKind::Div(a, b))
}

pub fn neg(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Neg(e))
}

pub fn pow(base: Expr, exp: Expr) -> Expr {
    Expr::from_kind(ExprKind::Pow(base, exp))
}

pub fn sin(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Sin(e))
}

pub fn cos(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Cos(e))
}

pub fn tan(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Tan(e))
}

pub fn atan(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Atan(e))
}

pub fn cot(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Cot(e))
}

pub fn sec(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Sec(e))
}

pub fn csc(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Csc(e))
}

pub fn asin(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Asin(e))
}

pub fn acos(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Acos(e))
}

pub fn acot(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Acot(e))
}

pub fn asec(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Asec(e))
}

pub fn acsc(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Acsc(e))
}

pub fn sinh(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Sinh(e))
}

pub fn cosh(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Cosh(e))
}

pub fn tanh(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Tanh(e))
}

pub fn coth(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Coth(e))
}

pub fn sech(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Sech(e))
}

pub fn csch(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Csch(e))
}

pub fn asinh(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Asinh(e))
}

pub fn acosh(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Acosh(e))
}

pub fn atanh(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Atanh(e))
}

pub fn acoth(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Acoth(e))
}

pub fn asech(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Asech(e))
}

pub fn acsch(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Acsch(e))
}

pub fn exp(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Exp(e))
}

pub fn ln(e: Expr) -> Expr {
    Expr::from_kind(ExprKind::Ln(e))
}

impl Add for Expr {
    type Output = Expr;
    fn add(self, rhs: Expr) -> Expr {
        add(self, rhs)
    }
}

impl Sub for Expr {
    type Output = Expr;
    fn sub(self, rhs: Expr) -> Expr {
        sub(self, rhs)
    }
}

impl Mul for Expr {
    type Output = Expr;
    fn mul(self, rhs: Expr) -> Expr {
        mul(self, rhs)
    }
}

impl Div for Expr {
    type Output = Expr;
    fn div(self, rhs: Expr) -> Expr {
        div(self, rhs)
    }
}

impl Neg for Expr {
    type Output = Expr;
    fn neg(self) -> Expr {
        neg(self)
    }
}

impl From<Symbol> for Expr {
    fn from(s: Symbol) -> Self {
        Expr::var(s)
    }
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        crate::display::fmt_expr(self, f, 0)
    }
}
