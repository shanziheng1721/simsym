use crate::rational::Rational;
use crate::symbol::Symbol;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::rc::Rc;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ExprKind {
    Const(Rational),
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
    Atan(Expr),
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

    pub fn const_(c: Rational) -> Self {
        Self::from_kind(ExprKind::Const(c))
    }

    pub fn var(s: Symbol) -> Self {
        Self::from_kind(ExprKind::Var(s))
    }

    pub fn pow(self, exp: impl Into<Expr>) -> Self {
        crate::expr::pow(self, exp.into())
    }

    pub fn simplify(self) -> Self {
        crate::simplify::simplify(self)
    }

    pub fn diff(self, var: Symbol) -> Self {
        crate::calculus::diff::diff(self, var)
    }

    /// Derivative without simplifying (see [`crate::calculus::diff_without_simplify`]).
    pub fn diff_without_simplify(self, var: Symbol) -> Self {
        crate::calculus::diff::diff_without_simplify(self, var)
    }

    pub fn integrate(self, var: Symbol) -> Result<Self, crate::calculus::integrate::IntegrateError> {
        crate::calculus::integrate::integrate(self, var)
    }

    pub fn gradient(self, vars: &[Symbol]) -> Vec<Self> {
        crate::calculus::multivar::gradient(self, vars)
    }

    pub fn hessian(self, vars: &[Symbol]) -> Vec<Vec<Self>> {
        crate::calculus::multivar::hessian(self, vars)
    }

    pub fn eval(
        self,
        env: &[(Symbol, Rational)],
    ) -> Result<Rational, crate::eval::EvalError> {
        crate::eval::eval(self, env)
    }

    pub fn eval_f64(
        self,
        env: &[(Symbol, f64)],
    ) -> Result<f64, crate::eval::EvalError> {
        crate::eval::eval_f64(self, env)
    }

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
    Expr::const_(c)
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

impl From<i32> for Expr {
    fn from(n: i32) -> Self {
        const_(Rational::from(n))
    }
}

impl From<i64> for Expr {
    fn from(n: i64) -> Self {
        const_(Rational::from(n))
    }
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        crate::display::fmt_expr(self, f, 0)
    }
}
