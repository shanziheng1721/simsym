//! Derivatives of trigonometric and hyperbolic unary functions.

use crate::expr::{
    add, const_, cosh, cot, coth, csc, csch, mul, pow, sec, sech, sinh, sub, tan, tanh, Expr,
};
use crate::expr::ExprKind;
use crate::rational::Rational;
use crate::symbol::Symbol;

use super::diff_expr;

fn half() -> Expr {
    const_(Rational::new(1, 2))
}

fn sqrt_expr(e: Expr) -> Expr {
    pow(e, half())
}

fn one_plus_sq(e: &Expr) -> Expr {
    add(const_(Rational::one()), pow(e.clone(), const_(Rational::from(2))))
}

fn one_minus_sq(e: &Expr) -> Expr {
    sub(const_(Rational::one()), pow(e.clone(), const_(Rational::from(2))))
}

fn sq_minus_one(e: &Expr) -> Expr {
    sub(pow(e.clone(), const_(Rational::from(2))), const_(Rational::one()))
}

pub(crate) fn diff_trig_kind(kind: &ExprKind, var: Symbol) -> Option<Expr> {
    let d = |e: &Expr| diff_expr(e, var);
    Some(match kind {
        ExprKind::Cot(e) => -d(e) * pow(csc(e.clone()), const_(Rational::from(2))),
        ExprKind::Sec(e) => d(e) * mul(sec(e.clone()), tan(e.clone())),
        ExprKind::Csc(e) => -d(e) * mul(csc(e.clone()), cot(e.clone())),
        ExprKind::Asin(e) => d(e) / sqrt_expr(one_minus_sq(e)),
        ExprKind::Acos(e) => -d(e) / sqrt_expr(one_minus_sq(e)),
        ExprKind::Acot(e) => -d(e) / one_plus_sq(e),
        ExprKind::Asec(e) => d(e) / mul(e.clone(), sqrt_expr(sq_minus_one(e))),
        ExprKind::Acsc(e) => -d(e) / mul(e.clone(), sqrt_expr(sq_minus_one(e))),
        ExprKind::Sinh(e) => d(e) * cosh(e.clone()),
        ExprKind::Cosh(e) => d(e) * sinh(e.clone()),
        ExprKind::Tanh(e) => d(e) * pow(sech(e.clone()), const_(Rational::from(2))),
        ExprKind::Coth(e) => -d(e) * pow(csch(e.clone()), const_(Rational::from(2))),
        ExprKind::Sech(e) => -d(e) * mul(sech(e.clone()), tanh(e.clone())),
        ExprKind::Csch(e) => -d(e) * mul(csch(e.clone()), coth(e.clone())),
        ExprKind::Asinh(e) => d(e) / sqrt_expr(one_plus_sq(e)),
        ExprKind::Acosh(e) => d(e) / sqrt_expr(sq_minus_one(e)),
        ExprKind::Atanh(e) => d(e) / one_minus_sq(e),
        ExprKind::Acoth(e) => d(e) / one_minus_sq(e),
        ExprKind::Asech(e) => -d(e) / mul(e.clone(), sqrt_expr(one_minus_sq(e))),
        ExprKind::Acsch(e) => -d(e) / mul(e.clone(), sqrt_expr(one_plus_sq(e))),
        _ => return None,
    })
}
