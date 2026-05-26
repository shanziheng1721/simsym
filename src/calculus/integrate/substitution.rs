use crate::expr::{div, mul, pow, Expr, ExprKind};
use crate::rational::Rational;
use crate::symbol::Symbol;

use super::util::{as_const, chain_rule_scale, linear_factor};

/// Try ∫ f(u) · u' dx when `f` is a power or elementary function of `u`.
pub fn try_u_substitution(expr: &Expr, var: Symbol) -> Option<Expr> {
    match expr.kind() {
        ExprKind::Mul(f, g) => try_substitution_product(f, g, var)
            .or_else(|| try_substitution_product(g, f, var)),
        _ => None,
    }
}

fn try_substitution_product(f: &Expr, gp: &Expr, var: Symbol) -> Option<Expr> {
    let gp = gp.clone().simplify();

    if let ExprKind::Pow(base, exp) = f.kind() {
        if let Some(n) = as_const(exp).and_then(|r| r.as_integer()) {
            if n != -1 {
                let u = base.clone();
                let du = u.clone().diff(var).simplify();
                if let Some(scale) = linear_factor(&gp, &du) {
                    let new_n = n + 1;
                    return Some(div(
                        pow(u, crate::expr::const_(Rational::from(new_n))),
                        crate::expr::const_(Rational::from(new_n) * scale),
                    ));
                }
            }
        }
        if let Some(n) = as_const(exp).and_then(|r| r.as_integer()) {
            if n == -2 {
                if let ExprKind::Cos(inner) = base.kind() {
                    if let Some(scale) = chain_rule_scale(&gp, inner, var) {
                        return Some(mul(
                            crate::expr::const_(scale),
                            crate::expr::tan(inner.clone()),
                        ));
                    }
                }
            }
        }
    }

    match f.kind() {
        ExprKind::Sin(inner) => chain_elementary(
            inner,
            &gp,
            var,
            |u| -crate::expr::cos(u.clone()),
        ),
        ExprKind::Cos(inner) => chain_elementary(inner, &gp, var, |u| crate::expr::sin(u.clone())),
        ExprKind::Tan(inner) => chain_elementary(
            inner,
            &gp,
            var,
            |u| -crate::expr::ln(crate::expr::cos(u.clone())),
        ),
        ExprKind::Exp(inner) => chain_elementary(inner, &gp, var, |u| crate::expr::exp(u.clone())),
        ExprKind::Ln(inner) => chain_elementary(
            inner,
            &gp,
            var,
            |u| mul(u.clone(), crate::expr::ln(u.clone())) - u.clone(),
        ),
        ExprKind::Atan(inner) => chain_elementary(
            inner,
            &gp,
            var,
            |u| {
                mul(u.clone(), crate::expr::atan(u.clone()))
                    - div(
                        crate::expr::ln(
                            crate::expr::const_(Rational::one())
                                + pow(u.clone(), crate::expr::const_(Rational::from(2))),
                        ),
                        crate::expr::const_(Rational::from(2)),
                    )
            },
        ),
        _ => None,
    }
}

fn chain_elementary(
    inner: &Expr,
    gp: &Expr,
    var: Symbol,
    antideriv_in_u: fn(&Expr) -> Expr,
) -> Option<Expr> {
    let scale = chain_rule_scale(gp, inner, var)?;
    Some(mul(crate::expr::const_(scale), antideriv_in_u(inner)))
}
