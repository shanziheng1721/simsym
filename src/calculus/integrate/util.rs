use crate::expr::{add, mul, Expr, ExprKind};
use crate::rational::Rational;
use crate::symbol::Symbol;
use std::collections::BTreeMap;

use super::IntegrateError;

pub fn contains_var(expr: &Expr, var: Symbol) -> bool {
    match expr.kind() {
        ExprKind::Var(s) => *s == var,
        ExprKind::Const(_) => false,
        ExprKind::Add(a, b)
        | ExprKind::Sub(a, b)
        | ExprKind::Mul(a, b)
        | ExprKind::Div(a, b)
        | ExprKind::Pow(a, b) => contains_var(a, var) || contains_var(b, var),
        ExprKind::Neg(e)
        | ExprKind::Sin(e)
        | ExprKind::Cos(e)
        | ExprKind::Tan(e)
        | ExprKind::Exp(e)
        | ExprKind::Ln(e)
        | ExprKind::Atan(e) => contains_var(e, var),
    }
}

pub fn is_var(e: &Expr, var: Symbol) -> bool {
    matches!(e.kind(), ExprKind::Var(s) if *s == var)
}

pub fn is_euler_var(e: &Expr) -> bool {
    matches!(e.kind(), ExprKind::Var(s) if s.name() == "e")
}

pub fn is_exp_of_var(e: &Expr, var: Symbol) -> bool {
    exp_affine(e, var).is_some()
}

/// `e^{k*x + b}` as `(k, b)`.
pub fn exp_affine(e: &Expr, var: Symbol) -> Option<(Rational, Rational)> {
    match e.kind() {
        ExprKind::Exp(inner) => affine_form(inner, var),
        ExprKind::Pow(base, exp) if is_euler_var(base) => affine_form(exp, var),
        _ => None,
    }
}

/// `sin|cos(k*x + b)` as `(k, b)`.
pub fn trig_affine(e: &Expr, var: Symbol, want_sin: bool) -> Option<(Rational, Rational)> {
    let inner = match e.kind() {
        ExprKind::Sin(i) if want_sin => i,
        ExprKind::Cos(i) if !want_sin => i,
        _ => return None,
    };
    affine_form(inner, var)
}

/// If `a == scale * b` with rational `scale`, return `scale`.
pub fn linear_factor(a: &Expr, b: &Expr) -> Option<Rational> {
    if a == b {
        return Some(Rational::one());
    }
    match a.kind() {
        ExprKind::Mul(l, r) => {
            if let ExprKind::Const(c) = l.kind() {
                if r.clone().simplify() == *b {
                    return Some(*c);
                }
            }
            if let ExprKind::Const(c) = r.kind() {
                if l.clone().simplify() == *b {
                    return Some(*c);
                }
            }
        }
        _ => {}
    }
    None
}

/// `gp` is a constant multiple of `d(inner)/dx`.
pub fn chain_rule_scale(gp: &Expr, inner: &Expr, var: Symbol) -> Option<Rational> {
    let du = inner.clone().diff(var).simplify();
    linear_factor(gp, &du)
}

/// Inner expression `k * var + b` (or just `var` when k=1, b=0).
pub fn affine_form(e: &Expr, var: Symbol) -> Option<(Rational, Rational)> {
    match e.kind() {
        ExprKind::Const(c) => Some((Rational::zero(), *c)),
        ExprKind::Var(s) if *s == var => Some((Rational::one(), Rational::zero())),
        ExprKind::Mul(l, r) => {
            if let ExprKind::Const(k) = l.kind() {
                if is_var(r, var) {
                    return Some((*k, Rational::zero()));
                }
                if let Some((k2, b)) = affine_form(r, var) {
                    return Some((*k * k2, *k * b));
                }
            }
            if let ExprKind::Const(k) = r.kind() {
                if is_var(l, var) {
                    return Some((*k, Rational::zero()));
                }
                if let Some((k2, b)) = affine_form(l, var) {
                    return Some((*k * k2, *k * b));
                }
            }
            None
        }
        ExprKind::Add(l, r) => {
            let (k1, b1) = affine_form(l, var)?;
            let (k2, b2) = affine_form(r, var)?;
            if !k2.is_zero() {
                return None;
            }
            Some((k1, b1 + b2))
        }
        ExprKind::Sub(l, r) => {
            if let Some(c) = as_const(r) {
                let (k, b) = affine_form(l, var)?;
                return Some((k, b - c));
            }
            let (k1, b1) = affine_form(l, var)?;
            let (k2, b2) = affine_form(r, var)?;
            if !k2.is_zero() {
                return None;
            }
            Some((k1, b1 - b2))
        }
        ExprKind::Neg(inner) => {
            let (k, b) = affine_form(inner, var)?;
            Some((-k, -b))
        }
        _ => None,
    }
}

/// `x - c` or `x + c` as `(offset)` meaning variable `x + offset` (so `x - 2` → offset -2).
pub fn var_plus_const(e: &Expr, var: Symbol) -> Option<Rational> {
    match e.kind() {
        ExprKind::Var(s) if *s == var => Some(Rational::zero()),
        ExprKind::Add(l, r) => {
            if is_var(l, var) {
                if let ExprKind::Const(c) = r.kind() {
                    return Some(*c);
                }
            }
            if is_var(r, var) {
                if let ExprKind::Const(c) = l.kind() {
                    return Some(*c);
                }
            }
            None
        }
        ExprKind::Sub(l, r) => {
            if is_var(l, var) {
                if let ExprKind::Const(c) = r.kind() {
                    return Some(-*c);
                }
            }
            None
        }
        ExprKind::Neg(inner) => {
            if let ExprKind::Var(s) = inner.kind() {
                if *s == var {
                    return Some(Rational::zero());
                }
            }
            None
        }
        _ => None,
    }
}

pub fn as_const(e: &Expr) -> Option<Rational> {
    match e.kind() {
        ExprKind::Const(c) => Some(*c),
        _ => None,
    }
}

pub fn collect_polynomial_terms(
    expr: &Expr,
    var: Symbol,
) -> Result<BTreeMap<i64, Rational>, IntegrateError> {
    let mut map = BTreeMap::new();
    collect_poly_rec(expr, var, Rational::one(), &mut map)?;
    Ok(map)
}

fn collect_poly_rec(
    expr: &Expr,
    var: Symbol,
    coeff: Rational,
    out: &mut BTreeMap<i64, Rational>,
) -> Result<(), IntegrateError> {
    match expr.kind() {
        ExprKind::Const(c) => {
            if !coeff.is_zero() {
                *out.entry(0).or_insert(Rational::zero()) += coeff * *c;
            }
            Ok(())
        }
        ExprKind::Var(s) if *s == var => {
            *out.entry(1).or_insert(Rational::zero()) += coeff;
            Ok(())
        }
        ExprKind::Var(_) => Ok(()),
        ExprKind::Add(a, b) => {
            collect_poly_rec(a, var, coeff, out)?;
            collect_poly_rec(b, var, coeff, out)
        }
        ExprKind::Sub(a, b) => {
            collect_poly_rec(a, var, coeff, out)?;
            collect_poly_rec(b, var, -coeff, out)
        }
        ExprKind::Neg(e) => collect_poly_rec(e, var, -coeff, out),
        ExprKind::Mul(l, r) => {
            if let ExprKind::Const(c) = l.kind() {
                collect_poly_rec(r, var, coeff * *c, out)
            } else if let ExprKind::Const(c) = r.kind() {
                collect_poly_rec(l, var, coeff * *c, out)
            } else {
                Err(IntegrateError::NoRule)
            }
        }
        ExprKind::Pow(base, exp) => {
            if let Some(offset) = var_plus_const(base, var) {
                if offset.is_zero() {
                    if let ExprKind::Const(n) = exp.kind() {
                        if let Some(k) = n.as_integer() {
                            if k >= 0 {
                                *out.entry(k).or_insert(Rational::zero()) += coeff;
                                return Ok(());
                            }
                        }
                    }
                }
            }
            if matches!(base.kind(), ExprKind::Var(s) if *s == var) {
                if let ExprKind::Const(n) = exp.kind() {
                    if let Some(k) = n.as_integer() {
                        if k >= 0 {
                            *out.entry(k).or_insert(Rational::zero()) += coeff;
                            return Ok(());
                        }
                    }
                }
            }
            Err(IntegrateError::NoRule)
        }
        _ => Err(IntegrateError::NoRule),
    }
}

/// Build `k*var + b` as an expression.
pub fn affine_angle(k: Rational, b: Rational, var: Symbol) -> Expr {
    let x = Expr::var(var);
    if b.is_zero() {
        mul(crate::expr::const_(k), x)
    } else if k.is_one() {
        add(x, crate::expr::const_(b))
    } else {
        add(mul(crate::expr::const_(k), x), crate::expr::const_(b))
    }
}

/// Scale chain rule factor 1/k for ∫ f(kx+b) dx.
pub fn affine_outer_integral(k: Rational, inner_antideriv: Expr) -> Expr {
    if k.is_one() {
        inner_antideriv
    } else {
        inner_antideriv / crate::expr::const_(k)
    }
}
