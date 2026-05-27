use crate::expr::{add, const_, div, mul, sub, Expr, ExprKind};
use crate::rational::Rational;
use crate::symbol::Symbol;

use super::util::{as_const, contains_var, is_var, var_plus_const};
use super::{integrate, integrate_expr, IntegrateError};

pub fn integrate_div(f: Expr, g: Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if !contains_var(&g, var) {
        return Ok(integrate_expr(f, var)? / g);
    }

    // 1/(x-a)/(x-b) encoded as Div(Div(1, x-a), x-b)
    if let ExprKind::Div(num, den) = f.kind() {
        if as_const(num).is_some_and(|c| c.is_one()) {
            let combined = mul(den.clone(), g.clone());
            if let Some(integ) = integrate_two_linear_factors(&const_(Rational::one()), &combined, var)
            {
                return Ok(integ);
            }
        }
    }

    if let Some(integ) = integrate_reciprocal_linear(&f, &g, var) {
        return Ok(integ);
    }

    if let Some(integ) = integrate_two_linear_factors(&f, &g, var) {
        return Ok(integ);
    }

    if super::util::collect_polynomial_terms(&f, var).is_ok()
        && super::util::collect_polynomial_terms(&g, var).is_ok()
    {
        if let Ok(integ) =
            super::rational_function::integrate_rational_function(f.clone(), g.clone(), var)
        {
            return Ok(integ);
        }
    }

    if let Some(integ) = integrate_polynomial_over_linear(f, g, var) {
        return Ok(integ);
    }

    Err(IntegrateError::NoRule)
}

/// ∫ c / (x - a) dx or ∫ c / x dx
fn integrate_reciprocal_linear(f: &Expr, g: &Expr, var: Symbol) -> Option<Expr> {
    let c = as_const(f)?;
    if let Some(offset) = var_plus_const(g, var) {
        let arg = if offset.is_zero() {
            Expr::var(var)
        } else {
            add(Expr::var(var), const_(offset))
        };
        return Some(mul(const_(c), crate::expr::ln(arg)));
    }
    if is_var(g, var) {
        return Some(mul(const_(c), crate::expr::ln(Expr::var(var))));
    }
    None
}

/// ∫ 1 / ((x-a)(x-b)) dx = (ln|x-a| - ln|x-b|) / (a-b)  (for distinct a,b)
fn integrate_two_linear_factors(f: &Expr, g: &Expr, var: Symbol) -> Option<Expr> {
    if !as_const(f)?.is_one() {
        return None;
    }
    let (a, b) = product_of_two_linear_factors(g, var)?;
    if a == b {
        return None;
    }
    let denom = a - b;
    let xa = add(Expr::var(var), const_(-a));
    let xb = add(Expr::var(var), const_(-b));
    Some(div(
        sub(crate::expr::ln(xa), crate::expr::ln(xb)),
        const_(denom),
    ))
}

fn product_of_two_linear_factors(g: &Expr, var: Symbol) -> Option<(Rational, Rational)> {
    let (l, r) = match g.kind() {
        ExprKind::Mul(a, b) => (a, b),
        _ => return None,
    };
    let root = |e: &Expr| -> Option<Rational> {
        var_plus_const(e, var).map(|offset| -offset)
    };
    Some((root(l)?, root(r)?))
}

/// ∫ (px + q) / (x - a) dx — polynomial long division then integrate
fn integrate_polynomial_over_linear(f: Expr, g: Expr, var: Symbol) -> Option<Expr> {
    let offset = var_plus_const(&g, var)?;
    let a = -offset;
    let num_terms = super::util::collect_polynomial_terms(&f, var).ok()?;
    if num_terms.is_empty() {
        return Some(const_(Rational::zero()));
    }
    let degree = *num_terms.keys().max()?;
    if degree != 1 {
        return None;
    }
    // (px+q)/(x-a) = p + (pa+q)/(x-a) when written as division; general: synthetic division
    let p = num_terms.get(&1).copied().unwrap_or(Rational::zero());
    let q = num_terms.get(&0).copied().unwrap_or(Rational::zero());
    let remainder = p * a + q;
    let arg = add(Expr::var(var), const_(-a));
    let integ = if remainder.is_zero() {
        mul(const_(p), Expr::var(var))
    } else {
        mul(const_(p), Expr::var(var)) + mul(const_(remainder), crate::expr::ln(arg))
    };
    integrate(integ, var).ok()
}
