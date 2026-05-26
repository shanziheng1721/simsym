use crate::expr::{add, const_, mul, pow, Expr, ExprKind};
use crate::rational::Rational;
use crate::symbol::Symbol;
use std::collections::BTreeMap;
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum IntegrateError {
    #[error("no integration rule applies to this expression")]
    NoRule,
    #[error("logarithm of negative expression is not supported in exact integration")]
    NonPositiveLogArgument,
}

pub fn integrate(expr: Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    let result = integrate_kind(expr.kind(), var)?;
    Ok(result.simplify())
}

fn integrate_kind(kind: &ExprKind, var: Symbol) -> Result<Expr, IntegrateError> {
    match kind {
        ExprKind::Const(c) => Ok(const_(*c) * Expr::var(var)),
        ExprKind::Var(s) => {
            if *s == var {
                Ok(pow(Expr::var(var), const_(Rational::from(2))) / const_(Rational::from(2)))
            } else {
                Ok(Expr::var(*s) * Expr::var(var))
            }
        }
        ExprKind::Add(a, b) => Ok(integrate(a.clone(), var)? + integrate(b.clone(), var)?),
        ExprKind::Sub(a, b) => Ok(integrate(a.clone(), var)? - integrate(b.clone(), var)?),
        ExprKind::Neg(e) => Ok(-integrate(e.clone(), var)?),
        ExprKind::Mul(f, g) => integrate_product(f.clone(), g.clone(), var),
        ExprKind::Div(f, g) => integrate_div(f.clone(), g.clone(), var),
        ExprKind::Pow(base, exp) => integrate_pow(base.clone(), exp.clone(), var),
        ExprKind::Sin(e) => integrate_sin(e, var),
        ExprKind::Cos(e) => integrate_cos(e, var),
        ExprKind::Exp(e) => integrate_exp(e, var),
        ExprKind::Tan(_) | ExprKind::Ln(_) => Err(IntegrateError::NoRule),
    }
}

fn integrate_sin(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        return Ok(-crate::expr::cos(e.clone()));
    }
    if linear_in_var(e, var) {
        let k = linear_coeff(e, var)?;
        return Ok(-crate::expr::cos(e.clone()) / const_(k));
    }
    if !contains_var(e, var) {
        return Ok(crate::expr::sin(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

fn integrate_cos(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        return Ok(crate::expr::sin(e.clone()));
    }
    if linear_in_var(e, var) {
        let k = linear_coeff(e, var)?;
        return Ok(crate::expr::sin(e.clone()) / const_(k));
    }
    if !contains_var(e, var) {
        return Ok(crate::expr::cos(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

fn integrate_exp(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        return Ok(crate::expr::exp(e.clone()));
    }
    if linear_in_var(e, var) {
        let k = linear_coeff(e, var)?;
        return Ok(crate::expr::exp(e.clone()) / const_(k));
    }
    if !contains_var(e, var) {
        return Ok(crate::expr::exp(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

fn integrate_product(f: Expr, g: Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if !contains_var(&f, var) {
        return Ok(f * integrate(g, var)?);
    }
    if !contains_var(&g, var) {
        return Ok(g * integrate(f, var)?);
    }
    if let Some(integ) = integrate_exp_times_polynomial(&f, &g, var) {
        return Ok(integ);
    }
    if let Some(integ) = integrate_exp_times_polynomial(&g, &f, var) {
        return Ok(integ);
    }
    // integration by parts for simple u * v' patterns: ∫ u v' = u v - ∫ u' v
    let fp = f.clone().diff(var);
    if !contains_var(&fp, var) {
        let g_int = integrate(g, var)?;
        return Ok(f * g_int.clone() - fp * integrate(g_int, var)?);
    }
    Err(IntegrateError::NoRule)
}

fn integrate_div(f: Expr, g: Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if !contains_var(&g, var) {
        return Ok(integrate(f, var)? / g);
    }
    // 1/x
    if matches!(f.kind(), ExprKind::Const(c) if c.is_one())
        && matches!(g.kind(), ExprKind::Var(s) if *s == var)
    {
        return Ok(crate::expr::ln(Expr::var(var)));
    }
    if matches!(f.kind(), ExprKind::Const(c) if !contains_var(&f, var))
        && matches!(g.kind(), ExprKind::Var(s) if *s == var)
    {
        if let ExprKind::Const(c) = f.kind() {
            return Ok(const_(*c) * crate::expr::ln(Expr::var(var)));
        }
    }
    Err(IntegrateError::NoRule)
}

/// ∫ e^x * P(x) dx when `exp_side` is `exp(x)` and `poly` is a polynomial in `var`.
fn integrate_exp_times_polynomial(
    exp_side: &Expr,
    poly: &Expr,
    var: Symbol,
) -> Option<Expr> {
    if !is_exp_of_var(exp_side, var) {
        return None;
    }
    let terms = collect_polynomial_terms(poly, var).ok()?;
    if terms.is_empty() {
        return Some(const_(Rational::zero()));
    }
    let mut poly_part = const_(Rational::zero());
    for (n, coeff) in terms.iter().rev() {
        if coeff.is_zero() {
            continue;
        }
        let p = exp_integrate_polynomial_factor(*n, var);
        poly_part = add(poly_part, mul(const_(*coeff), p));
    }
    Some(mul(
        crate::expr::exp(Expr::var(var)),
        poly_part,
    ))
}

fn collect_polynomial_terms(
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
        ExprKind::Pow(base, exp) if matches!(base.kind(), ExprKind::Var(s) if *s == var) => {
            if let ExprKind::Const(n) = exp.kind() {
                if let Some(k) = n.as_integer() {
                    if k >= 0 {
                        *out.entry(k).or_insert(Rational::zero()) += coeff;
                        return Ok(());
                    }
                }
            }
            Err(IntegrateError::NoRule)
        }
        _ => Err(IntegrateError::NoRule),
    }
}

/// Polynomial P_n(x) such that ∫ x^n e^x dx = e^x * P_n(x), in expanded form.
fn exp_integrate_polynomial_factor(n: i64, var: Symbol) -> Expr {
    let mut coeffs: BTreeMap<i64, Rational> = BTreeMap::new();
    coeffs.insert(0, Rational::one());
    for i in 1..=n {
        let mut next = BTreeMap::new();
        next.insert(i, Rational::one());
        for (k, c) in &coeffs {
            let slot = next.entry(*k).or_insert(Rational::zero());
            *slot = *slot - Rational::from(i) * *c;
        }
        coeffs = next;
    }
    let x = Expr::var(var);
    let mut sum = const_(Rational::zero());
    for (k, c) in coeffs {
        if c.is_zero() {
            continue;
        }
        let term = if k == 0 {
            const_(c)
        } else if k == 1 {
            mul(const_(c), x.clone())
        } else {
            mul(const_(c), pow(x.clone(), const_(Rational::from(k))))
        };
        sum = add(sum, term);
    }
    sum
}

fn integrate_pow(base: Expr, exp: Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if let ExprKind::Var(s) = base.kind() {
        if *s != var {
            return Ok(pow(base, exp) * Expr::var(var));
        }
        if let ExprKind::Const(n) = exp.kind() {
            if let Some(k) = n.as_integer() {
                if k == -1 {
                    return Ok(crate::expr::ln(Expr::var(var)));
                }
                if k < -1 || k == 0 {
                    return Err(IntegrateError::NoRule);
                }
                let new_exp = k + 1;
                return Ok(
                    pow(Expr::var(var), const_(Rational::from(new_exp)))
                        / const_(Rational::from(new_exp)),
                );
            }
        }
    }
    if !contains_var(&base, var) && !contains_var(&exp, var) {
        return Ok(pow(base, exp) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

fn contains_var(expr: &Expr, var: Symbol) -> bool {
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
        | ExprKind::Ln(e) => contains_var(e, var),
    }
}

fn is_var(e: &Expr, var: Symbol) -> bool {
    matches!(e.kind(), ExprKind::Var(s) if *s == var)
}

fn is_euler_var(e: &Expr) -> bool {
    matches!(e.kind(), ExprKind::Var(s) if s.name() == "e")
}

fn is_exp_of_var(e: &Expr, var: Symbol) -> bool {
    match e.kind() {
        ExprKind::Exp(inner) => is_var(inner, var),
        ExprKind::Pow(base, exp) => is_euler_var(base) && is_var(exp, var),
        _ => false,
    }
}

fn linear_in_var(e: &Expr, var: Symbol) -> bool {
    matches!(e.kind(), ExprKind::Mul(l, r) if (
        matches!(l.kind(), ExprKind::Const(_)) && matches!(r.kind(), ExprKind::Var(s) if *s == var)
    ) || (
        matches!(r.kind(), ExprKind::Const(_)) && matches!(l.kind(), ExprKind::Var(s) if *s == var)
    ))
}

fn linear_coeff(e: &Expr, var: Symbol) -> Result<Rational, IntegrateError> {
    match e.kind() {
        ExprKind::Mul(l, r) => {
            if matches!(r.kind(), ExprKind::Var(s) if *s == var) {
                if let ExprKind::Const(c) = l.kind() {
                    return Ok(*c);
                }
            }
            if matches!(l.kind(), ExprKind::Var(s) if *s == var) {
                if let ExprKind::Const(c) = r.kind() {
                    return Ok(*c);
                }
            }
        }
        _ => {}
    }
    Err(IntegrateError::NoRule)
}
