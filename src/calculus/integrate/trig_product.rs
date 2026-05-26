//! ∫ sin^m(x) cos^n(x) and affine variants via peeling / reduction.

use crate::expr::{const_, mul, pow, Expr, ExprKind};
use crate::rational::Rational;
use crate::symbol::Symbol;

use super::util::{affine_form, affine_outer_integral, as_const, contains_var};
use super::IntegrateError;

/// Parse `sin(inner)^m * cos(inner)^n` (factors may appear in any order).
pub fn try_sin_cos_power_product(expr: &Expr, var: Symbol) -> Option<Expr> {
    let (inner, m, n) = parse_sin_cos_powers(expr)?;
    if !contains_var(&inner, var) {
        return None;
    }
    let (k, _) = affine_form(&inner, var)?;
    if k.is_zero() {
        return None;
    }
    let result = integrate_sin_m_cos_n(&inner, m, n, var).ok()?;
    Some(affine_outer_integral(k, result))
}

fn parse_sin_cos_powers(expr: &Expr) -> Option<(Expr, i64, i64)> {
    let mut sin_m = 0i64;
    let mut cos_n = 0i64;
    let mut inner: Option<Expr> = None;
    let mut rest = vec![expr.clone()];

    while let Some(f) = rest.pop() {
        match f.kind() {
            ExprKind::Mul(l, r) => {
                rest.push(l.clone());
                rest.push(r.clone());
            }
            ExprKind::Pow(base, exp) => {
                let p = as_const(exp)?.as_integer()?;
                if p < 0 {
                    return None;
                }
                match base.kind() {
                    ExprKind::Sin(i) => {
                        inner = merge_inner(inner.take(), i.clone());
                        sin_m += p;
                    }
                    ExprKind::Cos(i) => {
                        inner = merge_inner(inner.take(), i.clone());
                        cos_n += p;
                    }
                    _ => return None,
                }
            }
            ExprKind::Sin(i) => {
                inner = merge_inner(inner.take(), i.clone());
                sin_m += 1;
            }
            ExprKind::Cos(i) => {
                inner = merge_inner(inner.take(), i.clone());
                cos_n += 1;
            }
            ExprKind::Const(c) if c.is_one() => {}
            _ => return None,
        }
    }
    if sin_m == 0 && cos_n == 0 {
        return None;
    }
    let inner = inner?;
    Some((inner, sin_m, cos_n))
}

fn merge_inner(a: Option<Expr>, b: Expr) -> Option<Expr> {
    match a {
        None => Some(b),
        Some(prev) if prev == b => Some(prev),
        _ => None,
    }
}

fn integrate_sin_m_cos_n(inner: &Expr, m: i64, n: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    if m < 0 || n < 0 {
        return Err(IntegrateError::NoRule);
    }
    if m == 0 && n == 0 {
        return Ok(Expr::var(var));
    }
    if m == 0 {
        return super::power::integrate_cos_power_direct(inner, n, var);
    }
    if n == 0 {
        return super::power::integrate_sin_power_direct(inner, m, var);
    }
    if m % 2 == 1 {
        return sin_odd_peel(inner, m, n, var);
    }
    if n % 2 == 1 {
        return cos_odd_peel(inner, m, n, var);
    }
    both_even_reduction(inner, m, n, var)
}

/// ∫ sin^m cos^n = -sin^{m-1} cos^{n+1}/(n+1) + (m-1)/(n+1) ∫ sin^{m-2} cos^{n+2}
fn sin_odd_peel(inner: &Expr, m: i64, n: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    let sin_u = crate::expr::sin(inner.clone());
    let cos_u = crate::expr::cos(inner.clone());
    if m == 1 {
        return Ok(-pow(cos_u, const_(Rational::from(n + 1))) / const_(Rational::from(n + 1)));
    }
    let part = -mul(
        pow(sin_u, const_(Rational::from(m - 1))),
        pow(cos_u.clone(), const_(Rational::from(n + 1))),
    ) / const_(Rational::from(n + 1));
    let rest = integrate_sin_m_cos_n(inner, m - 2, n + 2, var)?;
    Ok(part + mul(const_(Rational::new(m - 1, n + 1)), rest))
}

/// ∫ sin^m cos^n = sin^{m+1} cos^{n-1}/(m+1) + (n-1)/(m+1) ∫ sin^{m+2} cos^{n-2}
fn cos_odd_peel(inner: &Expr, m: i64, n: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    let sin_u = crate::expr::sin(inner.clone());
    let cos_u = crate::expr::cos(inner.clone());
    if n == 1 {
        return Ok(pow(sin_u, const_(Rational::from(m + 1))) / const_(Rational::from(m + 1)));
    }
    let part = mul(
        pow(sin_u.clone(), const_(Rational::from(m + 1))),
        pow(cos_u, const_(Rational::from(n - 1))),
    ) / const_(Rational::from(m + 1));
    let rest = integrate_sin_m_cos_n(inner, m + 2, n - 2, var)?;
    Ok(part + mul(const_(Rational::new(n - 1, m + 1)), rest))
}

/// Both even: use sin² = 1 − cos² (or symmetric) to reduce total degree.
fn both_even_reduction(inner: &Expr, m: i64, n: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    if m >= 2 {
        let a = integrate_sin_m_cos_n(inner, m - 2, n, var)?;
        let b = integrate_sin_m_cos_n(inner, m - 2, n + 2, var)?;
        return Ok(a - b);
    }
    let a = integrate_sin_m_cos_n(inner, m, n - 2, var)?;
    let b = integrate_sin_m_cos_n(inner, m + 2, n - 2, var)?;
    Ok(a - b)
}
