use crate::expr::{add, const_, div, mul, pow, Expr, ExprKind};
use crate::rational::Rational;
use crate::symbol::Symbol;

use super::util::{as_const, contains_var, is_var, var_plus_const};
use super::IntegrateError;

pub fn integrate_pow(base: Expr, exp: Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if let Some(n) = as_const(&exp).and_then(|r| r.as_integer()) {
        if n >= 2 {
            if let ExprKind::Sin(inner) = base.kind() {
                if n % 2 == 0 {
                    return integrate_sin_even_power(inner, n, var);
                }
                return integrate_sin_odd_power(inner, n, var);
            }
            if let ExprKind::Cos(inner) = base.kind() {
                if n % 2 == 0 {
                    return integrate_cos_even_power(inner, n, var);
                }
                return integrate_cos_odd_power(inner, n, var);
            }
            if let ExprKind::Tan(inner) = base.kind() {
                if is_var(inner, var) || super::util::affine_form(inner, var).is_some() {
                    return integrate_tan_power(inner, n, var);
                }
            }
        }
        if n < -1 && n % 2 == 0 {
            if let ExprKind::Cos(inner) = base.kind() {
                if let Some((k, _)) = super::util::affine_form(inner, var) {
                    let m = (-n) as i64;
                    let inner_int = integrate_sec_power(inner, m, var)?;
                    return Ok(super::util::affine_outer_integral(k, inner_int));
                }
            }
        }
    }

    if let ExprKind::Var(s) = base.kind() {
        if *s != var {
            return Ok(pow(base, exp) * Expr::var(var));
        }
        if let Some(k) = as_const(&exp).and_then(|r| r.as_integer()) {
            return integrate_var_power(k, var);
        }
    }

    if let Some(offset) = var_plus_const(&base, var) {
        if let Some(k) = as_const(&exp).and_then(|r| r.as_integer()) {
            if k == -1 {
                let shifted = if offset.is_zero() {
                    Expr::var(var)
                } else {
                    add(Expr::var(var), const_(offset))
                };
                return Ok(crate::expr::ln(shifted));
            }
            if k < -1 {
                return Err(IntegrateError::NoRule);
            }
            if k == 0 {
                return Ok(Expr::var(var));
            }
            let shifted = add(Expr::var(var), const_(offset));
            let new_exp = k + 1;
            return Ok(pow(shifted, const_(Rational::from(new_exp))) / const_(Rational::from(new_exp)));
        }
    }

    if !contains_var(&base, var) && !contains_var(&exp, var) {
        return Ok(pow(base, exp) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub(crate) fn integrate_var_power(k: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    if k == -1 {
        return Ok(crate::expr::ln(Expr::var(var)));
    }
    if k < -1 || k == 0 {
        return Err(IntegrateError::NoRule);
    }
    let new_exp = k + 1;
    Ok(pow(Expr::var(var), const_(Rational::from(new_exp))) / const_(Rational::from(new_exp)))
}

/// ∫ sin(x)^n dx = -sin^(n-1)(x)cos(x)/n + (n-1)/n ∫ sin^(n-2)(x) dx
fn integrate_sin_even_power(inner: &Expr, n: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    if n == 2 {
        let x = inner.clone();
        let two_x = mul(const_(Rational::from(2)), x.clone());
        return Ok(
            div(Expr::var(var), const_(Rational::from(2)))
                - div(crate::expr::sin(two_x), const_(Rational::from(4))),
        );
    }
    if n < 4 || n % 2 != 0 {
        return Err(IntegrateError::NoRule);
    }
    sin_power_reduction(inner, n, var, true)
}

fn integrate_sin_odd_power(inner: &Expr, n: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    if n == 1 {
        return super::elementary::integrate_sin(inner, var);
    }
    if n < 3 || n % 2 == 0 {
        return Err(IntegrateError::NoRule);
    }
    sin_power_reduction(inner, n, var, true)
}

fn integrate_cos_even_power(inner: &Expr, n: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    if n == 2 {
        let x = inner.clone();
        let two_x = mul(const_(Rational::from(2)), x.clone());
        return Ok(
            div(Expr::var(var), const_(Rational::from(2)))
                + div(crate::expr::sin(two_x), const_(Rational::from(4))),
        );
    }
    if n < 4 || n % 2 != 0 {
        return Err(IntegrateError::NoRule);
    }
    sin_power_reduction(inner, n, var, false)
}

fn integrate_cos_odd_power(inner: &Expr, n: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    if n == 1 {
        return super::elementary::integrate_cos(inner, var);
    }
    if n < 3 || n % 2 == 0 {
        return Err(IntegrateError::NoRule);
    }
    sin_power_reduction(inner, n, var, false)
}

fn sin_power_reduction(
    inner: &Expr,
    n: i64,
    var: Symbol,
    for_sin: bool,
) -> Result<Expr, IntegrateError> {
    let x = inner.clone();
    let trig = if for_sin {
        crate::expr::sin(x.clone())
    } else {
        crate::expr::cos(x.clone())
    };
    let partner = if for_sin {
        crate::expr::cos(x.clone())
    } else {
        crate::expr::sin(x.clone())
    };
    let sign = if for_sin { -1 } else { 1 };
    let part = mul(
        const_(Rational::new(sign, n)),
        mul(pow(trig, const_(Rational::from(n - 1))), partner),
    );
    let rest = if for_sin {
        integrate_sin_power_n(&x, n - 2, var)
    } else {
        integrate_cos_power_n(&x, n - 2, var)
    }?;
    Ok(part + mul(const_(Rational::new(n - 1, n)), rest))
}

pub(crate) fn integrate_sin_power_direct(
    inner: &Expr,
    n: i64,
    var: Symbol,
) -> Result<Expr, IntegrateError> {
    integrate_sin_power_n(inner, n, var)
}

pub(crate) fn integrate_cos_power_direct(
    inner: &Expr,
    n: i64,
    var: Symbol,
) -> Result<Expr, IntegrateError> {
    integrate_cos_power_n(inner, n, var)
}

/// ∫ sec^m(x) dx via ∫ cos(x)^{-m} dx (m ≥ 2).
pub(crate) fn integrate_sec_power(inner: &Expr, m: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    integrate_cos_negative_power(inner, m, var)
}

fn integrate_cos_negative_power(
    inner: &Expr,
    m: i64,
    var: Symbol,
) -> Result<Expr, IntegrateError> {
    if m < 2 {
        return Err(IntegrateError::NoRule);
    }
    if m == 2 {
        return Ok(crate::expr::tan(inner.clone()));
    }
    let sin_u = crate::expr::sin(inner.clone());
    let cos_u = crate::expr::cos(inner.clone());
    let part = mul(sin_u, pow(cos_u.clone(), const_(Rational::from(1 - m))))
        / const_(Rational::from(m - 1));
    let rest = integrate_cos_negative_power(inner, m - 2, var)?;
    Ok(part + mul(const_(Rational::new(m - 2, m - 1)), rest))
}

fn integrate_sin_power_n(inner: &Expr, n: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    match n {
        1 => super::elementary::integrate_sin(inner, var),
        2 => integrate_sin_even_power(inner, 2, var),
        _ if n % 2 == 0 => integrate_sin_even_power(inner, n, var),
        _ => integrate_sin_odd_power(inner, n, var),
    }
}

fn integrate_cos_power_n(inner: &Expr, n: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    match n {
        1 => super::elementary::integrate_cos(inner, var),
        2 => integrate_cos_even_power(inner, 2, var),
        _ if n % 2 == 0 => integrate_cos_even_power(inner, n, var),
        _ => integrate_cos_odd_power(inner, n, var),
    }
}

/// ∫ tan^n(x) dx = tan^(n-1)/(n-1) − ∫ tan^(n-2) dx
fn integrate_tan_power(inner: &Expr, n: i64, var: Symbol) -> Result<Expr, IntegrateError> {
    if n < 0 {
        return Err(IntegrateError::NoRule);
    }
    if n == 0 {
        return Ok(Expr::var(var));
    }
    if n == 1 {
        return super::elementary::integrate_tan(inner, var);
    }
    let tan_x = crate::expr::tan(inner.clone());
    let part = pow(tan_x, const_(Rational::from(n - 1))) / const_(Rational::from(n - 1));
    let rest = integrate_tan_power(inner, n - 2, var)?;
    Ok(part - rest)
}
