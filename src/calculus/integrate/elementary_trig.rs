//! Elementary integrals for extended trigonometric and hyperbolic functions.

use crate::expr::{
    acos, acosh, acot, acoth, acsc, acsch, add, asin, asinh, asec, asech, atan, atanh, const_,
    cosh, cot, csc,
    div, ln, mul, neg, pow, sec, sin, sinh, sub, tan, tanh, Expr,
};
use crate::rational::Rational;
use crate::symbol::Symbol;

use super::util::{affine_form, affine_outer_integral, contains_var, is_var};
use super::IntegrateError;

fn half() -> Expr {
    const_(Rational::new(1, 2))
}

fn sqrt_expr(e: Expr) -> Expr {
    pow(e, half())
}

fn integrate_affine(
    e: &Expr,
    var: Symbol,
    at_var: fn(Expr) -> Expr,
) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        return Ok(at_var(Expr::var(var)));
    }
    if let Some((k, _)) = affine_form(e, var) {
        return Ok(affine_outer_integral(k, at_var(e.clone())));
    }
    if !contains_var(e, var) {
        return Ok(at_var(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_cot(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    integrate_affine(e, var, |u| ln(sin(u)))
}

pub fn integrate_sec(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    integrate_affine(e, var, |u| ln(abs_sum(sec(u.clone()), tan(u))))
}

pub fn integrate_csc(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    integrate_affine(e, var, |u| neg(ln(abs_sum(csc(u.clone()), cot(u)))))
}

fn abs_sum(a: Expr, b: Expr) -> Expr {
    // ln(|a+b|) is represented as ln(a+b); sign not tracked in exact integration.
    add(a, b)
}

pub fn integrate_asin(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        let x = Expr::var(var);
        return Ok(
            mul(x.clone(), asin(x.clone()))
                + sqrt_expr(sub(const_(Rational::one()), pow(x, const_(Rational::from(2))))),
        );
    }
    if !contains_var(e, var) {
        return Ok(asin(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_acos(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        let x = Expr::var(var);
        return Ok(
            mul(x.clone(), acos(x.clone()))
                - sqrt_expr(sub(const_(Rational::one()), pow(x, const_(Rational::from(2))))),
        );
    }
    if !contains_var(e, var) {
        return Ok(acos(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_acot(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        let x = Expr::var(var);
        return Ok(
            mul(x.clone(), acot(x.clone()))
                + div(
                    ln(add(const_(Rational::one()), pow(x, const_(Rational::from(2))))),
                    const_(Rational::from(2)),
                ),
        );
    }
    if !contains_var(e, var) {
        return Ok(acot(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_sinh(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    integrate_affine(e, var, cosh)
}

pub fn integrate_cosh(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    integrate_affine(e, var, sinh)
}

pub fn integrate_tanh(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    integrate_affine(e, var, |u| ln(cosh(u)))
}

pub fn integrate_coth(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    integrate_affine(e, var, |u| ln(sinh(u)))
}

pub fn integrate_sech(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    integrate_affine(e, var, |u| atan(sinh(u)))
}

pub fn integrate_csch(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    integrate_affine(e, var, |u| ln(tanh(div(u, const_(Rational::from(2))))))
}

pub fn integrate_asinh(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        let x = Expr::var(var);
        return Ok(
            mul(x.clone(), asinh(x.clone()))
                + sqrt_expr(add(const_(Rational::one()), pow(x, const_(Rational::from(2))))),
        );
    }
    if !contains_var(e, var) {
        return Ok(asinh(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_acosh(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        let x = Expr::var(var);
        return Ok(
            mul(x.clone(), acosh(x.clone()))
                - sqrt_expr(sub(pow(x, const_(Rational::from(2))), const_(Rational::one()))),
        );
    }
    if !contains_var(e, var) {
        return Ok(acosh(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_atanh(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        let x = Expr::var(var);
        return Ok(
            mul(x.clone(), atanh(x.clone()))
                + div(
                    ln(sub(const_(Rational::one()), pow(x, const_(Rational::from(2))))),
                    const_(Rational::from(2)),
                ),
        );
    }
    if !contains_var(e, var) {
        return Ok(atanh(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_acoth(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        let x = Expr::var(var);
        return Ok(
            mul(x.clone(), acoth(x.clone()))
                + div(
                    ln(sub(pow(x, const_(Rational::from(2))), const_(Rational::one()))),
                    const_(Rational::from(2)),
                ),
        );
    }
    if !contains_var(e, var) {
        return Ok(acoth(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

/// ∫ asec(x) dx, ∫ asech(x) dx, ∫ acsc(x) dx, ∫ acsch(x) dx — constant inner only for now.
pub fn integrate_asec(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if !contains_var(e, var) {
        return Ok(asec(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_acsc(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if !contains_var(e, var) {
        return Ok(acsc(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_asech(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if !contains_var(e, var) {
        return Ok(asech(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_acsch(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if !contains_var(e, var) {
        return Ok(acsch(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}
