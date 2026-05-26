use crate::expr::{const_, div, mul, pow, Expr};
use crate::rational::Rational;
use crate::symbol::Symbol;

use super::util::{affine_form, affine_outer_integral, contains_var, is_var};
use super::IntegrateError;

pub fn integrate_sin(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        return Ok(-crate::expr::cos(e.clone()));
    }
    if let Some((k, _b)) = affine_form(e, var) {
        if _b.is_zero() {
            return Ok(affine_outer_integral(k, -crate::expr::cos(e.clone())));
        }
        // ∫ sin(kx+b) dx = -cos(kx+b)/k
        return Ok(affine_outer_integral(k, -crate::expr::cos(e.clone())));
    }
    if !contains_var(e, var) {
        return Ok(crate::expr::sin(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_cos(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        return Ok(crate::expr::sin(e.clone()));
    }
    if let Some((k, _b)) = affine_form(e, var) {
        return Ok(affine_outer_integral(k, crate::expr::sin(e.clone())));
    }
    if !contains_var(e, var) {
        return Ok(crate::expr::cos(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_exp(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        return Ok(crate::expr::exp(e.clone()));
    }
    if let Some((k, _b)) = affine_form(e, var) {
        return Ok(affine_outer_integral(k, crate::expr::exp(e.clone())));
    }
    if !contains_var(e, var) {
        return Ok(crate::expr::exp(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

pub fn integrate_tan(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        return Ok(-crate::expr::ln(crate::expr::cos(e.clone())));
    }
    if let Some((k, _b)) = affine_form(e, var) {
        return Ok(affine_outer_integral(
            k,
            -crate::expr::ln(crate::expr::cos(e.clone())),
        ));
    }
    if !contains_var(e, var) {
        return Ok(crate::expr::tan(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

/// ∫ ln(x) dx = x ln(x) - x; for ln(kx+b) use affine substitution.
pub fn integrate_ln(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        let x = Expr::var(var);
        return Ok(mul(x.clone(), crate::expr::ln(x.clone())) - x);
    }
    if let Some((k, _b)) = affine_form(e, var) {
        if k.is_zero() {
            return Err(IntegrateError::NoRule);
        }
        // ∫ ln(kx+b) dx = ((kx+b) ln(kx+b) - (kx+b)) / k
        let inner = e.clone();
        let term = mul(inner.clone(), crate::expr::ln(inner.clone())) - inner;
        return Ok(affine_outer_integral(k, term));
    }
    if !contains_var(e, var) {
        return Ok(crate::expr::ln(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}

/// ∫ atan(x) dx = x·atan(x) − ½ ln(1+x²).
pub fn integrate_atan(e: &Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if is_var(e, var) {
        let x = Expr::var(var);
        return Ok(
            mul(x.clone(), crate::expr::atan(x.clone()))
                - div(
                    crate::expr::ln(
                        const_(Rational::one())
                            + pow(x.clone(), const_(Rational::from(2))),
                    ),
                    const_(Rational::from(2)),
                ),
        );
    }
    if let Some((k, _b)) = affine_form(e, var) {
        if k.is_zero() {
            return Err(IntegrateError::NoRule);
        }
        let u = e.clone();
        let term = mul(u.clone(), crate::expr::atan(u.clone()))
            - div(
                crate::expr::ln(
                    const_(Rational::one()) + pow(u, const_(Rational::from(2))),
                ),
                const_(Rational::from(2)),
            );
        return Ok(affine_outer_integral(k, term));
    }
    if !contains_var(e, var) {
        return Ok(crate::expr::atan(e.clone()) * Expr::var(var));
    }
    Err(IntegrateError::NoRule)
}
