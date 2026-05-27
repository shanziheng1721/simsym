//! Symbolic integration (rule-based / heuristic).
//!
//! Algorithm stack (increasing strength):
//! 1. Linearity — sum, constant factor
//! 2. Elementary templates — polynomials, `1/x`, `sin`/`cos`/`exp`/`tan`/`ln` with affine inner `k*x+b`
//! 3. Products — `e^{ax+b}×sin/cos(cx+d)`, `sin^m cos^n`, `x^n ln x`, `sin(ax+b)cos(cx+d)`
//! 4. Powers — `(x+c)^n`, `sin^n`/`cos^n`/`tan^n`/`sec^n` reductions, `cos⁻²` → `tan`
//! 5. Rational — `c/(x-a)`, `1/((x-a)(x-b))`, `P/Q` partial fractions (`deg Q ≤ 6`), `atan` for `1/(x²+a²)`
//! 6. Substitution — `u^n u'`, elementary `f(u)·u'`
//! 7. Integration by parts — depth-limited, with direct integrate fallback
//!
//! Not implemented: full Risch decision procedure, algebraic extensions,
//! or RUBI-scale rule databases.

mod elementary;
mod power;
mod product;
mod rational;
mod rational_function;
mod substitution;
mod trig_product;
mod util;

use crate::expr::{Expr, ExprKind};
use crate::symbol::Symbol;

pub use util::contains_var;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum IntegrateError {
    #[error("no integration rule applies to this expression")]
    NoRule,
    #[error("logarithm of negative expression is not supported in exact integration")]
    NonPositiveLogArgument,
}

fn finalize_integral(expr: Expr) -> Expr {
    #[cfg(feature = "simplify")]
    {
        expr.simplify()
    }
    #[cfg(not(feature = "simplify"))]
    {
        expr
    }
}

pub fn integrate(expr: Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if let Some(result) = substitution::try_u_substitution(&expr, var) {
        return Ok(finalize_integral(result));
    }
    let result = integrate_kind(expr.kind(), var)?;
    Ok(finalize_integral(result))
}

fn integrate_kind(kind: &ExprKind, var: Symbol) -> Result<Expr, IntegrateError> {
    match kind {
        ExprKind::Const(c) => Ok(crate::expr::const_(*c) * Expr::var(var)),
        ExprKind::Var(s) => {
            if *s == var {
                Ok(power::integrate_var_power(1, var)?)
            } else {
                Ok(Expr::var(*s) * Expr::var(var))
            }
        }
        ExprKind::Add(a, b) => Ok(integrate(a.clone(), var)? + integrate(b.clone(), var)?),
        ExprKind::Sub(a, b) => Ok(integrate(a.clone(), var)? - integrate(b.clone(), var)?),
        ExprKind::Neg(e) => Ok(-integrate(e.clone(), var)?),
        ExprKind::Mul(f, g) => product::integrate_product(f.clone(), g.clone(), var),
        ExprKind::Div(f, g) => rational::integrate_div(f.clone(), g.clone(), var),
        ExprKind::Pow(base, exp) => power::integrate_pow(base.clone(), exp.clone(), var),
        ExprKind::Sin(e) => elementary::integrate_sin(e, var),
        ExprKind::Cos(e) => elementary::integrate_cos(e, var),
        ExprKind::Exp(e) => elementary::integrate_exp(e, var),
        ExprKind::Tan(e) => elementary::integrate_tan(e, var),
        ExprKind::Ln(e) => elementary::integrate_ln(e, var),
        ExprKind::Atan(e) => elementary::integrate_atan(e, var),
    }
}
