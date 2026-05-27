mod float;

use crate::expr::{Expr, ExprKind};
use crate::rational::{rat_pow_int, Rational, RationalPowError};
use crate::symbol::Symbol;

pub use float::eval_float;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum EvalError {
    #[error("division by zero")]
    DivisionByZero,
    #[error("undefined for given values: {0}")]
    Undefined(&'static str),
    #[error("symbol {0} is not bound")]
    UnboundSymbol(Symbol),
    #[error("overflow or unsupported rational operation")]
    Overflow,
    #[error(transparent)]
    RationalPow(#[from] RationalPowError),
    #[error("coefficient does not fit in exact i64 rational evaluation; enable `bigint` or use eval_f64")]
    CoefficientTooLarge,
}

pub fn eval(expr: &Expr, env: &[(Symbol, Rational)]) -> Result<Rational, EvalError> {
    match expr.kind() {
        ExprKind::Const(c) => c
            .try_as_rational()
            .ok_or(EvalError::CoefficientTooLarge),
        ExprKind::Var(s) => env
            .iter()
            .find(|(sym, _)| sym == s)
            .map(|(_, v)| *v)
            .ok_or(EvalError::UnboundSymbol(*s)),
        ExprKind::Add(a, b) => Ok(eval(a, env)? + eval(b, env)?),
        ExprKind::Sub(a, b) => Ok(eval(a, env)? - eval(b, env)?),
        ExprKind::Mul(a, b) => Ok(eval(a, env)? * eval(b, env)?),
        ExprKind::Div(a, b) => {
            let denom = eval(b, env)?;
            if denom.is_zero() {
                return Err(EvalError::DivisionByZero);
            }
            Ok(eval(a, env)? / denom)
        }
        ExprKind::Neg(e) => Ok(-eval(e, env)?),
        ExprKind::Pow(base, exp) => eval_pow(base, exp, env),
        ExprKind::Sin(_)
        | ExprKind::Cos(_)
        | ExprKind::Tan(_)
        | ExprKind::Cot(_)
        | ExprKind::Sec(_)
        | ExprKind::Csc(_)
        | ExprKind::Asin(_)
        | ExprKind::Acos(_)
        | ExprKind::Atan(_)
        | ExprKind::Acot(_)
        | ExprKind::Asec(_)
        | ExprKind::Acsc(_)
        | ExprKind::Sinh(_)
        | ExprKind::Cosh(_)
        | ExprKind::Tanh(_)
        | ExprKind::Coth(_)
        | ExprKind::Sech(_)
        | ExprKind::Csch(_)
        | ExprKind::Asinh(_)
        | ExprKind::Acosh(_)
        | ExprKind::Atanh(_)
        | ExprKind::Acoth(_)
        | ExprKind::Asech(_)
        | ExprKind::Acsch(_)
        | ExprKind::Exp(_)
        | ExprKind::Ln(_) => Err(EvalError::Undefined(
            "transcendental functions require eval_f32 or eval_f64",
        )),
    }
}

fn eval_pow(base: &Expr, exp: &Expr, env: &[(Symbol, Rational)]) -> Result<Rational, EvalError> {
    let b = eval(base, env)?;
    if let ExprKind::Const(c) = exp.kind() {
        if let Some(e) = c.try_as_rational() {
            if let Some(n) = e.as_integer() {
                return rat_pow_int(b, n).map_err(EvalError::from);
            }
        }
    }
    Err(EvalError::Undefined("non-integer exponent in exact eval"))
}

/// Evaluate numerically with `f64` bindings.
pub fn eval_f64(expr: &Expr, env: &[(Symbol, f64)]) -> Result<f64, EvalError> {
    eval_float(expr, env)
}

/// Evaluate numerically with `f32` bindings.
pub fn eval_f32(expr: &Expr, env: &[(Symbol, f32)]) -> Result<f32, EvalError> {
    eval_float(expr, env)
}
