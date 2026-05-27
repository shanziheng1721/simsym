use crate::expr::{Expr, ExprKind};
use crate::rational::{rat_pow_int, Rational, RationalPowError};
use crate::symbol::Symbol;
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
}

pub fn eval(expr: &Expr, env: &[(Symbol, Rational)]) -> Result<Rational, EvalError> {
    match expr.kind() {
        ExprKind::Const(c) => Ok(*c),
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
        | ExprKind::Atan(_)
        | ExprKind::Exp(_)
        | ExprKind::Ln(_) => Err(EvalError::Undefined(
            "transcendental functions require eval_f64",
        )),
    }
}

fn eval_pow(base: &Expr, exp: &Expr, env: &[(Symbol, Rational)]) -> Result<Rational, EvalError> {
    let b = eval(base, env)?;
    if let ExprKind::Const(e) = exp.kind() {
        if let Some(n) = e.as_integer() {
            return rat_pow_int(b, n).map_err(EvalError::from);
        }
    }
    Err(EvalError::Undefined("non-integer exponent in exact eval"))
}

pub fn eval_f64(expr: &Expr, env: &[(Symbol, f64)]) -> Result<f64, EvalError> {
    match expr.kind() {
        ExprKind::Const(c) => Ok((*c).to_f64()),
        ExprKind::Var(s) => env
            .iter()
            .find(|(sym, _)| sym == s)
            .map(|(_, v)| *v)
            .ok_or(EvalError::UnboundSymbol(*s)),
        ExprKind::Add(a, b) => Ok(eval_f64(a, env)? + eval_f64(b, env)?),
        ExprKind::Sub(a, b) => Ok(eval_f64(a, env)? - eval_f64(b, env)?),
        ExprKind::Mul(a, b) => Ok(eval_f64(a, env)? * eval_f64(b, env)?),
        ExprKind::Div(a, b) => {
            let d = eval_f64(b, env)?;
            if d == 0.0 {
                return Err(EvalError::DivisionByZero);
            }
            Ok(eval_f64(a, env)? / d)
        }
        ExprKind::Neg(e) => Ok(-eval_f64(e, env)?),
        ExprKind::Pow(base, exp) => Ok(eval_f64(base, env)?.powf(eval_f64(exp, env)?)),
        ExprKind::Sin(e) => Ok(eval_f64(e, env)?.sin()),
        ExprKind::Cos(e) => Ok(eval_f64(e, env)?.cos()),
        ExprKind::Tan(e) => Ok(eval_f64(e, env)?.tan()),
        ExprKind::Atan(e) => Ok(eval_f64(e, env)?.atan()),
        ExprKind::Exp(e) => Ok(eval_f64(e, env)?.exp()),
        ExprKind::Ln(e) => {
            let v = eval_f64(e, env)?;
            if v <= 0.0 {
                Err(EvalError::Undefined("ln of non-positive"))
            } else {
                Ok(v.ln())
            }
        }
    }
}
