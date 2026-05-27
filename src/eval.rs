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
        ExprKind::Cot(e) => Ok(1.0 / eval_f64(e, env)?.tan()),
        ExprKind::Sec(e) => Ok(1.0 / eval_f64(e, env)?.cos()),
        ExprKind::Csc(e) => Ok(1.0 / eval_f64(e, env)?.sin()),
        ExprKind::Asin(e) => Ok(eval_f64(e, env)?.asin()),
        ExprKind::Acos(e) => Ok(eval_f64(e, env)?.acos()),
        ExprKind::Atan(e) => Ok(eval_f64(e, env)?.atan()),
        ExprKind::Acot(e) => Ok(std::f64::consts::FRAC_PI_2 - eval_f64(e, env)?.atan()),
        ExprKind::Asec(e) => Ok((1.0 / eval_f64(e, env)?).acos()),
        ExprKind::Acsc(e) => Ok((1.0 / eval_f64(e, env)?).asin()),
        ExprKind::Sinh(e) => Ok(eval_f64(e, env)?.sinh()),
        ExprKind::Cosh(e) => Ok(eval_f64(e, env)?.cosh()),
        ExprKind::Tanh(e) => Ok(eval_f64(e, env)?.tanh()),
        ExprKind::Coth(e) => Ok(1.0 / eval_f64(e, env)?.tanh()),
        ExprKind::Sech(e) => Ok(1.0 / eval_f64(e, env)?.cosh()),
        ExprKind::Csch(e) => Ok(1.0 / eval_f64(e, env)?.sinh()),
        ExprKind::Asinh(e) => Ok(eval_f64(e, env)?.asinh()),
        ExprKind::Acosh(e) => Ok(eval_f64(e, env)?.acosh()),
        ExprKind::Atanh(e) => Ok(eval_f64(e, env)?.atanh()),
        ExprKind::Acoth(e) => Ok(acoth_f64(eval_f64(e, env)?)),
        ExprKind::Asech(e) => Ok(asech_f64(eval_f64(e, env)?)),
        ExprKind::Acsch(e) => Ok(acsch_f64(eval_f64(e, env)?)),
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

fn acoth_f64(x: f64) -> f64 {
    0.5 * ((x + 1.0) / (x - 1.0)).ln()
}

fn asech_f64(x: f64) -> f64 {
    (1.0 / x).acosh()
}

fn acsch_f64(x: f64) -> f64 {
    (1.0 / x).asinh()
}
