//! Generic floating-point evaluation (`f32` / `f64`) via [`num_traits::Float`].

use crate::expr::{Expr, ExprKind};
use crate::symbol::Symbol;
use num_traits::float::Float;
use num_traits::FromPrimitive;

use super::EvalError;

pub fn eval_float<F>(expr: &Expr, env: &[(Symbol, F)]) -> Result<F, EvalError>
where
    F: Float + FromPrimitive,
{
    match expr.kind() {
        ExprKind::Const(c) => Ok(F::from_f64(c.to_f64()).expect("finite constant")),
        ExprKind::Var(s) => env
            .iter()
            .find(|(sym, _)| sym == s)
            .map(|(_, v)| *v)
            .ok_or(EvalError::UnboundSymbol(*s)),
        ExprKind::Add(a, b) => Ok(eval_float(a, env)? + eval_float(b, env)?),
        ExprKind::Sub(a, b) => Ok(eval_float(a, env)? - eval_float(b, env)?),
        ExprKind::Mul(a, b) => Ok(eval_float(a, env)? * eval_float(b, env)?),
        ExprKind::Div(a, b) => {
            let d = eval_float(b, env)?;
            if d.is_zero() {
                return Err(EvalError::DivisionByZero);
            }
            Ok(eval_float(a, env)? / d)
        }
        ExprKind::Neg(e) => Ok(-eval_float(e, env)?),
        ExprKind::Pow(base, exp) => Ok(eval_float(base, env)?.powf(eval_float(exp, env)?)),
        ExprKind::Sin(e) => Ok(eval_float(e, env)?.sin()),
        ExprKind::Cos(e) => Ok(eval_float(e, env)?.cos()),
        ExprKind::Tan(e) => Ok(eval_float(e, env)?.tan()),
        ExprKind::Cot(e) => Ok(F::one() / eval_float(e, env)?.tan()),
        ExprKind::Sec(e) => Ok(F::one() / eval_float(e, env)?.cos()),
        ExprKind::Csc(e) => Ok(F::one() / eval_float(e, env)?.sin()),
        ExprKind::Asin(e) => Ok(eval_float(e, env)?.asin()),
        ExprKind::Acos(e) => Ok(eval_float(e, env)?.acos()),
        ExprKind::Atan(e) => Ok(eval_float(e, env)?.atan()),
        ExprKind::Acot(e) => Ok(frac_pi_2::<F>() - eval_float(e, env)?.atan()),
        ExprKind::Asec(e) => Ok((F::one() / eval_float(e, env)?).acos()),
        ExprKind::Acsc(e) => Ok((F::one() / eval_float(e, env)?).asin()),
        ExprKind::Sinh(e) => Ok(eval_float(e, env)?.sinh()),
        ExprKind::Cosh(e) => Ok(eval_float(e, env)?.cosh()),
        ExprKind::Tanh(e) => Ok(eval_float(e, env)?.tanh()),
        ExprKind::Coth(e) => Ok(F::one() / eval_float(e, env)?.tanh()),
        ExprKind::Sech(e) => Ok(F::one() / eval_float(e, env)?.cosh()),
        ExprKind::Csch(e) => Ok(F::one() / eval_float(e, env)?.sinh()),
        ExprKind::Asinh(e) => Ok(eval_float(e, env)?.asinh()),
        ExprKind::Acosh(e) => Ok(eval_float(e, env)?.acosh()),
        ExprKind::Atanh(e) => Ok(eval_float(e, env)?.atanh()),
        ExprKind::Acoth(e) => Ok(acoth(eval_float(e, env)?)),
        ExprKind::Asech(e) => Ok(asech(eval_float(e, env)?)),
        ExprKind::Acsch(e) => Ok(acsch(eval_float(e, env)?)),
        ExprKind::Exp(e) => Ok(eval_float(e, env)?.exp()),
        ExprKind::Ln(e) => {
            let v = eval_float(e, env)?;
            if v <= F::zero() {
                Err(EvalError::Undefined("ln of non-positive"))
            } else {
                Ok(v.ln())
            }
        }
    }
}

fn frac_pi_2<F: Float + FromPrimitive>() -> F {
    F::from_f64(std::f64::consts::FRAC_PI_2).expect("pi/2")
}

fn acoth<F: Float + FromPrimitive>(x: F) -> F {
    F::from(0.5).unwrap() * ((x + F::one()) / (x - F::one())).ln()
}

fn asech<F: Float + FromPrimitive>(x: F) -> F {
    (F::one() / x).acosh()
}

fn acsch<F: Float + FromPrimitive>(x: F) -> F {
    (F::one() / x).asinh()
}
