use crate::calculus::integrate::{integrate, IntegrateError};
use crate::eval::EvalError;
use crate::expr::Expr;
use crate::rational::Rational;
use crate::symbol::Symbol;

#[derive(Debug, Clone, thiserror::Error)]
pub enum DefiniteIntegralError {
    #[error(transparent)]
    Integrate(#[from] IntegrateError),
    #[error(transparent)]
    Eval(#[from] EvalError),
    #[error("numeric integration failed: {0}")]
    Numeric(String),
}

#[derive(Debug, Clone, Copy)]
pub struct NumericOptions {
    pub tol: f64,
    pub max_depth: u32,
}

impl Default for NumericOptions {
    fn default() -> Self {
        Self {
            tol: 1e-10,
            max_depth: 20,
        }
    }
}

pub fn integrate_numeric(
    expr: &Expr,
    var: Symbol,
    a: f64,
    b: f64,
    env: &[(Symbol, f64)],
    opts: NumericOptions,
) -> Result<f64, DefiniteIntegralError> {
    adaptive_simpson(expr, var, a, b, env, opts.tol, opts.max_depth)
        .map_err(DefiniteIntegralError::Numeric)
}

pub fn integrate_definite(
    expr: Expr,
    var: Symbol,
    a: Rational,
    b: Rational,
    env: &[(Symbol, Rational)],
) -> Result<Rational, DefiniteIntegralError> {
    match integrate(expr.clone(), var) {
        Ok(antiderivative) => {
            let fa = substitute_and_eval(antiderivative.clone(), var, a, env)?;
            let fb = substitute_and_eval(antiderivative, var, b, env)?;
            Ok(fb - fa)
        }
        Err(IntegrateError::NoRule) => {
            let f64_env: Vec<(Symbol, f64)> = env.iter().map(|(s, v)| (*s, v.to_f64())).collect();
            let val = integrate_numeric(
                &expr,
                var,
                a.to_f64(),
                b.to_f64(),
                &f64_env,
                NumericOptions::default(),
            )?;
            float_to_rational_approx(val)
                .ok_or_else(|| DefiniteIntegralError::Numeric("could not convert numeric result to rational".into()))
        }
        Err(e) => Err(e.into()),
    }
}

fn float_to_rational_approx(v: f64) -> Option<Rational> {
    if v.fract() == 0.0 && v.is_finite() && v >= i64::MIN as f64 && v <= i64::MAX as f64 {
        return Some(Rational::from_integer(v as i64));
    }
    const DEN: i64 = 1_000_000_000;
    let n = (v * DEN as f64).round() as i64;
    Some(Rational::new(n, DEN))
}

fn substitute_and_eval(
    expr: Expr,
    var: Symbol,
    at: Rational,
    env: &[(Symbol, Rational)],
) -> Result<Rational, EvalError> {
    let mut full: Vec<(Symbol, Rational)> = env.to_vec();
    if let Some(slot) = full.iter_mut().find(|(s, _)| *s == var) {
        slot.1 = at;
    } else {
        full.push((var, at));
    }
    expr.eval(&full)
}

fn adaptive_simpson(
    expr: &Expr,
    var: Symbol,
    a: f64,
    b: f64,
    env: &[(Symbol, f64)],
    tol: f64,
    depth: u32,
) -> Result<f64, String> {
    if depth == 0 {
        return Err("max recursion depth".into());
    }
    let c = (a + b) / 2.0;
    let fa = eval_at(expr, var, a, env).map_err(|e| e.to_string())?;
    let fb = eval_at(expr, var, b, env).map_err(|e| e.to_string())?;
    let fc = eval_at(expr, var, c, env).map_err(|e| e.to_string())?;
    let h = b - a;
    let whole = h / 6.0 * (fa + 4.0 * fc + fb);

    let d = (a + c) / 2.0;
    let e_mid = (c + b) / 2.0;
    let fd = eval_at(expr, var, d, env).map_err(|e| e.to_string())?;
    let fe = eval_at(expr, var, e_mid, env).map_err(|e| e.to_string())?;
    let left = h / 12.0 * (fa + 4.0 * fd + fc);
    let right = h / 12.0 * (fc + 4.0 * fe + fb);
    let delta = left + right - whole;

    if delta.abs() <= 15.0 * tol {
        return Ok(left + right + delta / 15.0);
    }
    let left_val = adaptive_simpson(expr, var, a, c, env, tol / 2.0, depth - 1)?;
    let right_val = adaptive_simpson(expr, var, c, b, env, tol / 2.0, depth - 1)?;
    Ok(left_val + right_val)
}

fn eval_at(expr: &Expr, var: Symbol, t: f64, env: &[(Symbol, f64)]) -> Result<f64, EvalError> {
    let mut full: Vec<(Symbol, f64)> = env.to_vec();
    if let Some(slot) = full.iter_mut().find(|(s, _)| *s == var) {
        slot.1 = t;
    } else {
        full.push((var, t));
    }
    expr.clone().eval_f64(&full)
}
