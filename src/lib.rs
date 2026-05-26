//! simsym — a simple symbolic computation library.
//!
//! ```rust
//! use simsym::prelude::*;
//!
//! let x = symbol("x");
//! let f = x.pow(2) + rational(2, 1) * x;
//! let df = f.clone().diff(x);
//! let val = f.eval(&[(x, rational(1, 2))]).unwrap();
//! assert_eq!(val, rational(5, 4));
//! ```

pub mod calculus;
pub mod display;
pub mod eval;
pub mod expr;
mod ops_ext;
pub mod poly;
pub mod rational;
pub mod simplify;
pub mod symbol;

#[cfg(feature = "bigint")]
pub mod rational_big;

pub use calculus::{
    gradient, hessian, integrate_definite, integrate_numeric, DefiniteIntegralError,
    IntegrateError, NumericOptions,
};
pub use eval::EvalError;
pub use expr::Expr;
pub use rational::{rational, rational_from_i32, Rational};
pub use symbol::{symbol, Symbol};

pub use simsym_macros::expr;

/// Elementary function constructors.
pub fn sin(e: impl Into<Expr>) -> Expr {
    expr::sin(e.into())
}
pub fn cos(e: impl Into<Expr>) -> Expr {
    expr::cos(e.into())
}
pub fn tan(e: impl Into<Expr>) -> Expr {
    expr::tan(e.into())
}
pub fn exp(e: impl Into<Expr>) -> Expr {
    expr::exp(e.into())
}
pub fn ln(e: impl Into<Expr>) -> Expr {
    expr::ln(e.into())
}

pub mod prelude {
    pub use crate::{
        cos, exp, expr, gradient, hessian, integrate_definite, integrate_numeric, ln, rational,
        rational_from_i32, sin, symbol, tan, DefiniteIntegralError, EvalError, Expr,
        IntegrateError, NumericOptions, Rational, Symbol,
    };
}

#[cfg(test)]
mod tests {
    use super::prelude::*;

    #[test]
    fn diff_power() {
        let x = symbol("x");
        let f = x.pow(2);
        let df = f.diff(x).simplify();
        let expected = rational(2, 1) * x.to_expr();
        assert_eq!(df.simplify().to_string(), expected.simplify().to_string());
    }

    #[test]
    fn eval_exact() {
        let x = symbol("x");
        let f = x.pow(2) + rational(2, 1) * x;
        let v = f.eval(&[(x, rational(1, 2))]).unwrap();
        assert_eq!(v, rational(5, 4));
    }

    #[test]
    fn integrate_polynomial() {
        let x = symbol("x");
        let f = x.pow(2);
        let antideriv = f.integrate(x).unwrap().simplify();
        let expected = x.pow(3) / rational(3, 1);
        assert_eq!(antideriv.to_string(), expected.simplify().to_string());
    }

    #[test]
    fn integrate_exp_manual() {
        let x = symbol("x");
        let f = exp(x) * (x.pow(3) + rational(2, 1) * x);
        assert!(f.integrate(x).is_ok());
    }

    #[test]
    fn diff_exp_times_polynomial() {
        let x = symbol("x");
        let f = exp(x) * (x.pow(3) + rational(2, 1) * x);
        let df = f.clone().diff(x).simplify();
        let x0 = rational(1, 2);
        let f_val = f.clone().eval_f64(&[(x, x0.to_f64())]).unwrap();
        let df_val = df.eval_f64(&[(x, x0.to_f64())]).unwrap();
        let h = 1e-6;
        let numeric = (f.clone().eval_f64(&[(x, x0.to_f64() + h)]).unwrap()
            - f.eval_f64(&[(x, x0.to_f64() - h)]).unwrap())
            / (2.0 * h);
        assert!((df_val - numeric).abs() < 1e-4, "f'={df_val} numeric={numeric}");
        assert!((df_val - f_val).abs() > 0.0);
    }

    #[test]
    fn integrate_exp_times_polynomial() {
        let x = symbol("x");
        let f = exp(x) * (x.pow(3) + rational(2, 1) * x);
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        let recovered = antideriv.diff(x).simplify();
        let x0 = rational(1, 3);
        let env = [(x, x0.to_f64())];
        let f_val = f.eval_f64(&env).unwrap();
        let recovered_val = recovered.eval_f64(&env).unwrap();
        assert!(
            (f_val - recovered_val).abs() < 1e-6,
            "d/dx(F) should match f: {recovered_val} vs {f_val}"
        );
    }

    #[test]
    fn polynomial_normal_form_orders_terms() {
        let x = symbol("x");
        let messy = x.pow(3) - rational(3, 1) * x.pow(2) + rational(6, 1) * x + (-rational(6, 1))
            + rational(2, 1) * x
            + (-rational(2, 1));
        let neat =
            x.pow(3) - rational(3, 1) * x.pow(2) + rational(8, 1) * x + (-rational(8, 1));
        assert_eq!(messy.simplify().to_string(), neat.simplify().to_string());
    }

    #[test]
    fn integrate_cancels_coefficients_in_quotient() {
        let x = symbol("x");
        let f = x.pow(3) + rational(2, 1) * x;
        let antideriv = f.integrate(x).unwrap().simplify();
        assert_eq!(
            antideriv.to_string(),
            (x.pow(4) / rational(4, 1) + x.pow(2)).simplify().to_string()
        );
    }

    #[test]
    fn sin_diff_and_integrate() {
        let x = symbol("x");
        let f = sin(x);
        let df = f.clone().diff(x).simplify();
        assert_eq!(df.to_string(), cos(x).simplify().to_string());
        let antideriv = f.integrate(x).unwrap().simplify();
        assert_eq!(antideriv.to_string(), (-cos(x)).simplify().to_string());
    }

    #[test]
    fn gradient_two_vars() {
        let x = symbol("x");
        let y = symbol("y");
        let f = x * y;
        let g = f.gradient(&[x, y]);
        assert_eq!(g[0].clone().simplify().to_string(), y.to_string());
        assert_eq!(g[1].clone().simplify().to_string(), x.to_string());
    }

    #[test]
    fn hessian_bilinear() {
        let x = symbol("x");
        let y = symbol("y");
        let f = x * y;
        let h = f.hessian(&[x, y]);
        assert_eq!(h[0][0].clone().simplify().to_string(), "0");
        assert_eq!(h[0][1].clone().simplify().to_string(), "1");
        assert_eq!(h[1][0].clone().simplify().to_string(), "1");
        assert_eq!(h[1][1].clone().simplify().to_string(), "0");
    }

    #[test]
    fn numeric_integral_sin() {
        let x = symbol("x");
        let f = sin(x);
        let pi = std::f64::consts::PI;
        let val = integrate_numeric(&f, x, 0.0, pi, &[], NumericOptions::default()).unwrap();
        assert!((val - 2.0).abs() < 1e-6);
    }

}
