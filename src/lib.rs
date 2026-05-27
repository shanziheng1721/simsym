//! simsym — a simple symbolic computation library.
//!
//! ## Cargo features
//!
//! - **`simplify`** (default): algebraic simplification via [`Expr::simplify`]
//! - **`diff`** (default): symbolic differentiation, [`Expr::gradient`], [`Expr::hessian`]
//! - **`integrate`** (default, implies `diff`): symbolic integration and [`Expr::integrate_definite`]
//!
//! Build a minimal core (expression AST, rationals, evaluation only):
//!
//! ```bash
//! cargo build --no-default-features
//! ```
//!
//! ```rust
//! use simsym::prelude::*;
//! let x = symbol("x");
//! let f = x.pow(2) + rational(2, 1) * x;
//! let val = f.eval(&[(x, rational(1, 2))]).unwrap();
//! assert_eq!(val, rational(5, 4));
//! ```

pub mod display;
pub mod eval;
pub mod expr;
mod ops_ext;
pub mod poly;
pub mod rational;
pub mod symbol;

#[cfg(feature = "simplify")]
pub mod simplify;
#[cfg(not(feature = "simplify"))]
pub mod simplify_nop;
#[cfg(not(feature = "simplify"))]
pub use simplify_nop as simplify;

#[cfg(any(feature = "diff", feature = "integrate"))]
pub mod calculus;

#[cfg(feature = "bigint")]
pub mod rational_big;

pub use eval::EvalError;
pub use expr::Expr;
pub use rational::{rational, rational_from_i32, Rational};
pub use symbol::{symbol, Symbol};

#[cfg(any(feature = "diff", feature = "integrate"))]
pub use calculus::{integrate_numeric, DefiniteIntegralError, NumericOptions};

#[cfg(feature = "integrate")]
pub use calculus::{integrate, integrate_definite, IntegrateError};

#[cfg(feature = "diff")]
pub use calculus::{diff, diff_without_simplify, gradient, hessian};

pub use simsym_macros::expr;

macro_rules! unary_fn {
    ($($name:ident),* $(,)?) => {
        $(pub fn $name(e: impl Into<Expr>) -> Expr {
            expr::$name(e.into())
        })*
    };
}

unary_fn!(
    sin, cos, tan, cot, sec, csc, asin, acos, atan, acot, asec, acsc, sinh, cosh, tanh, coth,
    sech, csch, asinh, acosh, atanh, acoth, asech, acsch, exp, ln,
);

pub mod prelude {
    pub use crate::{
        acos, acosh, acot, acoth, acsc, acsch, asec, asech, asin, asinh, atan, atanh, cos, cosh,
        cot, coth, csc, csch, exp, expr, ln, rational, rational_from_i32, sec, sech, sin, sinh,
        symbol, tan, tanh, EvalError, Expr, Rational, Symbol,
    };
    #[cfg(any(feature = "diff", feature = "integrate"))]
    pub use crate::{DefiniteIntegralError, NumericOptions, integrate_numeric};
    #[cfg(feature = "integrate")]
    pub use crate::{integrate, integrate_definite, IntegrateError};
    #[cfg(feature = "diff")]
    pub use crate::{diff_without_simplify, gradient, hessian};
}

#[cfg(test)]
mod tests {
    use super::prelude::*;

    #[test]
    fn eval_exact() {
        let x = symbol("x");
        let f = x.pow(2) + rational(2, 1) * x;
        let v = f.eval(&[(x, rational(1, 2))]).unwrap();
        assert_eq!(v, rational(5, 4));
    }

    #[test]
    fn extended_trig_eval_f64() {
        let x = symbol("x");
        let env = [(x, 0.5)];
        assert!((cot(x).eval_f64(&env).unwrap() - 1.0 / 0.5f64.tan()).abs() < 1e-10);
        assert!((sinh(x).eval_f64(&env).unwrap() - 0.5f64.sinh()).abs() < 1e-10);
        assert_eq!(cot(x).to_string(), "cot(x)");
        assert_eq!(sinh(x).to_string(), "sinh(x)");
    }

    #[test]
    #[cfg(all(feature = "diff", feature = "simplify"))]
    fn diff_power() {
        let x = symbol("x");
        let f = x.pow(2);
        let df = f.diff(x).simplify();
        let expected = rational(2, 1) * x.to_expr();
        assert_eq!(df.simplify().to_string(), expected.simplify().to_string());
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff"))]
    fn diff_without_simplify_matches_numeric() {
        let x = symbol("x");
        let f = sin(x).pow(3) * cos(x).pow(2);
        let fx = f.clone().integrate(x).unwrap();
        let df = fx.diff_without_simplify(x);
        let v = f.eval_f64(&[(x, 0.5)]).unwrap();
        let rv = df.eval_f64(&[(x, 0.5)]).unwrap();
        assert!((v - rv).abs() < 1e-5);
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "simplify"))]
    fn integrate_polynomial() {
        let x = symbol("x");
        let f = x.pow(2);
        let antideriv = f.integrate(x).unwrap().simplify();
        let expected = x.pow(3) / rational(3, 1);
        assert_eq!(antideriv.to_string(), expected.simplify().to_string());
    }

    #[test]
    #[cfg(feature = "integrate")]
    fn integrate_exp_manual() {
        let x = symbol("x");
        let f = exp(x) * (x.pow(3) + rational(2, 1) * x);
        assert!(f.integrate(x).is_ok());
    }

    #[test]
    #[cfg(all(feature = "diff", feature = "simplify"))]
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
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
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
    #[cfg(feature = "simplify")]
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
    #[cfg(all(feature = "integrate", feature = "simplify"))]
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
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn sin_diff_and_integrate() {
        let x = symbol("x");
        let f = sin(x);
        let df = f.clone().diff(x).simplify();
        assert_eq!(df.to_string(), cos(x).simplify().to_string());
        let antideriv = f.integrate(x).unwrap().simplify();
        assert_eq!(antideriv.to_string(), (-cos(x)).simplify().to_string());
    }

    #[test]
    #[cfg(all(feature = "diff", feature = "simplify"))]
    fn gradient_two_vars() {
        let x = symbol("x");
        let y = symbol("y");
        let f = x * y;
        let g = f.gradient(&[x, y]);
        assert_eq!(g[0].clone().simplify().to_string(), y.to_string());
        assert_eq!(g[1].clone().simplify().to_string(), x.to_string());
    }

    #[test]
    #[cfg(all(feature = "diff", feature = "simplify"))]
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
    #[cfg(all(feature = "integrate", feature = "simplify"))]
    fn integrate_tan_and_ln() {
        let x = symbol("x");
        let tan_int = crate::expr::tan(x.to_expr()).integrate(x).unwrap().simplify();
        assert_eq!(
            tan_int.to_string(),
            (-crate::expr::ln(crate::expr::cos(x.to_expr()))).simplify().to_string()
        );
        let ln_int = crate::expr::ln(x.to_expr()).integrate(x).unwrap().simplify();
        let x_expr = x.to_expr();
        let expected_ln = x_expr.clone() * crate::expr::ln(x_expr.clone()) - x_expr;
        assert_eq!(ln_int.to_string(), expected_ln.simplify().to_string());
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_sin_squared() {
        let x = symbol("x");
        let f = sin(x).pow(2);
        let f_int = f.clone().integrate(x).unwrap().simplify();
        let recovered = f_int.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.3)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.3)]).unwrap();
        assert!((v - rv).abs() < 1e-6);
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_sin_cos_product() {
        let x = symbol("x");
        let f = sin(x) * cos(x);
        let f_int = f.clone().integrate(x).unwrap().simplify();
        let recovered = f_int.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.4)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.4)]).unwrap();
        assert!((v - rv).abs() < 1e-5);
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "simplify"))]
    fn integrate_partial_fractions() {
        use crate::expr::const_;
        let x = symbol("x");
        let f = const_(Rational::one())
            / (x.to_expr() - const_(rational(1, 1)))
            / (x.to_expr() - const_(rational(2, 1)));
        let antideriv = f.integrate(x).unwrap().simplify();
        assert!(
            antideriv.to_string().contains("ln"),
            "expected ln terms, got {antideriv}"
        );
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_one_over_x_squared_plus_one() {
        use crate::expr::const_;
        let x = symbol("x");
        let f = const_(Rational::one()) / (x.to_expr().pow(2) + rational(1, 1));
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        assert!(
            antideriv.to_string().contains("atan"),
            "expected atan, got {antideriv}"
        );
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.5)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.5)]).unwrap();
        assert!((v - rv).abs() < 1e-6, "{rv} vs {v}");
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_rational_polynomial_quotient() {
        use crate::expr::const_;
        let x = symbol("x");
        let f = (x.to_expr().pow(2) + rational(1, 1))
            / (x.to_expr() - const_(rational(1, 1)));
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 2.0)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 2.0)]).unwrap();
        assert!((v - rv).abs() < 1e-5, "{rv} vs {v}");
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_sin_fourth_power() {
        let x = symbol("x");
        let f = sin(x).pow(4);
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.25)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.25)]).unwrap();
        assert!((v - rv).abs() < 1e-5, "{rv} vs {v}");
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_atan_x() {
        let x = symbol("x");
        let f = crate::expr::atan(x.to_expr());
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        assert!(antideriv.to_string().contains("atan"));
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 1.0)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 1.0)]).unwrap();
        assert!((v - rv).abs() < 1e-5);
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_exp_sin() {
        let x = symbol("x");
        let f = exp(x) * sin(x);
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.7)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.7)]).unwrap();
        assert!((v - rv).abs() < 1e-5, "{rv} vs {v}");
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_x_squared_ln() {
        let x = symbol("x");
        let f = x.to_expr().pow(2) * crate::expr::ln(x.to_expr());
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 2.0)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 2.0)]).unwrap();
        assert!((v - rv).abs() < 1e-5);
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_tan_cubed() {
        let x = symbol("x");
        let f = crate::expr::tan(x.to_expr()).pow(3);
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.3)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.3)]).unwrap();
        assert!((v - rv).abs() < 1e-5);
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_one_over_x_squared_minus_one() {
        use crate::expr::const_;
        let x = symbol("x");
        let f = const_(Rational::one()) / (x.to_expr().pow(2) - const_(rational(1, 1)));
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        assert!(antideriv.to_string().contains("ln"));
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 2.0)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 2.0)]).unwrap();
        assert!((v - rv).abs() < 1e-5);
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_exp_sin_with_phase() {
        let x = symbol("x");
        let f = exp(rational(2, 1) * x.to_expr() + rational(1, 1))
            * sin(rational(3, 1) * x.to_expr() + rational(1, 4));
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.2)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.2)]).unwrap();
        assert!((v - rv).abs() < 1e-5, "{rv} vs {v}");
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_sin_squared_cos_squared() {
        let x = symbol("x");
        let f = sin(x).pow(2) * cos(x).pow(2);
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.4)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.4)]).unwrap();
        assert!((v - rv).abs() < 1e-5);
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_sin_cubed_cos_squared() {
        let x = symbol("x");
        let f = sin(x).pow(3) * cos(x).pow(2);
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.35)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.35)]).unwrap();
        assert!((v - rv).abs() < 1e-5);
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_sec_fourth() {
        let x = symbol("x");
        let f = cos(x).pow(-4);
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.25)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.25)]).unwrap();
        assert!((v - rv).abs() < 1e-5, "{rv} vs {v}");
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_sin_cos_with_phase() {
        let x = symbol("x");
        let f = sin(x.to_expr() + rational(1, 4)) * cos(rational(2, 1) * x.to_expr());
        let antideriv = f.clone().integrate(x).unwrap().simplify();
        let recovered = antideriv.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.5)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.5)]).unwrap();
        assert!((v - rv).abs() < 1e-5);
    }

    #[test]
    #[cfg(all(feature = "integrate", feature = "diff", feature = "simplify"))]
    fn integrate_x_exp_by_parts() {
        let x = symbol("x");
        let f = x.to_expr() * exp(x);
        let f_int = f.clone().integrate(x).unwrap().simplify();
        let recovered = f_int.diff(x).simplify();
        let v = f.eval_f64(&[(x, 1.0)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 1.0)]).unwrap();
        assert!((v - rv).abs() < 1e-6);
    }

    #[test]
    #[cfg(feature = "integrate")]
    fn numeric_integral_sin() {
        let x = symbol("x");
        let f = sin(x);
        let pi = std::f64::consts::PI;
        let val = integrate_numeric(&f, x, 0.0, pi, &[], NumericOptions::default()).unwrap();
        assert!((val - 2.0).abs() < 1e-6);
    }
}
