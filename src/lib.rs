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
    diff_without_simplify, gradient, hessian, integrate_definite, integrate_numeric,
    DefiniteIntegralError, IntegrateError, NumericOptions,
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
pub fn atan(e: impl Into<Expr>) -> Expr {
    expr::atan(e.into())
}

pub mod prelude {
    pub use crate::{
        atan, cos, diff_without_simplify, exp, expr, gradient, hessian, integrate_definite,
        integrate_numeric, ln, rational, rational_from_i32, sin, symbol, tan,
        DefiniteIntegralError, EvalError, Expr, IntegrateError, NumericOptions, Rational, Symbol,
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
    fn integrate_tan_and_ln() {
        let x = symbol("x");
        let tan_int = crate::expr::tan(x.to_expr()).integrate(x).unwrap().simplify();
        assert_eq!(
            tan_int.to_string(),
            (-crate::expr::ln(crate::expr::cos(x.to_expr()))).simplify().to_string()
        );
        let ln_int = crate::expr::ln(x.to_expr()).integrate(x).unwrap().simplify();
        let x_expr = x.to_expr();
        let expected_ln =
            x_expr.clone() * crate::expr::ln(x_expr.clone()) - x_expr;
        assert_eq!(ln_int.to_string(), expected_ln.simplify().to_string());
    }

    #[test]
    fn integrate_sin_squared() {
        let x = symbol("x");
        let f = sin(x).pow(2);
        let F = f.clone().integrate(x).unwrap().simplify();
        let recovered = F.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.3)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.3)]).unwrap();
        assert!((v - rv).abs() < 1e-6);
    }

    #[test]
    fn integrate_sin_cos_product() {
        let x = symbol("x");
        let f = sin(x) * cos(x);
        let F = f.clone().integrate(x).unwrap().simplify();
        let recovered = F.diff(x).simplify();
        let v = f.eval_f64(&[(x, 0.4)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 0.4)]).unwrap();
        assert!((v - rv).abs() < 1e-5);
    }

    #[test]
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
    fn integrate_x_exp_by_parts() {
        let x = symbol("x");
        let f = x.to_expr() * exp(x);
        let F = f.clone().integrate(x).unwrap().simplify();
        let recovered = F.diff(x).simplify();
        let v = f.eval_f64(&[(x, 1.0)]).unwrap();
        let rv = recovered.eval_f64(&[(x, 1.0)]).unwrap();
        assert!((v - rv).abs() < 1e-6);
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
