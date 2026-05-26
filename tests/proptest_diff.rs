use proptest::prelude::*;
use simsym::prelude::*;

proptest! {
    #[test]
    fn diff_polynomial_matches_numeric(
        a in -5i32..5,
        b in -5i32..5,
        x0 in -2.0f64..2.0,
        h in 1e-6f64..1e-3,
    ) {
        let x = symbol("x");
        let f = rational(i64::from(a), 1) * x.pow(2) + rational(i64::from(b), 1) * x;
        let df = f.clone().diff(x).simplify();
        let env = [(x, x0)];
        let f_plus = f.clone().eval_f64(&[(x, x0 + h)])?;
        let f_minus = f.eval_f64(&[(x, x0 - h)])?;
        let numeric = (f_plus - f_minus) / (2.0 * h);
        let symbolic = df.eval_f64(&env)?;
        prop_assert!((numeric - symbolic).abs() < 1e-4);
    }
}
