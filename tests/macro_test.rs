use simsym::prelude::*;

#[test]
fn expr_macro_euler_exp() {
    let x = symbol("x");
    let e = expr!(e^x);
    let expected = exp(x);
    assert_eq!(e.simplify().to_string(), expected.simplify().to_string());
}

#[test]
fn expr_macro_integrate_exp_poly() {
    let x = symbol("x");
    let f_macro = expr!(e^x * (x^3 + 2*x));
    let f_manual = exp(x) * (x.pow(3) + rational(2, 1) * x);
    assert_eq!(
        f_macro.clone().simplify().to_string(),
        f_manual.simplify().to_string(),
        "macro and manual forms should match"
    );
    assert!(f_macro.clone().integrate(x).is_ok(), "{:?}", f_macro.integrate(x));
}

#[test]
fn expr_macro_matches_operators() {
    let x = symbol("x");
    let a = expr!(x^2 + 2*x + 1);
    let b = x.pow(2) + rational(2, 1) * x + rational(1, 1);
    assert_eq!(a.simplify().to_string(), b.simplify().to_string());
}
