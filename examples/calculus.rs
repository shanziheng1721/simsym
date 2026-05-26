use simsym::prelude::*;

fn check_integral(f: Expr, x: Symbol) -> Result<(), Box<dyn std::error::Error>> {
    let f = f.simplify();
    println!("=== {f} ===");
    match f.clone().integrate(x) {
        Ok(fx) => {
            let fx = fx.simplify();
            println!("∫f dx = {fx}");
            let check = fx.diff_without_simplify(x);
            if let Ok(v) = f.eval_f64(&[(x, 0.5)]) {
                if let Ok(rv) = check.eval_f64(&[(x, 0.5)]) {
                    println!("check f(0.5)={v}  d/dx(F)(0.5)={rv}");
                }
            }
        }
        Err(e) => println!("integrate: {e}"),
    }
    println!();
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let x = symbol("x");

    check_integral(expr!(e ^ x * (x ^ 3 + 2 * x)), x)?;
    check_integral(exp(x) * sin(x), x)?;
    check_integral(exp(expr!(2 * x + 1)) * sin(expr!(3 * x + 1 / 4)), x)?;
    check_integral(sin(x).pow(2) * cos(x).pow(2), x)?;
    check_integral(sin(x).pow(3) * cos(x).pow(2), x)?;
    check_integral(cos(x).pow(-4), x)?;
    check_integral(atan(x), x)?;
    check_integral(expr!(1 / (x ^ 2 + 1)), x)?;

    Ok(())
}
