use simsym::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let x = symbol("x");
    let f = sin(x.pow(2)); // no elementary antiderivative
    let a = 0.0;
    let b = 1.0;
    match f.clone().integrate(x) {
        Ok(fx) => println!("symbolic: {fx}"),
        Err(e) => println!("symbolic integrate failed: {e}"),
    }
    let num = integrate_numeric(&f, x, a, b, &[], NumericOptions::default())?;
    println!("∫_0^1 sin(x^2) dx ≈ {num}");
    Ok(())
}
