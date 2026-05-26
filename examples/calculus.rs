use simsym::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let x = symbol("x");
    let f = expr!(e ^ x * (x ^ 3 + 2 * x)).simplify();
    let df = f.clone().diff(x).simplify();
    println!("f  = {f}");
    println!("f' = {df}");

    let antideriv = f.clone().integrate(x)?.simplify();
    println!("∫f dx = {antideriv}");

    let val = f.eval_f64(&[(x, 2.0)])?;
    println!("f(2) = {val}");
    Ok(())
}
