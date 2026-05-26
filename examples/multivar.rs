use simsym::prelude::*;

fn main() {
    let x = symbol("x");
    let y = symbol("y");
    let f = expr!(x^2*y + sin(y));
    let grad = f.clone().gradient(&[x, y]);
    println!("f = {f}");
    for (i, g) in grad.iter().enumerate() {
        println!("∂f/∂var{i} = {}", g.clone().simplify());
    }
}
