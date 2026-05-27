use crate::calculus::diff::diff;
use crate::expr::Expr;
use crate::symbol::Symbol;

pub fn gradient(expr: Expr, vars: &[Symbol]) -> Vec<Expr> {
    vars.iter()
        .map(|&v| diff(expr.clone(), v))
        .collect()
}

pub fn hessian(expr: Expr, vars: &[Symbol]) -> Vec<Vec<Expr>> {
    let grad = gradient(expr, vars);
    grad.into_iter()
        .map(|g| vars.iter().map(|&v| diff(g.clone(), v)).collect())
        .collect()
}
