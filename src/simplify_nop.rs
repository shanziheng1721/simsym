//! Identity simplification when the `simplify` feature is disabled.

use crate::expr::Expr;

#[inline]
pub fn simplify(expr: Expr) -> Expr {
    expr
}
