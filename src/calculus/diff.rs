#[path = "trig_diff.rs"]
mod trig_diff;

use crate::expr::{const_, pow, Expr, ExprKind};
use crate::rational::Rational;
use crate::symbol::Symbol;

pub fn diff(expr: Expr, var: Symbol) -> Expr {
    let raw = diff_without_simplify(expr, var);
    #[cfg(feature = "simplify")]
    {
        raw.simplify()
    }
    #[cfg(not(feature = "simplify"))]
    {
        raw
    }
}

/// Symbolic derivative without calling [`Expr::simplify`] (no-op when the `simplify` feature is off).
///
/// Prefer this when the simplified form is expensive to compute but a raw
/// derivative tree is enough (e.g. numeric spot-checks via [`Expr::eval_f64`]).
pub fn diff_without_simplify(expr: Expr, var: Symbol) -> Expr {
    diff_expr(&expr, var)
}

/// Differentiate without simplify (internal building block).
pub(crate) fn diff_expr(e: &Expr, var: Symbol) -> Expr {
    diff_kind(e.kind(), var)
}

fn diff_kind(kind: &ExprKind, var: Symbol) -> Expr {
    match kind {
        ExprKind::Const(_) => const_(Rational::zero()),
        ExprKind::Var(s) => {
            if *s == var {
                const_(Rational::one())
            } else {
                const_(Rational::zero())
            }
        }
        ExprKind::Add(a, b) => diff_expr(a, var) + diff_expr(b, var),
        ExprKind::Sub(a, b) => diff_expr(a, var) - diff_expr(b, var),
        ExprKind::Neg(e) => -diff_expr(e, var),
        ExprKind::Mul(f, g) => {
            let fp = diff_expr(f, var);
            let gp = diff_expr(g, var);
            fp * g.clone() + f.clone() * gp
        }
        ExprKind::Div(f, g) => {
            let fp = diff_expr(f, var);
            let gp = diff_expr(g, var);
            let g_c = g.clone();
            (fp * g_c.clone() - f.clone() * gp) / pow(g_c, const_(Rational::from(2)))
        }
        ExprKind::Pow(base, exp) => diff_pow(base, exp, var),
        ExprKind::Sin(e) => diff_expr(e, var) * crate::expr::cos(e.clone()),
        ExprKind::Cos(e) => -diff_expr(e, var) * crate::expr::sin(e.clone()),
        ExprKind::Tan(e) => {
            diff_expr(e, var) / pow(crate::expr::cos(e.clone()), const_(Rational::from(2)))
        }
        ExprKind::Exp(e) => diff_expr(e, var) * crate::expr::exp(e.clone()),
        ExprKind::Ln(e) => diff_expr(e, var) / e.clone(),
        ExprKind::Atan(e) => {
            diff_expr(e, var) / (const_(Rational::one()) + pow(e.clone(), const_(Rational::from(2))))
        }
        k => trig_diff::diff_trig_kind(k, var)
            .expect("extended trigonometric kinds are handled in trig_diff"),
    }
}

fn diff_pow(base: &Expr, exp: &Expr, var: Symbol) -> Expr {
    if let ExprKind::Var(s) = base.kind() {
        if s.name() == "e" {
            if is_var_expr(exp, var) {
                return diff_expr(&crate::expr::exp(exp.clone()), var);
            }
        }
    }
    if let ExprKind::Const(n) = exp.kind() {
        if let Some(n) = n.try_as_rational() {
            if let Some(k) = n.as_integer() {
            if k == 0 {
                return const_(Rational::zero());
            }
            if k == 1 {
                return diff_expr(base, var);
            }
            return const_(Rational::from(k))
                * pow(base.clone(), const_(Rational::from(k - 1)))
                * diff_expr(base, var);
            }
        }
    }
    let up = diff_expr(base, var);
    let vp = diff_expr(exp, var);
    let term1 = exp.clone() * pow(base.clone(), exp.clone() - const_(Rational::one())) * up;
    let u_ln = base.clone();
    let term2 = pow(base.clone(), exp.clone()) * crate::expr::ln(u_ln) * vp;
    term1 + term2
}

fn is_var_expr(e: &Expr, var: Symbol) -> bool {
    matches!(e.kind(), ExprKind::Var(s) if *s == var)
}
