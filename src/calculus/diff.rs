use crate::expr::{const_, pow, Expr, ExprKind};
use crate::rational::Rational;
use crate::symbol::Symbol;

pub fn diff(expr: Expr, var: Symbol) -> Expr {
    diff_without_simplify(expr, var).simplify()
}

/// Symbolic derivative without calling [`crate::simplify::simplify`].
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
            fp.clone() * g.clone() + f.clone() * gp
        }
        ExprKind::Div(f, g) => {
            let fp = diff_expr(f, var);
            let gp = diff_expr(g, var);
            (fp * g.clone() - f.clone() * gp) / pow(g.clone(), const_(Rational::from(2)))
        }
        ExprKind::Pow(base, exp) => diff_pow(base.clone(), exp.clone(), var),
        ExprKind::Sin(e) => diff_expr(e, var) * crate::expr::cos(e.clone()),
        ExprKind::Cos(e) => -diff_expr(e, var) * crate::expr::sin(e.clone()),
        ExprKind::Tan(e) => {
            let u = e.clone();
            diff_expr(&u, var) / pow(crate::expr::cos(u), const_(Rational::from(2)))
        }
        ExprKind::Exp(e) => diff_expr(e, var) * crate::expr::exp(e.clone()),
        ExprKind::Ln(e) => diff_expr(e, var) / e.clone(),
        ExprKind::Atan(e) => diff_expr(e, var) / (const_(Rational::one()) + pow(e.clone(), const_(Rational::from(2)))),
    }
}

fn diff_pow(base: Expr, exp: Expr, var: Symbol) -> Expr {
    if let ExprKind::Var(s) = base.kind() {
        if s.name() == "e" {
            if is_var_expr(&exp, var) {
                return diff_expr(&crate::expr::exp(exp), var);
            }
        }
    }
    if let ExprKind::Const(n) = exp.kind() {
        if let Some(k) = n.as_integer() {
            if k == 0 {
                return const_(Rational::zero());
            }
            if k == 1 {
                return diff_expr(&base, var);
            }
            return const_(Rational::from(k))
                * pow(base.clone(), const_(Rational::from(k - 1)))
                * diff_expr(&base, var);
        }
    }
    let u = base;
    let v = exp;
    let up = diff_expr(&u, var);
    let vp = diff_expr(&v, var);
    let term1 = v.clone() * pow(u.clone(), v.clone() - const_(Rational::one())) * up;
    let term2 = pow(u.clone(), v) * crate::expr::ln(u) * vp;
    term1 + term2
}

fn is_var_expr(e: &Expr, var: Symbol) -> bool {
    matches!(e.kind(), ExprKind::Var(s) if *s == var)
}
