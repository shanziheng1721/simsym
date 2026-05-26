use crate::expr::{const_, pow, Expr, ExprKind};
use crate::rational::Rational;
use crate::symbol::Symbol;
pub fn diff(expr: Expr, var: Symbol) -> Expr {
    let d = diff_kind(expr.kind(), var);
    d.simplify()
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
        ExprKind::Add(a, b) => diff(a.clone(), var) + diff(b.clone(), var),
        ExprKind::Sub(a, b) => diff(a.clone(), var) - diff(b.clone(), var),
        ExprKind::Neg(e) => -diff(e.clone(), var),
        ExprKind::Mul(f, g) => {
            let fp = diff(f.clone(), var);
            let gp = diff(g.clone(), var);
            fp.clone() * g.clone() + f.clone() * gp
        }
        ExprKind::Div(f, g) => {
            let fp = diff(f.clone(), var);
            let gp = diff(g.clone(), var);
            (fp * g.clone() - f.clone() * gp) / pow(g.clone(), const_(Rational::from(2)))
        }
        ExprKind::Pow(base, exp) => diff_pow(base.clone(), exp.clone(), var),
        ExprKind::Sin(e) => diff(e.clone(), var) * crate::expr::cos(e.clone()),
        ExprKind::Cos(e) => -diff(e.clone(), var) * crate::expr::sin(e.clone()),
        ExprKind::Tan(e) => {
            let u = e.clone();
            diff(u.clone(), var) / pow(crate::expr::cos(u), const_(Rational::from(2)))
        }
        ExprKind::Exp(e) => diff(e.clone(), var) * crate::expr::exp(e.clone()),
        ExprKind::Ln(e) => diff(e.clone(), var) / e.clone(),
    }
}

fn diff_pow(base: Expr, exp: Expr, var: Symbol) -> Expr {
    // Treat symbol `e` raised to x as exp(x) for differentiation
    if let ExprKind::Var(s) = base.kind() {
        if s.name() == "e" {
            if is_var_expr(&exp, var) {
                return crate::expr::exp(exp).diff(var);
            }
        }
    }
    if let ExprKind::Const(n) = exp.kind() {
        if let Some(k) = n.as_integer() {
            if k == 0 {
                return const_(Rational::zero());
            }
            if k == 1 {
                return diff(base.clone(), var);
            }
            return const_(Rational::from(k))
                * pow(base.clone(), const_(Rational::from(k - 1)))
                * diff(base, var);
        }
    }
    // general u^v: v*u^(v-1)*u' + u^v*ln(u)*v'
    let u = base;
    let v = exp;
    let up = diff(u.clone(), var);
    let vp = diff(v.clone(), var);
    let term1 = v.clone() * pow(u.clone(), v.clone() - const_(Rational::one())) * up;
    let term2 = pow(u.clone(), v) * crate::expr::ln(u) * vp;
    term1 + term2
}

fn is_var_expr(e: &Expr, var: Symbol) -> bool {
    matches!(e.kind(), ExprKind::Var(s) if *s == var)
}
