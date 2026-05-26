use crate::expr::{Expr, ExprKind};
use crate::poly::try_polynomial_normal_form;
use crate::rational::Rational;
use std::fmt;

pub fn fmt_expr(e: &Expr, f: &mut fmt::Formatter<'_>, parent_prec: u8) -> fmt::Result {
    if let Some(poly) = try_polynomial_normal_form(e.clone()) {
        return fmt_polynomial(&poly, f, parent_prec);
    }
    fmt_expr_inner(e, f, parent_prec)
}

fn fmt_polynomial(e: &Expr, f: &mut fmt::Formatter<'_>, parent_prec: u8) -> fmt::Result {
    let prec = 1u8;
    let needs_paren = prec < parent_prec;
    if needs_paren {
        write!(f, "(")?;
    }
    fmt_polynomial_descending(e, f)?;
    if needs_paren {
        write!(f, ")")?;
    }
    Ok(())
}

/// Format an already-normalized polynomial (descending powers, merged coefficients).
fn fmt_polynomial_descending(e: &Expr, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match e.kind() {
        ExprKind::Const(c) => write!(f, "{c}")?,
        ExprKind::Add(a, b) => {
            fmt_polynomial_descending(a, f)?;
            write!(f, " + ")?;
            fmt_polynomial_term(b, f, true)?;
        }
        ExprKind::Sub(a, b) => {
            fmt_polynomial_descending(a, f)?;
            write!(f, " - ")?;
            fmt_polynomial_term(b, f, false)?;
        }
        _ => fmt_expr_inner(e, f, 0)?,
    }
    Ok(())
}

fn fmt_polynomial_term(e: &Expr, f: &mut fmt::Formatter<'_>, force_add_sign: bool) -> fmt::Result {
    match e.kind() {
        ExprKind::Sub(a, b) => {
            fmt_polynomial_descending(a, f)?;
            write!(f, " - ")?;
            fmt_polynomial_term(b, f, false)
        }
        ExprKind::Mul(l, r) => {
            if let ExprKind::Const(c) = l.kind() {
                if c.is_one() {
                    return fmt_expr_inner(r, f, 10);
                }
                if *c == -Rational::one() {
                    write!(f, "-")?;
                    return fmt_expr_inner(r, f, 10);
                }
            }
            fmt_expr_inner(e, f, 10)
        }
        ExprKind::Neg(inner) => {
            write!(f, "-")?;
            fmt_expr_inner(inner, f, 10)
        }
        _ => {
            if force_add_sign {
                fmt_expr_inner(e, f, 10)
            } else {
                fmt_polynomial_descending(e, f)
            }
        }
    }
}

fn fmt_expr_inner(e: &Expr, f: &mut fmt::Formatter<'_>, parent_prec: u8) -> fmt::Result {
    let prec = precedence(e.kind());
    let needs_paren = prec < parent_prec;
    if needs_paren {
        write!(f, "(")?;
    }
    match e.kind() {
        ExprKind::Const(c) => write!(f, "{c}")?,
        ExprKind::Var(s) => write!(f, "{s}")?,
        ExprKind::Add(a, b) => {
            fmt_expr_inner(a, f, prec)?;
            write!(f, " + ")?;
            fmt_expr_inner(b, f, prec + 1)?;
        }
        ExprKind::Sub(a, b) => {
            fmt_expr_inner(a, f, prec)?;
            write!(f, " - ")?;
            fmt_expr_inner(b, f, prec + 1)?;
        }
        ExprKind::Mul(a, b) => {
            fmt_expr_inner(a, f, prec)?;
            write!(f, "*")?;
            fmt_expr_inner(b, f, prec)?;
        }
        ExprKind::Div(a, b) => {
            fmt_expr_inner(a, f, prec)?;
            write!(f, "/")?;
            fmt_expr_inner(b, f, prec)?;
        }
        ExprKind::Neg(inner) => {
            write!(f, "-")?;
            fmt_expr_inner(inner, f, prec)?;
        }
        ExprKind::Pow(base, exp) => {
            fmt_expr_inner(base, f, prec + 1)?;
            write!(f, "^")?;
            fmt_expr_inner(exp, f, prec + 1)?;
        }
        ExprKind::Sin(e) => write!(f, "sin({e})")?,
        ExprKind::Cos(e) => write!(f, "cos({e})")?,
        ExprKind::Tan(e) => write!(f, "tan({e})")?,
        ExprKind::Exp(inner) => {
            if matches!(inner.kind(), ExprKind::Var(_)) {
                write!(f, "e^")?;
                fmt_expr_inner(inner, f, prec + 1)?;
            } else {
                write!(f, "exp(")?;
                fmt_expr_inner(inner, f, 0)?;
                write!(f, ")")?;
            }
        }
        ExprKind::Ln(e) => write!(f, "ln({e})")?,
    }
    if needs_paren {
        write!(f, ")")?;
    }
    Ok(())
}

fn precedence(kind: &ExprKind) -> u8 {
    match kind {
        ExprKind::Add(..) | ExprKind::Sub(..) => 1,
        ExprKind::Neg(_) => 2,
        ExprKind::Mul(..) | ExprKind::Div(..) => 3,
        ExprKind::Pow(_, _) => 4,
        ExprKind::Exp(_) => 4,
        _ => 10,
    }
}
