//! Univariate polynomial normal form: combine like terms, order by descending degree.

use crate::expr::{add, const_, mul, neg, pow, sub, Expr, ExprKind};
use crate::rational::Rational;
use crate::symbol::Symbol;
use std::collections::BTreeMap;

/// If `expr` is a polynomial in a single variable, return it in standard form (highest degree first).
pub fn polynomial_normal_form(expr: Expr) -> Expr {
    try_polynomial_normal_form(expr.clone()).unwrap_or(expr)
}

pub fn try_polynomial_normal_form(expr: Expr) -> Option<Expr> {
    let mut summands = Vec::new();
    flatten_summands(&expr, Rational::one(), &mut summands);
    build_polynomial(&summands)
}

fn flatten_summands(expr: &Expr, sign: Rational, out: &mut Vec<Expr>) {
    if sign.is_zero() {
        return;
    }
    match expr.kind() {
        ExprKind::Add(a, b) => {
            flatten_summands(a, sign, out);
            flatten_summands(b, sign, out);
        }
        ExprKind::Sub(a, b) => {
            flatten_summands(a, sign, out);
            flatten_summands(b, -sign, out);
        }
        ExprKind::Neg(e) => flatten_summands(e, -sign, out),
        _ => out.push(scale_expr(sign, expr)),
    }
}

fn scale_expr(sign: Rational, e: &Expr) -> Expr {
    if sign.is_one() {
        e.clone()
    } else if sign == -Rational::one() {
        neg(e.clone())
    } else {
        mul(const_(sign), e.clone())
    }
}

fn build_polynomial(summands: &[Expr]) -> Option<Expr> {
    let mut by_degree: BTreeMap<i64, Rational> = BTreeMap::new();
    let mut var: Option<Symbol> = None;

    for term in summands {
        let (v_opt, degree, coeff) = parse_univariate_term(term)?;
        if let Some(v) = v_opt {
            if let Some(existing) = var {
                if existing != v {
                    return None;
                }
            } else {
                var = Some(v);
            }
        } else if degree != 0 {
            return None;
        }
        *by_degree.entry(degree).or_insert(Rational::zero()) += coeff;
    }

    let var = match var {
        Some(v) => v,
        None => {
            let c = by_degree.get(&0).copied().unwrap_or(Rational::zero());
            return Some(const_(c));
        }
    };
    let mut degrees: Vec<i64> = by_degree
        .iter()
        .filter(|(_, c)| !c.is_zero())
        .map(|(&d, _)| d)
        .collect();
    if degrees.is_empty() {
        return Some(const_(Rational::zero()));
    }
    degrees.sort_by(|a, b| b.cmp(a));

    let mut acc = term_for_degree(var, degrees[0], by_degree[&degrees[0]])?;
    for &d in &degrees[1..] {
        let c = by_degree[&d];
        let t = term_for_degree(var, d, if c.is_positive() { c } else { -c })?;
        acc = if c.is_positive() {
            add(acc, t)
        } else {
            sub(acc, t)
        };
    }
    Some(acc)
}

fn parse_univariate_term(e: &Expr) -> Option<(Option<Symbol>, i64, Rational)> {
    match e.kind() {
        ExprKind::Neg(inner) => {
            let (v, d, c) = parse_univariate_term(inner)?;
            Some((v, d, -c))
        }
        _ => {
            let (coeff, mono) = coeff_monomial(e)?;
            match mono.kind() {
                ExprKind::Const(_) => Some((None, 0, coeff)),
                ExprKind::Var(s) => Some((Some(*s), 1, coeff)),
                ExprKind::Pow(base, exp) => {
                    let s = match base.kind() {
                        ExprKind::Var(s) => *s,
                        _ => return None,
                    };
                    let n = as_const(exp)?.as_integer()?;
                    if n < 0 {
                        return None;
                    }
                    Some((Some(s), n, coeff))
                }
                _ => None,
            }
        }
    }
}

fn term_for_degree(var: Symbol, degree: i64, coeff: Rational) -> Option<Expr> {
    if coeff.is_zero() {
        return Some(const_(Rational::zero()));
    }
    if degree == 0 {
        return Some(const_(coeff));
    }
    let x = Expr::var(var);
    let mono = if degree == 1 {
        x
    } else {
        pow(x, const_(Rational::from(degree)))
    };
    if coeff.is_one() {
        Some(mono)
    } else {
        Some(mul(const_(coeff), mono))
    }
}

fn as_const(e: &Expr) -> Option<Rational> {
    match e.kind() {
        ExprKind::Const(c) => c.try_as_rational(),
        _ => None,
    }
}

fn coeff_monomial(e: &Expr) -> Option<(Rational, Expr)> {
    match e.kind() {
        ExprKind::Const(c) => c
            .try_as_rational()
            .map(|r| (r, const_(Rational::one()))),
        ExprKind::Var(_) => Some((Rational::one(), e.clone())),
        ExprKind::Mul(l, r) => {
            if let Some(c) = as_const(l) {
                let (_, m) = coeff_monomial(r)?;
                return Some((c, m));
            }
            if let Some(c) = as_const(r) {
                let (_, m) = coeff_monomial(l)?;
                return Some((c, m));
            }
            None
        }
        ExprKind::Pow(base, exp) => {
            let n = as_const(exp)?;
            if let Some(k) = n.as_integer() {
                if k >= 0 {
                    return Some((Rational::one(), pow(base.clone(), const_(n))));
                }
            }
            None
        }
        ExprKind::Div(num, den) => {
            let cd = as_const(den)?;
            if cd.is_zero() {
                return None;
            }
            let (cn, m) = coeff_monomial(num)?;
            Some((cn / cd, m))
        }
        _ => None,
    }
}
