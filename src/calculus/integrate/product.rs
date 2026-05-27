use crate::expr::{add, const_, div, mul, pow, sub, Expr, ExprKind};
use crate::rational::Rational;
use crate::symbol::Symbol;
use std::collections::BTreeMap;

use super::util::{
    affine_angle, collect_polynomial_terms, contains_var, exp_affine, is_exp_of_var, is_var,
    trig_affine,
};
use super::trig_product::try_sin_cos_power_product;
use super::{integrate, integrate_expr, IntegrateError};

const PARTS_MAX_DEPTH: u8 = 4;

pub fn integrate_product(f: Expr, g: Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    if !contains_var(&f, var) {
        return Ok(f * integrate_expr(g, var)?);
    }
    if !contains_var(&g, var) {
        return Ok(g * integrate_expr(f, var)?);
    }

    if let Some(integ) = integrate_exp_times_polynomial(&f, &g, var) {
        return Ok(integ);
    }
    if let Some(integ) = integrate_exp_times_polynomial(&g, &f, var) {
        return Ok(integ);
    }

    if let Some(integ) = try_sin_cos_power_product(&(f.clone() * g.clone()).simplify(), var) {
        return Ok(integ);
    }

    if let Some(integ) = integrate_sin_cos_product(&f, &g, var) {
        return Ok(integ);
    }
    if let Some(integ) = integrate_sin_cos_product(&g, &f, var) {
        return Ok(integ);
    }

    if let Some(integ) = integrate_exp_times_trig(&f, &g, var) {
        return Ok(integ);
    }
    if let Some(integ) = integrate_exp_times_trig(&g, &f, var) {
        return Ok(integ);
    }

    if let Some(integ) = integrate_polynomial_times_ln(&f, &g, var) {
        return Ok(integ);
    }
    if let Some(integ) = integrate_polynomial_times_ln(&g, &f, var) {
        return Ok(integ);
    }

    if let Some(integ) = try_integration_by_parts(&f, &g, var, PARTS_MAX_DEPTH) {
        return Ok(integ);
    }
    if let Some(integ) = try_integration_by_parts(&g, &f, var, PARTS_MAX_DEPTH) {
        return Ok(integ);
    }

    Err(IntegrateError::NoRule)
}

/// ∫ sin(ax+b) cos(cx+d) dx via product-to-sum.
fn integrate_sin_cos_product(f: &Expr, g: &Expr, var: Symbol) -> Option<Expr> {
    let (ka, ba) = trig_affine(f, var, true)?;
    let (kc, bc) = trig_affine(g, var, false)?;
    let sum_k = ka + kc;
    let diff_k = ka - kc;
    let sum_b = ba + bc;
    let diff_b = ba - bc;
    if diff_k.is_zero() {
        let angle = affine_angle(Rational::from(2) * ka, ba + bc, var);
        let mut sum = -crate::expr::cos(angle) / const_(Rational::from(4) * ka);
        if diff_b != Rational::zero() {
            sum = sum
                + mul(
                    Expr::var(var),
                    crate::expr::sin(crate::expr::const_(diff_b)),
                ) / const_(Rational::from(2));
        }
        return Some(sum);
    }
    if sum_k.is_zero() {
        return None;
    }
    let sum_angle = affine_angle(sum_k, sum_b, var);
    let diff_angle = affine_angle(diff_k, diff_b, var);
    Some(
        -crate::expr::cos(sum_angle) / const_(Rational::from(2) * sum_k)
            - crate::expr::cos(diff_angle) / const_(Rational::from(2) * diff_k),
    )
}

fn try_integration_by_parts(u: &Expr, dv: &Expr, var: Symbol, depth: u8) -> Option<Expr> {
    if depth == 0 {
        return None;
    }
    let v = integrate(dv.clone(), var).ok()?;
    let du = u.clone().diff(var).simplify();
    if !contains_var(&du, var) {
        let second = if matches!(du.kind(), ExprKind::Const(c) if c.is_zero()) {
            const_(Rational::zero())
        } else {
            integrate(v.clone(), var).ok()? * du
        };
        return Some(u.clone() * v - second);
    }
    if let Ok(rest) = integrate(v.clone() * du.clone(), var) {
        return Some(u.clone() * v - rest);
    }
    if depth <= 1 {
        return None;
    }
    let inner = try_integration_by_parts(&v, &du, var, depth - 1)?;
    Some(u.clone() * v - inner)
}

/// ∫ e^{ax+b} sin(cx+d) dx and ∫ e^{ax+b} cos(cx+d) dx.
fn integrate_exp_times_trig(f: &Expr, g: &Expr, var: Symbol) -> Option<Expr> {
    exp_trig_pair(f, g, var)
}

fn exp_trig_pair(exp_side: &Expr, trig_side: &Expr, var: Symbol) -> Option<Expr> {
    let (a, eb) = exp_affine(exp_side, var)?;
    let (c, ed, want_sin) = if let Some(tc) = trig_affine(trig_side, var, true) {
        (tc.0, tc.1, true)
    } else {
        let (c, ed) = trig_affine(trig_side, var, false)?;
        (c, ed, false)
    };
    exp_trig_template(a, eb, c, ed, want_sin, var)
}

/// ∫ e^{ax+b} sin(cx+d) or cos(cx+d) dx = e^{ax+b}(a·trig −/+ c·partner) / (a²+c²).
fn exp_trig_template(
    a: Rational,
    eb: Rational,
    c: Rational,
    ed: Rational,
    want_sin: bool,
    var: Symbol,
) -> Option<Expr> {
    let denom = a * a + c * c;
    if denom.is_zero() {
        return None;
    }
    let exp_part = crate::expr::exp(affine_angle(a, eb, var));
    let angle = affine_angle(c, ed, var);
    if want_sin {
        Some(div(
            mul(
                exp_part,
                sub(
                    mul(const_(a), crate::expr::sin(angle.clone())),
                    mul(const_(c), crate::expr::cos(angle)),
                ),
            ),
            const_(denom),
        ))
    } else {
        Some(div(
            mul(
                exp_part,
                add(
                    mul(const_(a), crate::expr::cos(angle.clone())),
                    mul(const_(c), crate::expr::sin(angle)),
                ),
            ),
            const_(denom),
        ))
    }
}

/// ∫ x^n ln(x) dx for polynomial `f` and ln(var).
fn integrate_polynomial_times_ln(f: &Expr, g: &Expr, var: Symbol) -> Option<Expr> {
    match g.kind() {
        ExprKind::Ln(inner) if is_var(inner, var) => {}
        _ => return None,
    }
    let terms = collect_polynomial_terms(f, var).ok()?;
    let mut sum = const_(Rational::zero());
    for (&n, coeff) in terms.iter().rev() {
        if coeff.is_zero() {
            continue;
        }
        if n < 0 {
            return None;
        }
        let np = n + 1;
        let denom = Rational::from(np) * Rational::from(np);
        let x = Expr::var(var);
        let xp = if np == 1 {
            x.clone()
        } else {
            pow(x.clone(), const_(Rational::from(np)))
        };
        let term = mul(
            const_(*coeff),
            mul(
                xp,
                sub(
                    mul(const_(Rational::from(np)), crate::expr::ln(x.clone())),
                    const_(Rational::one()),
                ),
            ),
        ) / const_(denom);
        sum = add(sum, term);
    }
    Some(sum)
}

fn integrate_exp_times_polynomial(
    exp_side: &Expr,
    poly: &Expr,
    var: Symbol,
) -> Option<Expr> {
    if !is_exp_of_var(exp_side, var) {
        return None;
    }
    let terms = collect_polynomial_terms(poly, var).ok()?;
    if terms.is_empty() {
        return Some(const_(Rational::zero()));
    }
    let mut poly_part = const_(Rational::zero());
    for (n, coeff) in terms.iter().rev() {
        if coeff.is_zero() {
            continue;
        }
        let p = exp_integrate_polynomial_factor(*n, var);
        poly_part = add(poly_part, mul(const_(*coeff), p));
    }
    Some(mul(crate::expr::exp(Expr::var(var)), poly_part))
}

fn exp_integrate_polynomial_factor(n: i64, var: Symbol) -> Expr {
    let mut coeffs: BTreeMap<i64, Rational> = BTreeMap::new();
    coeffs.insert(0, Rational::one());
    for i in 1..=n {
        let mut next = BTreeMap::new();
        next.insert(i, Rational::one());
        for (k, c) in &coeffs {
            let slot = next.entry(*k).or_insert(Rational::zero());
            *slot = *slot - Rational::from(i) * *c;
        }
        coeffs = next;
    }
    let x = Expr::var(var);
    let mut sum = const_(Rational::zero());
    for (k, c) in coeffs {
        if c.is_zero() {
            continue;
        }
        let term = if k == 0 {
            const_(c)
        } else if k == 1 {
            mul(const_(c), x.clone())
        } else {
            mul(const_(c), pow(x.clone(), const_(Rational::from(k))))
        };
        sum = add(sum, term);
    }
    sum
}
