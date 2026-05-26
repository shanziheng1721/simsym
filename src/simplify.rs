use crate::expr::{add, const_, div, mul, neg, pow, sub, Expr, ExprKind};
use crate::rational::Rational;
use std::collections::BTreeMap;

pub fn simplify(expr: Expr) -> Expr {
    let e = simplify_once(expr);
    let e = simplify_factor_exp_terms(e, 0);
    normalize_polynomial_shapes(e)
}

/// Apply polynomial normal form to sums and to `e^x * poly` factors.
fn normalize_polynomial_shapes(e: Expr) -> Expr {
    match e.kind().clone() {
        ExprKind::Mul(l, r) => {
            let l = normalize_polynomial_shapes(l);
            let r = normalize_polynomial_shapes(r);
            if let ExprKind::Exp(_) = l.kind() {
                if let Some(inner) = crate::poly::try_polynomial_normal_form(r.clone()) {
                    return mul(l, inner);
                }
            }
            if let ExprKind::Exp(_) = r.kind() {
                if let Some(inner) = crate::poly::try_polynomial_normal_form(l.clone()) {
                    return mul(r, inner);
                }
            }
            mul(l, r)
        }
        ExprKind::Add(..) | ExprKind::Sub(..) | ExprKind::Neg(_) => {
            crate::poly::try_polynomial_normal_form(e.clone()).unwrap_or(e)
        }
        ExprKind::Exp(inner) => crate::expr::exp(normalize_polynomial_shapes(inner)),
        ExprKind::Sin(inner) => crate::expr::sin(normalize_polynomial_shapes(inner)),
        ExprKind::Cos(inner) => crate::expr::cos(normalize_polynomial_shapes(inner)),
        ExprKind::Tan(inner) => crate::expr::tan(normalize_polynomial_shapes(inner)),
        ExprKind::Ln(inner) => crate::expr::ln(normalize_polynomial_shapes(inner)),
        _ => e,
    }
}

const SIMPLIFY_MAX_DEPTH: u32 = 64;

fn simplify_once(expr: Expr) -> Expr {
    match expr.kind().clone() {
        ExprKind::Const(c) => const_(c),
        ExprKind::Var(s) => Expr::var(s),
        ExprKind::Add(a, b) => simplify_add(simplify_once(a), simplify_once(b)),
        ExprKind::Sub(a, b) => simplify_sub(simplify_once(a), simplify_once(b)),
        ExprKind::Mul(a, b) => simplify_mul(simplify_once(a), simplify_once(b)),
        ExprKind::Div(a, b) => simplify_div(simplify_once(a), simplify_once(b)),
        ExprKind::Neg(e) => simplify_neg(simplify_once(e)),
        ExprKind::Pow(base, exp) => simplify_pow(simplify_once(base), simplify_once(exp)),
        ExprKind::Sin(e) => crate::expr::sin(simplify_trig_arg(simplify_once(e))),
        ExprKind::Cos(e) => crate::expr::cos(simplify_trig_arg(simplify_once(e))),
        ExprKind::Tan(e) => crate::expr::tan(simplify_trig_arg(simplify_once(e))),
        ExprKind::Exp(e) => crate::expr::exp(simplify_once(e)),
        ExprKind::Ln(e) => simplify_ln(simplify_once(e)),
    }
}

/// `e^x*a + e^x*b` → `e^x*(a+b)`
fn simplify_factor_exp_terms(e: Expr, depth: u32) -> Expr {
    if depth > SIMPLIFY_MAX_DEPTH {
        return e;
    }
    let d = depth + 1;
    match e.kind().clone() {
        ExprKind::Add(a, b) => {
            let a = simplify_factor_exp_terms(simplify_once(a), d);
            let b = simplify_factor_exp_terms(simplify_once(b), d);
            if let Some(factored) = try_factor_exp_add(&a, &b) {
                return simplify_factor_exp_terms(simplify_once(factored), d);
            }
            add(a, b)
        }
        ExprKind::Sub(a, b) => sub(
            simplify_factor_exp_terms(simplify_once(a), d),
            simplify_factor_exp_terms(simplify_once(b), d),
        ),
        ExprKind::Mul(l, r) => mul(
            simplify_factor_exp_terms(simplify_once(l), d),
            simplify_factor_exp_terms(simplify_once(r), d),
        ),
        ExprKind::Neg(inner) => neg(simplify_factor_exp_terms(simplify_once(inner), d)),
        _ => e,
    }
}

fn try_factor_exp_add(a: &Expr, b: &Expr) -> Option<Expr> {
    let (ea, ra) = as_exp_times(a)?;
    let (eb, rb) = as_exp_times(b)?;
    if ea != eb {
        return None;
    }
    let sum = crate::poly::polynomial_normal_form(add(ra, rb));
    Some(mul(crate::expr::exp(ea), sum))
}

fn as_exp_times(e: &Expr) -> Option<(Expr, Expr)> {
    match e.kind() {
        ExprKind::Mul(l, r) => {
            if let ExprKind::Exp(inner) = l.kind() {
                return Some((inner.clone(), r.clone()));
            }
            if let ExprKind::Exp(inner) = r.kind() {
                return Some((inner.clone(), l.clone()));
            }
        }
        ExprKind::Exp(inner) => {
            return Some((inner.clone(), const_(Rational::one())));
        }
        _ => {}
    }
    None
}

fn simplify_ln(e: Expr) -> Expr {
    if let ExprKind::Exp(inner) = e.kind() {
        return simplify_once(inner.clone());
    }
    crate::expr::ln(e)
}

fn simplify_trig_arg(e: Expr) -> Expr {
    if let ExprKind::Const(c) = e.kind() {
        if c.is_zero() {
            return const_(Rational::zero());
        }
        if c.is_one() {
            // sin(1), cos(1) etc. — keep symbolic
        }
    }
    e
}

fn simplify_add(a: Expr, b: Expr) -> Expr {
    let sum = add(simplify_once(a.clone()), simplify_once(b.clone()));
    crate::poly::try_polynomial_normal_form(sum).unwrap_or_else(|| {
        let mut terms = Vec::new();
        flatten_add_terms(&a, &mut terms);
        flatten_add_terms(&b, &mut terms);
        let mut acc = const_(Rational::zero());
        for t in terms {
            acc = simplify_add_step(acc, t);
        }
        acc
    })
}

fn simplify_add_step(a: Expr, b: Expr) -> Expr {
    if is_zero(&a) {
        return b;
    }
    if is_zero(&b) {
        return a;
    }
    if let (Some(ca), Some(cb)) = (as_const(&a), as_const(&b)) {
        return const_(ca + cb);
    }
    if a == b {
        return mul(const_(Rational::from(2)), a);
    }
    try_merge_polynomial_terms(&a, &b)
        .unwrap_or_else(|| add(a, b))
}

fn flatten_add_terms(e: &Expr, out: &mut Vec<Expr>) {
    match e.kind() {
        ExprKind::Add(l, r) => {
            flatten_add_terms(l, out);
            flatten_add_terms(r, out);
        }
        ExprKind::Sub(l, r) => {
            flatten_add_terms(l, out);
            flatten_add_terms(r, out);
            // Sub handled as add(a, neg(b)) at upper level if needed
        }
        _ => out.push(e.clone()),
    }
}

fn simplify_sub(a: Expr, b: Expr) -> Expr {
    if is_zero(&b) {
        return a;
    }
    if a == b {
        return const_(Rational::zero());
    }
    if let (Some(ca), Some(cb)) = (as_const(&a), as_const(&b)) {
        return const_(ca - cb);
    }
    sub(a, b)
}

fn simplify_mul(a: Expr, b: Expr) -> Expr {
    if is_zero(&a) || is_zero(&b) {
        return const_(Rational::zero());
    }
    if is_one(&a) {
        return simplify_once(b);
    }
    if is_one(&b) {
        return simplify_once(a);
    }
    if let (Some(ca), Some(cb)) = (as_const(&a), as_const(&b)) {
        return const_(ca * cb);
    }
    if let Some(scale) = scale_monomial_product(&a, &b) {
        return scale;
    }
    if let Some(c) = as_const(&a) {
        if let ExprKind::Mul(l, r) = b.kind() {
            if let Some(d) = as_const(l) {
                return simplify_once(distribute_const_mul(c * d, r));
            }
            if let Some(d) = as_const(r) {
                return simplify_once(distribute_const_mul(c * d, l));
            }
        }
        return simplify_once(distribute_const_mul(c, &b));
    }
    if let Some(c) = as_const(&b) {
        if let ExprKind::Mul(l, r) = a.kind() {
            if let Some(d) = as_const(l) {
                return simplify_once(distribute_const_mul(c * d, r));
            }
            if let Some(d) = as_const(r) {
                return simplify_once(distribute_const_mul(c * d, l));
            }
        }
        return simplify_once(distribute_const_mul(c, &a));
    }
    mul(simplify_once(a), simplify_once(b))
}

fn distribute_const_mul(c: Rational, expr: &Expr) -> Expr {
    match expr.kind() {
        ExprKind::Add(l, r) => add(mul(const_(c), l.clone()), mul(const_(c), r.clone())),
        ExprKind::Sub(l, r) => sub(mul(const_(c), l.clone()), mul(const_(c), r.clone())),
        ExprKind::Neg(e) => mul(const_(c), neg(e.clone())),
        _ => mul(const_(c), expr.clone()),
    }
}

/// Combine `c * (d * m)` or `(c * m) * d` into a single scaled monomial.
fn scale_monomial_product(a: &Expr, b: &Expr) -> Option<Expr> {
    if let Some(ca) = as_const(a) {
        if let Some((cb, m)) = coeff_monomial(b) {
            return Some(scaled_monomial(ca * cb, m));
        }
    }
    if let Some(cb) = as_const(b) {
        if let Some((ca, m)) = coeff_monomial(a) {
            return Some(scaled_monomial(ca * cb, m));
        }
    }
    None
}

fn scaled_monomial(coeff: Rational, m: Expr) -> Expr {
    if coeff.is_zero() {
        const_(Rational::zero())
    } else if coeff.is_one() {
        simplify_once(m)
    } else {
        mul(const_(coeff), simplify_once(m))
    }
}

fn simplify_div(a: Expr, b: Expr) -> Expr {
    if is_one(&b) {
        return simplify_once(a);
    }
    if is_zero(&b) {
        return div(a, b);
    }
    if let Some(cb) = as_const(&b) {
        if cb.is_zero() {
            return div(a, b);
        }
        if let Some((ca, m)) = coeff_monomial(&a) {
            let coeff = ca / cb;
            if coeff.is_zero() {
                return const_(Rational::zero());
            }
            if coeff.is_one() {
                return simplify_once(m);
            }
            return simplify_mul(const_(coeff), m);
        }
    }
    if let (Some(ca), Some(cb)) = (as_const(&a), as_const(&b)) {
        return const_(ca / cb);
    }
    div(simplify_once(a), simplify_once(b))
}

fn simplify_neg(e: Expr) -> Expr {
    if let Some(c) = as_const(&e) {
        return const_(-c);
    }
    if let ExprKind::Neg(inner) = e.kind() {
        return inner.clone();
    }
    neg(e)
}

fn simplify_pow(base: Expr, exp: Expr) -> Expr {
    if is_euler_base(&base) {
        return crate::expr::exp(simplify_once(exp));
    }
    if is_one(&exp) {
        return base;
    }
    if is_zero(&exp) {
        return const_(Rational::one());
    }
    if let (Some(b), Some(e)) = (as_const(&base), as_const(&exp)) {
        if let Some(n) = e.as_integer() {
            if let Ok(v) = crate::rational::rat_pow_int(b, n) {
                return const_(v);
            }
        }
    }
    if let Some(e) = as_const(&exp) {
        if e.is_zero() {
            return const_(Rational::one());
        }
        if e.is_one() {
            return base;
        }
    }
    pow(base, exp)
}

fn is_zero(e: &Expr) -> bool {
    matches!(e.kind(), ExprKind::Const(c) if c.is_zero())
}

fn is_one(e: &Expr) -> bool {
    matches!(e.kind(), ExprKind::Const(c) if c.is_one())
}

fn is_euler_base(e: &Expr) -> bool {
    matches!(e.kind(), ExprKind::Var(s) if s.name() == "e")
}

fn as_const(e: &Expr) -> Option<Rational> {
    match e.kind() {
        ExprKind::Const(c) => Some(*c),
        _ => None,
    }
}

/// Merge `a + b` when both are scalar multiples of the same monomial (single-var polynomial terms).
fn try_merge_polynomial_terms(a: &Expr, b: &Expr) -> Option<Expr> {
    let (ca, ma) = coeff_monomial(a)?;
    let (cb, mb) = coeff_monomial(b)?;
    if ma != mb {
        return None;
    }
    let sum = ca + cb;
    if sum.is_zero() {
        return Some(const_(Rational::zero()));
    }
    Some(mul(const_(sum), ma))
}

fn coeff_monomial(e: &Expr) -> Option<(Rational, Expr)> {
    match e.kind() {
        ExprKind::Const(c) => Some((*c, const_(Rational::one()))),
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

/// Flatten a sum into (coeff, monomial) pairs for one variable — used in tests / future poly normal form.
#[allow(dead_code)]
fn collect_terms(expr: &Expr, var: crate::symbol::Symbol) -> BTreeMap<i64, Rational> {
    let mut map = BTreeMap::new();
    collect_terms_rec(expr, var, Rational::one(), &mut map);
    map
}

fn collect_terms_rec(
    expr: &Expr,
    var: crate::symbol::Symbol,
    coeff: Rational,
    out: &mut BTreeMap<i64, Rational>,
) {
    match expr.kind() {
        ExprKind::Const(_) => {}
        ExprKind::Var(s) if *s == var => {
            *out.entry(1).or_insert(Rational::zero()) += coeff;
        }
        ExprKind::Add(a, b) | ExprKind::Sub(a, b) => {
            collect_terms_rec(a, var, coeff, out);
            let cb = if matches!(expr.kind(), ExprKind::Sub(_, _)) {
                -coeff
            } else {
                coeff
            };
            collect_terms_rec(b, var, cb, out);
        }
        ExprKind::Mul(l, r) => {
            if let Some(c) = as_const(l) {
                collect_terms_rec(r, var, coeff * c, out);
            } else if let Some(c) = as_const(r) {
                collect_terms_rec(l, var, coeff * c, out);
            }
        }
        ExprKind::Pow(base, exp) if matches!(base.kind(), ExprKind::Var(s) if *s == var) => {
            if let Some(n) = as_const(exp).and_then(|r| r.as_integer()) {
                *out.entry(n).or_insert(Rational::zero()) += coeff;
            }
        }
        _ => {}
    }
}
