use crate::expr::{add, const_, div, mul, neg, pow, sub, Expr, ExprKind};
use crate::rational::Rational;
use std::cell::Cell;
use std::collections::BTreeMap;

thread_local! {
    static SIMPLIFY_DEPTH: Cell<u32> = const { Cell::new(0) };
}
const SIMPLIFY_RECURSION_LIMIT: u32 = 256;

pub fn simplify(expr: Expr) -> Expr {
    let e = simplify_once(expr);
    let e = simplify_factor_exp_terms(e, 0);
    normalize_polynomial_shapes(e)
}

/// Apply polynomial normal form to sums and to `e^x * poly` factors.
fn normalize_polynomial_shapes(e: Expr) -> Expr {
    match e.into_kind() {
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
        kind @ (ExprKind::Add(..) | ExprKind::Sub(..) | ExprKind::Neg(_)) => {
            let e = Expr::from_kind(kind);
            crate::poly::try_polynomial_normal_form(e.clone()).unwrap_or(e)
        }
        ExprKind::Exp(inner) => crate::expr::exp(normalize_polynomial_shapes(inner)),
        ExprKind::Sin(inner) => crate::expr::sin(normalize_polynomial_shapes(inner)),
        ExprKind::Cos(inner) => crate::expr::cos(normalize_polynomial_shapes(inner)),
        ExprKind::Tan(inner) => crate::expr::tan(normalize_polynomial_shapes(inner)),
        ExprKind::Cot(inner) => crate::expr::cot(normalize_polynomial_shapes(inner)),
        ExprKind::Sec(inner) => crate::expr::sec(normalize_polynomial_shapes(inner)),
        ExprKind::Csc(inner) => crate::expr::csc(normalize_polynomial_shapes(inner)),
        ExprKind::Asin(inner) => crate::expr::asin(normalize_polynomial_shapes(inner)),
        ExprKind::Acos(inner) => crate::expr::acos(normalize_polynomial_shapes(inner)),
        ExprKind::Atan(inner) => crate::expr::atan(normalize_polynomial_shapes(inner)),
        ExprKind::Acot(inner) => crate::expr::acot(normalize_polynomial_shapes(inner)),
        ExprKind::Asec(inner) => crate::expr::asec(normalize_polynomial_shapes(inner)),
        ExprKind::Acsc(inner) => crate::expr::acsc(normalize_polynomial_shapes(inner)),
        ExprKind::Sinh(inner) => crate::expr::sinh(normalize_polynomial_shapes(inner)),
        ExprKind::Cosh(inner) => crate::expr::cosh(normalize_polynomial_shapes(inner)),
        ExprKind::Tanh(inner) => crate::expr::tanh(normalize_polynomial_shapes(inner)),
        ExprKind::Coth(inner) => crate::expr::coth(normalize_polynomial_shapes(inner)),
        ExprKind::Sech(inner) => crate::expr::sech(normalize_polynomial_shapes(inner)),
        ExprKind::Csch(inner) => crate::expr::csch(normalize_polynomial_shapes(inner)),
        ExprKind::Asinh(inner) => crate::expr::asinh(normalize_polynomial_shapes(inner)),
        ExprKind::Acosh(inner) => crate::expr::acosh(normalize_polynomial_shapes(inner)),
        ExprKind::Atanh(inner) => crate::expr::atanh(normalize_polynomial_shapes(inner)),
        ExprKind::Acoth(inner) => crate::expr::acoth(normalize_polynomial_shapes(inner)),
        ExprKind::Asech(inner) => crate::expr::asech(normalize_polynomial_shapes(inner)),
        ExprKind::Acsch(inner) => crate::expr::acsch(normalize_polynomial_shapes(inner)),
        ExprKind::Ln(inner) => crate::expr::ln(normalize_polynomial_shapes(inner)),
        k => Expr::from_kind(k),
    }
}

const SIMPLIFY_MAX_DEPTH: u32 = 64;

fn simplify_once(expr: Expr) -> Expr {
    let depth = SIMPLIFY_DEPTH.with(|d| d.get());
    if depth > SIMPLIFY_RECURSION_LIMIT {
        return expr;
    }
    SIMPLIFY_DEPTH.with(|d| d.set(depth + 1));
    let out = simplify_once_inner(expr);
    SIMPLIFY_DEPTH.with(|d| d.set(depth));
    out
}

fn simplify_once_inner(expr: Expr) -> Expr {
    match expr.into_kind() {
        ExprKind::Const(c) => crate::constant::constant(c.clone()),
        ExprKind::Var(s) => Expr::var(s),
        ExprKind::Add(a, b) => simplify_add(simplify_once_inner(a), simplify_once_inner(b)),
        ExprKind::Sub(a, b) => simplify_sub(simplify_once_inner(a), simplify_once_inner(b)),
        ExprKind::Mul(a, b) => simplify_mul(simplify_once_inner(a), simplify_once_inner(b)),
        ExprKind::Div(a, b) => simplify_div(simplify_once_inner(a), simplify_once_inner(b)),
        ExprKind::Neg(e) => simplify_neg(simplify_once_inner(e)),
        ExprKind::Pow(base, exp) => simplify_pow(simplify_once_inner(base), simplify_once_inner(exp)),
        ExprKind::Sin(e) => crate::expr::sin(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Cos(e) => crate::expr::cos(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Tan(e) => crate::expr::tan(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Cot(e) => crate::expr::cot(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Sec(e) => crate::expr::sec(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Csc(e) => crate::expr::csc(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Asin(e) => crate::expr::asin(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Acos(e) => crate::expr::acos(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Atan(e) => crate::expr::atan(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Acot(e) => crate::expr::acot(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Asec(e) => crate::expr::asec(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Acsc(e) => crate::expr::acsc(simplify_trig_arg(simplify_once_inner(e))),
        ExprKind::Sinh(e) => crate::expr::sinh(simplify_once_inner(e)),
        ExprKind::Cosh(e) => crate::expr::cosh(simplify_once_inner(e)),
        ExprKind::Tanh(e) => crate::expr::tanh(simplify_once_inner(e)),
        ExprKind::Coth(e) => crate::expr::coth(simplify_once_inner(e)),
        ExprKind::Sech(e) => crate::expr::sech(simplify_once_inner(e)),
        ExprKind::Csch(e) => crate::expr::csch(simplify_once_inner(e)),
        ExprKind::Asinh(e) => crate::expr::asinh(simplify_once_inner(e)),
        ExprKind::Acosh(e) => crate::expr::acosh(simplify_once_inner(e)),
        ExprKind::Atanh(e) => crate::expr::atanh(simplify_once_inner(e)),
        ExprKind::Acoth(e) => crate::expr::acoth(simplify_once_inner(e)),
        ExprKind::Asech(e) => crate::expr::asech(simplify_once_inner(e)),
        ExprKind::Acsch(e) => crate::expr::acsch(simplify_once_inner(e)),
        ExprKind::Exp(e) => crate::expr::exp(simplify_once_inner(e)),
        ExprKind::Ln(e) => simplify_ln(simplify_once_inner(e)),
    }
}

/// `e^x*a + e^x*b` → `e^x*(a+b)`
fn simplify_factor_exp_terms(e: Expr, depth: u32) -> Expr {
    if depth > SIMPLIFY_MAX_DEPTH {
        return e;
    }
    let d = depth + 1;
    match e.into_kind() {
        ExprKind::Add(a, b) => {
            let a = simplify_factor_exp_terms(simplify_once(a), d);
            let b = simplify_factor_exp_terms(simplify_once(b), d);
            if let Some(factored) = try_factor_exp_add(&a, &b) {
                let unchanged = add(a.clone(), b.clone());
                if factored != unchanged {
                    return simplify_factor_exp_terms(simplify_once(factored), d);
                }
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
        k => Expr::from_kind(k),
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
    let a = simplify_once_inner(a);
    let b = simplify_once_inner(b);
    if let Some(p) = crate::poly::try_polynomial_normal_form(add(a.clone(), b.clone())) {
        return p;
    }
    simplify_add_step(a, b)
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
        ExprKind::Const(c) => c.try_as_rational(),
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
