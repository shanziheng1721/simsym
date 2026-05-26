//! Integration of rational functions P(x)/Q(x): polynomial division and partial fractions.

use crate::expr::{add, const_, div, mul, pow, sub, Expr};
use crate::rational::Rational;
use crate::symbol::Symbol;
use std::collections::BTreeMap;

use super::util::collect_polynomial_terms;
use super::IntegrateError;

type Poly = BTreeMap<i64, Rational>;

pub fn integrate_rational_function(f: Expr, g: Expr, var: Symbol) -> Result<Expr, IntegrateError> {
    let num = collect_polynomial_terms(&f, var)?;
    let den = collect_polynomial_terms(&g, var)?;
    if den.is_empty() || poly_degree(&den) == 0 {
        return Err(IntegrateError::NoRule);
    }
    if poly_degree(&den) > 6 {
        return Err(IntegrateError::NoRule);
    }
    let (quot, rem) = poly_div(num, &den)?;
    let mut result = integrate_polynomial(&quot, var)?;
    if !poly_is_zero(&rem) {
        if poly_degree(&rem) >= poly_degree(&den) {
            return Err(IntegrateError::NoRule);
        }
        result = result + integrate_partial_fractions(&rem, &den, var)?;
    }
    Ok(result)
}

fn integrate_partial_fractions(num: &Poly, den: &Poly, var: Symbol) -> Result<Expr, IntegrateError> {
    let terms = partial_fraction_decompose(num, den)?;
    let mut sum = const_(Rational::zero());
    for t in terms {
        sum = sum + integrate_pf_term(&t, var)?;
    }
    Ok(sum)
}

#[derive(Clone)]
struct PfTerm {
    numer: Poly,
    quad_b: Option<Rational>,
    quad_c: Option<Rational>,
    root: Option<Rational>,
    power: u32,
}

fn integrate_pf_term(t: &PfTerm, var: Symbol) -> Result<Expr, IntegrateError> {
    let x = Expr::var(var);
    if let Some(a) = t.root {
        let xa = sub(x.clone(), const_(a));
        let c = t.numer.get(&0).copied().unwrap_or(Rational::zero());
        let k = t.power as i64;
        if k == 1 {
            return Ok(mul(const_(c), crate::expr::ln(xa)));
        }
        return Ok(mul(const_(c), pow(xa, const_(Rational::from(1 - k))))
            / const_(Rational::from(1 - k)));
    }
    let b = t.quad_b.ok_or(IntegrateError::NoRule)?;
    let c = t.quad_c.ok_or(IntegrateError::NoRule)?;
    let a_lin = t.numer.get(&1).copied().unwrap_or(Rational::zero());
    let b_const = t.numer.get(&0).copied().unwrap_or(Rational::zero());
    if t.power != 1 {
        return Err(IntegrateError::NoRule);
    }
    let disc = b * b - Rational::from(4) * c;
    if disc >= Rational::zero() {
        return Err(IntegrateError::NoRule);
    }
    let half_b = b / Rational::from(2);
    let d = c - half_b * half_b;
    if d <= Rational::zero() {
        return Err(IntegrateError::NoRule);
    }
    let s = rational_sqrt(d).ok_or(IntegrateError::NoRule)?;
    let u = sub(x.clone(), const_(half_b));
    let mut sum = const_(Rational::zero());
    if !a_lin.is_zero() {
        let quad = add(
            pow(x.clone(), const_(Rational::from(2))),
            add(mul(const_(b), x.clone()), const_(c)),
        );
        sum = sum + mul(const_(a_lin / Rational::from(2)), crate::expr::ln(quad));
    }
    let atan_coeff = b_const - a_lin * half_b;
    if !atan_coeff.is_zero() {
        sum = sum + mul(const_(atan_coeff / s), crate::expr::atan(div(u, const_(s))));
    }
    Ok(sum)
}

fn partial_fraction_decompose(num: &Poly, den: &Poly) -> Result<Vec<PfTerm>, IntegrateError> {
    let factors = factor_polynomial(den)?;
    let mut terms = Vec::new();
    let mut n = num.clone();
    for (fac, mult) in factors {
        match fac {
            Factor::Linear(root) => {
                for k in (1..=mult).rev() {
                    let cofactor = build_cofactor(den, &Factor::Linear(root), k)?;
                    let coeff = poly_eval(&n, root) / poly_eval(&cofactor, root);
                    if !coeff.is_zero() {
                        let mut numer = BTreeMap::new();
                        numer.insert(0, coeff);
                        terms.push(PfTerm {
                            numer,
                            quad_b: None,
                            quad_c: None,
                            root: Some(root),
                            power: k,
                        });
                    }
                    let term_poly = poly_scale(&cofactor, coeff);
                    n = poly_sub(n, term_poly);
                    if k > 1 {
                        let lin = linear_poly(root);
                        let (n2, _) = poly_div(n, &lin)?;
                        n = n2;
                    }
                }
            }
            Factor::Quadratic { b, c } => {
                if mult != 1 {
                    return Err(IntegrateError::NoRule);
                }
                let cofactor = build_cofactor(den, &Factor::Quadratic { b, c }, 1)?;
                let q = quadratic_poly(b, c);
                let (a_lin, b_const) = solve_linear_over_quadratic(&n, &q, &cofactor)?;
                if !a_lin.is_zero() || !b_const.is_zero() {
                    let mut numer = BTreeMap::new();
                    if !b_const.is_zero() {
                        numer.insert(0, b_const);
                    }
                    if !a_lin.is_zero() {
                        numer.insert(1, a_lin);
                    }
                    terms.push(PfTerm {
                        numer,
                        quad_b: Some(b),
                        quad_c: Some(c),
                        root: None,
                        power: 1,
                    });
                }
                let term_poly = poly_add(
                    poly_scale(&q, a_lin),
                    poly_scale(&cofactor, b_const),
                );
                n = poly_sub(n, term_poly);
            }
        }
    }
    if !poly_is_zero(&n) {
        return Err(IntegrateError::NoRule);
    }
    Ok(terms)
}

fn solve_linear_over_quadratic(
    num: &Poly,
    quad: &Poly,
    cofactor: &Poly,
) -> Result<(Rational, Rational), IntegrateError> {
    if poly_degree(num) == 0 {
        let b = poly_eval(num, Rational::zero()) / poly_eval(cofactor, Rational::zero());
        return Ok((Rational::zero(), b));
    }
    let n0 = poly_eval(num, Rational::zero());
    let n1 = poly_coeff(num, 1);
    let q0 = poly_eval(quad, Rational::zero());
    let q1 = poly_coeff(quad, 1);
    let c0 = poly_eval(cofactor, Rational::zero());
    let c1 = poly_coeff(cofactor, 1);
    let det = q1 * c0 - q0 * c1;
    if det.is_zero() {
        if cofactor.len() == 1 {
            let a = n1 / q1;
            let b = (n0 - a * q0) / c0;
            return Ok((a, b));
        }
        return Err(IntegrateError::NoRule);
    }
    let a = (n1 * c0 - n0 * c1) / det;
    let b = (n0 * q1 - n1 * q0) / det;
    Ok((a, b))
}

fn build_cofactor(den: &Poly, fac: &Factor, power: u32) -> Result<Poly, IntegrateError> {
    let fpoly = factor_to_poly(fac, power);
    let (cof, _) = poly_div(den.clone(), &fpoly)?;
    Ok(cof)
}

#[derive(Clone)]
enum Factor {
    Linear(Rational),
    Quadratic { b: Rational, c: Rational },
}

fn factor_to_poly(f: &Factor, mult: u32) -> Poly {
    let base = match f {
        Factor::Linear(a) => linear_poly(*a),
        Factor::Quadratic { b, c } => quadratic_poly(*b, *c),
    };
    if mult == 1 {
        base
    } else {
        poly_pow(&base, mult as i64)
    }
}

fn quadratic_poly(b: Rational, c: Rational) -> Poly {
    BTreeMap::from([(0, c), (1, b), (2, Rational::one())])
}

fn factor_polynomial(f: &Poly) -> Result<Vec<(Factor, u32)>, IntegrateError> {
    if poly_degree(f) <= 0 {
        return Ok(vec![]);
    }
    let mut f = f.clone();
    let mut result = Vec::new();
    for root in rational_roots(&f) {
        let mult = multiplicity(&f, root);
        result.push((Factor::Linear(root), mult));
        let lin = linear_poly(root);
        for _ in 0..mult {
            let (q, _) = poly_div(f, &lin)?;
            f = q;
        }
    }
    loop {
        let deg = poly_degree(&f);
        if deg == 0 {
            break;
        }
        if deg == 2 {
            let b = f.get(&1).copied().unwrap_or(Rational::zero());
            let c = f.get(&0).copied().unwrap_or(Rational::zero());
            let disc = b * b - Rational::from(4) * c;
            if disc < Rational::zero() {
                result.push((Factor::Quadratic { b, c }, 1));
            } else if let Some((r1, r2)) = quadratic_roots(b, c) {
                result.push((Factor::Linear(r1), 1));
                result.push((Factor::Linear(r2), 1));
            } else {
                return Err(IntegrateError::NoRule);
            }
            break;
        }
        let roots = rational_roots(&f);
        if roots.is_empty() {
            return Err(IntegrateError::NoRule);
        }
        for root in roots {
            let mult = multiplicity(&f, root);
            result.push((Factor::Linear(root), mult));
            let lin = linear_poly(root);
            for _ in 0..mult {
                let (q, _) = poly_div(f, &lin)?;
                f = q;
            }
        }
    }
    Ok(result)
}

fn integrate_polynomial(p: &Poly, var: Symbol) -> Result<Expr, IntegrateError> {
    let mut sum = const_(Rational::zero());
    for (&deg, coeff) in p.iter().rev() {
        if coeff.is_zero() {
            continue;
        }
        if deg == 0 {
            sum = sum + mul(const_(*coeff), Expr::var(var));
        } else {
            let np = deg + 1;
            sum = sum
                + mul(const_(*coeff), pow(Expr::var(var), const_(Rational::from(np))))
                    / const_(Rational::from(np));
        }
    }
    Ok(sum)
}

fn multiplicity(f: &Poly, root: Rational) -> u32 {
    let mut count = 0u32;
    let mut cur = f.clone();
    let lin = linear_poly(root);
    loop {
        let (_, r) = match poly_div(cur.clone(), &lin) {
            Ok(v) => v,
            Err(_) => break,
        };
        if poly_is_zero(&r) {
            count += 1;
            let (q, _) = poly_div(cur, &lin).unwrap();
            cur = q;
        } else {
            break;
        }
    }
    count.max(1)
}

fn rational_roots(f: &Poly) -> Vec<Rational> {
    let degree = poly_degree(f);
    if degree == 0 {
        return vec![];
    }
    let lc = *f.get(&degree).unwrap_or(&Rational::one());
    let ac = *f.get(&0).unwrap_or(&Rational::zero());
    let mut roots = Vec::new();
    for p in divisors_i64(lc.numer().abs()) {
        for q in divisors_i64(ac.denom().abs().max(1)) {
            for lp in divisors_i64(lc.denom().abs().max(1)) {
                for sign_p in [1i64, -1i64] {
                    for sign_q in [1i64, -1i64] {
                        let r = Rational::new(sign_p * p * sign_q, lp * q);
                        if poly_eval(f, r).is_zero() && !roots.contains(&r) {
                            roots.push(r);
                        }
                    }
                }
            }
        }
    }
    roots
}

fn quadratic_roots(b: Rational, c: Rational) -> Option<(Rational, Rational)> {
    let disc = b * b - Rational::from(4) * c;
    let s = rational_sqrt(disc)?;
    let two = Rational::from(2);
    Some(((-b + s) / two, (-b - s) / two))
}

fn rational_sqrt(r: Rational) -> Option<Rational> {
    if r < Rational::zero() {
        return None;
    }
    let n = isqrt_i64(r.numer())?;
    let d = isqrt_i64(r.denom())?;
    Some(Rational::new(n, d))
}

fn isqrt_i64(n: i64) -> Option<i64> {
    if n < 0 {
        return None;
    }
    if n == 0 {
        return Some(0);
    }
    let mut x = n;
    let mut y = (x + 1) / 2;
    while y < x {
        x = y;
        y = (x + n / x) / 2;
    }
    if x * x == n { Some(x) } else { None }
}

fn divisors_i64(n: i64) -> Vec<i64> {
    let n = n.abs().max(1);
    let mut out = Vec::new();
    let mut i = 1i64;
    while i * i <= n {
        if n % i == 0 {
            out.push(i);
            if i != n / i {
                out.push(n / i);
            }
        }
        i += 1;
    }
    out
}

fn poly_eval(p: &Poly, x: Rational) -> Rational {
    let mut acc = Rational::zero();
    for (&deg, c) in p {
        let mut pow_x = Rational::one();
        for _ in 0..deg {
            pow_x = pow_x * x;
        }
        acc = acc + *c * pow_x;
    }
    acc
}

fn poly_coeff(p: &Poly, deg: i64) -> Rational {
    p.get(&deg).copied().unwrap_or(Rational::zero())
}

fn linear_poly(root: Rational) -> Poly {
    BTreeMap::from([(0, -root), (1, Rational::one())])
}

fn poly_degree(p: &Poly) -> i64 {
    p.keys().max().copied().unwrap_or(0)
}

fn poly_is_zero(p: &Poly) -> bool {
    p.values().all(|c| c.is_zero())
}

fn poly_pow(base: &Poly, exp: i64) -> Poly {
    if exp == 0 {
        return BTreeMap::from([(0, Rational::one())]);
    }
    let mut r = base.clone();
    for _ in 1..exp {
        r = poly_mul(&r, base);
    }
    r
}

fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    let mut out = BTreeMap::new();
    for (&da, ca) in a {
        for (&db, cb) in b {
            *out.entry(da + db).or_insert(Rational::zero()) += *ca * *cb;
        }
    }
    normalize(&mut out);
    out
}

fn poly_add(a: Poly, b: Poly) -> Poly {
    let mut out = a;
    for (d, c) in b {
        *out.entry(d).or_insert(Rational::zero()) += c;
    }
    normalize(&mut out);
    out
}

fn poly_sub(a: Poly, b: Poly) -> Poly {
    let mut out = a;
    for (d, c) in b {
        *out.entry(d).or_insert(Rational::zero()) -= c;
    }
    normalize(&mut out);
    out
}

fn poly_scale(p: &Poly, s: Rational) -> Poly {
    p.iter().map(|(&d, c)| (d, *c * s)).collect()
}

fn poly_div(mut num: Poly, den: &Poly) -> Result<(Poly, Poly), IntegrateError> {
    let den_deg = poly_degree(den);
    let den_lc = *den.get(&den_deg).ok_or(IntegrateError::NoRule)?;
    if den_lc.is_zero() {
        return Err(IntegrateError::NoRule);
    }
    let mut quot = BTreeMap::new();
    normalize(&mut num);
    while poly_degree(&num) >= den_deg && !poly_is_zero(&num) {
        let nd = poly_degree(&num);
        let nc = *num.get(&nd).ok_or(IntegrateError::NoRule)?;
        let shift = nd - den_deg;
        let q = nc / den_lc;
        *quot.entry(shift).or_insert(Rational::zero()) += q;
        for (&d, &dc) in den {
            *num.entry(d + shift).or_insert(Rational::zero()) -= q * dc;
        }
        normalize(&mut num);
    }
    normalize(&mut quot);
    Ok((quot, num))
}

fn normalize(p: &mut Poly) {
    p.retain(|_, c| !c.is_zero());
}
