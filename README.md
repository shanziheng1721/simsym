# simsym

A simple symbolic computation library for Rust, with exact rational arithmetic.

## Quick start

```rust
use simsym::prelude::*;

let x = symbol("x");
let y = symbol("y");

// Operators
let f = x.pow(2) + rational(2, 1) * x * y + sin(y);

// Macro DSL
let g = expr!(x^2 + 2*x*y + sin(y));

let df_dx = f.clone().diff(x);
let grad = f.gradient(&[x, y]);

let exact = f.eval(&[(x, rational(1, 2)), (y, rational(1, 3))])?;
let approx = f.eval_f64(&[(x, 0.5), (y, 1.0 / 3.0)])?;

let F = f.integrate(x)?; // symbolic when a rule applies
let area = f.integrate_definite(x, rational(0, 1), rational(1, 1), &[])?;

// Numeric fallback
let num = integrate_numeric(&sin(x), x, 0.0, std::f64::consts::PI, &[], NumericOptions::default())?;
```

## Features

- Expression AST with `+ - * / ^`, full circular & hyperbolic trig (`sin`…`acsch`), plus `exp`, `ln`
- Exact `Rational` constants (`i64` ratios)
- Simplification (algebraic folding, some like-term merging)
- Partial derivatives and gradients / Hessians
- Symbolic integration for polynomials and common elementary forms
- Definite integrals (symbolic antiderivative or adaptive Simpson)
- `expr!` procedural macro

## Limitations

Symbolic integration is **rule-based** (polynomials, affine trig/exp, `e^x×P`, parts, `sin^n`/`cos^n` reduction, partial fractions for low-degree `P/Q` including `atan` terms, `sin·cos` products, etc.) — not a full Risch algorithm. See [docs/INTEGRATION.md](docs/INTEGRATION.md) for the algorithm stack. When `integrate` returns `IntegrateError::NoRule`, use `integrate_numeric` or `integrate_definite` (which falls back automatically).

Transcendental functions are not supported in exact `eval`; use `eval_f64`.

## Examples

```bash
cargo run --example calculus
cargo run --example multivar
cargo run --example numeric_integral
```

## Optional features

Cargo features (enabled by default):

| Feature | Enables |
|---------|---------|
| `simplify` | [`Expr::simplify`](https://docs.rs/simsym/latest/simsym/struct.Expr.html#method.simplify) |
| `diff` | [`Expr::diff`](https://docs.rs/simsym/latest/simsym/struct.Expr.html#method.diff), `gradient`, `hessian` |
| `integrate` | [`Expr::integrate`](https://docs.rs/simsym/latest/simsym/struct.Expr.html#method.integrate), `integrate_definite` (implies `diff`) |

Minimal build (AST + rationals + evaluation only):

```bash
cargo build --no-default-features
```

Other optional features:

- `serde` — serialize rationals as `(numer, denom)`
- `bigint` — `BigRational` alias for wider coefficients (see `simsym::rational_big`)
