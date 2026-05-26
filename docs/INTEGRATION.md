# Symbolic integration in simsym

## How production CAS systems integrate

| Approach | Idea | Completeness | Typical use |
|----------|------|--------------|-------------|
| **Risch algorithm** | Decision procedure in differential fields (Bronstein) | Complete for elementary integrals in theory; implementations often incomplete on algebraic extensions | FriCAS, partial Maple/Mathematica |
| **Risch–Norman** | Ansatz: guess antiderivative template, solve for coefficients | Fast heuristic, incomplete | SymPy `heurisch`, many systems as first pass |
| **Rule databases (RUBI)** | ~6700 pattern rules | High coverage for definite/univariate integrals, not a decision procedure | Mathematica packages |
| **Rule + fallback** | Polynomial/rational + parts + substitution, then numeric | Practical | **simsym** (this crate) |

simsym intentionally stays in the last row: a **small, maintainable rule stack** with **numeric quadrature** as fallback (`integrate_numeric`, `integrate_definite`).

## simsym evaluation order

1. **U-substitution** (limited): `∫ u^n u' dx`
2. **Structural dispatch** on `ExprKind`
3. **Simplify** the result

### Rules implemented

- **Polynomials** in `x` and `(x+c)^n`, `n ≠ -1`
- **Logarithmic**: `∫ ln(x) dx`, `∫ ln(kx+b) dx`
- **Trigonometric**: `sin`, `cos`, `tan` with affine inner `kx+b`
- **Inverse trig**: `∫ atan(x) dx`, `∫ atan(kx+b) dx`
- **Exponential**: `e^x`, `e^{kx+b}`, `e^x × polynomial` (template / Risch-Norman style)
- **Products**: `sin(ax+b)cos(cx+d)` via product-to-sum; `e^{ax+b} sin(cx+d)` / `e^{ax+b} cos(cx+d)`; `sin^m cos^n` (same affine inner); `x^n ln(x)`
- **Powers**: `sin^n`, `cos^n`, `tan^n`, `sec^n` (`cos^{-n}`, even `n ≥ 2`) via standard reductions
- **Substitution**: `u^n u'`, and `f(u)·u'` for `sin/cos/tan/exp/ln/atan` (including `cos⁻² u' → tan u`)
- **Rational (heuristic)**:
  - `c/(x-a)`, `1/((x-a)(x-b))`, `(px+q)/(x-a)`
  - **Polynomial `P/Q`** with `deg(Q) ≤ 6`: long division + partial fractions
  - **`∫ 1/(x²±a²)`** and irreducible quadratics → `atan` / `ln` terms
- **`atan`** in the AST (differentiation, `eval_f64`, display, integration)
- **Integration by parts** (depth ≤ 4), with `∫ v du` direct attempt

### Not implemented

- Full **Risch** decision procedure
- **Hermite–Rothstein–Trager** in full generality (only the partial-fraction / `atan` slice above for low-degree denominators)
- Arbitrary-degree denominators or algebraic extensions
- **Non-elementary** integrals (e.g. `∫ e^{x²} dx`) — use numeric integration

## References

- Manuel Bronstein, *Symbolic Integration I: Transcendental Functions*
- Robert Risch (1969), decision procedure for elementary integration
- Geddes et al., survey of integration in computer algebra
- RUBI rule-based integration (Rich, 2018)
