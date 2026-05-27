//! Operator implementations between `Symbol`/`Rational` and `Expr`.

use crate::expr::{add, const_, div, mul, Expr};
use crate::rational::Rational;
use crate::symbol::Symbol;
use std::ops::{Add, Div, Mul, Neg, Sub};

macro_rules! impl_int_expr_mul_add {
    ($($ty:ty),* $(,)?) => {
        $(
            impl Mul<Expr> for $ty {
                type Output = Expr;
                fn mul(self, rhs: Expr) -> Expr {
                    const_(Rational::try_from(self).expect("integer fits in rational"))
                        * rhs
                }
            }
            impl Add<Expr> for $ty {
                type Output = Expr;
                fn add(self, rhs: Expr) -> Expr {
                    const_(Rational::try_from(self).expect("integer fits in rational"))
                        + rhs
                }
            }
            impl Mul<$ty> for Expr {
                type Output = Expr;
                fn mul(self, rhs: $ty) -> Expr {
                    self * const_(Rational::try_from(rhs).expect("integer fits in rational"))
                }
            }
            impl Add<$ty> for Expr {
                type Output = Expr;
                fn add(self, rhs: $ty) -> Expr {
                    self + const_(Rational::try_from(rhs).expect("integer fits in rational"))
                }
            }
            impl Mul<Symbol> for $ty {
                type Output = Expr;
                fn mul(self, rhs: Symbol) -> Expr {
                    const_(Rational::try_from(self).expect("integer fits in rational")) * rhs
                }
            }
        )*
    };
}

impl_int_expr_mul_add!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize,
);

impl Symbol {
    pub fn to_expr(self) -> Expr {
        Expr::from(self)
    }

    pub fn pow(self, exp: impl Into<Expr>) -> Expr {
        self.to_expr().pow(exp)
    }
}

impl Add<Expr> for Symbol {
    type Output = Expr;
    fn add(self, rhs: Expr) -> Expr {
        self.to_expr() + rhs
    }
}

impl Add<Symbol> for Symbol {
    type Output = Expr;
    fn add(self, rhs: Symbol) -> Expr {
        self.to_expr() + rhs.to_expr()
    }
}

impl Add<Rational> for Symbol {
    type Output = Expr;
    fn add(self, rhs: Rational) -> Expr {
        self.to_expr() + rhs
    }
}

impl Mul<Expr> for Symbol {
    type Output = Expr;
    fn mul(self, rhs: Expr) -> Expr {
        self.to_expr() * rhs
    }
}

impl Mul<Symbol> for Symbol {
    type Output = Expr;
    fn mul(self, rhs: Symbol) -> Expr {
        self.to_expr() * rhs.to_expr()
    }
}

impl Mul<Rational> for Symbol {
    type Output = Expr;
    fn mul(self, rhs: Rational) -> Expr {
        self.to_expr() * rhs
    }
}

impl Mul<Symbol> for Rational {
    type Output = Expr;
    fn mul(self, rhs: Symbol) -> Expr {
        const_(self) * rhs
    }
}

impl Sub<Expr> for Symbol {
    type Output = Expr;
    fn sub(self, rhs: Expr) -> Expr {
        self.to_expr() - rhs
    }
}

impl Neg for Symbol {
    type Output = Expr;
    fn neg(self) -> Expr {
        -self.to_expr()
    }
}

impl Add<Symbol> for Rational {
    type Output = Expr;
    fn add(self, rhs: Symbol) -> Expr {
        const_(self) + rhs
    }
}

impl Add<Symbol> for Expr {
    type Output = Expr;
    fn add(self, rhs: Symbol) -> Expr {
        self + rhs.to_expr()
    }
}

impl Mul<Symbol> for Expr {
    type Output = Expr;
    fn mul(self, rhs: Symbol) -> Expr {
        self * rhs.to_expr()
    }
}

impl Add<Rational> for Expr {
    type Output = Expr;
    fn add(self, rhs: Rational) -> Expr {
        add(self, const_(rhs))
    }
}

impl Mul<Rational> for Expr {
    type Output = Expr;
    fn mul(self, rhs: Rational) -> Expr {
        mul(self, const_(rhs))
    }
}

impl Mul<Expr> for Rational {
    type Output = Expr;
    fn mul(self, rhs: Expr) -> Expr {
        const_(self) * rhs
    }
}

impl Div<Rational> for Expr {
    type Output = Expr;
    fn div(self, rhs: Rational) -> Expr {
        div(self, const_(rhs))
    }
}
