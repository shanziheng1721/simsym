//! Wide coefficients behind the `bigint` feature.

#[cfg(feature = "bigint")]
pub use num_rational::Ratio;
#[cfg(feature = "bigint")]
pub use num_bigint::BigInt;

/// Rational with arbitrary-size integers (`Ratio<BigInt>`).
#[cfg(feature = "bigint")]
pub type BigRational = Ratio<BigInt>;

#[cfg(feature = "bigint")]
pub fn big_rational(num: impl Into<BigInt>, den: impl Into<BigInt>) -> BigRational {
    Ratio::new(num.into(), den.into())
}
