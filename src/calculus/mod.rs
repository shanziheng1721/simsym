//! Symbolic and numeric calculus (enabled via Cargo features).

#[cfg(feature = "diff")]
pub mod diff;

#[cfg(feature = "integrate")]
pub mod integrate;

#[cfg(feature = "diff")]
pub mod multivar;

pub mod numeric;

#[cfg(feature = "diff")]
pub use diff::{diff, diff_without_simplify};

#[cfg(feature = "integrate")]
pub use integrate::{integrate, IntegrateError};

#[cfg(feature = "diff")]
pub use multivar::{gradient, hessian};

pub use numeric::{integrate_numeric, DefiniteIntegralError, NumericOptions};

#[cfg(feature = "integrate")]
pub use numeric::integrate_definite;
