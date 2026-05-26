pub mod diff;
pub mod integrate;
pub mod multivar;
pub mod numeric;

pub use diff::diff;
pub use integrate::{integrate, IntegrateError};
pub use multivar::{gradient, hessian};
pub use numeric::{integrate_definite, integrate_numeric, DefiniteIntegralError, NumericOptions};
