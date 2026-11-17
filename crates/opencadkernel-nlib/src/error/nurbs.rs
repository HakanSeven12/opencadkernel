//! Why a NURBS routine fails (Solidcore's `error/nurbs.rs`, the `NurbsFailure` part).

use super::SolidError;

/// Why a NURBS routine fails.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum NurbsFailure {
    /// The curve's control point, knot and degree counts do not agree (`n + p + 1 != m`).
    CurveDefinition,
    /// The surface's control net, knot and degree counts do not agree in a direction.
    SurfaceDefinition,
    /// The degree is not allowed here (too low, too high, or not one the routine
    /// handles).
    Degree,
    /// The degree `given` is above the largest `max` the routines handle.
    DegreeTooHigh {
        /// The degree asked for.
        given: usize,
        /// The largest degree.
        max: usize,
    },
    /// The knot vector is not valid here.
    Knots,
    /// The knot vector decreases at knot `index` (it is smaller than the knot before it).
    KnotsDecrease {
        /// The index of the knot that is smaller than its predecessor.
        index: usize,
    },
    /// A weight is not valid here.
    Weight,
    /// The weight of control point `index` is `weight`, outside the range `[min, max]`
    /// the routines accept.
    WeightOutOfRange {
        /// The control point.
        index: usize,
        /// Its weight.
        weight: f64,
        /// The smallest weight accepted.
        min: f64,
        /// The largest weight accepted.
        max: f64,
    },
    /// A parameter is outside the domain.
    Parameter,
    /// The parameter is outside the curve or surface's domain `[lo, hi]`.
    ParameterOutside {
        /// The parameter.
        parameter: f64,
        /// The start of the domain.
        lo: f64,
        /// The end of the domain.
        hi: f64,
    },
    /// An index is out of range.
    Index,
    /// The index `index` is not below the length `len`.
    IndexOutOfRange {
        /// The index.
        index: usize,
        /// The number of items.
        len: usize,
    },
    /// A derivative is zero where a direction is needed.
    ZeroDerivative,
    /// A buffer supplied by the caller is too small for the result.
    Storage,
    /// A numerical failure (a division by zero, an overflow, a value that is not finite).
    Numerical,
    /// An iteration did not converge.
    NotConverged,
    /// A routine it calls rejected its arguments, or the direction is not one the
    /// routine handles.
    CallFailed,
    /// The input data is not valid.
    Input,
    /// The memory for the result cannot be allocated.
    Memory,
    /// The system of equations has no unique solution (it is singular).
    Singular,
    /// The geometry does not allow the operation (a degenerate or non-planar input).
    Geometry,
    /// The output has no row for even the point itself, so no derivative order can be
    /// evaluated into it.
    NoDerivativeRows,
    /// The derivative order `given` is above the largest `max` the routines evaluate.
    DerivativeOrder {
        /// The order asked for.
        given: usize,
        /// The largest order.
        max: usize,
    },
}

impl NurbsFailure {
    /// The `NL_*_ERR` number the C++ stores for the failure (a comparison detail).
    #[must_use]
    pub const fn nlib_number(&self) -> Option<i32> {
        Some(match self {
            Self::CurveDefinition => 1,
            Self::SurfaceDefinition => 2,
            Self::Degree | Self::DegreeTooHigh { .. } => 3,
            Self::Knots | Self::KnotsDecrease { .. } => 4,
            Self::Weight | Self::WeightOutOfRange { .. } => 5,
            Self::Parameter | Self::ParameterOutside { .. } => 6,
            Self::Index | Self::IndexOutOfRange { .. } => 7,
            Self::ZeroDerivative => 8,
            Self::Storage => 9,
            Self::Numerical => 11,
            Self::NotConverged => 12,
            Self::CallFailed => 13,
            Self::Input => 14,
            Self::Memory => 15,
            Self::Singular => 16,
            Self::Geometry => 17,
            Self::NoDerivativeRows | Self::DerivativeOrder { .. } => 20,
        })
    }
}

/// The `NL_*_ERR` number of a NURBS routine's failure, `None` for any other error.
#[must_use]
pub fn nlib_number(error: &SolidError) -> Option<i32> {
    error.nurbs_failure().and_then(NurbsFailure::nlib_number)
}


impl core::fmt::Display for NurbsFailure {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::CurveDefinition => {
                f.write_str("the curve's control point, knot and degree counts do not agree")
            }
            Self::SurfaceDefinition => {
                f.write_str("the surface's control net, knot and degree counts do not agree")
            }
            Self::Degree => f.write_str("the degree is not allowed here"),
            Self::DegreeTooHigh { given, max } => {
                write!(f, "degree {given} is above the largest degree, {max}")
            }
            Self::Knots => f.write_str("the knot vector is not valid here"),
            Self::KnotsDecrease { index } => {
                write!(f, "the knot vector decreases at knot {index}")
            }
            Self::Weight => f.write_str("a weight is not valid here"),
            Self::WeightOutOfRange {
                index,
                weight,
                min,
                max,
            } => write!(
                f,
                "the weight of control point {index} is {weight}, outside [{min}, {max}]"
            ),
            Self::Parameter => f.write_str("a parameter is outside the domain"),
            Self::ParameterOutside { parameter, lo, hi } => {
                write!(
                    f,
                    "parameter {parameter} is outside the domain [{lo}, {hi}]"
                )
            }
            Self::Index => f.write_str("an index is out of range"),
            Self::IndexOutOfRange { index, len } => {
                write!(f, "index {index} is out of range for {len} items")
            }
            Self::ZeroDerivative => f.write_str("a derivative is zero where a direction is needed"),
            Self::Storage => f.write_str("a buffer supplied by the caller is too small"),
            Self::Numerical => f.write_str("a numerical failure (a division by zero or overflow)"),
            Self::NotConverged => f.write_str("the iteration did not converge"),
            Self::CallFailed => f.write_str(
                "a routine it calls rejected its arguments, or the direction is not handled",
            ),
            Self::Input => f.write_str("the input data is not valid"),
            Self::Memory => f.write_str("the memory for the result cannot be allocated"),
            Self::Singular => f.write_str("the system of equations is singular"),
            Self::Geometry => f.write_str("the geometry does not allow the operation"),
            Self::NoDerivativeRows => f.write_str("the output has no row for the point itself"),
            Self::DerivativeOrder { given, max } => {
                write!(
                    f,
                    "derivative order {given} is above the largest order, {max}"
                )
            }
        }
    }
}
