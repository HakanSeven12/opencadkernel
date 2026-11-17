//! Shim for Solidcore's `math`: the constants and helpers NLib uses.

pub use crate::geom::vector2d::{Point2, Vec2};
pub use crate::geom::vector3d::{Point3, Vec3};

pub const PI: f64 = core::f64::consts::PI;
pub const HALF_PI: f64 = core::f64::consts::FRAC_PI_2;
pub const TWO_PI: f64 = core::f64::consts::TAU;
pub const EFF_ZERO: f64 = 1.0e-12;
pub const EFF_ZERO_SQ: f64 = 1.0e-24;
pub const EFF_ZERO_SQRT: f64 = 1.0e-6;
pub const EFF_ZERO_PARAM: f64 = 1.2e-8;
pub const EFF_ZERO_DEG: f64 = 1.0e-4;
pub const UNDEF_DOUBLE: f64 = 1.23456e-61;
pub const BIG_DOUBLE: f64 = 1.0e20;

#[must_use]
pub fn max2(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}

#[must_use]
pub fn min2(a: f64, b: f64) -> f64 {
    if a < b {
        a
    } else {
        b
    }
}
