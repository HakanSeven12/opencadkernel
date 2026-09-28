//! 2D curve algebra.
//!
//! This layer answers the questions the rest of the kernel asks constantly:
//! where do two curves meet, what does this loop offset to, is this point
//! inside, which fragments chain into a closed boundary. Editing commands
//! sit directly on it, and the B-rep layer reaches down into it whenever it
//! has to work in a face's parameter space.
//!
//! Everything here obeys the coordinate policy in [`frame`].

pub mod angle;
pub mod arclength;
pub mod area;
pub mod arrangement;
pub mod band;
mod boundary_recovery;
pub mod clip;
pub mod centerline;
pub mod containment;
pub mod construct;
#[cfg(feature = "brep")]
pub(crate) mod constrained;
pub mod constraint_inference;
pub mod cross;
pub mod curve;
pub mod deviation;
pub mod dimension;
pub mod fillet;
pub mod frame;
pub mod gradient;
pub mod intersect;
pub mod nurbs;
#[cfg(feature = "offset")]
pub mod offset;
pub mod polyline;
pub mod polyline_edit;
pub mod snap;
pub mod tessellate;
pub mod triangulate;
pub mod transform;
pub mod vec;

pub use angle::{angle_within_arc, arc_parameter, arc_span, normalize_angle};
pub use arrangement::{bounded_faces, segment_crossing, signed_area, SegmentCrossing};
pub use band::{
    polyline_band_boundary, BandBoundaryEdge, BandStationPiece, PolylineBandBoundary,
};
pub use boundary_recovery::refine_spline_boundary;
pub use clip::{break_spans, inside_pieces, inside_spans, trim_spans};
pub use centerline::{centerline_between, CenterLineGeometry};
pub use containment::{closest_point, contains, distance_to, nearest_of, Closest};
pub use constraint_inference::{
    infer_constraints, infer_constraints_with_settings, Endpoint as ConstraintEndpoint,
    InferenceKind, InferenceSettings, InferredConstraint, ParametricPrimitive,
};
pub use dimension::linear_dimension_leader;
pub use construct::{
    arc_from_endpoints_angle, arc_from_endpoints_radius, arc_from_sagitta,
    arc_from_start_tangent, arc_sweep_from_chord, arc_through_points, bounded_arc,
};
pub use cross::{intersect, Crossing};
pub use curve::{Arc, Circle, Curve, EllipseArc, Extent, Line, Ray, XLine};
pub use nurbs::{NurbsCurve, Parameterization};
pub use frame::Frame;
pub use fillet::{fillet_between_rays, fillets_between, Fillet};
pub use gradient::{gradient_frame, GradientFrame};
pub use intersect::{
    circle_circle_angles, circle_circle_points, ellipse_circle, line_circle, line_ellipse,
    line_line,
};
#[cfg(feature = "offset")]
pub use offset::offset_polyline;
pub use polyline::{
    BulgeArc, Polyline, PolylineRange, PolylineRangeSegment, PolylineVertex, RectangleFrame,
};
pub use polyline_edit::{
    move_polyline_segment_parallel, polyline_arc_radius_from_point,
    polyline_segment_parallel_offset, resize_polyline_arc,
};
pub use snap::{
    characteristic_points, nearest_to, perpendicular_from, tangent_from, SnapKind, SnapPoint,
};
pub use tessellate::{arc, ellipse_arc, lerp, DEFAULT_SEGMENTS_PER_RADIAN};
pub use transform::Transform;
pub use triangulate::{
    nesting_depths as ring_nesting_depths, polygon as triangulate, polygon_frame,
    rings as triangulate_rings, PolygonFrame,
};
pub use vec::Vec2;

/// An ellipse in the plane.
///
/// Grouped rather than passed as loose scalars because the four fields are
/// only meaningful together, and because a parameter is not interpretable
/// without the axis it is measured from — see [`Ellipse::point_at`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipse {
    /// Centre in world coordinates.
    pub centre: [f64; 2],
    /// Half-length along [`major_axis`](Self::major_axis).
    pub major_radius: f64,
    /// Half-length across it.
    pub minor_radius: f64,
    /// Unit vector along the major axis. The minor direction is this turned a
    /// quarter turn counter-clockwise.
    pub major_axis: [f64; 2],
}

impl Ellipse {
    /// The unit minor direction: [`major_axis`](Self::major_axis) rotated a
    /// quarter turn counter-clockwise.
    pub fn minor_axis(&self) -> [f64; 2] {
        [-self.major_axis[1], self.major_axis[0]]
    }

    /// The point at parameter `t`.
    ///
    /// `t` is the ellipse's own parameter, not an angle measured at the
    /// centre: the two agree only when the radii are equal.
    pub fn point_at(&self, t: f64) -> [f64; 2] {
        let along = self.major_radius * t.cos();
        let across = self.minor_radius * t.sin();
        let minor = self.minor_axis();
        [
            self.centre[0] + along * self.major_axis[0] + across * minor[0],
            self.centre[1] + along * self.major_axis[1] + across * minor[1],
        ]
    }

    /// Whether the radii are large enough to describe a curve at all.
    pub(crate) fn is_degenerate(&self) -> bool {
        self.major_radius.abs() < 1e-20 || self.minor_radius.abs() < 1e-20
    }
}

/// Distance below which two positions are treated as one.
///
/// Carried explicitly rather than read from a global, because a drawing in
/// millimetres and a drawing in survey feet do not want the same value, and
/// because a boolean needs to widen it locally without disturbing anything
/// else.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tolerance {
    linear: f64,
}

impl Tolerance {
    /// A tolerance of `linear` world units.
    ///
    /// # Panics
    ///
    /// If `linear` is not finite and positive. A non-positive tolerance makes
    /// every comparison in the kernel meaningless, so it is worth refusing at
    /// the boundary rather than producing degenerate topology later.
    pub fn new(linear: f64) -> Self {
        assert!(
            linear.is_finite() && linear > 0.0,
            "tolerance must be finite and positive, got {linear}"
        );
        Self { linear }
    }

    /// The linear tolerance in world units.
    pub const fn linear(&self) -> f64 {
        self.linear
    }

    /// Whether `a` and `b` are within tolerance of each other.
    pub fn are_equal(&self, a: f64, b: f64) -> bool {
        (a - b).abs() <= self.linear
    }
}

impl Default for Tolerance {
    /// 1e-7, matching the point-equality tolerance ACIS records in its
    /// header, so geometry surviving a round trip is judged by the same
    /// yardstick at both ends.
    fn default() -> Self {
        Self::new(1e-7)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equality_respects_the_configured_tolerance() {
        let tol = Tolerance::new(1e-3);
        assert!(tol.are_equal(1.0, 1.0005));
        assert!(!tol.are_equal(1.0, 1.002));
    }

    #[test]
    #[should_panic(expected = "finite and positive")]
    fn zero_tolerance_is_refused() {
        Tolerance::new(0.0);
    }
}

mod arc_fit;
pub use arc_fit::{fit_arc_chain, ArcFitVertex};
mod bounds;
pub use bounds::analytic_curve_bounds;
