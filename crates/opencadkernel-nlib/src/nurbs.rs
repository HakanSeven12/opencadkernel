//! NURBS basis functions, knot vectors, and B-spline/NURBS curve and surface
//! evaluation, refinement, fitting and approximation (Solidcore's `nurbs` module,
//! NLib part only).

pub mod basis_adv;
pub mod basis_basic;
pub mod bezier;
pub mod crv_adv;
pub mod crv_approx;
pub mod crv_basic;
pub mod crv_common;
pub mod crv_deviation;
pub mod crv_fit;
pub mod crv_geom;
pub mod crv_nurb;
pub mod crv_offset;
pub mod crv_proportion;
pub mod crv_tool;
pub(crate) mod fit_arrays;
pub mod funcs_adv;
pub mod funcs_basic;
pub mod geom;
pub(crate) mod index;
pub mod knots;
pub mod math;
pub mod srf_adv;
pub mod srf_basic;
pub mod srf_common;
pub mod srf_deviation;
pub mod srf_fit;
pub mod srf_match;
pub mod srf_nurb;
pub mod srf_offset;
pub mod srf_seam;
pub mod srf_tool;
