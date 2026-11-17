# opencadkernel-nlib

NURBS curve and surface routines: bases and knot vectors, evaluation and
derivatives, curvature, interpolation and least-squares fitting, knot
removal and degree change within a tolerance, extension, offsets, deviation
bounds, edge matching, Coons/skinned/n-sided fills.

`src/nurbs/` is the NLib part of the `nurbs` module of
[Solidcore](https://github.com/tangency-labs/solidcore) (`crates/solidcore`,
commit `e261b7d`), itself a Rust port of NVIDIA's SMLib (USD BRep v1.0.1),
taken under the Apache License 2.0 (see `LICENSE` and `NOTICE`). The files
taken as they are:

- `src/nurbs/`: `geom, math, knots, basis_basic, basis_adv, bezier, index,
  fit_arrays, funcs_basic, funcs_adv, crv_basic, crv_nurb, crv_geom,
  crv_tool/, crv_adv/, crv_fit/, crv_approx, crv_offset, crv_common,
  crv_deviation, crv_proportion, srf_basic, srf_nurb, srf_tool/, srf_fit,
  srf_common, srf_adv/, srf_offset, srf_match/, srf_deviation, srf_seam`
- `src/limits.rs`, `src/geom/polynomial.rs`, `src/geom/vector2d.rs`,
  `src/geom/vector3d.rs`
- `src/error/nurbs.rs` holds `NurbsFailure` and its `Display` from Solidcore's
  `error/nurbs.rs`.

The rest (`lib.rs`, `nurbs.rs`, `error.rs`, `math.rs`, `kernel.rs`,
`geom.rs`, `surface.rs`) is a shim standing in for the parts of Solidcore
those files reach for: its error type, constants, cancellation token,
`reverse_knots`, `above`, `ExtendRefusal`. The only edits to the taken files
are two tests gated off with `#[cfg(any())]` (`srf_adv/p4.rs`
`the_upright_rate_is_the_frames_derivative`, `srf_adv/p6.rs`
`only_a_fixed_profile_end_skins_at_the_path_parameter`): they need
Solidcore's B-spline wrappers, which are not taken.

To take a later Solidcore, copy the same files over and rebuild; extend the
shim where it no longer covers what they use.

opencadkernel reaches it through `brep::nlib` (behind the `brep` feature).
