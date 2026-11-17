# Third-party code

## Solidcore

[Solidcore](https://github.com/tangency-labs/solidcore) is a Rust port of
NVIDIA's SMLib (USD BRep v1.0.1). Copyright (c) 2026 The solidcore
contributors; the upstream work is Copyright (c) 2020-2026 NVIDIA
CORPORATION & AFFILIATES. Licensed under the Apache License, Version 2.0.
"SMLib" and "NVIDIA" identify the origin of the code only.

### NLib NURBS layer (`crates/opencadkernel-nlib`)

Taken from `crates/solidcore` at commit `e261b7d`, with its `LICENSE` and
`NOTICE` alongside (see the crate's `README.md` for the exact file list and
the shim around it):

- `src/nurbs/{geom, math, knots, basis_basic, basis_adv, bezier, index,
  fit_arrays, funcs_basic, funcs_adv, crv_basic, crv_nurb, crv_geom,
  crv_tool/, crv_adv/, crv_fit/, crv_approx, crv_offset, crv_common,
  crv_deviation, crv_proportion, srf_basic, srf_nurb, srf_tool/, srf_fit,
  srf_common, srf_adv/, srf_offset, srf_match/, srf_deviation, srf_seam}.rs`
- `src/limits.rs`, `src/geom/{polynomial, vector2d, vector3d}.rs`
- `src/error/nurbs.rs` (`NurbsFailure`)

Ported onto the kernel's own types:

- `src/space/nurbs_tools.rs`: the surface curvatures from
  `src/surface/surface/p2.rs` (`evaluate_geometric`); the rest wraps the
  vendored routines.
- `src/brep/nurbs_form.rs`: the quarter-grid scheme of exact NURBS forms
  from `src/eval/form.rs`.
- `src/brep/geometry.rs` (torus ray hits): the closed-form quartic of
  `src/geom/polynomial.rs`.

### Mass properties, validity, draft and certified offsets

Ported from `crates/solidcore` at commit `e261b7d` onto the kernel's own
types (algorithms and constants adapted, not copied file for file):

- `src/brep/mass.rs` (`integrated_mass_properties`, `kronrod`): the
  Gauss–Kronrod 7/15 tables, the adaptive quadrature and the trim integral
  (volume integrals as boundary fluxes, Green's theorem into iso-line
  integrals split at knots) of `src/mass/integrator/properties.rs`
  (`kronrod_integrate_with`, `trim_integrand`).
- `src/brep/audit.rs`: the validity tests of `docs/validity-and-tolerances.md`
  §4 and `src/topology/brep/validation/assert_valid/` (edge ends on its
  vertices, edge on its faces, pcurves, loop nesting).
- `src/brep/edit.rs` (`draft_faces`): `plane_draft` and `cone_draft` of
  `src/offset/draft.rs`.
- `src/space/bernstein.rs`: `Ivl` of `src/blend/fillet_section/interval.rs`
  and `Bernstein2` of `src/offset/construction/bernstein.rs`.
- `src/brep/certified_offset.rs`: `certified_offset`, `certify_regular`,
  `largest_wall` and their forms from `src/offset/construction/sheets.rs`.
- `src/brep/wall_room.rs`: `shell_wall_room` of `src/offset/wall_room.rs`.
