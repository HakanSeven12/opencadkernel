# Third-party code

## Solidcore

[Solidcore](https://github.com/tangency-labs/solidcore) is a Rust port of
NVIDIA's SMLib (USD BRep v1.0.1). Copyright (c) 2026 The solidcore
contributors; the upstream work is Copyright (c) 2020-2026 NVIDIA
CORPORATION & AFFILIATES. Licensed under the Apache License, Version 2.0.
"SMLib" and "NVIDIA" identify the origin of the code only.

### Booleans

Ideas and algorithms from `crates/solidcore` at commit `e261b7d`, rewritten
onto the kernel's own types:

- `src/brep/boolean.rs` (`local_side`): the radial-sector classification of
  a face at an edge it shares with the other body, from
  `src/query/edge.rs` (`radial_sector_classify_between`) as driven by
  `src/boolean/merge/p4.rs` (`classify_edge_faces`), reduced to the smooth
  sector; within ten degrees of tangency it leaves the face to ray casting.
- `src/brep/boolean.rs` (`combine_all`): one Boolean against a list of tools
  combined into a single operand first, after `src/boolean/merge/p3.rs`
  (`boolean_lists`, `merge_breps`).
- `src/brep/imprint.rs` (`tangent_cap`): a crossing curve bounding a cap
  within the tolerance of both surfaces is a touch, from
  `src/topology/face/p3.rs` (`cap_centre_within_tolerance`).
- `src/brep/ssi.rs` (`apart`): proving two faces apart on the convex hulls of
  their surfaces' Bézier pieces, nearest pair first with de Casteljau
  subdivision, after `src/query/distance_bound.rs` and `src/query/contact.rs`.
