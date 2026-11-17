# Third-party code

## Solidcore

[Solidcore](https://github.com/tangency-labs/solidcore) is a Rust port of
NVIDIA's SMLib (USD BRep v1.0.1). Copyright (c) 2026 The solidcore
contributors; the upstream work is Copyright (c) 2020-2026 NVIDIA
CORPORATION & AFFILIATES. Licensed under the Apache License, Version 2.0.
"SMLib" and "NVIDIA" identify the origin of the code only.

### Sewing and trim-loop repair of imported bodies

Ported from `crates/solidcore` at commit `e261b7d` onto the kernel's own
types (algorithms adapted, not copied file for file):

- `src/brep/heal.rs` (`sew`, `sew_into_solid`): the stitching of
  `src/heal/stitch/p2.rs` (`Stitch::do_stitching`, `stitch_into_solid`) as
  `docs/validity-and-tolerances.md` §7.3 describes it: vertices sorted along
  (1, 1, 1) and glued within the tolerance, edges split where a vertex lands
  inside them, coincident edges glued, the tolerance doubled while free edges
  remain up to 1/200 of the box diagonal.
- `src/brep/heal.rs` (`fix_loops`): `order_curves_into_loops`,
  `fix_trim_loop_order` and `fix_loop_orientation` of
  `src/topology/trimming_tools.rs`.
