# Third-party code

## Solidcore

[Solidcore](https://github.com/tangency-labs/solidcore) is a Rust port of
NVIDIA's SMLib (USD BRep v1.0.1). Copyright (c) 2026 The solidcore
contributors; the upstream work is Copyright (c) 2020-2026 NVIDIA
CORPORATION & AFFILIATES. Licensed under the Apache License, Version 2.0.
"SMLib" and "NVIDIA" identify the origin of the code only.

### Blends: corners and sizes that fit

Ported from `crates/solidcore` at commit `e261b7d` onto the kernel's own
types (algorithms adapted, not copied file for file):

- `src/brep/blend.rs` (`fillet_fit`, `chamfer_fit`, `band_room`): the room
  a round's band has along a straight edge between two planes, a radius `r`
  touching each face `r / tan(phi / 2)` from the edge, read at stations of
  the edge against the face's other boundary edges, from
  `src/blend/band_room.rs` (`band_room`, `EdgeBand::room_at`), without its
  interval certification.
- `src/brep/blend.rs` (`corner_blends`): the ball's centre at a corner of
  any three planes by Cramer's rule, from `planes_meet` of
  `src/blend/blend_corner/ball.rs`; the corner itself stays the kernel's
  analytic sphere rather than Solidcore's rational triangle patch.
