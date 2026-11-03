# opencadkernel-ssi

Surface/surface intersection of rational B-spline patches: transversal,
tangent and singular intersections, boundaries, coincident regions,
near-degenerate pairs, seams and poles.

This crate is the `ssi` crate of [Solidcore](https://github.com/tangency-labs/solidcore)
(`crates/ssi`, commit `845bdb4`), taken as it is under its MIT OR Apache-2.0
license (see `LICENSE-MIT` and `LICENSE-APACHE`). Its method and its record on
[jyviko's SSI corpus](https://github.com/jyviko/ssi-corpus) are described in
the upstream README.

opencadkernel reaches it through `brep::ssi`, which turns the kernel's
surfaces into patches and the results back into the kernel's parameters.
Keep changes here few and upstreamable, so a later Solidcore commit can be
taken over the same way.

## Score on the corpus

```sh
cargo build --release -p opencadkernel-ssi
SSI_ENGINE=ssi-rs SSI_RUNNER=target/release/ssi_corpus_runner \
    python3 -I <ssi-corpus>/scoring/score.py fs 1e-6
```
