# Third-party material

## Solidcore

Solidcore, github.com/NVIDIA/solidcore, Copyright NVIDIA Corporation,
Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0), at commit
`e261b7d`:

- `examples/kernel_cases.rs` reads Solidcore's kernel-case corpus
  (`crates/solidcore/tests/kernel_cases`) and its portable B-rep record
  (`crates/solidcore-record/src/brep.rs`, whose documented parameterisations
  the reader follows). The corpus is not part of this repository:
  `scripts/fetch-kernel-cases.sh` unpacks it from a local Solidcore clone at
  the pinned commit.
