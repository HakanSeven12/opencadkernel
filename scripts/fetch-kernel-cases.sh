#!/bin/sh
# Unpacks Solidcore's kernel-case corpus (Apache-2.0, github.com/NVIDIA/solidcore,
# crates/solidcore/tests/kernel_cases) at a pinned commit, for examples/kernel_cases.rs.
# The corpus is not vendored here: it is read from a local Solidcore clone.
#
#   scripts/fetch-kernel-cases.sh <solidcore clone> [out dir, default target/kernel_cases]
set -eu
PIN=e261b7d
src=${1:?usage: fetch-kernel-cases.sh <solidcore clone> [out dir]}
out=${2:-target/kernel_cases}
git -C "$src" cat-file -e "$PIN^{commit}" || {
    echo "commit $PIN is not in $src: fetch it first" >&2
    exit 1
}
rm -rf "$out"
mkdir -p "$out"
git -C "$src" archive "$PIN" crates/solidcore/tests/kernel_cases |
    tar -x -C "$out" --strip-components=4
find "$out" -name '*.jpg' -delete
find "$out" -name '*.bodies.json.gz' -exec gzip -d -f {} +
echo "$(find "$out" -mindepth 2 -name '*.json' ! -name '*.bodies.json' | wc -l) cases in $out"
