//! Solidcore's NLib layer: NURBS bases, knots, curve and surface evaluation, fitting,
//! approximation, knot removal, degree change, extension, offsets, deviation, matching
//! and filling.
//!
//! `nurbs/`, `limits.rs` and `geom/vector{2,3}d.rs` are Solidcore's files taken as they
//! are (see `README.md`). The other modules are a small shim standing in for the parts of
//! Solidcore those files reach outside NLib for (its error type, constants, cancellation
//! token), so a later Solidcore can be copied over the same files.

#![forbid(unsafe_code)]
#![allow(
    dead_code,
    unused_imports,
    unused_variables,
    unused_mut,
    unused_assignments,
    unused_qualifications,
    unreachable_patterns,
    rustdoc::broken_intra_doc_links,
    rustdoc::private_intra_doc_links,
    clippy::all
)]

pub mod error;
pub mod geom;
pub mod kernel;
pub mod limits;
pub mod math;
pub mod nurbs;
pub mod surface;

pub use error::{Result, SolidError};
