//! Shim for Solidcore's `error`: the subset of `SolidError` NLib builds and reads.

use core::fmt;

pub mod nurbs;

pub use nurbs::NurbsFailure;

/// Status the comparison hooks report for a failed NURBS routine.
pub const SM_ERR: i32 = 1001;

/// Result of a NLib routine.
pub type Result<T> = core::result::Result<T, SolidError>;

/// The class of an error.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ErrorKind {
    InvalidInput,
    Degenerate,
    NotConverged,
    Unsupported,
    Failed,
    Cancelled,
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

/// Why a [`Failure`] happened.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Reason {
    /// A NURBS routine failed.
    Nurbs(NurbsFailure),
    /// A query could not compute what it was asked.
    Query(&'static str),
    /// A geometry utility failed.
    Geom(&'static str),
}

/// A typed failure of an operation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Failure {
    pub operation: &'static str,
    pub reason: Reason,
    pub status: i32,
}

impl Failure {
    #[must_use]
    pub const fn new(operation: &'static str, failure: NurbsFailure) -> Self {
        Self {
            operation,
            reason: Reason::Nurbs(failure),
            status: SM_ERR,
        }
    }

    #[must_use]
    pub const fn into_error(self) -> SolidError {
        SolidError::Failed(self)
    }

    #[must_use]
    pub const fn kind(&self) -> ErrorKind {
        ErrorKind::Failed
    }
}

/// An error of a NLib routine.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum SolidError {
    /// An argument was out of range, malformed or inconsistent.
    InvalidInput {
        operation: &'static str,
        reason: &'static str,
    },
    /// The geometry is degenerate for the operation.
    Degenerate {
        operation: &'static str,
        entity: Option<&'static str>,
    },
    /// An iterative solver did not converge.
    NotConverged {
        operation: &'static str,
        iterations: u32,
    },
    /// The input is valid but not supported.
    Unsupported {
        operation: &'static str,
        reason: &'static str,
    },
    /// Cancelled through a [`crate::kernel::CancelToken`].
    Cancelled { operation: &'static str },
    /// An extension refused.
    ExtendRefused {
        operation: &'static str,
        refusal: crate::surface::extend::ExtendRefusal,
        distance: Option<f64>,
    },
    /// A typed failure.
    Failed(Failure),
}

impl SolidError {
    #[must_use]
    pub const fn invalid_input(operation: &'static str, reason: &'static str) -> Self {
        Self::InvalidInput { operation, reason }
    }

    #[must_use]
    pub const fn nurbs(operation: &'static str, failure: NurbsFailure) -> Self {
        Self::Failed(Failure::new(operation, failure))
    }

    #[must_use]
    pub const fn query_failed(operation: &'static str, what: &'static str) -> Self {
        Self::Failed(Failure {
            operation,
            reason: Reason::Query(what),
            status: SM_ERR,
        })
    }

    #[must_use]
    pub const fn geom_failed(operation: &'static str, what: &'static str) -> Self {
        Self::Failed(Failure {
            operation,
            reason: Reason::Geom(what),
            status: SM_ERR,
        })
    }

    /// The NURBS failure behind this error, if it is one.
    #[must_use]
    pub fn nurbs_failure(&self) -> Option<&NurbsFailure> {
        match self {
            Self::Failed(Failure {
                reason: Reason::Nurbs(n),
                ..
            }) => Some(n),
            _ => None,
        }
    }

    /// True for the failure of a NURBS routine.
    #[must_use]
    pub fn is_nurbs_failure(&self) -> bool {
        self.nurbs_failure().is_some()
    }

    /// True for a port-side invariant failure (none here).
    #[must_use]
    pub const fn is_port_failure(&self) -> bool {
        false
    }

    /// True for an error that answers with an upstream status.
    #[must_use]
    pub const fn is_status_answer(&self) -> bool {
        self.sm_status().is_some()
    }

    /// The upstream status of a typed failure.
    #[must_use]
    pub const fn sm_status(&self) -> Option<i32> {
        match self {
            Self::Failed(f) => Some(f.status),
            _ => None,
        }
    }

    #[must_use]
    pub const fn kind(&self) -> ErrorKind {
        match self {
            Self::InvalidInput { .. } => ErrorKind::InvalidInput,
            Self::Degenerate { .. } => ErrorKind::Degenerate,
            Self::NotConverged { .. } => ErrorKind::NotConverged,
            Self::Unsupported { .. } | Self::ExtendRefused { .. } => ErrorKind::Unsupported,
            Self::Cancelled { .. } => ErrorKind::Cancelled,
            Self::Failed(f) => f.kind(),
        }
    }

    #[must_use]
    pub const fn operation(&self) -> &'static str {
        match self {
            Self::InvalidInput { operation, .. }
            | Self::Degenerate { operation, .. }
            | Self::NotConverged { operation, .. }
            | Self::Unsupported { operation, .. }
            | Self::Cancelled { operation }
            | Self::ExtendRefused { operation, .. } => operation,
            Self::Failed(f) => f.operation,
        }
    }
}

impl fmt::Display for SolidError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput { operation, reason } => {
                write!(f, "{operation}: invalid input: {reason}")
            }
            Self::Degenerate { operation, .. } => write!(f, "{operation}: degenerate geometry"),
            Self::NotConverged {
                operation,
                iterations,
            } => write!(f, "{operation}: did not converge after {iterations} iterations"),
            Self::Unsupported { operation, reason } => {
                write!(f, "{operation}: unsupported: {reason}")
            }
            Self::Cancelled { operation } => write!(f, "{operation}: cancelled"),
            Self::ExtendRefused {
                operation, refusal, ..
            } => write!(f, "{operation}: extension refused ({refusal:?})"),
            Self::Failed(x) => match x.reason {
                Reason::Nurbs(n) => write!(f, "{}: {n}", x.operation),
                Reason::Query(w) | Reason::Geom(w) => write!(f, "{}: {w}", x.operation),
            },
        }
    }
}

impl std::error::Error for SolidError {}
