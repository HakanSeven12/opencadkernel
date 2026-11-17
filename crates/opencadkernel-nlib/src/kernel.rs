//! Shim for Solidcore's `kernel`: cooperative cancellation and the opt-in upstream fixes
//! NLib reads.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use crate::error::{Result, SolidError};

/// A shared cancellation flag: cloning shares it.
#[derive(Debug, Clone, Default)]
pub struct CancelToken {
    flag: Arc<AtomicBool>,
}

impl CancelToken {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Signals cancellation to every clone of this token.
    pub fn cancel(&self) {
        self.flag.store(true, Ordering::Relaxed);
    }

    /// Clears the signal.
    pub fn reset(&self) {
        self.flag.store(false, Ordering::Relaxed);
    }

    #[must_use]
    pub fn is_cancelled(&self) -> bool {
        self.flag.load(Ordering::Relaxed)
    }

    /// `Err(SolidError::Cancelled)` once cancelled.
    pub fn check(&self, operation: &'static str) -> Result<()> {
        if self.is_cancelled() {
            return Err(SolidError::Cancelled { operation });
        }
        Ok(())
    }
}

/// Opt-in corrections of upstream behaviour (all off by default).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct UpstreamFixes {
    /// At a closed smooth surface's pole, compare the pole's tangent limits.
    pub closed_smooth_pole_tangent_limit: bool,
}
