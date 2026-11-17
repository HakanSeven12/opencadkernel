//! Shim for Solidcore's `geom`: its vector types (taken as they are) and the few helpers
//! NLib calls.

pub mod polynomial;
pub mod vector2d;
pub mod vector3d;

/// The tolerance tests the vector types use, at the default (unscaled) model size.
pub mod tol {
    use crate::math::{EFF_ZERO, EFF_ZERO_SQ};

    #[must_use]
    pub fn in_tol(val: f64, tol: f64) -> bool {
        val.abs() <= tol
    }

    #[must_use]
    pub fn is_zero(a: f64) -> bool {
        in_tol(a, EFF_ZERO)
    }

    #[must_use]
    pub fn is_zero_squared(a: f64) -> bool {
        in_tol(a, EFF_ZERO_SQ)
    }
}

pub mod certified {
    /// `a > b`, false for a NaN: a refusal written `!above(..)` refuses NaN too.
    #[must_use]
    pub fn above(a: f64, b: f64) -> bool {
        matches!(a.partial_cmp(&b), Some(core::cmp::Ordering::Greater))
    }
}

pub mod a_object {
    use crate::error::{Result, SolidError};

    /// Reverses knot vector `u[0..=m]` into `v`, keeping its domain.
    pub fn reverse_knots(u: &[f64], m: usize, v: &mut [f64]) -> Result<()> {
        const OP: &str = "geom::a_object::reverse_knots";
        let short = || SolidError::invalid_input(OP, "fewer than m + 1 knots");
        if u.len() <= m || v.len() <= m {
            return Err(short());
        }
        let u0 = *u.first().ok_or_else(short)?;
        let um = *u.get(m).ok_or_else(short)?;
        let c = u0 + um;
        let rev = |x: f64| {
            if x == um {
                u0
            } else if x == u0 {
                um
            } else {
                c - x
            }
        };
        for i in 0..=m / 2 {
            let a = *u.get(i).ok_or_else(short)?;
            let b = *u.get(m - i).ok_or_else(short)?;
            *v.get_mut(i).ok_or_else(short)? = rev(b);
            *v.get_mut(m - i).ok_or_else(short)? = rev(a);
        }
        Ok(())
    }
}
