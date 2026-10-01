//! Cutting each solid's faces where the other's pass through them.
//!
//! The step between finding the intersection curves and deciding what to
//! keep. Afterwards, every face of either solid is wholly inside the other or
//! wholly outside it — never partly both — so classifying a piece is one
//! question with one answer instead of a face that would have to be described
//! as "inside over here".
//!
//! Both solids are imprinted, not one. A union has to keep the outer part of
//! *each*, so each needs the other's curves cut into it.
//!
//! # Why it can fail, and why that is reported
//!
//! A face pair whose intersection has no closed form, a pair of coincident
//! faces, a cut this kernel cannot make: each leaves the imprint incomplete
//! in a way the caller cannot see by looking at the result. A boolean run on
//! a half-imprinted body produces a solid — one with a wall missing. So the
//! failures are returned, and a boolean refuses on them.

use super::bounds::{face_bounds, Aabb};
use super::geometry::{Curve3, Surface};
use super::intersect::{surfaces, Meeting};
use super::split::{split_edge, split_face};
use super::topology::{Body, FaceKey};
use crate::space::Vec3;

/// Why an imprint could not be completed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Snag {
    /// Two faces meet along something this kernel has no closed form for.
    NoClosedForm,
    /// Two faces lie on the same surface. Deciding what that leaves needs
    /// their overlap worked out in parameter space, which is its own
    /// operation.
    Coincident,
    /// The curves were found but a face could not be cut along one of them —
    /// a cut crossing the boundary more than twice, or closing inside the
    /// face.
    CutRefused,
}

/// What an imprint did.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Imprint {
    /// How many faces were cut in two, across both solids.
    pub cuts: usize,
    /// How many face pairs were found to meet at all.
    pub meetings: usize,
}

/// Cuts each body's faces along the curves it shares with the other.
///
/// Neither body's shape changes — an imprint only adds edges — so the result
/// still passes [`Body::validate`] and still has the same volume.
///
/// `tolerance` is passed through to the intersection and the cutting.
pub fn imprint(a: &mut Body, b: &mut Body, tolerance: f64) -> Result<Imprint, Snag> {
    let meetings = shared_curves(a, b, tolerance)?;
    let count = meetings.len();
    // Face pairs on the same two surfaces all report the same curve — every
    // part of a divided sphere meets a plane in one circle — and each copy
    // would be tried against every face again for nothing.
    let mut curves: Vec<Curve3> = Vec::new();
    for curve in meetings.iter().flat_map(|meeting| &meeting.curves) {
        if !curves.iter().any(|seen| same_curve(seen, curve, tolerance)) {
            curves.push(curve.clone());
        }
    }
    let mut cuts = cut_along(a, &curves, tolerance)?;
    cuts += cut_along(b, &curves, tolerance)?;
    // Stitching needs matching edge partitions on both bodies.
    align_edge_vertices(a, b, tolerance);
    align_edge_vertices(b, a, tolerance);
    Ok(Imprint {
        cuts,
        meetings: count,
    })
}

fn align_edge_vertices(source: &Body, target: &mut Body, tolerance: f64) {
    let points: Vec<[f64; 3]> = source.vertices.iter().map(|(_, vertex)| vertex.point).collect();
    for point in points {
        let edges: Vec<_> = target.edge_keys().collect();
        for key in edges {
            let Some(edge) = target.edges.get(key) else {
                continue;
            };
            if edge.coedges.is_empty() {
                continue;
            }
            let splittable = edge.coedges.iter().all(|coedge| {
                target.coedges.get(*coedge).is_some_and(|coedge| {
                    matches!(
                        coedge.pcurve.as_ref(),
                        None
                            | Some(&crate::geom2d::Curve::Line(_))
                            | Some(&crate::geom2d::Curve::Circle(_))
                            | Some(&crate::geom2d::Curve::Arc(_))
                            | Some(&crate::geom2d::Curve::Ellipse(_))
                            | Some(&crate::geom2d::Curve::Nurbs(_))
                    )
                })
            });
            if !splittable {
                continue;
            }
            let Some(curve) = target.curves.get(edge.curve) else {
                continue;
            };
            let mut parameter = curve.parameter_at(point);
            if matches!(curve, Curve3::Circle(_) | Curve3::Ellipse(_)) {
                let middle = 0.5 * (edge.start_parameter + edge.end_parameter);
                parameter += std::f64::consts::TAU
                    * ((middle - parameter) / std::f64::consts::TAU).round();
            }
            let span = edge.end_parameter - edge.start_parameter;
            if span == 0.0 {
                continue;
            }
            let across = (parameter - edge.start_parameter) / span;
            if !(1e-9..=1.0 - 1e-9).contains(&across)
                || Vec3::from(curve.point_at(parameter)).distance(Vec3::from(point)) > tolerance
            {
                continue;
            }
            split_edge(target, key, parameter);
        }
    }
}

/// Whether two curves trace the same points: one circle, or one line.
fn same_curve(one: &Curve3, other: &Curve3, tolerance: f64) -> bool {
    match (one, other) {
        (Curve3::Circle(_), Curve3::Circle(_)) => super::split::same_circle(one, other, tolerance),
        (Curve3::Line(one), Curve3::Line(other)) => {
            let direction = Vec3::from(one.direction);
            let Some(unit) = direction.normalize() else {
                return false;
            };
            let offset = Vec3::from(other.origin) - Vec3::from(one.origin);
            direction.is_parallel_to(Vec3::from(other.direction), tolerance)
                && (offset - unit * offset.dot(unit)).length() <= tolerance
        }
        _ => false,
    }
}

/// Whether two faces share any area rather than merely a plane.
///
/// Measured on their boxes, shrunk by the tolerance so a pair meeting along
/// an edge — which is what coplanar side walls of stacked solids do — reads
/// as apart rather than as an overlap with nothing in it.
///
/// A face that cannot be bounded is treated as overlapping, since "cannot
/// exclude" is the only safe answer a box test may give.
fn overlap(a: &Body, one: FaceKey, b: &Body, other: FaceKey, tolerance: f64) -> bool {
    let (Some(near), Some(far)) = (face_bounds(a, one), face_bounds(b, other)) else {
        return true;
    };
    shrunk(near, tolerance).overlaps(&shrunk(far, tolerance))
}

/// A box pulled in on every axis that has room to spare.
///
/// A face's box is flat in one direction by definition, and pulling that one
/// in turns it inside out — every coplanar pair then reads as apart, and a
/// wall that genuinely overlaps is passed over without a word.
fn shrunk(bounds: Aabb, by: f64) -> Aabb {
    let mut out = bounds;
    for axis in 0..3 {
        if bounds.max[axis] - bounds.min[axis] > 2.0 * by {
            out.min[axis] += by;
            out.max[axis] -= by;
        }
    }
    out
}

/// Whether two faces on the same surface cover the same region of it.
///
/// Compared by their corners rather than by area: two rings enclosing the
/// same area can be different shapes, and what matters here is that neither
/// face has any part the other does not.
pub fn same_ground(
    a: &Body,
    one: FaceKey,
    b: &Body,
    other: FaceKey,
    tolerance: f64,
) -> bool {
    let corners = |body: &Body, face: FaceKey| -> Vec<Vec3> {
        body.face_coedges(face)
            .iter()
            .filter_map(|coedge| body.coedge_vertices(*coedge))
            .filter_map(|(from, _)| Some(Vec3::from(body.vertices.get(from)?.point)))
            .collect()
    };
    let near = corners(a, one);
    let far = corners(b, other);
    if near.is_empty() || near.len() != far.len() {
        return false;
    }
    near.iter().all(|point| {
        far.iter()
            .any(|other| point.distance(*other) <= tolerance)
    }) && far.iter().all(|point| {
        near.iter()
            .any(|other| point.distance(*other) <= tolerance)
    })
}

/// Curves shared by a pair of faces.
struct Shared {
    curves: Vec<Curve3>,
}

/// The curves a face's own boundary runs along.
///
/// What cuts a partly shared wall into its shared and unshared parts. They
/// are used as whole curves rather than as the segments the edges cover, so
/// a boundary line may cut somewhere the edge itself does not reach — which
/// leaves an extra edge and never a different shape, since an imprint only
/// ever adds them.
fn boundary_curves(body: &Body, face: FaceKey) -> Vec<Curve3> {
    body.face_coedges(face)
        .iter()
        .filter_map(|coedge| {
            let edge = body.edges.get(body.coedges.get(*coedge)?.edge)?;
            body.curves.get(edge.curve).cloned()
        })
        .collect()
}

/// Every curve the two bodies' faces share.
fn shared_curves(a: &Body, b: &Body, tolerance: f64) -> Result<Vec<Shared>, Snag> {
    let near: Vec<(FaceKey, Option<Aabb>)> =
        a.face_keys().map(|key| (key, face_bounds(a, key))).collect();
    let far: Vec<(FaceKey, Option<Aabb>)> =
        b.face_keys().map(|key| (key, face_bounds(b, key))).collect();

    let mut out = Vec::new();
    let mut soups = None;
    let mut marched: Vec<(Surface, Surface, Vec<Vec<Vec3>>)> = Vec::new();
    for (one, one_box) in &near {
        for (other, other_box) in &far {
            // A box that could not be computed means "cannot exclude", so
            // the pair is tested rather than skipped.
            if let (Some(first), Some(second)) = (one_box, other_box) {
                if !first.grown(tolerance).overlaps(second) {
                    continue;
                }
            }
            let (Some(one_face), Some(other_face)) = (a.faces.get(*one), b.faces.get(*other))
            else {
                continue;
            };
            let (Some(one_surface), Some(other_surface)) = (
                a.surfaces.get(one_face.surface),
                b.surfaces.get(other_face.surface),
            ) else {
                continue;
            };
            match surfaces(one_surface, other_surface, tolerance) {
                Meeting::None | Meeting::Points(_) => {}
                Meeting::Curves(curves) => out.push(Shared { curves }),
                // Two faces on one surface. Where they cover exactly the same
                // ground there is nothing to imprint — the boolean decides
                // which copy of the shared wall survives from the two
                // normals.
                //
                // Where one covers more than the other, what separates the
                // shared part from the rest is the *other face's own
                // boundary*, already sitting in the shared plane. Cutting
                // each along it is the imprint: afterwards every piece is
                // either the whole of a shared wall or none of one, which is
                // the only case the boolean has an answer for.
                Meeting::Coincident => {
                    // Coplanar is not the same as overlapping: two boxes
                    // stacked face to face have four pairs of side walls on
                    // one plane each, meeting along a line and sharing no
                    // area at all. Only a genuine overlap has anything to
                    // decide.
                    if overlap(a, *one, b, *other, tolerance)
                        && !same_ground(a, *one, b, *other, tolerance)
                    {
                        let mut curves = boundary_curves(a, *one);
                        curves.extend(boundary_curves(b, *other));
                        if curves.is_empty() {
                            return Err(Snag::Coincident);
                        }
                        out.push(Shared { curves });
                    }
                }
                // No closed form for the pair. Their meshes may still show
                // the two faces never come near; where they do cross, the
                // crossing is traced onto both surfaces instead.
                Meeting::Unknown => {
                    let soups = soups.get_or_insert_with(|| {
                        (
                            super::near::FaceSoups::of(a, tolerance),
                            super::near::FaceSoups::of(b, tolerance),
                        )
                    });
                    if super::near::apart(&soups.0, *one, &soups.1, *other, tolerance) {
                        continue;
                    }
                    let known = match marched
                        .iter()
                        .position(|(s, t, _)| s == one_surface && t == other_surface)
                    {
                        Some(index) => index,
                        None => {
                            marched.push((
                                one_surface.clone(),
                                other_surface.clone(),
                                Vec::new(),
                            ));
                            marched.len() - 1
                        }
                    };
                    let curves = soups
                        .0
                        .triangles(*one)
                        .zip(soups.1.triangles(*other))
                        .and_then(|(first, second)| {
                            super::march::traced(
                                one_surface,
                                first,
                                other_surface,
                                second,
                                &mut marched[known].2,
                                tolerance,
                            )
                        })
                        .ok_or(Snag::NoClosedForm)?;
                    let curves = curves
                        .into_iter()
                        .flat_map(|curve| {
                            along_edges(a, &curve, tolerance)
                                .or_else(|| along_edges(b, &curve, tolerance))
                                .unwrap_or_else(|| vec![curve])
                        })
                        .collect();
                    out.push(Shared { curves });
                }
            }
        }
    }
    Ok(out)
}

/// The curves of a body's own edges that a traced curve runs along, if it
/// runs along them all the way. A shaft fitted in a hole meets the hole's
/// fillet where the hole already does, along an edge the file saved as a
/// spline: the traced curve is that edge again, a fit apart, and cutting
/// with both peels hairline slivers off between them. The edge's own curve
/// cuts the same and agrees with the vertices already there.
fn along_edges(body: &Body, curve: &Curve3, tolerance: f64) -> Option<Vec<Curve3>> {
    const SAMPLES: usize = 16;
    let Curve3::Nurbs(spline) = curve else {
        return None;
    };
    let (start, end) = spline.domain();
    let points: Vec<Vec3> = (0..=SAMPLES)
        .map(|index| {
            Vec3::from(curve.point_at(start + (end - start) * index as f64 / SAMPLES as f64))
        })
        .collect();
    let extent = points
        .iter()
        .map(|point| point.distance(points[0]))
        .fold(0.0, f64::max);
    let fit = tolerance.max(extent * 1e-3);
    let reach = |samples: &[Vec3]| {
        let low = samples.iter().fold(samples[0], |low, p| {
            Vec3::new(low.x.min(p.x), low.y.min(p.y), low.z.min(p.z))
        });
        let high = samples.iter().fold(samples[0], |high, p| {
            Vec3::new(high.x.max(p.x), high.y.max(p.y), high.z.max(p.z))
        });
        (low, high)
    };
    let (low, high) = reach(&points);
    // Only edges near the curve at all are worth the nearest-point search.
    let edges: Vec<&Curve3> = body
        .edges
        .iter()
        .filter_map(|(_, edge)| {
            let curve = body.curves.get(edge.curve)?;
            let samples: Vec<Vec3> = (0..=8)
                .map(|index| {
                    Vec3::from(curve.point_at(
                        edge.start_parameter
                            + (edge.end_parameter - edge.start_parameter) * index as f64 / 8.0,
                    ))
                })
                .collect();
            let (edge_low, edge_high) = reach(&samples);
            let margin = edge_low.distance(edge_high) * 0.25 + fit;
            let apart = edge_low.x > high.x + margin
                || edge_low.y > high.y + margin
                || edge_low.z > high.z + margin
                || edge_high.x < low.x - margin
                || edge_high.y < low.y - margin
                || edge_high.z < low.z - margin;
            (!apart).then_some(curve)
        })
        .collect();
    let mut used = vec![false; edges.len()];
    for point in &points {
        // Along the edge's whole curve, not just its span: a file keeps the
        // whole of an intersection loop as the curve of the part its face
        // still has, and the rest of the loop is the same curve again.
        let nearest = edges.iter().position(|edge| {
            Vec3::from(edge.point_at(edge.parameter_at(point.to_array()))).distance(*point) <= fit
        })?;
        used[nearest] = true;
    }
    let mut out: Vec<Curve3> = Vec::new();
    for (edge, used) in edges.into_iter().zip(used) {
        if used && !out.contains(edge) {
            out.push(edge.clone());
        }
    }
    Some(out)
}

/// Cuts every face of `body` that one of `curves` crosses.
///
/// A face split in two may still be crossed by the next curve, and by the
/// same one where a curve enters and leaves more than once, so the halves go
/// back into the list to be tried again.
fn cut_along(body: &mut Body, curves: &[Curve3], tolerance: f64) -> Result<usize, Snag> {
    let mut cuts = 0;
    for curve in curves {
        let mut pending: Vec<FaceKey> = body.face_keys().collect();
        // Each cut adds at most one face, so the work is bounded by however
        // many faces the curve can produce. The cap is a backstop against a
        // cut that somehow keeps splitting the same face rather than a limit
        // anything real reaches.
        let ceiling = body.faces.len() * 4 + 16;
        let mut done = 0;
        while let Some(face) = pending.pop() {
            done += 1;
            if done > ceiling {
                return Err(Snag::CutRefused);
            }
            if !body.faces.contains(face) {
                continue;
            }
            if let Some([kept, made]) = split_face(body, face, curve, tolerance) {
                cuts += 1;
                pending.push(kept);
                pending.push(made);
            }
        }
    }
    Ok(cuts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep::make::cuboid;

    const TOL: f64 = 1e-9;

    /// Two boxes overlapping in a corner.
    fn overlapping() -> (Body, Body) {
        (
            cuboid([0.0; 3], [10.0, 10.0, 10.0]).unwrap(),
            cuboid([5.0, 5.0, 5.0], [10.0, 10.0, 10.0]).unwrap(),
        )
    }

    #[test]
    fn imprinting_leaves_both_bodies_consistent() {
        let (mut a, mut b) = overlapping();
        imprint(&mut a, &mut b, TOL).expect("two boxes meet in planes");
        for (name, body) in [("a", &a), ("b", &b)] {
            let flaws = body.validate();
            assert!(flaws.is_empty(), "{name}: {flaws:?}");
            assert_eq!(body.euler_characteristic(), 2, "{name}");
        }
    }

    #[test]
    fn imprinting_adds_edges_and_faces_but_no_volume() {
        // An imprint only writes the shape down differently.
        let (mut a, mut b) = overlapping();
        let before = crate::brep::body_bounds(&a).unwrap();
        let faces = a.faces.len();
        let result = imprint(&mut a, &mut b, TOL).unwrap();
        assert!(result.cuts > 0, "the boxes overlap, so something is cut");
        assert!(a.faces.len() > faces);
        let after = crate::brep::body_bounds(&a).unwrap();
        assert_eq!(before, after, "the solid did not move or grow");
    }

    #[test]
    fn every_piece_is_wholly_in_or_wholly_out_afterwards() {
        // The property the whole step exists for. Before the imprint, the
        // faces of A that the corner of B passes through are partly inside
        // it; afterwards no face is.
        use crate::brep::{contains_point, Containment};
        let (mut a, mut b) = overlapping();
        imprint(&mut a, &mut b, TOL).unwrap();
        for face in a.face_keys() {
            let bounds = crate::brep::face_bounds(&a, face).unwrap();
            let mut seen = Vec::new();
            // The corners of the face's own box, nudged onto the face, are
            // enough to catch a face that straddles the boundary.
            for coedge in a.face_coedges(face) {
                let Some((from, _)) = a.coedge_vertices(coedge) else {
                    continue;
                };
                let point = a.vertices.get(from).unwrap().point;
                seen.push(contains_point(&b, point, 1e-6));
            }
            let _ = bounds;
            let inside = seen.iter().filter(|c| **c == Containment::Inside).count();
            let outside = seen.iter().filter(|c| **c == Containment::Outside).count();
            assert!(
                inside == 0 || outside == 0,
                "face {face:?} has corners both in and out: {seen:?}"
            );
        }
    }

    #[test]
    fn boxes_that_do_not_touch_are_left_alone() {
        let mut a = cuboid([0.0; 3], [1.0, 1.0, 1.0]).unwrap();
        let mut b = cuboid([50.0, 50.0, 50.0], [1.0, 1.0, 1.0]).unwrap();
        let faces = (a.faces.len(), b.faces.len());
        let result = imprint(&mut a, &mut b, TOL).unwrap();
        assert_eq!(result.cuts, 0);
        assert_eq!((a.faces.len(), b.faces.len()), faces);
    }

    #[test]
    fn faces_meeting_on_one_plane_need_nothing_when_they_cover_it_alike() {
        // Stacked exactly: the two that meet cover the same ground and the
        // four pairs of side walls share a plane without sharing any area.
        // Neither needs cutting, so there is nothing to imprint.
        let mut a = cuboid([0.0; 3], [10.0, 10.0, 10.0]).unwrap();
        let mut b = cuboid([0.0, 0.0, 10.0], [10.0, 10.0, 10.0]).unwrap();
        let result = imprint(&mut a, &mut b, TOL).expect("nothing to cut");
        assert_eq!(result.cuts, 0);
    }

    #[test]
    fn a_wall_shared_in_part_is_cut_where_the_sharing_stops() {
        // A smaller box on top: its bottom covers a corner of the other's
        // top. What separates the shared part from the rest is the small
        // box's own boundary, so cutting the big face along it leaves pieces
        // that are each wholly shared or wholly not — which is the only
        // shape the boolean has an answer for.
        let mut a = cuboid([0.0; 3], [10.0, 10.0, 10.0]).unwrap();
        let mut b = cuboid([0.0, 0.0, 10.0], [4.0, 4.0, 4.0]).unwrap();
        let before = (a.faces.len(), b.faces.len());
        let result = imprint(&mut a, &mut b, TOL).expect("a wall to cut");
        assert!(result.cuts > 0, "{result:?}");
        assert!(a.faces.len() > before.0, "the big face was divided");
        assert_eq!(b.faces.len(), before.1, "the small one had nothing to lose");
        // An imprint only adds edges, so both are still solids.
        assert!(a.validate().is_empty());
        assert!(b.validate().is_empty());
        assert_eq!(a.euler_characteristic(), 2);
        assert_eq!(b.euler_characteristic(), 2);

        // And the piece that matches the small box's footprint really is one
        // face now, rather than a corner of a larger one.
        let footprint = a.face_keys().filter(|face| {
            face_bounds(&a, *face).is_some_and(|box_| {
                (box_.min[2] - 10.0).abs() < TOL
                    && (box_.max[0] - 4.0).abs() < TOL
                    && (box_.max[1] - 4.0).abs() < TOL
                    && box_.min[0].abs() < TOL
                    && box_.min[1].abs() < TOL
            })
        });
        assert_eq!(footprint.count(), 1);
    }

    #[test]
    fn a_pair_with_no_closed_form_is_refused() {
        let mut a = cuboid([0.0; 3], [10.0, 10.0, 10.0]).unwrap();
        let mut b = cuboid([5.0; 3], [10.0, 10.0, 10.0]).unwrap();
        // Turn one of B's faces into a torus, which no pair here can meet.
        let face = b.face_keys().next().unwrap();
        let surface = b.faces.get(face).unwrap().surface;
        *b.surfaces.get_mut(surface).unwrap() =
            crate::brep::Surface::Torus(crate::brep::Torus {
                frame: crate::space::Plane::XY,
                major_radius: 4.0,
                minor_radius: 1.0,
            });
        assert_eq!(imprint(&mut a, &mut b, TOL), Err(Snag::NoClosedForm));
    }

    #[test]
    fn the_prefilter_does_not_lose_a_real_meeting() {
        // Boxes that share only an edge region: most face pairs are apart,
        // and the box test has to keep the few that are not.
        let (mut a, mut b) = overlapping();
        let result = imprint(&mut a, &mut b, TOL).unwrap();
        assert!(result.meetings > 0);
        assert!(
            result.meetings < 36,
            "the prefilter rejected nothing: {}",
            result.meetings
        );
    }

    #[test]
    fn imprinting_at_survey_coordinates_works_the_same() {
        let origin = [512_345.678, 4_512_345.678, 91.5];
        let mut a = cuboid(origin, [10.0, 10.0, 10.0]).unwrap();
        let mut b = cuboid(
            [origin[0] + 5.0, origin[1] + 5.0, origin[2] + 5.0],
            [10.0, 10.0, 10.0],
        )
        .unwrap();
        imprint(&mut a, &mut b, 1e-6).expect("the same two boxes, further out");
        assert!(a.validate().is_empty());
        assert!(b.validate().is_empty());
        assert!(a.worst_vertex_gap() < 1e-6);
    }
}
