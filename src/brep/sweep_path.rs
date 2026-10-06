//! Profile-preserving path sweeps, independent of application entities.
//!
//! Profile curves stay rational curves: a circle never becomes a polygon.
//! Undeformed straight and circular paths use the analytic sweep builders.
//! General paths use tolerance-controlled cubic transport patches, rational
//! in the profile parameter, with shared rims and rails between patches.

use super::geometry::{Curve3, Surface};
use super::nurbs_builder::{RationalCurve2, RationalCurve3};
use super::topology::{
    Body, Coedge, CoedgeKey, Edge, EdgeKey, Face, FaceKey, Loop, Lump, Shell, ShellKey, Vertex,
    VertexKey,
};
use super::Provenance;
use super::Placement;
use crate::geom2d::{Arc, Curve, Line, NurbsCurve, Transform};
use crate::space::{NurbsCurve3, NurbsSurface3, Plane, Vec3};
use std::f64::consts::{PI, TAU};

/// A connected path, expressed without projecting three-dimensional curves.
#[derive(Clone, Copy)]
pub enum SweepPath<'a> {
    Planar { plane: Plane, curves: &'a [Curve] },
    Polyline3d { points: &'a [[f64; 3]], closed: bool },
    Nurbs3(&'a NurbsCurve3),
}

/// Placement and deformation of a swept profile. Angles are radians.
#[derive(Clone, Copy, Debug)]
pub struct SweepOptions {
    /// Align the profile normal to the initial path tangent.
    pub align: bool,
    /// World-space point on the profile placed at the path start.
    pub base_point: Option<[f64; 3]>,
    /// Initial rotation about the profile normal.
    pub rotation: f64,
    /// Additional rotation accumulated over the path length.
    pub twist: f64,
    /// Final scale; scale changes linearly from one over the path length.
    pub scale: f64,
    /// Follow the path's curvature frame instead of minimal-twist transport.
    pub bank: bool,
    /// Omit caps, including when the supplied profile is closed.
    pub surface: bool,
}

impl Default for SweepOptions {
    fn default() -> Self {
        Self {
            align: true,
            base_point: None,
            rotation: 0.0,
            twist: 0.0,
            scale: 1.0,
            bank: false,
            surface: false,
        }
    }
}

/// Curve-length weighted common anchor for several independent profiles.
/// Relative profile positions are preserved when this same world-space
/// point is supplied as `base_point` for every member of a sweep selection.
pub fn sweep_profile_group_base(profiles: &[(Plane, Vec<Vec<Curve>>)]) -> Option<[f64; 3]> {
    let (first_plane, first_wires) = profiles.first()?;
    let origin = Vec3::from(sweep_profile_base(*first_plane, first_wires)?);
    let mut moment = Vec3::ZERO;
    let mut length = 0.0;
    for (plane, wires) in profiles {
        let base = Vec3::from(sweep_profile_base(*plane, wires)?);
        let mut profile_length = 0.0;
        for wire in wires {
            for curve in expanded(wire)? {
                profile_length += Piece::Planar(*plane, curve, true).length();
            }
        }
        moment = moment + (base - origin) * profile_length;
        length += profile_length;
    }
    if !length.is_finite() || length <= 1e-14 || !moment.is_finite() { return None; }
    Some((origin + moment / length).to_array())
}

/// Directed first point of a valid bounded path.
pub fn sweep_path_start(path: SweepPath<'_>) -> Option<[f64; 3]> {
    Some(path_pieces(path)?.first()?.point(0.0).to_array())
}

/// Unit tangent where a sweep path starts.
pub fn sweep_path_tangent(path: SweepPath<'_>) -> Option<[f64; 3]> {
    Some(path_pieces(path)?.first()?.tangent(0.0)?.to_array())
}

/// Rigid placement of the original profile onto the first sweep section.
/// Initial twist and scale are zero and one respectively; accumulated
/// deformation therefore does not affect this source-profile grip mapping.
pub fn sweep_profile_placement(
    plane: Plane, wires: &[Vec<Curve>], path: SweepPath<'_>, options: SweepOptions,
) -> Option<Placement> {
    let pieces = path_pieces(path)?;
    let first = pieces.first()?;
    let base = Vec3::from(options.base_point.or_else(|| sweep_profile_base(plane, wires))?);
    initial_placement(plane, base, first.point(0.0), first.tangent(0.0)?, options)
}

fn initial_placement(plane: Plane, base: Vec3, start: Vec3, tangent: Vec3, options: SweepOptions) -> Option<Placement> {
    if !base.is_finite() || !start.is_finite() || !options.rotation.is_finite() { return None; }
    let normal = Vec3::from(plane.normal()?);
    let source_x = Vec3::from(plane.x_axis).normalize()?;
    let source_y = normal.cross(source_x).normalize()?;
    let projected_x = source_x - tangent * source_x.dot(tangent);
    let (target_x, target_y, target_normal) = if projected_x.length() > 1e-10 {
        let x = projected_x.normalize()?;
        (x, tangent.cross(x).normalize()?, tangent)
    } else {
        // A path parallel to the profile's first axis has no projected
        // first direction. Keep the profile's cyclic frame in this case;
        // reversing that path must not turn the selected profile over.
        (source_y, normal, source_x)
    };
    let axis = if options.align { target_normal } else { normal };
    let map = |vector: Vec3| -> Option<Vec3> {
        let aligned = if options.align {
            target_x * vector.dot(source_x) + target_y * vector.dot(source_y) + target_normal * vector.dot(normal)
        } else { vector };
        Some(rotate(aligned, axis, options.rotation))
    };
    Some(Placement { x_axis: map(Vec3::X)?.to_array(), y_axis: map(Vec3::Y)?.to_array(),
        z_axis: map(Vec3::Z)?.to_array(), origin: (start - map(base)?).to_array() })
}

/// Default source-profile anchor. Closed conics use their centre and open
/// conics their middle point; other chains use
/// twenty equally spaced boundary samples including the directed endpoints.
/// The endpoint convention deliberately preserves a polyline's start vertex.
/// Multiple boundary loops contribute in proportion to their curve lengths.
pub fn sweep_profile_base(plane: Plane, wires: &[Vec<Curve>]) -> Option<[f64; 3]> {
    plane.normal()?;
    let origin = Vec3::from(plane.point_at(wires.first()?.first()?.point_at(0.0)));
    let mut total = 0.0;
    let mut moment = Vec3::ZERO;
    for wire in wires {
        let pieces = expanded(wire)?;
        let senses = chain_senses(&pieces)?;
        let path = pieces.iter().zip(senses).map(|(curve, forward)|
            Piece::Planar(plane, curve.clone(), forward)).collect::<Vec<_>>();
        let lengths = path.iter().map(Piece::length).collect::<Vec<_>>();
        let length = lengths.iter().sum::<f64>();
        if !length.is_finite() || length <= 1e-14 { return None; }
        let conic_anchor = conic_anchor(&pieces);
        let anchor = if let Some(point) = conic_anchor {
            Vec3::from(plane.point_at(point))
        } else {
            let mut samples = Vec3::ZERO;
            for sample in 0..20 {
                let mut distance = length * sample as f64 / 19.0;
                let mut index = 0;
                while index + 1 < path.len() && distance > lengths[index] {
                    distance -= lengths[index];
                    index += 1;
                }
                let parameter = if matches!(&pieces[index], Curve::Line(_) | Curve::Arc(_)) && plane.is_orthonormal() {
                    (distance / lengths[index]).clamp(0.0, 1.0)
                } else {
                    let (mut low, mut high) = (0.0, 1.0);
                    for _ in 0..40 {
                        let middle = (low + high) * 0.5;
                        if path[index].length_to(middle) < distance { low = middle; } else { high = middle; }
                    }
                    (low + high) * 0.5
                };
                samples = samples + (path[index].point(parameter) - origin);
            }
            origin + samples / 20.0
        };
        moment = moment + (anchor - origin) * length;
        total += length;
    }
    if !total.is_finite() || total <= 1e-14 || !moment.is_finite() { return None; }
    Some((origin + moment / total).to_array())
}

fn conic_anchor(pieces: &[Curve]) -> Option<[f64; 2]> {
    let center = match pieces.first()? {
        Curve::Arc(first) if pieces.iter().all(|piece| matches!(piece, Curve::Arc(arc)
            if near2(first.centre, arc.centre) && (first.radius - arc.radius).abs() <= first.radius.abs().max(1.0) * 1e-10)) => Some(first.centre),
        Curve::Ellipse(first) if pieces.iter().all(|piece| matches!(piece, Curve::Ellipse(arc)
            if first.ellipse == arc.ellipse)) => Some(first.ellipse.centre),
        _ => None,
    }?;
    let senses = chain_senses(pieces)?;
    if chain_closed(pieces, &senses) { return Some(center); }
    let spans = pieces.iter().map(|piece| match piece {
        Curve::Arc(arc) => arc.sweep(),
        Curve::Ellipse(arc) => arc.sweep(),
        _ => 0.0,
    }).collect::<Vec<_>>();
    let mut middle = spans.iter().sum::<f64>() * 0.5;
    let mut index = 0;
    while index + 1 < spans.len() && middle > spans[index] {
        middle -= spans[index];
        index += 1;
    }
    let parameter = middle / spans[index];
    Some(pieces[index].point_at(if senses[index] { parameter } else { 1.0 - parameter }))
}

/// Sweeps an outer profile and optional holes along a connected spatial path.
/// A single open wire produces a sheet even when `surface` is false.
///
/// Invalid/zero paths, a reversing cusp, non-positive scale, invalid profile
/// input, and incompatible deformation at a closed seam return `None`.
/// General-path fits are checked against 1e-7 of the local model extent;
/// excessive complexity is refused rather than silently reducing accuracy.
pub fn sweep_path(
    profile_plane: Plane,
    wires: &[Vec<Curve>],
    path: SweepPath<'_>,
    mut options: SweepOptions,
) -> Option<Body> {
    if !options.rotation.is_finite() || !options.twist.is_finite()
        || !options.scale.is_finite() || options.scale <= 1e-9
        || options.twist.abs() > TAU * 1024.0
    {
        return None;
    }
    let mut wires = prepare_wires(wires)?;
    let sheet = options.surface || !wires[0].closed;
    if !sheet && wires.iter().any(|wire| !wire.closed) {
        return None;
    }
    if wires.iter().any(|wire| !wire.closed) && wires.len() != 1 {
        return None;
    }
    let pieces = path_pieces(path)?;
    let start = pieces.first()?.point(0.0);
    let tangent = pieces.first()?.tangent(0.0)?;
    let end = pieces.last()?.point(1.0);
    let extent = pieces.iter().map(|piece| piece.length()).sum::<f64>();
    if !extent.is_finite() || extent <= 1e-12 {
        return None;
    }
    let closed = start.distance(end) <= extent.max(1.0) * 1e-9;
    // The reference modeler refuses to twist or scale along a path with a
    // corner (a mitred joint cannot change the section on both sides).
    let joints = pieces.windows(2).map(|pair| (&pair[0], &pair[1]))
        .chain(closed.then(|| (pieces.last().unwrap(), &pieces[0])));
    let mut cornered = false;
    for (before, after) in joints {
        cornered |= before.tangent(1.0)?.dot(after.tangent(0.0)?) < 1.0 - 1e-9;
    }
    if cornered && (options.twist.abs() > 1e-12 || (options.scale - 1.0).abs() > 1e-12) {
        return None;
    }
    let mut source_wires = wires.iter().map(|wire| wire.source.clone()).collect::<Vec<_>>();
    let base = Vec3::from(options.base_point
        .or_else(|| sweep_profile_base(profile_plane, &source_wires))?);
    if !base.is_finite() {
        return None;
    }
    // Keep local coordinates near the anchor. Otherwise a tiny rotation of a
    // profile far from the origin magnifies cancellation in transport fits.
    let uv = profile_plane.project(base.to_array())?;
    let recenter = Transform::translation([-uv[0], -uv[1]]);
    source_wires = source_wires.iter().map(|wire| wire.iter().map(|curve|
        curve.transformed(&recenter)).collect::<Option<Vec<_>>>()).collect::<Option<Vec<_>>>()?;
    wires = prepare_wires(&source_wires)?;
    let profile_plane = Plane::from_axes(profile_plane.point_at(uv), profile_plane.x_axis, profile_plane.y_axis);
    let placement = initial_placement(profile_plane, base, start, tangent, options)?;
    let x = Vec3::from(placement.vector(profile_plane.x_axis));
    let y = Vec3::from(placement.vector(profile_plane.y_axis));
    let aligned_normal = x.cross(y).normalize()?;
    let first = Frame { origin: Vec3::from(placement.point(profile_plane.origin)), x, y };
    let up = aligned_normal.dot(tangent);
    if !sheet && up.abs() <= 1e-8 {
        return None;
    }
    if up.abs() > 1.0 - 1e-10 && rotationally_invariant(&source_wires) {
        // Turning a complete concentric circle changes its seam parameter,
        // not its geometry. Avoid a helical spline skin for the same cone or
        // tube: this is an exact symmetry reduction, including closed paths.
        options.twist = 0.0;
        options.bank = false;
    }
    if closed && ((options.scale - 1.0).abs() > 1e-10
        || (options.twist / TAU - (options.twist / TAU).round()).abs() > 1e-9)
    {
        return None;
    }

    // Retain cylinders, planes, cones and tori in the ordinary cases.
    if options.twist.abs() <= 1e-12 && (options.scale - 1.0).abs() <= 1e-12
        && pieces.len() == 1
    {
        let plane = first.plane();
        let analytic = match &pieces[0] {
            Piece::Line(a, b) if sheet && wires.len() == 1 =>
                super::sweep::extrude_surface(plane, &source_wires[0], (*b - *a).to_array()),
            Piece::Line(a, b) if !sheet =>
                super::sweep::extrude_region(plane, &source_wires, (*b - *a).to_array()),
            Piece::Planar(path_plane, Curve::Arc(arc), forward) => {
                let pivot = path_plane.point_at(arc.centre);
                let axis = path_plane.normal()?;
                let angle = arc.sweep() * if *forward { 1.0 } else { -1.0 };
                if sheet {
                    super::sweep::revolve_surface_region(plane, &source_wires, pivot, axis, angle)
                } else {
                    super::sweep::revolve_region(plane, &source_wires, pivot, axis, angle)
                }
            }
            _ => None,
        };
        if let Some(body) = analytic {
            return Some(body);
        }
    }

    let radius = wires.iter().flat_map(|wire| &wire.curves)
        .flat_map(|curve| &curve.points)
        .map(|p| Vec3::from(profile_plane.point_at(*p)).distance(base))
        .fold(1.0_f64, f64::max) * options.scale.max(1.0);
    let patches = transported_patches(&pieces, first, options,
        extent, radius, closed)?;
    // A sheet need not sweep out any volume (for example an in-plane line
    // translated sideways), so only solid sections use the volume check.
    let outward = if sheet { up >= 0.0 } else { regular_transport(&wires, &patches)? };
    let body = build_body(&wires, &patches, sheet, closed, outward)?;
    if options.twist.abs() <= 1e-12 && (options.scale - 1.0).abs() <= 1e-12 {
        let mut exact = body.clone();
        if analytic_ruled_faces(&mut exact).is_some() && exact.validate().is_empty() {
            return Some(exact);
        }
    }
    Some(body)
}

/// Faces swept along straight and circular runs without twist or scaling
/// lie on planes, cylinders, cones, spheres and tori. Store those surfaces
/// exactly and let neighbouring faces on one surface share a face, as the
/// reference modeler does: a round profile gives one tube per run, and a
/// flat side shared by coplanar runs one face.
fn analytic_ruled_faces(body: &mut Body) -> Option<()> {
    let size = body.vertices.iter()
        .flat_map(|(_, vertex)| vertex.point)
        .map(f64::abs)
        .fold(1.0_f64, f64::max);
    let tolerance = size * 1e-9;
    for face in body.face_keys().collect::<Vec<_>>() {
        let (old, forward) = body.faces.get(face).map(|face| (face.surface, face.forward))?;
        let Some(Surface::Nurbs(surface)) = body.surfaces.get(old) else { continue };
        let found = ruled_analytic(surface, tolerance).or_else(|| revolved_analytic(surface, tolerance));
        let Some((analytic, agrees)) = found else { continue };
        let surface = body.surfaces.insert(analytic);
        let target = body.faces.get_mut(face)?;
        target.surface = surface;
        target.forward = forward == agrees;
        body.surfaces.remove(old);
        for coedge in body.face_coedges(face) {
            body.coedges.get_mut(coedge)?.pcurve = None;
        }
    }
    loop {
        let joinable = body.edge_keys().find(|edge| joins_one_surface(body, *edge, tolerance));
        let Some(edge) = joinable else { break };
        dissolve_edge(body, edge)?;
    }
    Some(())
}

/// A plane or circular cylinder that a ruled patch lies on exactly, and
/// whether that surface's natural normal agrees with the patch's.
fn ruled_analytic(surface: &NurbsSurface3, tolerance: f64) -> Option<(Surface, bool)> {
    let mut direction: Option<Vec3> = None;
    for (row, weights) in surface.control_points().iter().zip(surface.weights()) {
        let start = Vec3::from(*row.first()?);
        let along = (Vec3::from(*row.last()?) - start).normalize()?;
        match direction {
            Some(direction) if direction.cross(along).length() > 1e-9 || direction.dot(along) < 0.0 => return None,
            Some(_) => {}
            None => direction = Some(along),
        }
        for (point, weight) in row.iter().zip(weights) {
            let offset = Vec3::from(*point) - start;
            if (offset - along * offset.dot(along)).length() > tolerance
                || (weight - weights[0]).abs() > weights[0].abs() * 1e-12
            {
                return None;
            }
        }
    }
    let direction = direction?;
    let origin = Vec3::from(surface.point_at(0.0, 0.0));
    let flat = |point: Vec3| point - direction * (point - origin).dot(direction);
    let samples = (0..=8).map(|step| flat(Vec3::from(surface.point_at(step as f64 / 8.0, 0.0))))
        .collect::<Vec<_>>();
    let natural = Vec3::from(surface.normal_at(0.5, 0.5)?);
    let middle = Vec3::from(surface.point_at(0.5, 0.5));
    let chord = samples.iter().map(|point| *point - origin).find(|chord| chord.length() > tolerance)?;
    let normal = direction.cross(chord).normalize()?;
    let mut controls = surface.control_points().iter().flatten().map(|point| Vec3::from(*point));
    if controls.all(|point| (point - origin).dot(normal).abs() <= tolerance) {
        let plane = Plane::orthonormal(origin.to_array(), direction.to_array(), normal.to_array())?;
        return Some((Surface::Plane(plane), normal.dot(natural) > 0.0));
    }
    // The circle through the first, middle and last section points.
    let (a, b, c) = (samples[0], samples[4], samples[8]);
    let (ab, ac) = (b - a, c - a);
    let axis = ab.cross(ac);
    if axis.length() <= tolerance * tolerance.max(1.0) {
        return None;
    }
    let centre = a + (axis.cross(ab) * ac.dot(ac) + ac.cross(axis) * ab.dot(ab)) * (0.5 / axis.dot(axis));
    let radius = a.distance(centre);
    if axis.cross(direction).length() > axis.length() * 1e-9
        || samples.iter().any(|point| (point.distance(centre) - radius).abs() > tolerance)
    {
        return None;
    }
    let x_axis = (a - centre) * (1.0 / radius);
    let base = Plane::orthonormal(centre.to_array(), x_axis.to_array(), direction.to_array())?;
    let outward = flat(middle) - centre;
    Some((Surface::Cylinder(super::geometry::Cylinder { base, radius }), outward.dot(natural) > 0.0))
}

/// The plane, cylinder, cone, sphere or torus that a patch turned exactly
/// about one axis lies on, and whether its natural normal agrees with the
/// patch's.
fn revolved_analytic(surface: &NurbsSurface3, tolerance: f64) -> Option<(Surface, bool)> {
    // Each section point runs round a circle about the common axis.
    let mut axis: Option<(Vec3, Vec3)> = None;
    let mut meridian = Vec::new();
    for step in 0..=8 {
        let u = step as f64 / 8.0;
        let ring = [0.0, 0.25, 0.5, 0.75, 1.0].map(|v| Vec3::from(surface.point_at(u, v)));
        let (a, b, c) = (ring[0], ring[2], ring[4]);
        let (ab, ac) = (b - a, c - a);
        let normal = ab.cross(ac);
        if normal.length() <= tolerance * tolerance.max(1.0) {
            // A point on the axis stays put.
            if ring.iter().any(|point| point.distance(a) > tolerance) { return None; }
            meridian.push(a);
            continue;
        }
        let centre = a + (normal.cross(ab) * ac.dot(ac) + ac.cross(normal) * ab.dot(ab)) * (0.5 / normal.dot(normal));
        let radius = a.distance(centre);
        if ring.iter().any(|point| (point.distance(centre) - radius).abs() > tolerance) { return None; }
        let direction = normal.normalize()?;
        match axis {
            Some((origin, along)) => {
                let offset = centre - origin;
                if along.cross(direction).length() > 1e-9
                    || (offset - along * offset.dot(along)).length() > tolerance
                {
                    return None;
                }
            }
            None => axis = Some((centre, direction)),
        }
        meridian.push(a);
    }
    let (origin, along) = axis?;
    // Section points as (distance from the axis, height along it).
    let polar = |point: Vec3| {
        let offset = point - origin;
        let height = offset.dot(along);
        ((offset - along * height).length(), height)
    };
    let section = meridian.iter().map(|point| polar(*point)).collect::<Vec<_>>();
    let radial = (meridian[0] - origin) - along * (meridian[0] - origin).dot(along);
    let x_axis = radial.normalize().or_else(|| along.cross(Vec3::X).normalize()).or_else(|| along.cross(Vec3::Y).normalize())?;
    let natural = Vec3::from(surface.normal_at(0.5, 0.5)?);
    let middle = Vec3::from(surface.point_at(0.5, 0.5));
    let (middle_radius, middle_height) = polar(middle);
    let middle_radial = ((middle - origin) - along * middle_height).normalize()?;
    let (r0, h0) = section[0];
    let (r1, h1) = section[8];
    let frame_at = |height: f64| Plane::orthonormal((origin + along * height).to_array(), x_axis.to_array(), along.to_array());
    if section.iter().all(|(radius, _)| (radius - r0).abs() <= tolerance) && r0 > tolerance {
        let base = frame_at(h0)?;
        let surface = Surface::Cylinder(super::geometry::Cylinder { base, radius: r0 });
        return Some((surface, middle_radial.dot(natural) > 0.0));
    }
    if section.iter().all(|(_, height)| (height - h0).abs() <= tolerance) {
        let plane = frame_at(h0)?;
        return Some((Surface::Plane(plane), along.dot(natural) > 0.0));
    }
    // A straight meridian: a cone.
    let chord = (r1 - r0, h1 - h0);
    let length = chord.0.hypot(chord.1);
    let off_line = |(radius, height): (f64, f64)| ((radius - r0) * chord.1 - (height - h0) * chord.0).abs() / length;
    if length > tolerance && section.iter().all(|point| off_line(*point) <= tolerance) {
        let slope = chord.0 / chord.1;
        let half_angle = (-slope).atan();
        let (radius, height) = if r0 >= r1 { (r0, h0) } else { (r1, h1) };
        let base = frame_at(height)?;
        let outward = middle_radial + along * half_angle.tan();
        let surface = Surface::Cone(super::geometry::Cone { base, radius, half_angle });
        return Some((surface, outward.dot(natural) > 0.0));
    }
    // A circular meridian: a torus, or a sphere about a centre on the axis.
    let (a, b, c) = (section[0], section[4], section[8]);
    let determinant = 2.0 * (a.0 * (b.1 - c.1) + b.0 * (c.1 - a.1) + c.0 * (a.1 - b.1));
    if determinant.abs() <= tolerance * tolerance.max(1.0) { return None; }
    let square = |(x, y): (f64, f64)| x * x + y * y;
    let centre = (
        (square(a) * (b.1 - c.1) + square(b) * (c.1 - a.1) + square(c) * (a.1 - b.1)) / determinant,
        (square(a) * (c.0 - b.0) + square(b) * (a.0 - c.0) + square(c) * (b.0 - a.0)) / determinant,
    );
    let minor = (a.0 - centre.0).hypot(a.1 - centre.1);
    if section.iter().any(|point| ((point.0 - centre.0).hypot(point.1 - centre.1) - minor).abs() > tolerance) {
        return None;
    }
    let outward = (middle_radial * (middle_radius - centre.0) + along * (middle_height - centre.1)).normalize()?;
    let agrees = outward.dot(natural) > 0.0;
    let frame = frame_at(centre.1)?;
    if centre.0.abs() <= tolerance {
        return Some((Surface::Sphere(super::geometry::Sphere { frame, radius: minor }), agrees));
    }
    let torus = super::geometry::Torus { frame, major_radius: centre.0, minor_radius: minor };
    Some((Surface::Torus(torus), agrees))
}

/// Whether an edge separates two faces of one analytic surface, both facing
/// the same way, so that it can be removed.
fn joins_one_surface(body: &Body, edge: EdgeKey, tolerance: f64) -> bool {
    let Some(edge) = body.edges.get(edge) else { return false };
    let [first, second] = edge.coedges[..] else { return false };
    let side = |coedge| {
        let ring = body.coedges.get(coedge)?.owner;
        let face = body.loops.get(ring)?.owner;
        let value = body.faces.get(face)?;
        Some((face, ring, value.forward, body.surfaces.get(value.surface)?))
    };
    let (Some((face_a, loop_a, forward_a, surface_a)), Some((face_b, loop_b, forward_b, surface_b))) =
        (side(first), side(second))
    else {
        return false;
    };
    if face_a == face_b && loop_a != loop_b {
        return false;
    }
    match (surface_a, surface_b) {
        (Surface::Plane(a), Surface::Plane(b)) => {
            let (Some(normal_a), Some(normal_b)) = (a.normal(), b.normal()) else { return false };
            let (normal_a, normal_b) = (Vec3::from(normal_a), Vec3::from(normal_b));
            let facing = if forward_a == forward_b { 1.0 } else { -1.0 };
            normal_a.dot(normal_b) * facing > 1.0 - 1e-12
                && (Vec3::from(b.origin) - Vec3::from(a.origin)).dot(normal_a).abs() <= tolerance
        }
        (Surface::Cylinder(a), Surface::Cylinder(b)) => {
            let (Some(axis_a), Some(axis_b)) = (a.base.normal(), b.base.normal()) else { return false };
            let (axis_a, axis_b) = (Vec3::from(axis_a), Vec3::from(axis_b));
            let offset = Vec3::from(b.base.origin) - Vec3::from(a.base.origin);
            forward_a == forward_b
                && (a.radius - b.radius).abs() <= tolerance
                && axis_a.cross(axis_b).length() <= 1e-9
                && (offset - axis_a * offset.dot(axis_a)).length() <= tolerance
        }
        (Surface::Cone(a), Surface::Cone(b)) => {
            let (Some(axis_a), Some(axis_b)) = (a.base.normal(), b.base.normal()) else { return false };
            let (axis_a, axis_b) = (Vec3::from(axis_a), Vec3::from(axis_b));
            let offset = Vec3::from(b.base.origin) - Vec3::from(a.base.origin);
            let radius_a_at_b = a.radius - offset.dot(axis_a) * a.half_angle.tan();
            forward_a == forward_b
                && axis_a.dot(axis_b) > 1.0 - 1e-12
                && (a.half_angle - b.half_angle).abs() <= 1e-9
                && (offset - axis_a * offset.dot(axis_a)).length() <= tolerance
                && (radius_a_at_b - b.radius).abs() <= tolerance
        }
        (Surface::Sphere(a), Surface::Sphere(b)) => {
            forward_a == forward_b
                && (a.radius - b.radius).abs() <= tolerance
                && Vec3::from(a.frame.origin).distance(Vec3::from(b.frame.origin)) <= tolerance
        }
        (Surface::Torus(a), Surface::Torus(b)) => {
            let (Some(axis_a), Some(axis_b)) = (a.frame.normal(), b.frame.normal()) else { return false };
            forward_a == forward_b
                && Vec3::from(axis_a).cross(Vec3::from(axis_b)).length() <= 1e-9
                && Vec3::from(a.frame.origin).distance(Vec3::from(b.frame.origin)) <= tolerance
                && (a.major_radius - b.major_radius).abs() <= tolerance
                && (a.minor_radius - b.minor_radius).abs() <= tolerance
        }
        _ => false,
    }
}

/// Removes an edge between two faces of one surface: their loops are joined,
/// or, when both uses already belong to one loop, that loop is split in two.
fn dissolve_edge(body: &mut Body, edge: EdgeKey) -> Option<()> {
    let [first, second] = body.edges.get(edge)?.coedges[..] else { return None };
    let loop_a = body.coedges.get(first)?.owner;
    let loop_b = body.coedges.get(second)?.owner;
    let face_a = body.loops.get(loop_a)?.owner;
    let face_b = body.loops.get(loop_b)?.owner;
    // A ring rotated to start just after `coedge`, without it.
    let after = |ring: &[CoedgeKey], coedge: CoedgeKey| -> Option<Vec<CoedgeKey>> {
        let at = ring.iter().position(|key| *key == coedge)?;
        Some(ring[at + 1..].iter().chain(&ring[..at]).copied().collect())
    };
    if loop_a == loop_b {
        let ring = after(&body.loops.get(loop_a)?.coedges, first)?;
        let split = ring.iter().position(|key| *key == second)?;
        let (inner, outer) = (ring[..split].to_vec(), ring[split + 1..].to_vec());
        // Two uses running straight back on each other leave one ring, or
        // none at all when they were the whole of it.
        let (keep, moved) = if loop_area(body, face_a, &inner)? > loop_area(body, face_a, &outer)? {
            (inner, outer)
        } else {
            (outer, inner)
        };
        if keep.is_empty() {
            body.loops.remove(loop_a);
            body.faces.get_mut(face_a)?.loops.retain(|key| *key != loop_a);
        } else {
            body.loops.get_mut(loop_a)?.coedges = keep;
        }
        if !moved.is_empty() {
            let ring = body.loops.insert(Loop { coedges: moved.clone(), owner: face_a, provenance: Provenance::Synthesized });
            for coedge in moved {
                body.coedges.get_mut(coedge)?.owner = ring;
            }
            body.faces.get_mut(face_a)?.loops.push(ring);
        }
    } else {
        let mut ring = after(&body.loops.get(loop_a)?.coedges, first)?;
        let tail = after(&body.loops.get(loop_b)?.coedges, second)?;
        for coedge in &tail {
            body.coedges.get_mut(*coedge)?.owner = loop_a;
        }
        ring.extend(tail);
        body.loops.get_mut(loop_a)?.coedges = ring;
        body.loops.remove(loop_b);
        let others = body.faces.get(face_b)?.loops.iter().copied().filter(|key| *key != loop_b).collect::<Vec<_>>();
        for other in others {
            body.loops.get_mut(other)?.owner = face_a;
            body.faces.get_mut(face_a)?.loops.push(other);
        }
        let face = body.faces.remove(face_b)?;
        body.shells.get_mut(face.owner)?.faces.retain(|key| *key != face_b);
        body.surfaces.remove(face.surface);
    }
    body.coedges.remove(first);
    body.coedges.remove(second);
    let removed = body.edges.remove(edge)?;
    if !body.edges.iter().any(|(_, edge)| edge.curve == removed.curve) {
        body.curves.remove(removed.curve);
    }
    for vertex in [removed.start, removed.end] {
        if !body.edges.iter().any(|(_, edge)| edge.start == vertex || edge.end == vertex) {
            body.vertices.remove(vertex);
        }
    }
    Some(())
}

/// The area a loop encloses across its face's axis, from the polygon of its
/// coedge start points; it only orders the two loops a split leaves.
fn loop_area(body: &Body, face: FaceKey, ring: &[CoedgeKey]) -> Option<f64> {
    let normal = match body.surfaces.get(body.faces.get(face)?.surface)? {
        Surface::Plane(plane) => Vec3::from(plane.normal()?),
        Surface::Cylinder(cylinder) => Vec3::from(cylinder.base.normal()?),
        _ => return Some(0.0),
    };
    let points = ring.iter()
        .map(|coedge| {
            let (start, _) = body.coedge_vertices(*coedge)?;
            Some(Vec3::from(body.vertices.get(start)?.point))
        })
        .collect::<Option<Vec<_>>>()?;
    let area = (0..points.len()).fold(Vec3::ZERO, |sum, index| {
        sum + points[index].cross(points[(index + 1) % points.len()])
    });
    Some(area.dot(normal).abs() * 0.5)
}

fn rotationally_invariant(wires: &[Vec<Curve>]) -> bool {
    wires.iter().all(|wire| {
        let mut radius: Option<f64> = None;
        let mut angle = 0.0;
        for piece in wire {
            let Curve::Arc(arc) = piece else { return false; };
            if arc.centre[0].hypot(arc.centre[1]) > arc.radius.abs().max(1.0) * 1e-10 { return false; }
            if radius.is_some_and(|value| (value - arc.radius).abs() > value.abs().max(1.0) * 1e-10) { return false; }
            radius = Some(arc.radius);
            angle += arc.sweep();
        }
        (angle - TAU).abs() <= 1e-9
    })
}

const GAUSS: [(f64, f64); 5] = [
    (-0.906179845938664, 0.236926885056189),
    (-0.538469310105683, 0.478628670499366),
    (0.0, 0.568888888888889),
    (0.538469310105683, 0.478628670499366),
    (0.906179845938664, 0.236926885056189),
];

struct Wire {
    source: Vec<Curve>,
    curves: Vec<RationalCurve2>,
    closed: bool,
}

fn expanded(curves: &[Curve]) -> Option<Vec<Curve>> {
    let mut result = Vec::new();
    for curve in curves {
        match curve {
            Curve::Polyline(polyline) => {
                for (index, segment) in curve.segments().into_iter().enumerate() {
                    if near2(polyline.vertices[index].position, segment.point_at(0.0)) {
                        result.push(segment);
                    } else {
                        // Clockwise bulges are represented by a reversed
                        // CCW arc. Preserve the entity's directed start even
                        // when it contains only this one segment.
                        let rational = RationalCurve2::from_curve(&segment)?.reversed();
                        result.push(Curve::Nurbs(NurbsCurve::new_strict(rational.degree,
                            rational.points, rational.knots, rational.weights)?));
                    }
                }
            }
            Curve::Circle(circle) => {
                result.push(Curve::Arc(Arc { centre: circle.centre, radius: circle.radius,
                    start_angle: 0.0, end_angle: TAU }));
            }
            Curve::Ray(_) | Curve::XLine(_) => return None,
            _ => result.push(curve.clone()),
        }
    }
    (!result.is_empty()).then_some(result)
}

fn chain_senses(pieces: &[Curve]) -> Option<Vec<bool>> {
    [true, false].into_iter().find_map(|first| {
        let mut senses = vec![first];
        let mut point = pieces.first()?.point_at(if first { 1.0 } else { 0.0 });
        for piece in &pieces[1..] {
            let forward = near2(point, piece.point_at(0.0));
            if !forward && !near2(point, piece.point_at(1.0)) { return None; }
            senses.push(forward);
            point = piece.point_at(if forward { 1.0 } else { 0.0 });
        }
        Some(senses)
    })
}

fn near2(a: [f64; 2], b: [f64; 2]) -> bool {
    (a[0] - b[0]).hypot(a[1] - b[1]) <= 1e-8
        + a.into_iter().chain(b).map(f64::abs).fold(1.0_f64, f64::max) * f64::EPSILON * 64.0
}

fn chain_closed(pieces: &[Curve], senses: &[bool]) -> bool {
    near2(pieces[0].point_at(if senses[0] { 0.0 } else { 1.0 }),
        pieces.last().unwrap().point_at(if *senses.last().unwrap() { 1.0 } else { 0.0 }))
}

fn prepare_wires(wires: &[Vec<Curve>]) -> Option<Vec<Wire>> {
    if wires.is_empty() { return None; }
    wires.iter().enumerate().map(|(index, source)| {
        let source = expanded(source)?;
        let senses = chain_senses(&source)?;
        let closed = chain_closed(&source, &senses);
        let origin = source[0].point_at(0.0);
        let shift = Transform::translation([-origin[0], -origin[1]]);
        let area = source.iter().zip(&senses)
            .map(|(piece, forward)| Some(piece.transformed(&shift)?.enclosed_area()
                * if *forward { 1.0 } else { -1.0 }))
            .collect::<Option<Vec<_>>>()?.into_iter().sum::<f64>();
        if closed && (!area.is_finite() || area.abs() <= 1e-14) { return None; }
        let mut curves = source.iter().zip(senses).map(|(piece, forward)| {
            let curve = RationalCurve2::from_curve(piece)?;
            Some(if forward { curve } else { curve.reversed() })
        }).collect::<Option<Vec<_>>>()?;
        if closed && (area > 0.0) != (index == 0) {
            curves = curves.into_iter().rev().map(|curve| curve.reversed()).collect();
        }
        Some(Wire { source, curves, closed })
    }).collect()
}

#[derive(Clone)]
enum Piece {
    Line(Vec3, Vec3),
    Planar(Plane, Curve, bool),
    Spline(NurbsCurve3, f64, f64),
}

impl Piece {
    fn point(&self, t: f64) -> Vec3 {
        Vec3::from(match self {
            Self::Line(a, b) => a.lerp(*b, t).to_array(),
            Self::Planar(plane, curve, forward) =>
                plane.point_at(curve.point_at(if *forward { t } else { 1.0 - t })),
            Self::Spline(curve, a, b) => curve.point_at_knot(a + (b - a) * t),
        })
    }

    fn tangent(&self, t: f64) -> Option<Vec3> {
        let direction = match self {
            Self::Line(a, b) => *b - *a,
            Self::Planar(plane, curve, forward) => {
                Vec3::from(plane.vector_at(curve.tangent_at(if *forward { t } else { 1.0 - t })))
                    * if *forward { 1.0 } else { -1.0 }
            }
            Self::Spline(_, _, _) => self.spline_derivative(t),
        };
        direction.is_finite().then_some(())?;
        if direction.length() > 1e-12 { return direction.normalize(); }
        let delta = self.point((t + 1e-6).min(1.0)) - self.point((t - 1e-6).max(0.0));
        delta.normalize()
    }

    fn speed(&self, t: f64) -> f64 {
        match self {
            Self::Line(a, b) => a.distance(*b),
            Self::Planar(plane, curve, forward) => Vec3::from(plane.vector_at(
                curve.tangent_at(if *forward { t } else { 1.0 - t }))).length(),
            Self::Spline(_, _, _) => self.spline_derivative(t).length(),
        }
    }

    fn spline_derivative(&self, t: f64) -> Vec3 {
        // Stay inside this knot span. A global central difference straddles
        // repeated knots and rounds a deliberately sharp polyline corner.
        let a = (t - 1e-5).max(0.0);
        let b = (t + 1e-5).min(1.0);
        (self.point(b) - self.point(a)) / (b - a)
    }

    fn length_to(&self, end: f64) -> f64 {
        if let Self::Line(a, b) = self { return a.distance(*b) * end; }
        if let Self::Planar(plane, curve, _) = self {
            if let Curve::Line(line) = curve {
                return Vec3::from(plane.vector_at(line.direction())).length() * end;
            }
            if matches!(curve, Curve::Arc(_)) && plane.is_orthonormal() {
                return curve.length() * end;
            }
        }
        (0..16).map(|panel| GAUSS.into_iter().map(|(node, weight)| {
            let t = end * (panel as f64 + 0.5 + node * 0.5) / 16.0;
            self.speed(t) * weight * end / 32.0
        }).sum::<f64>()).sum()
    }

    fn length(&self) -> f64 { self.length_to(1.0) }
}

fn path_pieces(path: SweepPath<'_>) -> Option<Vec<Piece>> {
    let pieces = match path {
        SweepPath::Planar { plane, curves } => {
            plane.normal()?;
            let curves = expanded(curves)?;
            let senses = chain_senses(&curves)?;
            curves.into_iter().zip(senses).map(|(curve, forward)| match curve {
                Curve::Line(line) => Piece::Line(
                    Vec3::from(plane.point_at(if forward { line.start } else { line.end })),
                    Vec3::from(plane.point_at(if forward { line.end } else { line.start }))),
                _ => Piece::Planar(plane, curve, forward),
            }).collect::<Vec<_>>()
        }
        SweepPath::Polyline3d { points, closed } => {
            if points.len() < 2 || !points.iter().flatten().all(|v| v.is_finite()) { return None; }
            let mut pieces = points.windows(2).map(|pair|
                Piece::Line(Vec3::from(pair[0]), Vec3::from(pair[1]))).collect::<Vec<_>>();
            if closed && Vec3::from(points[0]).distance(Vec3::from(*points.last()?)) > 1e-9 {
                pieces.push(Piece::Line(Vec3::from(*points.last()?), Vec3::from(points[0])));
            }
            pieces
        }
        SweepPath::Nurbs3(curve) => {
            let (a, b) = curve.domain();
            let mut knots = curve.knots().iter().copied().filter(|k| *k >= a && *k <= b).collect::<Vec<_>>();
            knots.dedup_by(|a, b| (*a - *b).abs() <= 1e-14);
            knots.windows(2).map(|span| Piece::Spline(curve.clone(), span[0], span[1])).collect()
        }
    };
    if pieces.is_empty() || pieces.iter().any(|piece| !piece.length().is_finite() || piece.length() <= 1e-12) {
        return None;
    }
    Some(pieces)
}

#[derive(Clone, Copy)]
struct Frame { origin: Vec3, x: Vec3, y: Vec3 }

impl Frame {
    fn plane(self) -> Plane { Plane::from_axes(self.origin.to_array(), self.x.to_array(), self.y.to_array()) }
    fn point(self, p: [f64; 2]) -> Vec3 { self.origin + self.x * p[0] + self.y * p[1] }
    fn plus(self, other: Self) -> Self { Self { origin: self.origin + other.origin, x: self.x + other.x, y: self.y + other.y } }
    fn minus(self, other: Self) -> Self { self.plus(other.times(-1.0)) }
    fn times(self, t: f64) -> Self { Self { origin: self.origin * t, x: self.x * t, y: self.y * t } }
    fn lerp(self, other: Self, t: f64) -> Self { self.plus(other.minus(self).times(t)) }
}

fn rotate(value: Vec3, axis: Vec3, angle: f64) -> Vec3 {
    let (sin, cos) = angle.sin_cos();
    value * cos + axis.cross(value) * sin + axis * axis.dot(value) * (1.0 - cos)
}

fn transport(value: Vec3, from: Vec3, to: Vec3) -> Option<Vec3> {
    let cross = from.cross(to);
    let sin = cross.length();
    let cos = from.dot(to).clamp(-1.0, 1.0);
    if sin <= 1e-12 {
        if cos >= 0.0 { return Some(value); }
        let axis = if from.x.abs() < 0.8 { from.cross(Vec3::X) } else { from.cross(Vec3::Y) }.normalize()?;
        return Some(rotate(value, axis, PI));
    }
    Some(rotate(value, cross / sin, sin.atan2(cos)))
}

fn moved(frame: Frame, from: Vec3, to: Vec3, old_point: Vec3, point: Vec3) -> Option<Frame> {
    Some(Frame { origin: point + transport(frame.origin - old_point, from, to)?,
        x: transport(frame.x, from, to)?, y: transport(frame.y, from, to)? })
}

fn curved_transport(frame: Frame, piece: &Piece, a: f64, b: f64, bank: bool) -> Option<Frame> {
    let from = piece.tangent(a)?;
    let to = piece.tangent(b)?;
    let point = piece.point(b);
    let mut result = moved(frame, from, to, piece.point(a), point)?;
    if bank {
        let binormal = |t: f64| -> Option<Vec3> {
            let start = piece.tangent((t - 1e-3).max(0.0))?;
            let end = piece.tangent((t + 1e-3).min(1.0))?;
            let cross = start.cross(end);
            (cross.length() > 1e-9).then(|| cross.normalize()).flatten()
        };
        if let (Some(previous), Some(next)) = (binormal(a), binormal(b)) {
            let previous = transport(previous, from, to)?;
            let roll = to.dot(previous.cross(next)).atan2(previous.dot(next));
            result = Frame { origin: point + rotate(result.origin - point, to, roll),
                x: rotate(result.x, to, roll), y: rotate(result.y, to, roll) };
        }
    }
    Some(result)
}

fn miter(frame: Frame, point: Vec3, tangent: Vec3, other: Vec3) -> Option<Frame> {
    let normal = (tangent + other).normalize()?;
    let dot = normal.dot(tangent);
    if dot <= 1e-6 { return None; }
    let cut = |v: Vec3| v - tangent * (v.dot(normal) / dot);
    Some(Frame { origin: point + cut(frame.origin - point), x: cut(frame.x), y: cut(frame.y) })
}

fn twist_frame(frame: Frame, point: Vec3, angle: f64, scale: f64) -> Option<Frame> {
    let normal = frame.x.cross(frame.y).normalize()?;
    Some(Frame { origin: point + rotate(frame.origin - point, normal, angle) * scale,
        x: rotate(frame.x, normal, angle) * scale, y: rotate(frame.y, normal, angle) * scale })
}

/// Cubic Bezier control frames of one path patch, with their weights: one
/// for transport fits, rational for an exact turn about a circular run.
#[derive(Clone, Copy)]
struct Patch { frames: [Frame; 4], weights: [f64; 4] }

impl std::ops::Index<usize> for Patch {
    type Output = Frame;
    fn index(&self, index: usize) -> &Frame { &self.frames[index] }
}

impl std::ops::IndexMut<usize> for Patch {
    fn index_mut(&mut self, index: usize) -> &mut Frame { &mut self.frames[index] }
}

/// The circle a planar piece runs round: centre, unit axis (turning the
/// piece positively) and the angle it turns through.
fn circular_run(piece: &Piece) -> Option<(Vec3, Vec3, f64)> {
    let Piece::Planar(plane, curve, _) = piece else { return None };
    if !matches!(curve, Curve::Arc(_) | Curve::Nurbs(_)) { return None; }
    let (a, b, c) = (piece.point(0.0), piece.point(0.5), piece.point(1.0));
    let (ab, ac) = (b - a, c - a);
    let normal = ab.cross(ac);
    if normal.length() <= 1e-12 * ab.length().max(ac.length()).max(1.0).powi(2) { return None; }
    let centre = a + (normal.cross(ab) * ac.dot(ac) + ac.cross(normal) * ab.dot(ab)) * (0.5 / normal.dot(normal));
    let radius = a.distance(centre);
    let tolerance = radius.max(1.0) * 1e-10;
    if (0..=16).any(|step| (piece.point(step as f64 / 16.0).distance(centre) - radius).abs() > tolerance) {
        return None;
    }
    let plane_normal = Vec3::from(plane.normal()?);
    let axis = if (a - centre).cross(piece.tangent(0.0)?).dot(plane_normal) >= 0.0 { plane_normal } else { -plane_normal };
    let (from, to) = (a - centre, c - centre);
    let mut angle = axis.dot(from.cross(to)).atan2(from.dot(to));
    // The middle point decides whether the run goes the long way round.
    if angle <= 0.0 || axis.dot((b - centre).cross(to)) < 0.0 { angle += TAU; }
    if angle > TAU { angle -= TAU; }
    (angle > 1e-9).then_some((centre, axis, angle))
}

/// Exact rational patches turning `start` about the axis through `centre`,
/// each a quarter turn at most.
fn turning_patches(start: Frame, centre: Vec3, axis: Vec3, angle: f64, result: &mut Vec<Patch>) {
    let spans = (angle / (PI / 2.0) - 1e-9).ceil().max(1.0) as usize;
    let step = angle / spans as f64;
    let turned = |frame: Frame, by: f64| Frame {
        origin: centre + rotate(frame.origin - centre, axis, by),
        x: rotate(frame.x, axis, by),
        y: rotate(frame.y, axis, by),
    };
    // The middle control of a rational quadratic arc lies 1/cos(step/2) out
    // from the axis; along the axis nothing moves.
    let cosine = (step * 0.5).cos();
    let spread = |value: Vec3| axis * axis.dot(value) + (value - axis * axis.dot(value)) * (1.0 / cosine);
    for span in 0..spans {
        let first = turned(start, step * span as f64);
        let last = turned(start, step * (span + 1) as f64);
        let middle = turned(start, step * (span as f64 + 0.5));
        let corner = Frame {
            origin: centre + spread(middle.origin - centre),
            x: spread(middle.x),
            y: spread(middle.y),
        };
        // Degree elevation of weights (1, cos, 1).
        let inner = 1.0 + 2.0 * cosine;
        result.push(Patch {
            frames: [
                first,
                first.plus(corner.times(2.0 * cosine)).times(1.0 / inner),
                corner.times(2.0 * cosine).plus(last).times(1.0 / inner),
                last,
            ],
            weights: [1.0, inner / 3.0, inner / 3.0, 1.0],
        });
    }
}

fn transported_patches(
    pieces: &[Piece], first: Frame, options: SweepOptions,
    total: f64, radius: f64, closed: bool,
) -> Option<Vec<Patch>> {
    let mut frame = first;
    let mut previous_point = pieces[0].point(0.0);
    let mut previous_tangent = pieces[0].tangent(0.0)?;
    // Transport through small angular increments, also on nonplanar splines.
    let mut walks = Vec::new();
    for piece in pieces {
        let tangent = piece.tangent(0.0)?;
        if previous_tangent.dot(tangent) <= -1.0 + 1e-10 { return None; }
        frame = moved(frame, previous_tangent, tangent, previous_point, piece.point(0.0))?;
        let mut parameters = vec![0.0];
        divide_path(piece, 0.0, 1.0, 0, &mut parameters)?;
        let mut frames = vec![frame];
        for span in parameters.windows(2) {
            frame = curved_transport(frame, piece, span[0], span[1], options.bank)?;
            frames.push(frame);
        }
        previous_point = piece.point(1.0);
        previous_tangent = piece.tangent(1.0)?;
        walks.push((parameters, frames));
    }
    let closure_roll = if closed {
        let last = moved(frame, previous_tangent, pieces[0].tangent(0.0)?, previous_point, pieces[0].point(0.0))?;
        let axis = pieces[0].tangent(0.0)?;
        let a = (last.x - axis * last.x.dot(axis)).normalize()?;
        let b = (first.x - axis * first.x.dot(axis)).normalize()?;
        axis.dot(a.cross(b)).atan2(a.dot(b))
    } else { 0.0 };
    let mut patches = Vec::new();
    let mut travelled = 0.0;
    for (index, piece) in pieces.iter().enumerate() {
        let length = piece.length();
        let (parameters, frames) = &walks[index];
        let evaluate_raw = |t: f64| -> Option<Frame> {
            let slot = parameters.partition_point(|parameter| *parameter <= t).saturating_sub(1).min(parameters.len() - 2);
            let raw = curved_transport(frames[slot], piece, parameters[slot], t, options.bank)?;
            let fraction = (travelled + piece.length_to(t)) / total;
            twist_frame(raw, piece.point(t), options.twist * fraction + closure_roll * fraction,
                1.0 + (options.scale - 1.0) * fraction)
        };
        let raw_start = evaluate_raw(0.0)?;
        let raw_end = evaluate_raw(1.0)?;
        let previous = if index > 0 { Some(&pieces[index - 1]) } else if closed { pieces.last() } else { None };
        let next = if index + 1 < pieces.len() { Some(&pieces[index + 1]) } else if closed { pieces.first() } else { None };
        let cut_start = if let Some(previous) = previous {
            miter(raw_start, piece.point(0.0), piece.tangent(0.0)?, previous.tangent(1.0)?)?
        } else { raw_start };
        let cut_end = if let Some(next) = next {
            miter(raw_end, piece.point(1.0), piece.tangent(1.0)?, next.tangent(0.0)?)?
        } else { raw_end };
        let delta_start = cut_start.minus(raw_start);
        let delta_end = cut_end.minus(raw_end);
        let evaluate = |t: f64| -> Option<Frame> {
            Some(evaluate_raw(t)?.plus(delta_start.lerp(delta_end, t)))
        };
        let subdivisions = ((options.twist.abs() + closure_roll.abs()) * length / total / 0.1).ceil().max(1.0) as usize;
        if subdivisions > 8192 { return None; }
        let mut cuts = parameters.clone();
        cuts.extend((1..subdivisions).map(|i| i as f64 / subdivisions as f64));
        cuts.sort_by(f64::total_cmp);
        cuts.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        // An untwisted, unscaled run round a circle that meets its neighbours
        // without a corner is an exact turn: keep it exact rather than fitted.
        let delta = delta_start.origin.length() + delta_start.x.length() + delta_start.y.length()
            + delta_end.origin.length() + delta_end.x.length() + delta_end.y.length();
        let exact = (options.twist.abs() + closure_roll.abs() <= 1e-12
            && (options.scale - 1.0).abs() <= 1e-12 && delta <= 1e-12)
            .then(|| circular_run(piece)).flatten();
        if let Some((centre, axis, angle)) = exact {
            turning_patches(evaluate(0.0)?, centre, axis, angle, &mut patches);
        } else {
            for span in cuts.windows(2) {
                fit_patch(&evaluate, span[0], span[1], radius, total.max(radius) * 1e-7, 0, &mut patches)?;
            }
        }
        travelled += length;
    }
    // Shared topology requires exactly the same section on both sides of a
    // path corner. Minimal transport plus the bisector cut gives that map;
    // reject banking singularities rather than hiding a discontinuity.
    for index in 1..patches.len() {
        if frame_error(patches[index - 1][3], patches[index][0], radius) > total.max(radius) * 1e-6 { return None; }
        patches[index][0] = patches[index - 1][3];
    }
    if closed {
        let start = patches.first()?[0];
        if frame_error(patches.last()?[3], start, radius) > total.max(radius) * 1e-6 { return None; }
        patches.last_mut()?[3] = start;
    }
    Some(patches)
}

fn divide_path(piece: &Piece, a: f64, b: f64, depth: usize, result: &mut Vec<f64>) -> Option<()> {
    let middle = (a + b) * 0.5;
    let ta = piece.tangent(a)?;
    let tb = piece.tangent(b)?;
    let tm = piece.tangent(middle)?;
    let chord = piece.point(a).distance(piece.point(b));
    let broken = piece.point(a).distance(piece.point(middle)) + piece.point(middle).distance(piece.point(b));
    if ta.dot(tm) < 0.99875 || tm.dot(tb) < 0.99875 || broken - chord > broken.max(1.0) * 0.0001 {
        if depth >= 16 || result.len() >= 8192 { return None; }
        divide_path(piece, a, middle, depth + 1, result)?;
        divide_path(piece, middle, b, depth + 1, result)?;
    } else {
        result.push(b);
    }
    Some(())
}

fn bezier(p: &Patch, t: f64) -> Frame {
    let s = 1.0 - t;
    let basis = [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t];
    let weight = (0..4).map(|i| basis[i] * p.weights[i]).sum::<f64>();
    (0..4).fold(p[0].times(0.0), |sum, i| sum.plus(p[i].times(basis[i] * p.weights[i] / weight)))
}

fn bezier_derivative(p: &Patch, t: f64) -> Frame {
    let s = 1.0 - t;
    let basis = [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t];
    let slope = [-3.0 * s * s, 3.0 * s * s - 6.0 * s * t, 6.0 * s * t - 3.0 * t * t, 3.0 * t * t];
    let weight = (0..4).map(|i| basis[i] * p.weights[i]).sum::<f64>();
    let weight_slope = (0..4).map(|i| slope[i] * p.weights[i]).sum::<f64>();
    let point = bezier(p, t);
    (0..4).fold(point.times(-weight_slope / weight), |sum, i| {
        sum.plus(p[i].times(slope[i] * p.weights[i] / weight))
    })
}

/// A locally folded transport is not a regular swept body, even when its
/// shared-edge topology happens to be manifold. Detect its signed volume
/// Jacobian before handing a singular skin to adaptive face meshing. This
/// checks local regularity, not distant intersections between separate runs.
/// A consistently negative Jacobian remains regular for an off-path anchor;
/// return that orientation so the complete shell can be turned accordingly.
fn regular_transport(wires: &[Wire], patches: &[Patch]) -> Option<bool> {
    let mut points = Vec::new();
    for curve in wires.iter().flat_map(|wire| &wire.curves) {
        let spline = NurbsCurve::new_strict(curve.degree, curve.points.clone(),
            curve.knots.clone(), curve.weights.clone())?;
        let mut cuts = vec![0.0, 1.0];
        cuts.extend(curve.knots.iter().copied().filter(|knot| *knot > 0.0 && *knot < 1.0));
        cuts.sort_by(f64::total_cmp);
        cuts.dedup_by(|a, b| (*a - *b).abs() <= 1e-12);
        let samples = (curve.degree * 4).max(8);
        for span in cuts.windows(2) {
            for sample in 0..=samples {
                points.push(spline.point_at(span[0] + (span[1] - span[0]) * sample as f64 / samples as f64));
            }
        }
    }
    let mut orientation = None;
    for patch in patches {
        for at in [0.0, 0.125, 0.25, 0.5, 0.75, 0.875, 1.0] {
            let frame = bezier(patch, at);
            let derivative = bezier_derivative(patch, at);
            let normal = frame.x.cross(frame.y);
            let normal_length = normal.length();
            if !normal_length.is_finite() || normal_length <= 0.0 { return None; }
            for point in &points {
                let velocity = derivative.point(*point);
                let jacobian = normal.dot(velocity);
                let tolerance = normal_length * velocity.length() * 1e-10;
                if !jacobian.is_finite() || !tolerance.is_finite() || jacobian.abs() <= tolerance {
                    return None;
                }
                let forward = jacobian > 0.0;
                if orientation.is_some_and(|previous| previous != forward) { return None; }
                orientation = Some(forward);
            }
        }
    }
    orientation
}

fn frame_error(a: Frame, b: Frame, radius: f64) -> f64 {
    a.origin.distance(b.origin) + radius * (a.x.distance(b.x) + a.y.distance(b.y))
}

fn fit_patch(
    evaluate: &impl Fn(f64) -> Option<Frame>, a: f64, b: f64,
    radius: f64, tolerance: f64, depth: usize, result: &mut Vec<Patch>,
) -> Option<()> {
    let p0 = evaluate(a)?;
    let p3 = evaluate(b)?;
    let q1 = evaluate(a + (b - a) / 3.0)?;
    let q2 = evaluate(a + (b - a) * 2.0 / 3.0)?;
    let c = q1.times(27.0).minus(p0.times(8.0)).minus(p3);
    let d = q2.times(27.0).minus(p0).minus(p3.times(8.0));
    let patch = Patch { frames: [p0, c.times(2.0).minus(d).times(1.0 / 18.0),
        d.times(2.0).minus(c).times(1.0 / 18.0), p3], weights: [1.0; 4] };
    let mut error = 0.0_f64;
    for t in [0.125, 0.25, 0.5, 0.75, 0.875] {
        error = error.max(frame_error(bezier(&patch, t), evaluate(a + (b - a) * t)?, radius));
    }
    if !error.is_finite() { return None; }
    if error > tolerance {
        if depth >= 14 || result.len() >= 8192 { return None; }
        let middle = (a + b) * 0.5;
        fit_patch(evaluate, a, middle, radius, tolerance, depth + 1, result)?;
        fit_patch(evaluate, middle, b, radius, tolerance, depth + 1, result)?;
    } else {
        result.push(patch);
    }
    Some(())
}

fn build_body(wires: &[Wire], patches: &[Patch], sheet: bool, closed_path: bool, outward: bool) -> Option<Body> {
    if patches.is_empty() { return None; }
    let mut body = Body::new();
    let lump = body.lumps.insert(Lump { shells: Vec::new(), provenance: Provenance::Synthesized });
    let shell = body.shells.insert(Shell { faces: Vec::new(), owner: lump, provenance: Provenance::Synthesized });
    let mut cap_rims = Vec::new();
    for wire in wires {
        let mut coords = wire.curves.iter().map(|curve| rational_point(curve, 0.0)).collect::<Option<Vec<_>>>()?;
        if !wire.closed { coords.push(rational_point(wire.curves.last()?, 1.0)?); }
        let mut vertices: Vec<Vec<VertexKey>> = Vec::new();
        let mut rims: Vec<Vec<EdgeKey>> = Vec::new();
        for station in 0..=patches.len() {
            if closed_path && station == patches.len() {
                vertices.push(vertices[0].clone());
                rims.push(rims[0].clone());
                continue;
            }
            let frame = if station == 0 { patches[0][0] } else { patches[station - 1][3] };
            let ring = coords.iter().map(|p| body.vertices.insert(Vertex {
                point: frame.point(*p).to_array(), provenance: Provenance::Synthesized,
            })).collect::<Vec<_>>();
            let edges = wire.curves.iter().enumerate().map(|(index, curve)| {
                let lifted = curve.lifted(&frame.plane());
                add_curve_edge(&mut body, &lifted, ring[index], ring[(index + 1) % ring.len()])
            }).collect::<Option<Vec<_>>>()?;
            vertices.push(ring);
            rims.push(edges);
        }
        for (band, patch) in patches.iter().enumerate() {
            let rails = coords.iter().enumerate().map(|(index, point)| {
                let curve = RationalCurve3 { degree: 3, knots: cubic_knots(),
                    points: patch.frames.iter().map(|frame| frame.point(*point).to_array()).collect(),
                    weights: patch.weights.to_vec() };
                add_curve_edge(&mut body, &curve, vertices[band][index], vertices[band + 1][index])
            }).collect::<Option<Vec<_>>>()?;
            for (index, curve) in wire.curves.iter().enumerate() {
                let next = (index + 1) % coords.len();
                let points = curve.points.iter().map(|point|
                    patch.frames.iter().map(|frame| frame.point(*point).to_array()).collect()).collect();
                let weights = curve.weights.iter()
                    .map(|weight| patch.weights.iter().map(|along| weight * along).collect())
                    .collect();
                let surface = NurbsSurface3::new_strict(curve.degree, 3, points,
                    curve.knots.clone(), cubic_knots(), weights)?;
                let surface = body.surfaces.insert(Surface::Nurbs(surface));
                let mut circuit = vec![(rims[band][index], true), (rails[next], true),
                    (rims[band + 1][index], false), (rails[index], false)];
                let mut pcurves = vec![([0.0, 0.0], [1.0, 0.0]), ([1.0, 0.0], [1.0, 1.0]),
                    ([1.0, 1.0], [0.0, 1.0]), ([0.0, 1.0], [0.0, 0.0])];
                if !outward {
                    circuit = circuit.into_iter().rev().map(|(edge, forward)| (edge, !forward)).collect();
                    pcurves = pcurves.into_iter().rev().map(|(a, b)| (b, a)).collect();
                }
                let face = add_face(&mut body, shell, surface, outward);
                add_loop(&mut body, face, &circuit, Some(&pcurves))?;
            }
        }
        cap_rims.push((rims[0].clone(), rims.last()?.clone()));
    }
    if !sheet && !closed_path {
        for end in [false, true] {
            let frame = if end { patches.last()?[3] } else { patches[0][0] };
            let surface = body.surfaces.insert(Surface::Plane(frame.plane()));
            let forward = if end { outward } else { !outward };
            let face = add_face(&mut body, shell, surface, forward);
            for (wire, (start_edges, end_edges)) in wires.iter().zip(&cap_rims) {
                let edges = if end { end_edges } else { start_edges };
                let mut circuit = edges.iter().map(|edge| (*edge, true)).collect::<Vec<_>>();
                if !forward { circuit = circuit.into_iter().rev().map(|(edge, _)| (edge, false)).collect(); }
                let ring = add_loop(&mut body, face, &circuit, None)?;
                let curves = if forward { wire.curves.clone() } else {
                    wire.curves.iter().rev().map(RationalCurve2::reversed).collect()
                };
                for (coedge, curve) in body.loops.get(ring)?.coedges.clone().into_iter().zip(curves) {
                    body.coedges.get_mut(coedge)?.pcurve = Some(Curve::Nurbs(NurbsCurve::new_strict(
                        curve.degree, curve.points, curve.knots, curve.weights)?));
                }
            }
        }
    }
    body.lumps.get_mut(lump)?.shells = vec![shell];
    body.roots = vec![lump];
    body.validate().is_empty().then_some(body)
}

fn rational_point(curve: &RationalCurve2, parameter: f64) -> Option<[f64; 2]> {
    // Unclamped/periodic splines do not end at their outer control points.
    // Evaluate the actual boundary so rims and swept rails share endpoints.
    Some(NurbsCurve::new_strict(curve.degree, curve.points.clone(), curve.knots.clone(),
        curve.weights.clone())?.point_at(parameter))
}

fn cubic_knots() -> Vec<f64> { vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0] }

fn add_curve_edge(body: &mut Body, source: &RationalCurve3, start: VertexKey, end: VertexKey) -> Option<EdgeKey> {
    let curve = body.curves.insert(Curve3::Nurbs(source.curve()?));
    Some(body.edges.insert(Edge { curve, start_parameter: 0.0, end_parameter: 1.0, start, end,
        coedges: Vec::new(), provenance: Provenance::Synthesized }))
}

fn add_face(body: &mut Body, shell: ShellKey, surface: super::SurfaceKey, forward: bool) -> FaceKey {
    let face = body.faces.insert(Face { surface, forward, loops: Vec::new(), owner: shell, provenance: Provenance::Synthesized });
    body.shells.get_mut(shell).unwrap().faces.push(face);
    face
}

fn add_loop(body: &mut Body, face: FaceKey, circuit: &[(EdgeKey, bool)],
    pcurves: Option<&[([f64; 2], [f64; 2])]>,
) -> Option<super::LoopKey> {
    let ring = body.loops.insert(Loop { coedges: Vec::new(), owner: face, provenance: Provenance::Synthesized });
    for (index, (edge, forward)) in circuit.iter().enumerate() {
        let pcurve = pcurves.map(|curves| Curve::Line(Line { start: curves[index].0, end: curves[index].1 }));
        let coedge = body.coedges.insert(Coedge { edge: *edge, forward: *forward, pcurve, owner: ring,
            provenance: Provenance::Synthesized });
        body.edges.get_mut(*edge)?.coedges.push(coedge);
        body.loops.get_mut(ring)?.coedges.push(coedge);
    }
    body.faces.get_mut(face)?.loops.push(ring);
    Some(ring)
}
