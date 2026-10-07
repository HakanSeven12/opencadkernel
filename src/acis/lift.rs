//! Reading an ACIS document into kernel topology.
//!
//! Unsupported records remain attached through provenance and appear in
//! [`Loss`].

use cadcodec::entities::acis::types::{
    SatBody, SatCoedge, SatConeSurface, SatDocument, SatEdge, SatEllipseCurve, SatFace, SatIntCurve,
    SatLoop, SatLump, SatPCurve, SatPlaneSurface, SatPoint, SatPointer, SatRecord, SatShell,
    SatSphereSurface, SatSplineSurface, SatStraightCurve, SatTorusSurface, SatVertex, Sense,
};
use crate::brep::{
    Body, Circle3, Coedge, Cone, Curve3, CurveKey, Cylinder, Edge, EdgeKey, EllipticCone, Face,
    Line3, Loop, Lump, Provenance, Shell, SourceRef, Sphere, Surface, SurfaceKey, Torus, Vertex,
    VertexKey,
};
use crate::geom2d::{Curve as Curve2, NurbsCurve};
use crate::space::{NurbsCurve3, NurbsSurface3, Plane, Vec3};
use std::collections::HashMap;
use std::f64::consts::TAU;

/// What a lift could not represent.
///
/// Empty means the whole document is in the kernel's own terms and can be
/// edited freely. Anything listed is a node whose record is carried through
/// verbatim; an edit that touches one loses whatever the kernel does not
/// model about it.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Loss {
    /// Records whose surface kind has no kernel equivalent, by index.
    pub surfaces: Vec<usize>,
    /// Records whose curve kind has no kernel equivalent, by index.
    pub curves: Vec<usize>,
    /// Records the pointer graph named but that are missing or malformed.
    pub broken: Vec<usize>,
}

impl Loss {
    /// Whether the document lifted completely.
    pub fn is_empty(&self) -> bool {
        self.surfaces.is_empty() && self.curves.is_empty() && self.broken.is_empty()
    }
}

/// Every body in the document, in the order they appear.
pub fn lift(document: &SatDocument) -> (Vec<Body>, Loss) {
    let mut loss = Loss::default();
    // The record index each body came from, so its own provenance is set:
    // `bodies()` hands back views without saying which record each was.
    let indices: Vec<Option<u32>> = document
        .records_of_type("body")
        .iter()
        .map(|record| index_of(record))
        .collect();
    let bodies = document
        .bodies()
        .into_iter()
        .enumerate()
        .filter_map(|(order, body)| {
            lift_one(document, &body, indices.get(order).copied().flatten(), &mut loss)
        })
        .collect();
    (bodies, loss)
}

/// One body, named by the record index of its `body` record.
pub fn lift_body(document: &SatDocument, record: usize) -> Option<(Body, Loss)> {
    let source = document.record(record)?;
    let body = SatBody::from_record(source)?;
    let mut loss = Loss::default();
    let lifted = lift_one(document, &body, index_of(source), &mut loss)?;
    Some((lifted, loss))
}

/// Everything already built, so a record shared by several nodes becomes one
/// node rather than several copies.
#[derive(Default)]
struct Seen {
    vertices: HashMap<u32, VertexKey>,
    edges: HashMap<u32, EdgeKey>,
    curves: HashMap<u32, CurveKey>,
    surfaces: HashMap<u32, SurfaceKey>,
}

fn lift_one(
    document: &SatDocument,
    source: &SatBody<'_>,
    at: Option<u32>,
    loss: &mut Loss,
) -> Option<Body> {
    let mut body = Body::new();
    body.provenance = match at {
        Some(index) => Provenance::Clean(SourceRef::new(index)),
        None => Provenance::Synthesized,
    };
    let mut seen = Seen::default();

    // Lumps, shells, faces and loops are each a linked list rather than an
    // array: the record holds the first, and every one holds the next.
    let mut lump_pointer = source.lump();
    while let Some(record) = resolve(document, lump_pointer) {
        let Some(source_lump) = SatLump::from_record(record) else {
            note_broken(loss, record);
            break;
        };
        let lump = body.lumps.insert(Lump {
            shells: Vec::new(),
            provenance: clean(record),
        });
        body.roots.push(lump);

        let mut shell_pointer = source_lump.shell();
        while let Some(record) = resolve(document, shell_pointer) {
            let Some(source_shell) = SatShell::from_record(record) else {
                note_broken(loss, record);
                break;
            };
            let shell = body.shells.insert(Shell {
                faces: Vec::new(),
                owner: lump,
                provenance: clean(record),
            });
            body.lumps.get_mut(lump)?.shells.push(shell);

            let mut face_pointer = source_shell.face();
            while let Some(record) = resolve(document, face_pointer) {
                let Some(source_face) = SatFace::from_record(record) else {
                    note_broken(loss, record);
                    break;
                };
                lift_face(document, &mut body, &mut seen, loss, &source_face, record, shell);
                face_pointer = source_face.next_face();
            }
            shell_pointer = source_shell.next_shell();
        }
        lump_pointer = source_lump.next_lump();
    }

    (!body.roots.is_empty()).then_some(body)
}

fn lift_face(
    document: &SatDocument,
    body: &mut Body,
    seen: &mut Seen,
    loss: &mut Loss,
    source: &SatFace<'_>,
    record: &SatRecord,
    shell: crate::brep::ShellKey,
) -> Option<()> {
    let surface_record = resolve(document, source.surface())?;
    let reversed_v = analytic_surface_reversed(surface_record);
    let surface = surface_of(document, body, seen, loss, source.surface())?;
    let face = body.faces.insert(Face {
        surface,
        forward: (source.sense() == Sense::Forward) != reversed_v,
        loops: Vec::new(),
        owner: shell,
        provenance: clean(record),
    });
    body.shells.get_mut(shell)?.faces.push(face);

    let mut loop_pointer = source.first_loop();
    while let Some(record) = resolve(document, loop_pointer) {
        let Some(source_loop) = SatLoop::from_record(record) else {
            note_broken(loss, record);
            break;
        };
        let ring = body.loops.insert(Loop {
            coedges: Vec::new(),
            owner: face,
            provenance: clean(record),
        });
        body.faces.get_mut(face)?.loops.push(ring);

        // A loop's coedges are a ring joined by next pointers, so the walk
        // stops when it comes back to where it started rather than at a null.
        let first = source_loop.first_coedge();
        let mut pointer = first;
        let mut coedges = Vec::new();
        loop {
            let Some(record) = resolve(document, pointer) else {
                break;
            };
            let Some(source_coedge) = SatCoedge::from_record(record) else {
                note_broken(loss, record);
                break;
            };
            let pcurve = read_pcurve(document, source_coedge.pcurve(), reversed_v);
            if let Some(edge) = edge_of(
                document,
                body,
                seen,
                loss,
                &source_coedge,
                surface,
                pcurve.as_ref(),
            ) {
                let edge_forward = resolve(document, source_coedge.edge())
                    .and_then(SatEdge::from_record)
                    .is_some_and(|source_edge| {
                        let Some(edge) = body.edges.get(edge) else {
                            return false;
                        };
                        if edge.start == edge.end {
                            return source_edge.sense() == Sense::Forward;
                        }
                        let source_start = resolve(document, source_edge.start_vertex())
                            .and_then(index_of);
                        let kernel_start = body
                            .vertices
                            .get(edge.start)
                            .and_then(|vertex| vertex.provenance.source())
                            .map(|source| source.index() as u32);
                        source_start.is_some() && source_start == kernel_start
                    });
                let forward = (source_coedge.sense() == Sense::Forward) == edge_forward;
                let coedge = body.coedges.insert(Coedge {
                    edge,
                    forward,
                    pcurve,
                    owner: ring,
                    provenance: clean(record),
                });
                body.edges.get_mut(edge)?.coedges.push(coedge);
                coedges.push(coedge);
            }
            pointer = source_coedge.next();
            if pointer == first || pointer.is_null() {
                break;
            }
            // A malformed ring that never returns would spin here. The loop
            // cannot be longer than the document.
            if coedges.len() > document.record_count() {
                note_broken(loss, record);
                break;
            }
        }
        body.loops.get_mut(ring)?.coedges = coedges;
        loop_pointer = source_loop.next_loop();
    }
    Some(())
}

fn edge_of(
    document: &SatDocument,
    body: &mut Body,
    seen: &mut Seen,
    loss: &mut Loss,
    source_coedge: &SatCoedge<'_>,
    surface: SurfaceKey,
    pcurve: Option<&Curve2>,
) -> Option<EdgeKey> {
    let pointer = source_coedge.edge();
    let record = resolve(document, pointer)?;
    let index = index_of(record)?;
    if let Some(key) = seen.edges.get(&index) {
        return Some(*key);
    }
    let source = SatEdge::from_record(record).or_else(|| {
        note_broken(loss, record);
        None
    })?;
    let fallback = pcurve
        .and_then(|pcurve| surface_curve(body.surfaces.get(surface)?, pcurve))
        .or_else(|| partner_surface_curve(document, source_coedge));
    let curve = curve_of(document, body, seen, loss, source.curve(), fallback)?;
    let source_start = vertex_of(document, body, seen, loss, source.start_vertex())?;
    let source_end = vertex_of(document, body, seen, loss, source.end_vertex())?;
    let edge_forward = source.sense() == Sense::Forward;
    let (start, end) = if edge_forward {
        (source_start, source_end)
    } else {
        (source_end, source_start)
    };
    // Open curves have an unambiguous span, so their vertices resolve stale
    // stored parameters. Closed curves keep the stored choice of arc.
    let ends = (body.vertices.get(start)?.point, body.vertices.get(end)?.point);
    let apart = Vec3::from(ends.0).distance(Vec3::from(ends.1)) > 1e-9;
    let mut stored_low = source.start_param();
    let mut stored_high = source.end_param();
    if stored_high < stored_low {
        std::mem::swap(&mut stored_low, &mut stored_high);
    }
    let (low, high) = match body.curves.get(curve) {
        Some(shape @ Curve3::Line(_)) if apart => {
            (shape.parameter_at(ends.0), shape.parameter_at(ends.1))
        }
        Some(shape @ Curve3::Nurbs(nurbs))
            if start == end
                && nurbs.periodicity()
                && stored_high - stored_low >= (nurbs.domain().1 - nurbs.domain().0) * (1.0 - 1e-6) =>
        {
            let period = nurbs.domain().1 - nurbs.domain().0;
            let projected = shape.parameter_at(ends.0);
            let stored_gap = Vec3::from(shape.point_at(stored_low)).distance(Vec3::from(ends.0));
            let projected_gap = Vec3::from(shape.point_at(projected)).distance(Vec3::from(ends.0));
            let low = if projected_gap < stored_gap {
                projected
            } else {
                stored_low
            };
            (low, low + period)
        }
        Some(shape @ Curve3::Nurbs(curve)) if apart && !curve.periodicity() => {
            (shape.parameter_at(ends.0), shape.parameter_at(ends.1))
        }
        Some(shape @ (Curve3::Circle(_) | Curve3::Ellipse(_))) if apart => {
            let low = periodic_near(shape.parameter_at(ends.0), stored_low, TAU);
            let mut high = periodic_near(shape.parameter_at(ends.1), stored_high, TAU);
            while high <= low {
                high += TAU;
            }
            (low, high)
        }
        _ => (stored_low, stored_high),
    };
    let key = body.edges.insert(Edge {
        curve,
        start_parameter: low,
        end_parameter: high,
        start,
        end,
        coedges: Vec::new(),
        provenance: clean(record),
    });
    seen.edges.insert(index, key);
    Some(key)
}

fn periodic_near(value: f64, reference: f64, period: f64) -> f64 {
    value + period * ((reference - value) / period).round()
}

fn vertex_of(
    document: &SatDocument,
    body: &mut Body,
    seen: &mut Seen,
    loss: &mut Loss,
    pointer: SatPointer,
) -> Option<VertexKey> {
    let record = resolve(document, pointer)?;
    let index = index_of(record)?;
    if let Some(key) = seen.vertices.get(&index) {
        return Some(*key);
    }
    let source = SatVertex::from_record(record).or_else(|| {
        note_broken(loss, record);
        None
    })?;
    let point_record = resolve(document, source.point())?;
    let point = SatPoint::from_record(point_record).or_else(|| {
        note_broken(loss, point_record);
        None
    })?;
    let (x, y, z) = point.position();
    let key = body.vertices.insert(Vertex {
        point: [x, y, z],
        provenance: clean(record),
    });
    seen.vertices.insert(index, key);
    Some(key)
}

fn curve_of(
    document: &SatDocument,
    body: &mut Body,
    seen: &mut Seen,
    loss: &mut Loss,
    pointer: SatPointer,
    fallback: Option<Curve3>,
) -> Option<CurveKey> {
    let record = resolve(document, pointer)?;
    let index = index_of(record)?;
    if let Some(key) = seen.curves.get(&index) {
        return Some(*key);
    }
    let curve = read_curve(document, record).or(fallback).or_else(|| {
        // Not a kind the kernel has; the edge still exists and its record is
        // carried through, but nothing here can evaluate it.
        loss.curves.push(index as usize);
        None
    })?;
    let key = body.curves.insert(curve);
    seen.curves.insert(index, key);
    Some(key)
}

fn surface_curve(surface: &Surface, pcurve: &Curve2) -> Option<Curve3> {
    let Surface::Nurbs(surface) = surface else {
        return None;
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    let domain = [[u0, u1], [v0, v1]];
    let side = pcurve.rectangle_side(domain)?;
    let fixed = side / 2;
    Some(Curve3::Nurbs(
        surface.isocurve(fixed, domain[fixed][side % 2])?,
    ))
}

fn partner_surface_curve(
    document: &SatDocument,
    source: &SatCoedge<'_>,
) -> Option<Curve3> {
    let partner = SatCoedge::from_record(resolve(document, source.partner())?)?;
    let owner_loop = SatLoop::from_record(resolve(document, partner.owner_loop())?)?;
    let face = SatFace::from_record(resolve(document, owner_loop.face())?)?;
    let surface_record = resolve(document, face.surface())?;
    let reversed_v = analytic_surface_reversed(surface_record);
    let surface = read_surface(document, surface_record)?;
    let pcurve = read_pcurve(document, partner.pcurve(), reversed_v)?;
    surface_curve(&surface, &pcurve)
}

fn read_curve(document: &SatDocument, record: &SatRecord) -> Option<Curve3> {
    if let Some(line) = SatStraightCurve::from_record(record) {
        let (x, y, z) = line.root_point();
        let (dx, dy, dz) = line.direction();
        return Some(Curve3::Line(Line3 {
            origin: [x, y, z],
            direction: [dx, dy, dz],
        }));
    }
    if let Some(ellipse) = SatEllipseCurve::from_record(record) {
        let (cx, cy, cz) = ellipse.center();
        let (nx, ny, nz) = ellipse.normal();
        let (mx, my, mz) = ellipse.major_axis();
        let radius = Vec3::new(mx, my, mz).length();
        let plane = Plane::orthonormal([cx, cy, cz], [mx, my, mz], [nx, ny, nz])?;
        // A ratio of one is a circle, which is the overwhelming majority; the
        // rest is an ellipse and the kernel keeps it as one.
        return Some(if (ellipse.ratio() - 1.0).abs() < 1e-12 {
            Curve3::Circle(Circle3 { plane, radius })
        } else {
            Curve3::Ellipse(crate::brep::Ellipse3 {
                plane,
                major_radius: radius,
                minor_radius: radius * ellipse.ratio(),
            })
        });
    }
    if let Some(spline) = SatIntCurve::from_record(record) {
        let closed = spline.is_closed_in(document);
        let (degree, knots, controls) = spline.bspline_in(document)?;
        let mut points = Vec::with_capacity(controls.len());
        let mut weights = Vec::with_capacity(controls.len());
        for control in controls {
            let weight = control[3];
            if !weight.is_finite() || weight <= 0.0 {
                return None;
            }
            points.push([
                control[0] / weight,
                control[1] / weight,
                control[2] / weight,
            ]);
            weights.push(weight);
        }
        return Some(Curve3::Nurbs(NurbsCurve3::new_strict(
            degree, points, knots, weights,
        )?
        .with_periodicity(closed)));
    }
    None
}

fn read_pcurve(
    document: &SatDocument,
    pointer: SatPointer,
    reversed_v: bool,
) -> Option<Curve2> {
    let source = SatPCurve::from_record(resolve(document, pointer)?)?;
    let (degree, knots, controls) = source.bspline_in(document)?;
    let mut points = Vec::with_capacity(controls.len());
    let mut weights = Vec::with_capacity(controls.len());
    for control in controls {
        let weight = control[2];
        if !weight.is_finite() || weight <= 0.0 {
            return None;
        }
        points.push([control[0] / weight, control[1] / weight]);
        weights.push(weight);
    }
    let curve = Curve2::Nurbs(NurbsCurve::new_strict(degree, points, knots, weights)?);
    let curve = if reversed_v {
        curve.transformed(&crate::geom2d::Transform::scale(1.0, -1.0))?
    } else {
        curve
    };
    Some(curve)
}

fn analytic_surface_reversed(record: &SatRecord) -> bool {
    SatPlaneSurface::from_record(record)
        .map(|surface| surface.sense())
        .or_else(|| SatConeSurface::from_record(record).map(|surface| surface.sense()))
        .or_else(|| SatSphereSurface::from_record(record).map(|surface| surface.sense()))
        .or_else(|| SatTorusSurface::from_record(record).map(|surface| surface.sense()))
        .is_some_and(|sense| sense == Sense::Reversed)
}

fn surface_of(
    document: &SatDocument,
    body: &mut Body,
    seen: &mut Seen,
    loss: &mut Loss,
    pointer: SatPointer,
) -> Option<SurfaceKey> {
    let record = resolve(document, pointer)?;
    let index = index_of(record)?;
    if let Some(key) = seen.surfaces.get(&index) {
        return Some(*key);
    }
    let surface = read_surface(document, record).or_else(|| {
        loss.surfaces.push(index as usize);
        None
    })?;
    let key = body.surfaces.insert(surface);
    seen.surfaces.insert(index, key);
    Some(key)
}

fn read_surface(document: &SatDocument, record: &SatRecord) -> Option<Surface> {
    if let Some(plane) = SatPlaneSurface::from_record(record) {
        let (x, y, z) = plane.root_point();
        let (nx, ny, nz) = plane.normal();
        let (ux, uy, uz) = plane.u_direction();
        // The u direction is stored, so the frame comes from the file rather
        // than being invented — which is why the kernel's Plane takes axes.
        return Some(Surface::Plane(Plane::orthonormal(
            [x, y, z],
            [ux, uy, uz],
            [nx, ny, nz],
        )?));
    }
    if let Some(cone) = SatConeSurface::from_record(record) {
        let (cx, cy, cz) = cone.center();
        let (ax, ay, az) = cone.axis();
        let (mx, my, mz) = cone.major_axis();
        // The radius is the length of the major axis, not the `radius`
        // token — reading the token instead turns a disc into a ring.
        let radius = Vec3::new(mx, my, mz).length();
        let base = Plane::orthonormal([cx, cy, cz], [mx, my, mz], [ax, ay, az])?;
        // Only the ratio shapes the cone; a negative cosine flips the normal,
        // which `cone_points_inward` carries on the face.
        let (sine, cosine) = (cone.sin_half_angle(), cone.cos_half_angle());
        let (sine, cosine) = if cosine < 0.0 { (-sine, -cosine) } else { (sine, cosine) };
        let ratio = cone.ratio();
        let circular = (ratio - 1.0).abs() < 1e-12;
        return Some(if circular && sine.abs() < 1e-12 {
            Surface::Cylinder(Cylinder { base, radius })
        } else if circular {
            Surface::Cone(Cone {
                base,
                radius,
                half_angle: -sine.atan2(cosine),
            })
        } else {
            // An elliptical section, which is how native ACIS authors every
            // elliptical cylinder and cone: the one record the circular
            // kinds are the ratio-one special case of. Lifting it as a
            // circle would put the whole section several units off.
            Surface::EllipticCone(EllipticCone {
                base,
                radius,
                ratio,
                half_angle: if sine.abs() < 1e-12 {
                    0.0
                } else {
                    -sine.atan2(cosine)
                },
            })
        });
    }
    if let Some(sphere) = SatSphereSurface::from_record(record) {
        let (cx, cy, cz) = sphere.center();
        let (ux, uy, uz) = sphere.u_direction();
        let (px, py, pz) = sphere.pole();
        return Some(Surface::Sphere(Sphere {
            frame: Plane::orthonormal([cx, cy, cz], [ux, uy, uz], [px, py, pz])?,
            radius: sphere.radius(),
        }));
    }
    if let Some(torus) = SatTorusSurface::from_record(record) {
        let (cx, cy, cz) = torus.center();
        let (nx, ny, nz) = torus.normal();
        let (ux, uy, uz) = torus.u_direction();
        return Some(Surface::Torus(Torus {
            frame: Plane::orthonormal([cx, cy, cz], [ux, uy, uz], [nx, ny, nz])?,
            major_radius: torus.major_radius(),
            minor_radius: torus.minor_radius(),
        }));
    }
    if let Some(spline) = SatSplineSurface::from_record(record) {
        let reversed = spline.sense() == Sense::Reversed;
        let spline = spline.bspline(document)?;
        let closed = |value: &Option<String>| {
            value
                .as_deref()
                .is_some_and(|value| matches!(value, "closed" | "periodic"))
        };
        let (u_closed, v_closed) = (closed(&spline.u_closure), closed(&spline.v_closure));
        let mut points = vec![vec![[0.0; 3]; spline.control_count_v]; spline.control_count_u];
        let mut weights = vec![vec![1.0; spline.control_count_v]; spline.control_count_u];
        for v in 0..spline.control_count_v {
            for u in 0..spline.control_count_u {
                let control = spline.control_points[v * spline.control_count_u + u];
                let weight = control[3];
                if !weight.is_finite() || weight <= 0.0 {
                    return None;
                }
                points[u][v] = [
                    control[0] / weight,
                    control[1] / weight,
                    control[2] / weight,
                ];
                weights[u][v] = weight;
            }
        }
        return Some(Surface::Nurbs(NurbsSurface3::new_strict(
            spline.degree_u,
            spline.degree_v,
            points,
            spline.u_knots,
            spline.v_knots,
            weights,
        )?
        .with_periodicity(u_closed, v_closed)
        .with_v_reversed(reversed)));
    }
    None
}

fn resolve(document: &SatDocument, pointer: SatPointer) -> Option<&SatRecord> {
    (!pointer.is_null()).then(|| document.resolve(pointer)).flatten()
}

fn index_of(record: &SatRecord) -> Option<u32> {
    u32::try_from(record.index).ok()
}

fn clean(record: &SatRecord) -> Provenance {
    match index_of(record) {
        Some(index) => Provenance::Clean(SourceRef::new(index)),
        // A record with no usable index cannot be written back as itself, so
        // it is treated as something this kernel made up — which forces a
        // rebuild rather than a copy of a record it cannot find.
        None => Provenance::Synthesized,
    }
}

fn note_broken(loss: &mut Loss, record: &SatRecord) {
    if let Some(index) = index_of(record) {
        loss.broken.push(index as usize);
    }
}
