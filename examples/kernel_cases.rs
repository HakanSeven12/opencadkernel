//! Runner for Solidcore's kernel-case corpus (`crates/solidcore/tests/kernel_cases`
//! of github.com/NVIDIA/solidcore, Apache-2.0) through the kernel's own bodies
//! and operations.
//!
//! ```sh
//! scripts/fetch-kernel-cases.sh ~/solidcore            # unpacks into target/kernel_cases
//! cargo run --release --features brep --example kernel_cases -- target/kernel_cases
//! ```
//!
//! Each entry's operand bodies are read from their portable B-rep records
//! (`solidcore-record`'s `brep` form) into kernel [`Body`]s. Every operand is
//! checked with [`check`] and measured with [`mass_properties`]. The entry's
//! calls are then replayed where the kernel has the operation and its
//! parameters map cleanly (copies, transforms, booleans, equal-size rounds and
//! chamfers, shells, thickening), and the final bodies' volume is graded
//! against the entry's expected volume (its analytic one, else the reference
//! run's measured one) within the entry's tolerance, or 1e-3 relative.
//!
//! Each case runs in a child process so a hang or a crash costs one case:
//! `KERNEL_CASES_TIMEOUT_S` (60 by default) bounds it, `KERNEL_CASES_CASE`
//! keeps only cases whose path holds the substring, `KERNEL_CASES_AREA` only
//! the comma-separated areas.

#[path = "../crates/opencadkernel-ssi/src/bin/ssi_corpus_runner/json.rs"]
#[allow(dead_code)]
mod json;

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use json::Json;
use opencadkernel::brep::{
    chamfer_edges, check, combine, face_at_point, fillet_edges, mass_properties, shell, thicken,
    transform, Body, Circle3, Coedge, Cone, Curve3, Cylinder, Edge, EdgeKey, Ellipse3, Face, Line3,
    Loop, Lump, Operation, Placement, Provenance, Shell, Sphere, Surface, Torus, Vertex,
};
use opencadkernel::space::{NurbsCurve3, NurbsSurface3, Plane, Vec3};

// --- reading a record --------------------------------------------------------

fn nums(j: Option<&Json>) -> Option<Vec<f64>> {
    j?.nums()
}
fn v3(j: Option<&Json>) -> Option<[f64; 3]> {
    nums(j)?.try_into().ok()
}
fn num(j: Option<&Json>) -> Option<f64> {
    j?.num()
}
fn index(j: Option<&Json>) -> Option<usize> {
    num(j).map(|value| value as usize)
}
fn boolean(j: Option<&Json>) -> bool {
    matches!(j, Some(Json::Bool(true)))
}
fn list<'a>(j: Option<&'a Json>) -> &'a [Json] {
    j.and_then(Json::arr).unwrap_or(&[])
}
fn kind(j: &Json) -> &str {
    j.get("kind").and_then(Json::str).unwrap_or("")
}

fn frame(j: Option<&Json>) -> Option<Plane> {
    let j = j?;
    Some(Plane {
        origin: v3(j.get("origin"))?,
        x_axis: v3(j.get("x_axis"))?,
        y_axis: v3(j.get("y_axis"))?,
    })
}

/// Control points (homogeneous or not) as cartesian points and weights.
fn control(points: &[Json], rational: bool) -> Option<(Vec<[f64; 3]>, Option<Vec<f64>>)> {
    if !rational {
        return Some((points.iter().map(|p| v3(Some(p))).collect::<Option<_>>()?, None));
    }
    let mut out = Vec::with_capacity(points.len());
    let mut weights = Vec::with_capacity(points.len());
    for p in points {
        let [x, y, z, w]: [f64; 4] = p.nums()?.try_into().ok()?;
        out.push([x / w, y / w, z / w]);
        weights.push(w);
    }
    Some((out, Some(weights)))
}

fn spline_curve(j: &Json) -> Option<NurbsCurve3> {
    let rational = kind(j) == "rational_spline";
    let knots = nums(j.get("knots"))?;
    let points = list(j.get(if rational { "homogeneous" } else { "points" }));
    let degree = knots.len().checked_sub(points.len() + 1)?;
    let (points, weights) = control(points, rational)?;
    NurbsCurve3::new(degree, points, knots, weights)
}

fn spline_surface(j: &Json) -> Option<NurbsSurface3> {
    let rational = kind(j) == "rational_spline";
    let (u_knots, v_knots) = (nums(j.get("u_knots"))?, nums(j.get("v_knots"))?);
    let rows = list(j.get(if rational { "homogeneous" } else { "points" }));
    let mut net = Vec::with_capacity(rows.len());
    let mut weights = Vec::with_capacity(rows.len());
    for row in rows {
        let (points, w) = control(row.arr()?, rational)?;
        net.push(points);
        weights.push(w.unwrap_or_default());
    }
    let columns = net.first()?.len();
    let u_degree = u_knots.len().checked_sub(net.len() + 1)?;
    let v_degree = v_knots.len().checked_sub(columns + 1)?;
    NurbsSurface3::new(u_degree, v_degree, net, u_knots, v_knots, rational.then_some(weights))
}

/// A record curve in the kernel's terms. Conics are carried as themselves:
/// their NURBS form only sets which way the edge runs, which [`edge`] reads.
fn curve(j: &Json) -> Option<Curve3> {
    match kind(j) {
        "line" => {
            let (start, end) = (v3(j.get("start"))?, v3(j.get("end"))?);
            Some(Curve3::Line(Line3 {
                origin: start,
                direction: (Vec3::from(end) - Vec3::from(start)).to_array(),
            }))
        }
        "circle" => Some(Curve3::Circle(Circle3 {
            plane: frame(j.get("frame"))?,
            radius: num(j.get("radius"))?,
        })),
        "ellipse" => {
            let [major_radius, minor_radius]: [f64; 2] = nums(j.get("radii"))?.try_into().ok()?;
            Some(Curve3::Ellipse(Ellipse3 { plane: frame(j.get("frame"))?, major_radius, minor_radius }))
        }
        "polyline" => {
            let points: Vec<[f64; 3]> =
                list(j.get("points")).iter().map(|p| v3(Some(p))).collect::<Option<_>>()?;
            let last = points.len().checked_sub(1)? as f64;
            let mut knots = vec![0.0];
            knots.extend((0..points.len()).map(|i| i as f64));
            knots.push(last);
            NurbsCurve3::new(1, points, knots, None).map(Curve3::Nurbs)
        }
        "spline" | "rational_spline" => spline_curve(j).map(Curve3::Nurbs),
        // The traced curve itself is not the kernel's: its leader, which the
        // record says lies within `tolerance` of it, stands in.
        "surface_intersection" => curve(j.get("leader")?),
        _ => None,
    }
}

/// A record surface in the kernel's terms, on the record's own (u, v).
fn surface(j: &Json) -> Option<Surface> {
    match kind(j) {
        "plane" => Some(Surface::Plane(Plane {
            origin: v3(j.get("origin"))?,
            x_axis: v3(j.get("u_axis"))?,
            y_axis: v3(j.get("v_axis"))?,
        })),
        "plane_through" => {
            let origin = Vec3::from(v3(j.get("origin"))?);
            Some(Surface::Plane(Plane {
                origin: origin.to_array(),
                x_axis: (Vec3::from(v3(j.get("u_point"))?) - origin).to_array(),
                y_axis: (Vec3::from(v3(j.get("v_point"))?) - origin).to_array(),
            }))
        }
        "cylinder" => Some(Surface::Cylinder(Cylinder {
            base: frame(j.get("frame"))?,
            radius: num(j.get("radius"))?,
        })),
        // The record's radius grows along the axis at the half angle's
        // tangent; the kernel's shrinks.
        "cone" => Some(Surface::Cone(Cone {
            base: frame(j.get("frame"))?,
            radius: num(j.get("radius_at_origin"))?,
            half_angle: -num(j.get("half_angle"))?,
        })),
        "sphere" => Some(Surface::Sphere(Sphere {
            frame: frame(j.get("frame"))?,
            radius: num(j.get("radius"))?,
        })),
        // A tube starting off the outer equator only shifts `v`: the kernel
        // reads no pcurves from here, so the torus itself is the same.
        "torus" => Some(Surface::Torus(Torus {
            frame: frame(j.get("frame"))?,
            major_radius: num(j.get("major_radius"))?,
            minor_radius: num(j.get("minor_radius"))?,
        })),
        "spline" | "rational_spline" => spline_surface(j).map(Surface::Nurbs),
        _ => None,
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    Vec3::from(a).dot(Vec3::from(b))
}

/// The face's surface, and whether its natural normal agrees with the one
/// the record orients the face by: the NURBS form's, where it has one.
fn face_surface(face: &Json) -> Option<(Surface, bool)> {
    let named = face.get("surface")?;
    let Some(nurbs) = face.get("nurbs").and_then(spline_surface) else {
        return Some((surface(named)?, true));
    };
    let Some(analytic) = surface(named) else {
        // A ball patch, blend or revolution: its net is the face.
        return Some((Surface::Nurbs(nurbs), true));
    };
    let [[u0, u1], [v0, v1]]: [[f64; 2]; 2] = {
        let d = list(face.get("domain"));
        [nums(d.first())?.try_into().ok()?, nums(d.get(1))?.try_into().ok()?]
    };
    // A point well inside the domain, off any pole a corner might sit on.
    let (u, v) = (u0 + 0.37 * (u1 - u0), v0 + 0.41 * (v1 - v0));
    let reference = nurbs.normal_at_knot(u, v)?;
    let point = nurbs.point_at_knot(u, v);
    let (au, av) = analytic.parameters_at(point)?;
    let normal = analytic.normal_at(au, av)?;
    Some((analytic, dot(reference, normal) > 0.0))
}

/// The span of a conic edge on the kernel's angle, turning the conic's frame
/// over when the record's parameter runs the other way round.
fn conic_span(
    shape: &mut Curve3,
    nurbs: &NurbsCurve3,
    interval: [f64; 2],
    start: [f64; 3],
    end: [f64; 3],
) -> (f64, f64) {
    let tangent = nurbs.tangent_at_knot(interval[0]);
    let low = shape.parameter_at(nurbs.point_at_knot(interval[0]));
    if dot(shape.tangent_at(low), tangent) < 0.0 {
        let plane = match shape {
            Curve3::Circle(circle) => &mut circle.plane,
            Curve3::Ellipse(ellipse) => &mut ellipse.plane,
            _ => unreachable!(),
        };
        plane.y_axis = (-Vec3::from(plane.y_axis)).to_array();
    }
    let tau = std::f64::consts::TAU;
    let low = shape.parameter_at(start);
    let closed = Vec3::from(start).distance(Vec3::from(end)) <= 1e-9;
    if closed {
        return (low, low + tau);
    }
    let mut high = shape.parameter_at(end);
    while high <= low {
        high += tau;
    }
    while high - low > tau {
        high -= tau;
    }
    (low, high)
}

/// Reads a record into a body: every face shell of a material region becomes
/// part of that region's lump, every open shell a lump of its own. Pcurves
/// are left out; the kernel works from the edges' own curves.
fn lift(record: &Json) -> Result<Body, String> {
    let mut body = Body::new();
    body.provenance = Provenance::Synthesized;
    let vertices: Vec<_> = list(record.get("vertices"))
        .iter()
        .map(|v| {
            let point = v3(v.get("point")).ok_or("a vertex has no point")?;
            Ok(body.vertices.insert(Vertex { point, provenance: Provenance::Synthesized }))
        })
        .collect::<Result<_, String>>()?;
    let mut edges = Vec::new();
    for (at, e) in list(record.get("edges")).iter().enumerate() {
        let named = e.get("curve").ok_or("an edge has no curve")?;
        let mut shape = curve(named).ok_or_else(|| format!("edge {at}: no {} curve", kind(named)))?;
        let [start, end]: [usize; 2] = nums(e.get("vertices"))
            .and_then(|v| <[f64; 2]>::try_from(v).ok())
            .map(|v| [v[0] as usize, v[1] as usize])
            .ok_or("an edge has no vertices")?;
        let (start, end) = (
            *vertices.get(start).ok_or("an edge's vertex is missing")?,
            *vertices.get(end).ok_or("an edge's vertex is missing")?,
        );
        let interval: [f64; 2] =
            nums(e.get("interval")).and_then(|v| v.try_into().ok()).ok_or("an edge has no interval")?;
        let ends = (body.vertices.get(start).unwrap().point, body.vertices.get(end).unwrap().point);
        let (low, high) = match (&shape, e.get("nurbs").and_then(spline_curve)) {
            (Curve3::Circle(_) | Curve3::Ellipse(_), Some(nurbs)) => {
                conic_span(&mut shape, &nurbs, interval, ends.0, ends.1)
            }
            (Curve3::Circle(_) | Curve3::Ellipse(_), None) => (interval[0], interval[1]),
            (Curve3::Line(_), _) => (shape.parameter_at(ends.0), shape.parameter_at(ends.1)),
            _ if kind(named) == "surface_intersection" => {
                (shape.parameter_at(ends.0), shape.parameter_at(ends.1))
            }
            _ => (interval[0], interval[1]),
        };
        let curve = body.curves.insert(shape);
        edges.push(body.edges.insert(Edge {
            curve,
            start_parameter: low,
            end_parameter: high,
            start,
            end,
            coedges: Vec::new(),
            provenance: Provenance::Synthesized,
        }));
    }
    let faces = list(record.get("faces"));
    let shells = list(record.get("shells"));
    // Which lump each face shell goes to: the material regions first, then
    // the open shells the infinite region holds.
    let mut owner: Vec<Option<usize>> = vec![None; shells.len()];
    let mut lumps = Vec::new();
    let regions = list(record.get("regions"));
    for region in regions.iter().skip(1).filter(|r| !boolean(r.get("void"))) {
        let at = lumps.len();
        lumps.push(body.lumps.insert(Lump { shells: Vec::new(), provenance: Provenance::Synthesized }));
        for shell in list(region.get("shells")) {
            if let Some(slot) = shell.num().and_then(|s| owner.get_mut(s as usize)) {
                slot.get_or_insert(at);
            }
        }
    }
    for (at, s) in shells.iter().enumerate() {
        if kind(s) != "faces" || owner[at].is_some() {
            continue;
        }
        owner[at] = Some(lumps.len());
        lumps.push(body.lumps.insert(Lump { shells: Vec::new(), provenance: Provenance::Synthesized }));
    }
    body.roots = lumps.clone();
    for (at, s) in shells.iter().enumerate() {
        if kind(s) != "faces" {
            continue;
        }
        let Some(lump) = owner[at].map(|o| lumps[o]) else { continue };
        let shell = body.shells.insert(Shell {
            faces: Vec::new(),
            owner: lump,
            provenance: Provenance::Synthesized,
        });
        body.lumps.get_mut(lump).unwrap().shells.push(shell);
        for used in list(s.get("faces")) {
            let at = index(used.get("face")).ok_or("a shell's face has no index")?;
            let f = faces.get(at).ok_or("a shell's face is missing")?;
            let (geometry, agrees) =
                face_surface(f).ok_or_else(|| format!("face {at}: no {} surface", kind(f.get("surface").unwrap_or(&Json::Null))))?;
            // Loops run counter-clockwise about the surface's natural normal
            // in the record, about the face's outward normal in the kernel.
            let outward = boolean(used.get("forward"));
            let forward = outward == agrees;
            let surface = body.surfaces.insert(geometry);
            let face = body.faces.insert(Face {
                surface,
                forward,
                loops: Vec::new(),
                owner: shell,
                provenance: Provenance::Synthesized,
            });
            body.shells.get_mut(shell).unwrap().faces.push(face);
            for l in list(f.get("loops")) {
                if kind(l) != "edges" {
                    continue;
                }
                let ring = body.loops.insert(Loop {
                    coedges: Vec::new(),
                    owner: face,
                    provenance: Provenance::Synthesized,
                });
                body.faces.get_mut(face).unwrap().loops.push(ring);
                let mut uses: Vec<(EdgeKey, bool)> = list(l.get("edges"))
                    .iter()
                    .map(|u| {
                        let edge = *index(u.get("edge")).and_then(|e| edges.get(e)).ok_or("a loop's edge is missing")?;
                        Ok((edge, boolean(u.get("forward"))))
                    })
                    .collect::<Result<_, String>>()?;
                if !outward {
                    uses.reverse();
                    uses.iter_mut().for_each(|(_, f)| *f = !*f);
                }
                for (edge, forward) in uses {
                    let coedge = body.coedges.insert(Coedge {
                        edge,
                        forward,
                        pcurve: None,
                        owner: ring,
                        provenance: Provenance::Synthesized,
                    });
                    body.edges.get_mut(edge).unwrap().coedges.push(coedge);
                    body.loops.get_mut(ring).unwrap().coedges.push(coedge);
                }
            }
        }
    }
    if body.faces.is_empty() {
        return Err("no faces".into());
    }
    Ok(body)
}

// --- replaying a case ----------------------------------------------------------

/// The kernel edge a pick names: the one whose ends are the pick's and whose
/// samples pass nearest the picked point.
fn picked_edge(body: &Body, ends: &[Json], point: [f64; 3]) -> Option<EdgeKey> {
    let ends: Vec<Vec3> = ends.iter().filter_map(|p| v3(Some(p))).map(Vec3::from).collect();
    let size = opencadkernel::brep::body_bounds(body)
        .map(|b| Vec3::from(b.min).distance(Vec3::from(b.max)))
        .unwrap_or(1.0);
    let near = |a: Vec3, b: [f64; 3]| a.distance(Vec3::from(b)) <= 1e-6 * size.max(1.0);
    body.edges
        .iter()
        .filter(|(_, edge)| {
            let (a, b) = (body.vertices.get(edge.start).unwrap().point, body.vertices.get(edge.end).unwrap().point);
            ends.len() != 2
                || (near(ends[0], a) && near(ends[1], b))
                || (near(ends[0], b) && near(ends[1], a))
        })
        .filter_map(|(key, _)| {
            let samples = opencadkernel::brep::edge_points(body, key, size * 1e-4)?;
            let gap = samples
                .iter()
                .map(|s| Vec3::from(*s).distance(Vec3::from(point)))
                .fold(f64::INFINITY, f64::min);
            Some((key, gap))
        })
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(key, _)| key)
}

fn picked_edges(body: &Body, pick: Option<&Json>) -> Option<Vec<EdgeKey>> {
    let pick = pick?;
    let points = list(pick.get("points"));
    list(pick.get("picked"))
        .iter()
        .zip(points)
        .map(|(p, at)| picked_edge(body, list(p.get("at")), v3(Some(at))?))
        .collect()
}

enum Step {
    Made(Vec<Body>),
    Refused(String),
    Unported(String),
}

fn run_call(call: &Json, inputs: Vec<Body>) -> Step {
    let op = call.get("op").and_then(Json::str).unwrap_or("");
    let p = call.get("params").unwrap_or(&Json::Null);
    let tolerance = |bodies: &[Body]| {
        opencadkernel::brep::operation_tolerance(&bodies.iter().collect::<Vec<_>>())
    };
    let refused = |what: &str| Step::Refused(what.to_string());
    match op {
        "copy" | "tagged" => Step::Made(inputs),
        "transform" => {
            let m: Vec<Vec<f64>> = list(p.get("matrix")).iter().filter_map(Json::nums).collect();
            if m.len() != 4 || m.iter().any(|r| r.len() != 4) {
                return Step::Unported("transform matrix".into());
            }
            let place = Placement {
                x_axis: [m[0][0], m[0][1], m[0][2]],
                y_axis: [m[1][0], m[1][1], m[1][2]],
                z_axis: [m[2][0], m[2][1], m[2][2]],
                origin: [m[3][0], m[3][1], m[3][2]],
            };
            match inputs.iter().map(|b| transform(b, &place)).collect::<Option<Vec<_>>>() {
                Some(bodies) => Step::Made(bodies),
                None => refused("transform"),
            }
        }
        "boolean" | "union_all" => {
            let how = match p.get("op").and_then(Json::str).unwrap_or("union") {
                "union" => Operation::Union,
                "subtract" => Operation::Difference,
                "intersect" => Operation::Intersection,
                other => return Step::Unported(format!("boolean {other}")),
            };
            let tolerance = tolerance(&inputs);
            let mut bodies = inputs.into_iter();
            let Some(mut acc) = bodies.next() else { return refused("no operand") };
            for other in bodies {
                acc = match combine(acc, other, how, tolerance) {
                    Ok(body) => body,
                    Err(snag) => return Step::Refused(format!("{snag:?}")),
                };
            }
            Step::Made(vec![acc])
        }
        "blend_dropping" => {
            let blend = kind(p);
            let sizes = nums(p.get("sizes")).unwrap_or_default();
            if p.get("turn").is_some() || sizes.is_empty() || sizes.iter().any(|s| *s != sizes[0]) {
                return Step::Unported(format!("blend {blend} with a turn or several sizes"));
            }
            let [body] = <[Body; 1]>::try_from(inputs).ok().unwrap();
            let Some(edges) = picked_edges(&body, p.get("edges")) else {
                return refused("edge pick not found");
            };
            match blend {
                "round" => match fillet_edges(&body, &edges, sizes[0]) {
                    Ok(body) => Step::Made(vec![body]),
                    Err(e) => Step::Refused(format!("{e:?}")),
                },
                "chamfer" => {
                    let Some(face) = body
                        .edges
                        .get(edges[0])
                        .and_then(|edge| body.coedges.get(*edge.coedges.first()?))
                        .and_then(|coedge| body.loops.get(coedge.owner))
                        .map(|ring| ring.owner)
                    else {
                        return refused("edge without a face");
                    };
                    match chamfer_edges(&body, &edges, face, sizes[0], sizes[0]) {
                        Ok(body) => Step::Made(vec![body]),
                        Err(e) => Step::Refused(format!("{e:?}")),
                    }
                }
                other => Step::Unported(format!("blend {other}")),
            }
        }
        "shell" => {
            let [body] = <[Body; 1]>::try_from(inputs).ok().unwrap();
            let Some(thickness) = num(p.get("thickness")) else { return Step::Unported("shell".into()) };
            let faces: Option<Vec<_>> = list(p.get("open_faces").and_then(|f| f.get("points")))
                .iter()
                .map(|at| face_at_point(&body, v3(Some(at))?))
                .collect();
            let Some(faces) = faces else { return refused("face pick not found") };
            // The record's thickness is negative inward; the kernel's
            // distance is positive inward.
            match shell(&body, &faces, -thickness) {
                Ok(body) => Step::Made(vec![body]),
                Err(e) => Step::Refused(format!("{e:?}")),
            }
        }
        "thicken" if !boolean(p.get("symmetric")) => {
            let [body] = <[Body; 1]>::try_from(inputs).ok().unwrap();
            match num(p.get("thickness")).map(|t| thicken(&body, t)) {
                Some(Ok(body)) => Step::Made(vec![body]),
                Some(Err(e)) => Step::Refused(format!("{e:?}")),
                None => Step::Unported("thicken".into()),
            }
        }
        other => Step::Unported(other.to_string()),
    }
}

/// The expected volume as `[low, high]`.
fn expected_volume(case: &Json) -> Option<[f64; 2]> {
    let expected = case.get("expected");
    let from_rig = || {
        let rig = case.get("rig")?;
        (rig.get("status").and_then(Json::str) == Some("pass"))
            .then(|| num(rig.get("measured")?.get("volume_mm3")))
            .flatten()
    };
    let (value, tolerance) = match expected.and_then(|e| e.get("volume_mm3")) {
        Some(Json::Arr(range)) => {
            let [low, high]: [f64; 2] = range.iter().filter_map(Json::num).collect::<Vec<_>>().try_into().ok()?;
            return Some([low, high]);
        }
        Some(value) => (value.num()?, expected.and_then(|e| e.get("tolerance"))),
        None => (from_rig()?, None),
    };
    let slack = match tolerance {
        Some(t) if t.get("abs_mm").is_some() => num(t.get("abs_mm"))?,
        Some(t) => num(t.get("rel")).unwrap_or(1e-3) * value.abs(),
        None => 1e-3 * value.abs(),
    };
    Some([value - slack, value + slack])
}

/// Runs one case and prints `status \t detail \t operands lifted/checked/measured`.
fn run_case(path: &Path) -> Result<String, String> {
    let read = |p: &Path| std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()));
    let case = json::parse(&read(path)?)?;
    let records = match case.get("records").and_then(Json::str) {
        Some(name) => json::parse(&read(&path.with_file_name(name.trim_end_matches(".gz")))?)?,
        None => Json::Obj(Vec::new()),
    };
    let mut bodies: HashMap<String, Body> = HashMap::new();
    let (mut operands, mut checked, mut measured) = (0, 0, 0);
    let mut lift_error = None;
    if let Json::Obj(fields) = &records {
        for (name, record) in fields {
            operands += 1;
            match lift(record) {
                Ok(body) => {
                    checked += usize::from(check(&body).is_empty());
                    measured += usize::from(mass_properties(&body).is_some_and(|m| m.volume > 0.0));
                    bodies.insert(name.clone(), body);
                }
                Err(e) => {
                    lift_error.get_or_insert(format!("{name}: {e}"));
                }
            }
        }
    }
    let lifted = bodies.len();
    let stats = format!("{lifted}/{operands}\t{checked}\t{measured}");
    if let Some(e) = lift_error {
        return Ok(format!("lift-fail\t{e}\t{stats}"));
    }
    let calls = list(case.get("calls"));
    if calls.is_empty() || matches!(case.get("operation_expressible"), Some(Json::Bool(false))) {
        return Ok(format!("unsupported\tnot the kernel operations' alone\t{stats}"));
    }
    for call in calls {
        let inputs: Option<Vec<Body>> =
            list(call.get("inputs")).iter().map(|n| bodies.get(n.str()?).cloned()).collect();
        let Some(inputs) = inputs else {
            return Ok(format!("unported\tan operand of {} was not made\t{stats}", call.get("op").and_then(Json::str).unwrap_or("?")));
        };
        let reference_refused = call.get("error").is_some_and(|e| !matches!(e, Json::Null));
        let continued = boolean(call.get("continued"));
        match run_call(call, inputs) {
            Step::Unported(what) => return Ok(format!("unported\t{what}\t{stats}")),
            Step::Refused(why) if reference_refused => {
                if !continued {
                    return Ok(format!("pass\trefused as the reference did: {why}\t{stats}"));
                }
            }
            Step::Refused(why) => {
                return Ok(format!(
                    "fail\t{} refused: {why}\t{stats}",
                    call.get("op").and_then(Json::str).unwrap_or("?")
                ))
            }
            Step::Made(_) if reference_refused && !continued => {
                return Ok(format!("fail\tthe reference refused, the kernel did not\t{stats}"));
            }
            Step::Made(made) => {
                for (name, body) in list(call.get("outputs")).iter().zip(made) {
                    if let Some(name) = name.str() {
                        bodies.insert(name.to_string(), body);
                    }
                }
            }
        }
    }
    let Some([low, high]) = expected_volume(&case) else {
        return Ok(format!("ran\tno expected volume\t{stats}"));
    };
    let mut volume = 0.0;
    for result in list(case.get("result")) {
        let name = result.get("made").or_else(|| result.get("operand")).and_then(Json::str).unwrap_or("");
        let Some(body) = bodies.get(name) else {
            return Ok(format!("fail\tresult {name} missing\t{stats}"));
        };
        let Some(properties) = mass_properties(body) else {
            return Ok(format!("fail\tresult {name} not measurable\t{stats}"));
        };
        volume += properties.volume;
    }
    let status = if (low..=high).contains(&volume) { "pass" } else { "fail" };
    Ok(format!("{status}\tvolume {volume:.6} expected [{low:.6}, {high:.6}]\t{stats}"))
}

fn cases(root: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let mut stack = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else { continue };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else if dir != root
                && path.extension().is_some_and(|e| e == "json")
                && !path.to_string_lossy().ends_with(".bodies.json")
            {
                out.push(path);
            }
        }
    }
    out.sort();
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("--one") {
        let path = PathBuf::from(&args[2]);
        match run_case(&path) {
            Ok(line) => println!("{line}"),
            Err(e) => println!("error\t{e}\t0/0\t0\t0"),
        }
        return;
    }
    let Some(root) = args.get(1).map(PathBuf::from) else {
        eprintln!("usage: kernel_cases <kernel_cases dir>");
        std::process::exit(2);
    };
    let timeout = Duration::from_secs(
        std::env::var("KERNEL_CASES_TIMEOUT_S").ok().and_then(|s| s.parse().ok()).unwrap_or(60),
    );
    let only = std::env::var("KERNEL_CASES_CASE").ok();
    let areas: Option<Vec<String>> = std::env::var("KERNEL_CASES_AREA")
        .ok()
        .map(|a| a.split(',').map(str::to_string).collect());
    let exe = std::env::current_exe().expect("own path");
    let mut totals: std::collections::BTreeMap<String, [usize; 8]> = Default::default();
    let mut lift = [0usize; 4];
    for path in cases(&root) {
        let relative = path.strip_prefix(&root).unwrap_or(&path).to_string_lossy().to_string();
        let area = relative.split('/').next().unwrap_or("").to_string();
        if only.as_ref().is_some_and(|o| !relative.contains(o.as_str()))
            || areas.as_ref().is_some_and(|a| !a.contains(&area))
        {
            continue;
        }
        let started = Instant::now();
        let mut child = Command::new(&exe)
            .arg("--one")
            .arg(&path)
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .expect("spawn the case");
        let line = loop {
            match child.try_wait() {
                Ok(Some(_)) => {
                    let mut text = String::new();
                    use std::io::Read;
                    child.stdout.take().map(|mut s| s.read_to_string(&mut text));
                    break if text.trim().is_empty() {
                        "crash\tno output\t0/0\t0\t0".to_string()
                    } else {
                        text.trim().to_string()
                    };
                }
                Ok(None) if started.elapsed() > timeout => {
                    let _ = child.kill();
                    let _ = child.wait();
                    break "timeout\t\t0/0\t0\t0".to_string();
                }
                Ok(None) => std::thread::sleep(Duration::from_millis(20)),
                Err(e) => break format!("crash\t{e}\t0/0\t0\t0"),
            }
        };
        let fields: Vec<&str> = line.split('\t').collect();
        let status = fields.first().copied().unwrap_or("crash");
        let slot = ["pass", "fail", "unported", "unsupported", "lift-fail", "ran", "timeout"]
            .iter()
            .position(|s| *s == status)
            .unwrap_or(7);
        let total = totals.entry(area).or_default();
        total[slot] += 1;
        if let Some((done, all)) = fields.get(2).and_then(|f| f.split_once('/')) {
            lift[0] += all.parse::<usize>().unwrap_or(0);
            lift[1] += done.parse::<usize>().unwrap_or(0);
            lift[2] += fields.get(3).and_then(|f| f.parse::<usize>().ok()).unwrap_or(0);
            lift[3] += fields.get(4).and_then(|f| f.parse::<usize>().ok()).unwrap_or(0);
        }
        println!("{relative}\t{line}\t{:.1}s", started.elapsed().as_secs_f64());
    }
    println!("\narea\tpass\tfail\tunported\tunsupported\tlift-fail\tran\ttimeout\tcrash");
    let mut sum = [0usize; 8];
    for (area, t) in &totals {
        println!("{area}\t{}", t.map(|n| n.to_string()).join("\t"));
        for i in 0..8 {
            sum[i] += t[i];
        }
    }
    println!("all\t{}", sum.map(|n| n.to_string()).join("\t"));
    println!(
        "operand bodies: {} records, {} lifted, {} pass check, {} measured",
        lift[0], lift[1], lift[2], lift[3]
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A unit cube as the record writes one: six planes, twelve lines, the
    /// loops counter-clockwise about each plane's own normal and the faces
    /// turned outward by their shell use.
    #[test]
    fn a_unit_cube_record_lifts_to_a_closed_unit_volume() {
        let corner = |i: usize| [(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64];
        // Faces as corner rings, counter-clockwise seen from outside.
        let rings: [[usize; 4]; 6] = [
            [0, 2, 3, 1], // z = 0
            [4, 5, 7, 6], // z = 1
            [0, 1, 5, 4], // y = 0
            [2, 6, 7, 3], // y = 1
            [0, 4, 6, 2], // x = 0
            [1, 3, 7, 5], // x = 1
        ];
        let mut edges: Vec<[usize; 2]> = Vec::new();
        let mut text = String::from("{\"regions\":[{\"shells\":[0],\"void\":false},{\"shells\":[0],\"void\":false}],");
        let mut faces = Vec::new();
        for (f, ring) in rings.iter().enumerate() {
            let mut uses = Vec::new();
            for k in 0..4 {
                let (a, b) = (ring[k], ring[(k + 1) % 4]);
                let (low, high) = (a.min(b), a.max(b));
                let at = edges.iter().position(|e| *e == [low, high]).unwrap_or_else(|| {
                    edges.push([low, high]);
                    edges.len() - 1
                });
                uses.push(format!("{{\"edge\":{at},\"forward\":{}}}", a == low));
            }
            // Odd faces keep the ring about the plane's own normal; even ones
            // are written about the reversed normal and used backwards.
            let (o, u, v) = (corner(ring[0]), corner(ring[1]), corner(ring[3]));
            let axis = |p: [f64; 3]| [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
            let (mut u_axis, mut v_axis) = (axis(u), axis(v));
            let forward = f % 2 == 1;
            if !forward {
                std::mem::swap(&mut u_axis, &mut v_axis);
                uses.reverse();
                for use_ in &mut uses {
                    *use_ = if use_.contains("true") {
                        use_.replace("true", "false")
                    } else {
                        use_.replace("false", "true")
                    };
                }
            }
            faces.push((
                format!(
                    "{{\"surface\":{{\"kind\":\"plane\",\"origin\":{o:?},\"u_axis\":{u_axis:?},\"v_axis\":{v_axis:?}}},\"loops\":[{{\"kind\":\"edges\",\"edges\":[{}]}}]}}",
                    uses.join(",")
                ),
                forward,
            ));
        }
        let shell_faces: Vec<String> =
            (0..6).map(|f| format!("{{\"face\":{f},\"forward\":{}}}", faces[f].1)).collect();
        text += &format!("\"shells\":[{{\"kind\":\"faces\",\"faces\":[{}]}}],", shell_faces.join(","));
        text += &format!("\"faces\":[{}],", faces.iter().map(|f| f.0.clone()).collect::<Vec<_>>().join(","));
        let edge_text: Vec<String> = edges
            .iter()
            .map(|[a, b]| {
                format!(
                    "{{\"curve\":{{\"kind\":\"line\",\"start\":{:?},\"end\":{:?}}},\"interval\":[0,1],\"vertices\":[{a},{b}]}}",
                    corner(*a),
                    corner(*b)
                )
            })
            .collect();
        text += &format!("\"edges\":[{}],", edge_text.join(","));
        let vertex_text: Vec<String> = (0..8).map(|i| format!("{{\"point\":{:?}}}", corner(i))).collect();
        text += &format!("\"vertices\":[{}]}}", vertex_text.join(","));
        let body = lift(&json::parse(&text).unwrap()).unwrap();
        assert!(check(&body).is_empty(), "{:?}", check(&body));
        let volume = mass_properties(&body).unwrap().volume;
        assert!((volume - 1.0).abs() < 1e-9, "volume {volume}");
    }
}
