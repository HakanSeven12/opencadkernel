//! A ball rolled along a closed chain of edges between two faces of any
//! surfaces with an exact offset: the saddle where two cylinders join, a
//! ring where a cone meets a sphere. The ball's centre runs where the two
//! surfaces moved by the radius meet; it touches each face along a curve,
//! and the round between those two curves replaces the chain. The faces
//! keep their own surfaces and are trimmed back to the touching curves; an
//! edge meeting the chain — a seam crossing it — is followed to where it
//! meets the touching curve instead.

use super::shell::{least_norm_step, moved_surface};
use super::topology::{Coedge, CoedgeKey, Edge, EdgeKey, Face, FaceKey, Loop, Vertex, VertexKey};
use super::{Body, Curve3, Provenance, Surface};
use crate::space::spline::interpolate_periodic;
use crate::space::{NurbsCurve3, NurbsSurface3, Parameterization, Vec3};
use std::collections::{HashMap, HashSet};

/// Points the ball's path is followed through, round each closed chain.
const SAMPLES: usize = 160;

/// Rounds every closed chain the selected edges lie on. `None` where a
/// chain is open, its faces' surfaces have no exact offset, the ball cannot
/// be placed along it, or the rounded solid does not come out sound.
pub(super) fn fillet_rolling(body: &Body, selected: &[EdgeKey], radius: f64) -> Option<Body> {
    let tolerance = super::operation_tolerance(&[body]);
    let face_of = |body: &Body, coedge: CoedgeKey| -> Option<FaceKey> {
        Some(body.loops.get(body.coedges.get(coedge)?.owner)?.owner)
    };
    let pair_of = |body: &Body, edge: EdgeKey| -> Option<(FaceKey, FaceKey)> {
        let node = body.edges.get(edge)?;
        let [one, other] = node.coedges[..] else { return None };
        let (one, other) = (face_of(body, one)?, face_of(body, other)?);
        (one != other).then_some((one, other))
    };
    // Each pair of faces the selection runs between, with every edge they
    // share: one selected edge of a chain rounds the whole chain.
    let mut pairs: Vec<(FaceKey, FaceKey)> = Vec::new();
    for edge in selected {
        let (one, other) = pair_of(body, *edge)?;
        if !pairs.iter().any(|pair| *pair == (one, other) || *pair == (other, one)) {
            pairs.push((one, other));
        }
    }
    // Each chain, and which way its corner turns, read off the body as it
    // came: a round made already would be asked through its own spline.
    let mut chains = Vec::new();
    for (one, other) in pairs {
        let chain: Vec<EdgeKey> = body
            .edges
            .keys()
            .filter(|edge| {
                pair_of(body, *edge)
                    .is_some_and(|pair| pair == (one, other) || pair == (other, one))
            })
            .collect();
        // The pair may meet along more than one loop — a hole's two ends
        // between a wall and the bore — each its own chain.
        for chain in loops_of(body, &chain) {
            let concave = concave_chain(body, &chain, one, other, radius, tolerance)?;
            chains.push((chain, one, other, concave));
        }
    }
    let mut result = body.clone();
    for (chain, one, other, concave) in chains {
        round_chain(&mut result, &chain, one, other, radius, concave, tolerance)?;
    }
    if !result.validate().is_empty() {
        return None;
    }
    // A round adds what it fills into a concave chain and takes what it cuts
    // off a convex one, a sliver either way: a face turned inside out shows
    // as a large change even on a coarse mesh, which is all this asks.
    let volume = |body: &Body| -> Option<f64> {
        let mesh =
            super::mesh::tessellate(body, super::mesh::TessellationTolerance::new(0.3, tolerance));
        mesh.missing_faces.is_empty().then(|| mesh.mesh.mass_properties().map(|(v, _)| v))?
    };
    let (before, after) = (volume(body)?, volume(&result)?);
    ((after - before).abs() <= before.abs() * 0.1).then_some(result)
}

/// `edges` parted into the groups joined end to end.
fn loops_of(body: &Body, edges: &[EdgeKey]) -> Vec<Vec<EdgeKey>> {
    let mut left: Vec<EdgeKey> = edges.to_vec();
    let mut out = Vec::new();
    while let Some(seed) = left.pop() {
        let mut group = vec![seed];
        let mut grew = true;
        while grew {
            grew = false;
            let ends: HashSet<VertexKey> = group
                .iter()
                .filter_map(|edge| body.edges.get(*edge))
                .flat_map(|edge| [edge.start, edge.end])
                .collect();
            if let Some(at) = left.iter().position(|edge| {
                body.edges
                    .get(*edge)
                    .is_some_and(|node| ends.contains(&node.start) || ends.contains(&node.end))
            }) {
                group.push(left.remove(at));
                grew = true;
            }
        }
        out.push(group);
    }
    out
}

/// The chain's edges in order round it, each with whether it is walked
/// along its own direction. `None` unless they close into one loop.
fn ordered(body: &Body, chain: &[EdgeKey]) -> Option<Vec<(EdgeKey, bool)>> {
    let first = *chain.first()?;
    let mut out = vec![(first, true)];
    let mut at = body.edges.get(first)?.end;
    let start = body.edges.get(first)?.start;
    let mut used: HashSet<EdgeKey> = HashSet::from([first]);
    while at != start || out.len() < chain.len() {
        let next = chain.iter().copied().find_map(|edge| {
            let node = body.edges.get(edge)?;
            if used.contains(&edge) {
                return None;
            }
            if node.start == at {
                Some((edge, true))
            } else if node.end == at {
                Some((edge, false))
            } else {
                None
            }
        })?;
        used.insert(next.0);
        let node = body.edges.get(next.0)?;
        at = if next.1 { node.end } else { node.start };
        out.push(next);
        if out.len() > chain.len() {
            return None;
        }
    }
    (out.len() == chain.len()).then_some(out)
}

/// Concave where the solid takes in a point off one face's outside and the
/// other's inside — the union of the two sides — so the ball rolls outside
/// it and fills; convex, where the solid is only what is inside both, so it
/// rolls inside and cuts.
fn concave_chain(
    body: &Body,
    chain: &[EdgeKey],
    one: FaceKey,
    other: FaceKey,
    radius: f64,
    tolerance: f64,
) -> Option<bool> {
    let edge = body.edges.get(*chain.first()?)?;
    let curve = body.curves.get(edge.curve)?;
    let probe = Vec3::from(curve.point_at(0.5 * (edge.start_parameter + edge.end_parameter)));
    let outward = |face: FaceKey| -> Option<Vec3> {
        let node = body.faces.get(face)?;
        let surface = body.surfaces.get(node.surface)?;
        let (u, v) = surface.parameters_at(probe.to_array())?;
        let sense = if node.forward { 1.0 } else { -1.0 };
        Some(Vec3::from(surface.normal_at(u, v)?).normalize()? * sense)
    };
    let across = outward(one)? - outward(other)?;
    let at = probe + across * (radius * 0.01);
    let side = super::classify::contains_point(body, at.to_array(), tolerance);
    Some(side == super::Containment::Inside)
}

fn round_chain(
    body: &mut Body,
    chain: &[EdgeKey],
    one: FaceKey,
    other: FaceKey,
    radius: f64,
    concave: bool,
    tolerance: f64,
) -> Option<()> {
    let ring = ordered(body, chain)?;
    // Points along the chain, evenly by length — the round is fitted as
    // if they were, and a short edge given its own share bunched them —
    // and which way it runs there.
    let mut path: Vec<(Vec3, f64, EdgeKey, bool)> = Vec::new();
    for (edge, forward) in &ring {
        let node = body.edges.get(*edge)?;
        let span = node.end_parameter - node.start_parameter;
        let curve = body.curves.get(node.curve)?;
        for step in 0..64 {
            let fraction = step as f64 / 64.0;
            let fraction = if *forward { fraction } else { 1.0 - fraction };
            let t = node.start_parameter + span * fraction;
            path.push((Vec3::from(curve.point_at(t)), t, *edge, *forward));
        }
    }
    let mut lengths = vec![0.0];
    for index in 0..path.len() {
        let next = &path[(index + 1) % path.len()];
        lengths.push(lengths[index] + path[index].0.distance(next.0));
    }
    let total = *lengths.last()?;
    let mut along: Vec<(Vec3, Vec3)> = Vec::with_capacity(SAMPLES);
    let mut index = 0;
    for sample in 0..SAMPLES {
        let target = total * sample as f64 / SAMPLES as f64;
        while lengths[index + 1] < target {
            index += 1;
        }
        let (_, t, edge, forward) = path[index];
        let node = body.edges.get(edge)?;
        let curve = body.curves.get(node.curve)?;
        // Within one sampled step of its edge, the parameter in proportion.
        let next = &path[(index + 1) % path.len()];
        let piece = lengths[index + 1] - lengths[index];
        let share = if piece > 0.0 { (target - lengths[index]) / piece } else { 0.0 };
        let t = if next.2 == edge {
            t + (next.1 - t) * share
        } else {
            let end = if forward { node.end_parameter } else { node.start_parameter };
            t + (end - t) * share
        };
        let point = Vec3::from(curve.point_at(t));
        let step = (node.end_parameter - node.start_parameter).abs() * 1e-5;
        let ahead = Vec3::from(curve.point_at(t + step)) - Vec3::from(curve.point_at(t - step));
        let tangent = ahead.normalize()? * if forward { 1.0 } else { -1.0 };
        along.push((point, tangent));
    }
    let surface_of = |face: FaceKey| -> Option<(Surface, f64)> {
        let node = body.faces.get(face)?;
        Some((body.surfaces.get(node.surface)?.clone(), if node.forward { 1.0 } else { -1.0 }))
    };
    let ((first, first_sense), (second, second_sense)) = (surface_of(one)?, surface_of(other)?);
    let outward = |surface: &Surface, sense: f64, point: Vec3| -> Option<Vec3> {
        let (u, v) = surface.parameters_at(point.to_array())?;
        Some(Vec3::from(surface.normal_at(u, v)?).normalize()? * sense)
    };
    let side = if concave { 1.0 } else { -1.0 };
    let first_moved = moved_surface(&first, radius * side * first_sense)?;
    let second_moved = moved_surface(&second, radius * side * second_sense)?;
    // The ball's centre at each point, square to the chain there, and where
    // it touches each face.
    let mut centres = Vec::with_capacity(along.len());
    let mut touches = Vec::with_capacity(along.len());
    for (point, tangent) in &along {
        let lean =
            outward(&first, first_sense, *point)? + outward(&second, second_sense, *point)?;
        let mut centre = *point + lean.normalize()? * (radius * side);
        for _ in 0..40 {
            let gap = |surface: &Surface| -> Option<(Vec3, f64)> {
                let (u, v) = surface.parameters_at(centre.to_array())?;
                let normal = Vec3::from(surface.normal_at(u, v)?).normalize()?;
                Some((normal, surface.distance_to(centre.to_array())))
            };
            let square = (*tangent, (centre - *point).dot(*tangent));
            let rows = [gap(&first_moved)?, gap(&second_moved)?, square];
            if rows.iter().all(|(_, gap)| gap.abs() <= tolerance * 1e-3) {
                break;
            }
            centre = centre - least_norm_step(&rows)?;
        }
        let foot = |surface: &Surface| -> Option<Vec3> {
            let (u, v) = surface.parameters_at(centre.to_array())?;
            Some(Vec3::from(surface.point_at(u, v)))
        };
        let (near, far) = (foot(&first)?, foot(&second)?);
        if (near.distance(centre) - radius).abs() > radius * 1e-6
            || (far.distance(centre) - radius).abs() > radius * 1e-6
        {
            return None;
        }
        centres.push(centre);
        touches.push((near, far));
    }
    // The round: at each point the arc of the ball from one touch to the
    // other facing the chain, as a rational quadratic; the arcs joined
    // round the chain in homogeneous coordinates, so the weights follow.
    let mut rows: [Vec<[f64; 4]>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    for (centre, (near, far)) in centres.iter().zip(&touches) {
        let (a, b) = (*near - *centre, *far - *centre);
        let half = 0.5 * a.cross(b).length().atan2(a.dot(b));
        if half <= 1e-6 || half >= std::f64::consts::FRAC_PI_2 - 1e-6 {
            return None;
        }
        let weight = half.cos();
        let middle = *centre + (a + b).normalize()? * (radius / weight);
        let lift = |point: Vec3, weight: f64| {
            [point.x * weight, point.y * weight, point.z * weight, weight]
        };
        rows[0].push(lift(*near, 1.0));
        rows[1].push(lift(middle, weight));
        rows[2].push(lift(*far, 1.0));
    }
    let fitted: Vec<(Vec<[f64; 4]>, Vec<f64>)> = rows
        .iter()
        .map(|row| interpolate_periodic(row, Parameterization::Uniform))
        .collect::<Option<_>>()?;
    let knots = fitted[0].1.clone();
    let count = fitted[0].0.len();
    let mut net = Vec::with_capacity(count);
    let mut weights = Vec::with_capacity(count);
    for index in 0..count {
        let mut points = Vec::with_capacity(3);
        let mut row_weights = Vec::with_capacity(3);
        for (control, _) in &fitted {
            let [x, y, z, w] = control[index];
            points.push([x / w, y / w, z / w]);
            row_weights.push(w);
        }
        net.push(points);
        weights.push(row_weights);
    }
    let round = NurbsSurface3::new(
        3,
        2,
        net,
        knots,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        Some(weights),
    )?
    .with_periodicity(true, false);
    let touch_curve = |side: usize| -> Option<NurbsCurve3> {
        let points: Vec<[f64; 3]> = touches
            .iter()
            .map(|(near, far)| if side == 0 { near.to_array() } else { far.to_array() })
            .collect();
        NurbsCurve3::interpolate_periodic(&points, Parameterization::Uniform)
    };
    let (near_curve, far_curve) = (touch_curve(0)?, touch_curve(1)?);

    // Each face trimmed back to its touching curve.
    let shell = body.faces.get(one)?.owner;
    let near_edges = retrim(body, one, chain, &near_curve, tolerance)?;
    let far_edges = retrim(body, other, chain, &far_curve, tolerance)?;
    // The chain's own edges go, and its corners nothing else still uses.
    let mut corners = HashSet::new();
    for edge in chain {
        if let Some(node) = body.edges.remove(*edge) {
            corners.extend([node.start, node.end]);
            for coedge in node.coedges {
                body.coedges.remove(coedge);
            }
        }
    }
    let used: HashSet<VertexKey> =
        body.edges.iter().flat_map(|(_, edge)| [edge.start, edge.end]).collect();
    for vertex in corners {
        if !used.contains(&vertex) {
            body.vertices.remove(vertex);
        }
    }

    // The round's face, its outward side away from the material: toward
    // the ball's centre where it fills, away from it where it cuts.
    let surface = body.surfaces.insert(Surface::Nurbs(round.clone()));
    // Spaced evenly, sample `i` is at `u = i`.
    let index = centres.len() / 2;
    let normal = Vec3::from(round.normal_at_knot(index as f64, 0.5)?);
    let point = Vec3::from(round.point_at_knot(index as f64, 0.5));
    let toward_centre = (centres[index] - point).dot(normal) > 0.0;
    let forward = toward_centre == concave;
    let face = body.faces.insert(Face {
        surface,
        forward,
        loops: Vec::new(),
        owner: shell,
        provenance: Provenance::Synthesized,
    });
    for edges in [near_edges, far_edges] {
        let ring = body.loops.insert(Loop {
            coedges: Vec::new(),
            owner: face,
            provenance: Provenance::Synthesized,
        });
        // The other use of each edge, the opposite way, in order round it.
        let mut coedges = Vec::with_capacity(edges.len());
        let reversed = edges.first().is_some_and(|(_, forward)| *forward);
        let mut order = edges.clone();
        if reversed {
            order.reverse();
        }
        for (edge, forward) in order {
            let coedge = body.coedges.insert(Coedge {
                edge,
                forward: !forward,
                pcurve: None,
                owner: ring,
                provenance: Provenance::Synthesized,
            });
            body.edges.get_mut(edge)?.coedges.push(coedge);
            coedges.push(coedge);
        }
        body.loops.get_mut(ring)?.coedges = coedges;
        body.faces.get_mut(face)?.loops.push(ring);
    }
    body.shells.get_mut(shell)?.faces.push(face);
    Some(())
}

/// Rebuilds `face`'s loops with the chain's run through each replaced by
/// the touching curve between the same neighbours, each neighbour's end
/// moved to where it meets that curve. Returns the new edges in the order
/// the curve runs, each with the way the face walks it.
fn retrim(
    body: &mut Body,
    face: FaceKey,
    chain: &[EdgeKey],
    curve: &NurbsCurve3,
    tolerance: f64,
) -> Option<Vec<(EdgeKey, bool)>> {
    let (start, end) = curve.domain();
    let period = end - start;
    let curve_key = body.curves.insert(Curve3::Nurbs(curve.clone()));
    let chained: HashSet<EdgeKey> = chain.iter().copied().collect();
    // Where along the touching curve each corner the chain shares with a
    // neighbour edge goes: where the curve meets that edge's own curve,
    // the crossing nearest the corner.
    let mut moved: HashMap<VertexKey, (f64, VertexKey)> = HashMap::new();
    let mut made: Vec<(EdgeKey, bool, f64)> = Vec::new();
    let loops = body.faces.get(face)?.loops.clone();
    for ring in loops {
        let coedges = body.loops.get(ring)?.coedges.clone();
        let flags: Vec<bool> = coedges
            .iter()
            .map(|coedge| {
                body.coedges.get(*coedge).is_some_and(|node| chained.contains(&node.edge))
            })
            .collect();
        if !flags.iter().any(|flag| *flag) {
            continue;
        }
        let place: HashMap<CoedgeKey, usize> =
            coedges.iter().enumerate().map(|(index, coedge)| (*coedge, index)).collect();
        let on_chain = |coedge: &CoedgeKey| place.get(coedge).is_some_and(|index| flags[*index]);
        let walk = |body: &Body, coedge: CoedgeKey| -> Option<(VertexKey, VertexKey)> {
            body.coedge_vertices(coedge)
        };
        // Runs of chain coedges, each from the corner before to the one after.
        let count = coedges.len();
        let all = coedges.iter().all(on_chain);
        let begin = if all {
            0
        } else {
            (0..count).find(|index| {
                !on_chain(&coedges[*index]) && on_chain(&coedges[(index + 1) % count])
            })? + 1
        };
        let mut rebuilt: Vec<CoedgeKey> = Vec::new();
        let mut index = 0;
        while index < count {
            let at = (begin + index) % count;
            if !on_chain(&coedges[at]) {
                rebuilt.push(coedges[at]);
                index += 1;
                continue;
            }
            let mut run = Vec::new();
            while index < count && on_chain(&coedges[(begin + index) % count]) {
                run.push(coedges[(begin + index) % count]);
                index += 1;
            }
            // The way the face walks the run against the curve's own way.
            let (first_coedge, last_coedge) = (*run.first()?, *run.last()?);
            let probe = {
                let node = body.coedges.get(first_coedge)?;
                let edge = body.edges.get(node.edge)?;
                let source = body.curves.get(edge.curve)?;
                let t = 0.5 * (edge.start_parameter + edge.end_parameter);
                let sense = if node.forward { 1.0 } else { -1.0 };
                let tangent = Vec3::from(source.tangent_at(t)) * sense;
                (Vec3::from(source.point_at(t)), tangent)
            };
            let s = curve.parameter_at(probe.0.to_array());
            let along = Vec3::from(curve.derivative_at_knot(s)).dot(probe.1) > 0.0;
            let (from, to) = if all {
                let point = curve.point_at_knot(start);
                let vertex =
                    body.vertices.insert(Vertex { point, provenance: Provenance::Synthesized });
                ((start, vertex), (start + period, vertex))
            } else {
                let (a, _) = walk(body, first_coedge)?;
                let (_, b) = walk(body, last_coedge)?;
                let before = coedges[(begin + count + index - run.len() - 1) % count];
                let after = coedges[(begin + index) % count];
                let a = settle(body, &mut moved, curve, a, before, tolerance)?;
                let b = settle(body, &mut moved, curve, b, after, tolerance)?;
                (a, b)
            };
            // The stretch from the run's first corner to its last, the way
            // the face walks it.
            let (mut low, mut high, low_vertex, high_vertex) = if along {
                (from.0, to.0, from.1, to.1)
            } else {
                (to.0, from.0, to.1, from.1)
            };
            if all {
                (low, high) = (start, start + period);
            } else {
                high = low + (high - low).rem_euclid(period);
                if high - low <= period * 1e-9 {
                    high = low + period;
                }
            }
            let edge = body.edges.insert(Edge {
                curve: curve_key,
                start_parameter: low,
                end_parameter: high,
                start: low_vertex,
                end: high_vertex,
                coedges: Vec::new(),
                provenance: Provenance::Synthesized,
            });
            let coedge = body.coedges.insert(Coedge {
                edge,
                forward: along,
                pcurve: None,
                owner: ring,
                provenance: Provenance::Synthesized,
            });
            body.edges.get_mut(edge)?.coedges.push(coedge);
            rebuilt.push(coedge);
            made.push((edge, along, low));
        }
        body.loops.get_mut(ring)?.coedges = rebuilt;
    }
    made.sort_by(|a, b| a.2.total_cmp(&b.2));
    Some(made.into_iter().map(|(edge, along, _)| (edge, along)).collect())
}

/// Where corner `vertex` of the chain goes on the touching curve: where
/// the curve meets the curve of the edge `neighbour` walks, the meeting
/// nearest the corner. The neighbour's edge is moved onto it with a new
/// corner, made once for each old one.
fn settle(
    body: &mut Body,
    moved: &mut HashMap<VertexKey, (f64, VertexKey)>,
    curve: &NurbsCurve3,
    vertex: VertexKey,
    neighbour: CoedgeKey,
    tolerance: f64,
) -> Option<(f64, VertexKey)> {
    if let Some(known) = moved.get(&vertex) {
        return Some(*known);
    }
    let edge_key = body.coedges.get(neighbour)?.edge;
    let edge = body.edges.get(edge_key)?.clone();
    let other = body.curves.get(edge.curve)?.clone();
    let corner = Vec3::from(body.vertices.get(vertex)?.point);
    let (start, end) = curve.domain();
    let gap = |s: f64| {
        let point = curve.point_at_knot(s);
        let on = other.point_at(other.parameter_at(point));
        Vec3::from(point).distance(Vec3::from(on))
    };
    // Every place the touching curve comes to the neighbour, each settled,
    // the one nearest the corner kept.
    const STEPS: usize = 512;
    let at = |step: usize| start + (end - start) * step as f64 / STEPS as f64;
    let mut best: Option<(f64, f64)> = None;
    for step in 0..STEPS {
        let here = gap(at(step));
        let (before, after) = (gap(at((step + STEPS - 1) % STEPS)), gap(at(step + 1)));
        if !(here <= before && here <= after) {
            continue;
        }
        let reach = (end - start) / STEPS as f64;
        let (mut low, mut high) = (at(step) - reach, at(step) + reach);
        for _ in 0..80 {
            let third = (high - low) / 3.0;
            if gap(low + third) <= gap(high - third) {
                high -= third;
            } else {
                low += third;
            }
        }
        let s = 0.5 * (low + high);
        if gap(s) > tolerance * 10.0 {
            continue;
        }
        let distance = Vec3::from(curve.point_at_knot(s)).distance(corner);
        if best.is_none_or(|(_, kept)| distance < kept) {
            best = Some((s, distance));
        }
    }
    let (s, _) = best?;
    let point = curve.point_at_knot(s);
    let new_vertex =
        body.vertices.insert(Vertex { point, provenance: Provenance::Synthesized });
    let node = body.edges.get_mut(edge_key)?;
    let parameter = other.parameter_at(point);
    if node.start == vertex {
        node.start = new_vertex;
        node.start_parameter = parameter;
    }
    if node.end == vertex {
        node.end = new_vertex;
        node.end_parameter = parameter;
    }
    // A line's parameter runs from its origin; a closed curve's may need
    // a turn added to stay after the start.
    if node.end_parameter < node.start_parameter {
        if let Curve3::Circle(_) | Curve3::Ellipse(_) = other {
            node.end_parameter += std::f64::consts::TAU;
        }
    }
    moved.insert(vertex, (s, new_vertex));
    Some((s, new_vertex))
}
