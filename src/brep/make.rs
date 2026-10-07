//! Building a body from nothing.
//!
//! Two purposes. A modeller needs primitives to start from, and a kernel
//! needs a solid it built itself to test against ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â one lifted from a file
//! tests the reader as much as the topology, and when it fails there is no
//! saying which.
//!
//! Everything here comes out passing [`Body::validate`] with no flaws, which
//! is the property that makes it worth having: a builder that can produce an
//! inconsistent body is a builder that will.

use super::arena::Key;
use super::geometry::{
    Circle3, Cone, Curve3, Cylinder, Ellipse3, EllipticCone, Line3, Sphere, Surface, Torus,
};
use super::topology::{
    Body, Coedge, CoedgeKey, Edge, EdgeKey, Face, Loop, Lump, Shell, Vertex, VertexKey,
};
use super::Provenance;
use crate::geom2d::{Arc, Curve as Curve2, Line as Line2};
use crate::space::{Plane, Vec3};
use std::collections::{HashMap, VecDeque};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

/// A rectangular box with one corner at `origin` and the opposite at
/// `origin + size`.
///
/// Six planar faces, twelve edges, eight vertices ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â and every edge shared by
/// exactly two faces that run it opposite ways, which is what makes the
/// result a solid rather than six unrelated rectangles.
///
/// `None` for a size with a zero or negative component: a box of no thickness
/// has faces on top of each other, and no amount of care downstream recovers
/// from that.
pub fn cuboid(origin: [f64; 3], size: [f64; 3]) -> Option<Body> {
    // NaN needs naming: it compares false against everything, so a size test
    // alone would let it through and put every corner at nowhere.
    if size
        .iter()
        .any(|extent| extent.is_nan() || *extent <= 0.0 || !extent.is_finite())
    {
        return None;
    }
    let mut body = Body::new();
    let base = Vec3::from(origin);
    let (dx, dy, dz) = (size[0], size[1], size[2]);

    // The eight corners, indexed by bit: 1 = +x, 2 = +y, 4 = +z. That
    // numbering is what makes the face and edge tables below readable.
    let corners: Vec<VertexKey> = (0..8)
        .map(|bits| {
            let point = base
                + Vec3::new(
                    if bits & 1 != 0 { dx } else { 0.0 },
                    if bits & 2 != 0 { dy } else { 0.0 },
                    if bits & 4 != 0 { dz } else { 0.0 },
                );
            body.vertices.insert(Vertex {
                point: point.to_array(),
                provenance: Provenance::Synthesized,
            })
        })
        .collect();

    // The twelve edges as corner pairs, each running from the lower index to
    // the higher so a shared edge is found rather than duplicated.
    const EDGES: [(usize, usize); 12] = [
        (0, 1),
        (2, 3),
        (4, 5),
        (6, 7), // along x
        (0, 2),
        (1, 3),
        (4, 6),
        (5, 7), // along y
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7), // along z
    ];
    let edges: Vec<EdgeKey> = EDGES
        .iter()
        .map(|(from, to)| {
            let start = body.vertices.get(corners[*from])?.point;
            let end = body.vertices.get(corners[*to])?.point;
            let direction = Vec3::from(end) - Vec3::from(start);
            let curve = body.curves.insert(Curve3::Line(Line3 {
                origin: start,
                direction: direction.to_array(),
            }));
            Some(body.edges.insert(Edge {
                curve,
                // The curve is parameterised by the edge's own span, so the
                // edge runs 0 to 1 along it.
                start_parameter: 0.0,
                end_parameter: 1.0,
                start: corners[*from],
                end: corners[*to],
                coedges: Vec::new(),
                provenance: Provenance::Synthesized,
            }))
        })
        .collect::<Option<Vec<_>>>()?;

    let lump = body.lumps.insert(Lump {
        shells: Vec::new(),
        provenance: Provenance::Synthesized,
    });
    let shell = body.shells.insert(Shell {
        faces: Vec::new(),
        owner: lump,
        provenance: Provenance::Synthesized,
    });

    // Each face as the four corners of its outer loop, listed
    // counter-clockwise seen from outside the box ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â material on the left of
    // the way round, as ACIS and every other builder here have it. That
    // ordering is also what makes every edge come out traversed once each way.
    const FACES: [[usize; 4]; 6] = [
        [0, 4, 6, 2], // -x
        [1, 3, 7, 5], // +x
        [0, 1, 5, 4], // -y
        [2, 6, 7, 3], // +y
        [0, 2, 3, 1], // -z
        [4, 5, 7, 6], // +z
    ];
    for ring in FACES {
        let plane = face_plane(&body, &corners, ring)?;
        let surface = body.surfaces.insert(Surface::Plane(plane));
        let face = body.faces.insert(Face {
            surface,
            forward: true,
            loops: Vec::new(),
            owner: shell,
            provenance: Provenance::Synthesized,
        });
        let boundary = body.loops.insert(Loop {
            coedges: Vec::new(),
            owner: face,
            provenance: Provenance::Synthesized,
        });

        let mut coedges: Vec<CoedgeKey> = Vec::with_capacity(4);
        for index in 0..4 {
            let from = ring[index];
            let to = ring[(index + 1) % 4];
            let (position, forward) = find_edge(from, to)?;
            let edge = edges[position];
            let coedge = body.coedges.insert(Coedge {
                edge,
                forward,
                pcurve: None,
                owner: boundary,
                provenance: Provenance::Synthesized,
            });
            body.edges.get_mut(edge)?.coedges.push(coedge);
            coedges.push(coedge);
        }
        body.loops.get_mut(boundary)?.coedges = coedges;
        body.faces.get_mut(face)?.loops = vec![boundary];
        body.shells.get_mut(shell)?.faces.push(face);
    }

    body.lumps.get_mut(lump)?.shells = vec![shell];
    body.roots = vec![lump];
    Some(body)
}

/// Builds a closed planar-faced solid from an indexed polygon mesh.
///
/// Face winding is made consistent across every connected component and then
/// oriented so every loop runs counter-clockwise about its face's outward
/// normal, the way ACIS and the other builders wind them. Every undirected edge
/// must be shared by exactly two faces; open and non-manifold meshes are
/// rejected without producing a partial body.
pub fn faceted_solid(vertices: &[[f64; 3]], faces: &[Vec<usize>]) -> Option<Body> {
    if vertices.len() < 4
        || faces.len() < 4
        || vertices.iter().flatten().any(|value| !value.is_finite())
    {
        return None;
    }

    let mut faces = faces
        .iter()
        .map(|source| {
            let mut face = Vec::with_capacity(source.len());
            for &index in source {
                if index >= vertices.len() {
                    return None;
                }
                if face.last().copied() != Some(index) {
                    face.push(index);
                }
            }
            if face.first() == face.last() {
                face.pop();
            }
            (face.len() >= 3).then_some(face)
        })
        .collect::<Option<Vec<_>>>()?;

    type EdgeUse = (usize, bool);
    let edge_uses = |rings: &[Vec<usize>]| -> Option<HashMap<(usize, usize), Vec<EdgeUse>>> {
        let mut uses = HashMap::<(usize, usize), Vec<EdgeUse>>::new();
        for (face, ring) in rings.iter().enumerate() {
            for corner in 0..ring.len() {
                let from = ring[corner];
                let to = ring[(corner + 1) % ring.len()];
                if from == to {
                    return None;
                }
                let key = if from < to { (from, to) } else { (to, from) };
                uses.entry(key).or_default().push((face, from < to));
            }
        }
        uses.values()
            .all(|entries| entries.len() == 2)
            .then_some(uses)
    };

    let uses = edge_uses(&faces)?;
    let mut flips = vec![None; faces.len()];
    let mut components = Vec::<Vec<usize>>::new();
    for start in 0..faces.len() {
        if flips[start].is_some() {
            continue;
        }
        flips[start] = Some(false);
        let mut queue = VecDeque::from([start]);
        let mut component = Vec::new();
        while let Some(face) = queue.pop_front() {
            component.push(face);
            let flip = flips[face]?;
            let ring = &faces[face];
            for corner in 0..ring.len() {
                let from = ring[corner];
                let to = ring[(corner + 1) % ring.len()];
                let key = if from < to { (from, to) } else { (to, from) };
                let forward = from < to;
                for &(other, other_forward) in uses.get(&key)? {
                    if other == face {
                        continue;
                    }
                    let required = flip ^ (forward == other_forward);
                    match flips[other] {
                        Some(existing) if existing != required => return None,
                        Some(_) => {}
                        None => {
                            flips[other] = Some(required);
                            queue.push_back(other);
                        }
                    }
                }
            }
        }
        components.push(component);
    }
    // A disconnected or nested shell needs explicit lump/cavity containment.
    // Until that topology is represented, accept only one closed shell so a
    // cavity can never be silently converted into filled material.
    if components.len() != 1 {
        return None;
    }
    for (face, flip) in faces.iter_mut().zip(flips) {
        if flip? {
            face.reverse();
        }
    }

    for component in &components {
        let pivot = Vec3::from(vertices[faces[*component.first()?][0]]);
        let scale = component
            .iter()
            .flat_map(|&face| faces[face].iter().copied())
            .map(|vertex| (Vec3::from(vertices[vertex]) - pivot).length())
            .fold(0.0_f64, f64::max)
            .max(1e-12);
        let signed_volume = component
            .iter()
            .map(|&face| {
                let ring = &faces[face];
                let origin = Vec3::from(vertices[ring[0]]) - pivot;
                (1..ring.len() - 1)
                    .map(|corner| {
                        origin.dot(
                            (Vec3::from(vertices[ring[corner]]) - pivot)
                                .cross(Vec3::from(vertices[ring[corner + 1]]) - pivot),
                        ) / 6.0
                    })
                    .sum::<f64>()
            })
            .sum::<f64>();
        if !signed_volume.is_finite() || signed_volume.abs() <= scale.powi(3) * 1e-12 {
            return None;
        }
        // Counter-clockwise seen from outside encloses a positive volume.
        if signed_volume < 0.0 {
            for &face in component {
                faces[face].reverse();
            }
        }
    }

    let uses = edge_uses(&faces)?;
    if uses.values().any(|entries| entries[0].1 == entries[1].1) {
        return None;
    }

    let mut body = Body::new();
    let vertex_keys = vertices
        .iter()
        .map(|&point| {
            body.vertices.insert(Vertex {
                point,
                provenance: Provenance::Synthesized,
            })
        })
        .collect::<Vec<_>>();
    let mut edge_keys = HashMap::new();
    for &(from, to) in uses.keys() {
        let curve = body.curves.insert(Curve3::Line(Line3 {
            origin: vertices[from],
            direction: (Vec3::from(vertices[to]) - Vec3::from(vertices[from])).to_array(),
        }));
        let edge = body.edges.insert(Edge {
            curve,
            start_parameter: 0.0,
            end_parameter: 1.0,
            start: vertex_keys[from],
            end: vertex_keys[to],
            coedges: Vec::new(),
            provenance: Provenance::Synthesized,
        });
        edge_keys.insert((from, to), edge);
    }

    for component in components {
        let lump = body.lumps.insert(Lump {
            shells: Vec::new(),
            provenance: Provenance::Synthesized,
        });
        let shell = body.shells.insert(Shell {
            faces: Vec::new(),
            owner: lump,
            provenance: Provenance::Synthesized,
        });
        for face_index in component {
            let ring = &faces[face_index];
            let points = ring
                .iter()
                .map(|&index| vertices[index])
                .collect::<Vec<_>>();
            let outward = crate::space::polygon::normal(&points)?;
            let origin = Vec3::from(points[0]);
            let along = points[1..]
                .iter()
                .map(|point| Vec3::from(*point) - origin)
                .find(|vector| vector.length() > 1e-12)?;
            let plane = Plane::orthonormal(points[0], along.to_array(), outward)?;
            if points
                .iter()
                .any(|point| plane.distance_to(*point).is_none_or(|gap| gap.abs() > 1e-8))
            {
                return None;
            }
            let surface = body.surfaces.insert(Surface::Plane(plane));
            let face = body.faces.insert(Face {
                surface,
                forward: true,
                loops: Vec::new(),
                owner: shell,
                provenance: Provenance::Synthesized,
            });
            let boundary = body.loops.insert(Loop {
                coedges: Vec::new(),
                owner: face,
                provenance: Provenance::Synthesized,
            });
            let mut coedges = Vec::with_capacity(ring.len());
            for corner in 0..ring.len() {
                let from = ring[corner];
                let to = ring[(corner + 1) % ring.len()];
                let (key, forward) = if from < to {
                    ((from, to), true)
                } else {
                    ((to, from), false)
                };
                let edge = *edge_keys.get(&key)?;
                let coedge = body.coedges.insert(Coedge {
                    edge,
                    forward,
                    pcurve: None,
                    owner: boundary,
                    provenance: Provenance::Synthesized,
                });
                body.edges.get_mut(edge)?.coedges.push(coedge);
                coedges.push(coedge);
            }
            body.loops.get_mut(boundary)?.coedges = coedges;
            body.faces.get_mut(face)?.loops = vec![boundary];
            body.shells.get_mut(shell)?.faces.push(face);
        }
        body.lumps.get_mut(lump)?.shells = vec![shell];
        body.roots.push(lump);
    }

    body.validate().is_empty().then_some(body)
}

/// A right circular cylinder standing on `base`, `height` tall.
///
/// Three faces ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â the two discs and the wall ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â and the smallest topology that
/// describes them: two circular edges, one seam running between them, and two
/// vertices where the seam meets each rim.
///
/// The discs are bounded by a single coedge each, which is what a loop looks
/// like when its edge closes on itself. Requiring two would mean inventing a
/// split in a rim that has none.
///
/// `None` for a radius or height that is not positive.
pub fn cylinder(base: [f64; 3], radius: f64, height: f64) -> Option<Body> {
    if radius.is_nan() || height.is_nan() || radius <= 0.0 || height <= 0.0 {
        return None;
    }
    let mut body = Body::new();
    let bottom_plane = Plane::orthonormal(base, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0])?;
    let top_origin = (Vec3::from(base) + Vec3::Z * height).to_array();
    let top_plane = Plane::orthonormal(top_origin, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0])?;

    // The seam sits where the frame's own x axis meets each rim.
    let seam_low = bottom_plane.point_at([radius, 0.0]);
    let seam_high = top_plane.point_at([radius, 0.0]);
    let low = body.vertices.insert(Vertex {
        point: seam_low,
        provenance: Provenance::Synthesized,
    });
    let high = body.vertices.insert(Vertex {
        point: seam_high,
        provenance: Provenance::Synthesized,
    });

    let bottom_circle = body.curves.insert(Curve3::Circle(Circle3 {
        plane: bottom_plane,
        radius,
    }));
    let top_circle = body.curves.insert(Curve3::Circle(Circle3 {
        plane: top_plane,
        radius,
    }));
    let seam_line = body.curves.insert(Curve3::Line(Line3 {
        origin: seam_low,
        direction: (Vec3::from(seam_high) - Vec3::from(seam_low)).to_array(),
    }));

    let rim_low = body.edges.insert(Edge {
        curve: bottom_circle,
        start_parameter: 0.0,
        end_parameter: TAU,
        start: low,
        end: low,
        coedges: Vec::new(),
        provenance: Provenance::Synthesized,
    });
    let rim_high = body.edges.insert(Edge {
        curve: top_circle,
        start_parameter: 0.0,
        end_parameter: TAU,
        start: high,
        end: high,
        coedges: Vec::new(),
        provenance: Provenance::Synthesized,
    });
    let seam = body.edges.insert(Edge {
        curve: seam_line,
        start_parameter: 0.0,
        end_parameter: 1.0,
        start: low,
        end: high,
        coedges: Vec::new(),
        provenance: Provenance::Synthesized,
    });

    let lump = body.lumps.insert(Lump {
        shells: Vec::new(),
        provenance: Provenance::Synthesized,
    });
    let shell = body.shells.insert(Shell {
        faces: Vec::new(),
        owner: lump,
        provenance: Provenance::Synthesized,
    });

    // The bottom disc faces down, so its own plane ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â which faces up ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â is the
    // wrong way round and the face says so rather than a second plane being
    // made for it.
    let bottom = disc(&mut body, shell, bottom_plane, rim_low, false)?;
    let top = disc(&mut body, shell, top_plane, rim_high, true)?;
    let wall = body.surfaces.insert(Surface::Cylinder(Cylinder {
        base: bottom_plane,
        radius,
    }));
    let side = body.faces.insert(Face {
        surface: wall,
        forward: true,
        loops: Vec::new(),
        owner: shell,
        provenance: Provenance::Synthesized,
    });
    let ring = body.loops.insert(Loop {
        coedges: Vec::new(),
        owner: side,
        provenance: Provenance::Synthesized,
    });
    // Round the bottom rim, up the seam, back round the top rim, down again.
    let mut coedges = Vec::with_capacity(4);
    for (edge, forward) in [
        (rim_low, true),
        (seam, true),
        (rim_high, false),
        (seam, false),
    ] {
        let coedge = body.coedges.insert(Coedge {
            edge,
            forward,
            pcurve: None,
            owner: ring,
            provenance: Provenance::Synthesized,
        });
        body.edges.get_mut(edge)?.coedges.push(coedge);
        coedges.push(coedge);
    }
    body.loops.get_mut(ring)?.coedges = coedges;
    body.faces.get_mut(side)?.loops = vec![ring];
    body.shells.get_mut(shell)?.faces.push(side);

    let _ = (bottom, top);
    body.lumps.get_mut(lump)?.shells = vec![shell];
    body.roots = vec![lump];
    Some(body)
}

/// One end cap: a planar face bounded by the rim alone.
fn disc(
    body: &mut Body,
    shell: super::topology::ShellKey,
    plane: Plane,
    rim: EdgeKey,
    forward: bool,
) -> Option<super::topology::FaceKey> {
    let surface = body.surfaces.insert(Surface::Plane(plane));
    let face = body.faces.insert(Face {
        surface,
        forward,
        loops: Vec::new(),
        owner: shell,
        provenance: Provenance::Synthesized,
    });
    let ring = body.loops.insert(Loop {
        coedges: Vec::new(),
        owner: face,
        provenance: Provenance::Synthesized,
    });
    let coedge = body.coedges.insert(Coedge {
        edge: rim,
        // A cap's loop runs counter-clockwise seen from outside it, which for
        // the one facing down is clockwise seen from above ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â the rim
        // traversed backwards. So the sense follows the face's own, and the
        // wall takes the opposite of each: that is what makes a rim one
        // shared edge instead of two coincident ones.
        forward,
        pcurve: None,
        owner: ring,
        provenance: Provenance::Synthesized,
    });
    body.edges.get_mut(rim)?.coedges.push(coedge);
    body.loops.get_mut(ring)?.coedges = vec![coedge];
    body.faces.get_mut(face)?.loops = vec![ring];
    body.shells.get_mut(shell)?.faces.push(face);
    Some(face)
}

/// A sphere of `radius` about `centre`.
///
/// One face, as ACIS models it. A sphere has no edges of its own ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â it is
/// closed in both directions ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â so what bounds the face is a seam: a single
/// meridian, traversed once each way. The two poles are its vertices, and the
/// surface is singular there, which is a thing a B-rep says rather than
/// avoids.
///
/// `V ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€šÃ‚Â¹ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¾ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ E + F = 2 ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€šÃ‚Â¹ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¾ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ 1 + 1 = 2`, the same as any other closed shell.
pub fn sphere(centre: [f64; 3], radius: f64) -> Option<Body> {
    if radius.is_nan() || radius <= 0.0 {
        return None;
    }
    let mut body = Body::new();
    let frame = Plane::orthonormal(centre, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0])?;
    // The meridian the seam runs along: a great circle in the plane holding
    // the pole axis and the frame's own x direction, so the seam sits where
    // the surface's `u` reads zero.
    let meridian = Plane::orthonormal(centre, [1.0, 0.0, 0.0], [0.0, -1.0, 0.0])?;
    let south = add_vertex(&mut body, meridian.point_at([0.0, -radius]));
    let north = add_vertex(&mut body, meridian.point_at([0.0, radius]));

    let curve = body.curves.insert(Curve3::Circle(Circle3 {
        plane: meridian,
        radius,
    }));
    let seam = body.edges.insert(Edge {
        curve,
        start_parameter: -FRAC_PI_2,
        end_parameter: FRAC_PI_2,
        start: south,
        end: north,
        coedges: Vec::new(),
        provenance: Provenance::Synthesized,
    });

    let (lump, shell) = add_shell(&mut body);
    let surface = body
        .surfaces
        .insert(Surface::Sphere(Sphere { frame, radius }));
    // Up the seam and back down it: the loop closes on itself because the
    // face wraps the whole way round in between.
    close_shell(
        &mut body,
        lump,
        shell,
        surface,
        true,
        &[(seam, true), (seam, false)],
    )?;
    body.validate().is_empty().then_some(body)
}

/// The same sphere as [`sphere`], divided into eight triangular faces by
/// three perpendicular great circles.
///
/// A whole sphere has one face whose parameter space ends at two poles and
/// meets itself along a seam, and a cut that runs round a pole or across the
/// seam has no plain region there. Here each face has its own frame, with the
/// poles and seam well outside it, so every piece a boolean cuts from it is
/// an ordinary region of `(u, v)`. The dividing circles are tilted off the
/// world axes, so the planar cuts drawings mostly make do not run along them.
pub(crate) fn sphere_in_octants(centre: [f64; 3], radius: f64) -> Option<Body> {
    if !radius.is_finite() || radius <= 0.0 {
        return None;
    }
    let first = Vec3::new(0.9, 0.35, 0.26).normalize()?;
    let second = {
        let rough = Vec3::new(-0.31, 0.83, 0.17);
        (rough - first * rough.dot(first)).normalize()?
    };
    let axes = [first, second, first.cross(second)];
    let centre_vec = Vec3::from(centre);
    let mut body = Body::new();

    // The six places the circles cross, as (axis, positive side).
    let mut corners = HashMap::new();
    for axis in 0..3 {
        for positive in [true, false] {
            let direction = if positive { axes[axis] } else { -axes[axis] };
            let point = (centre_vec + direction * radius).to_array();
            corners.insert((axis, positive), add_vertex(&mut body, point));
        }
    }
    // Circle `k` is the one square to axis `k`; it starts on axis `k + 1` and
    // turns towards axis `k + 2`, a quarter at a time.
    let mut edges: HashMap<(VertexKey, VertexKey), EdgeKey> = HashMap::new();
    for axis in 0..3 {
        let (along, towards) = ((axis + 1) % 3, (axis + 2) % 3);
        let plane = Plane::orthonormal(centre, axes[along].to_array(), axes[axis].to_array())?;
        let curve = body.curves.insert(Curve3::Circle(Circle3 { plane, radius }));
        let round = [(along, true), (towards, true), (along, false), (towards, false)];
        for quarter in 0..4 {
            let start = corners[&round[quarter]];
            let end = corners[&round[(quarter + 1) % 4]];
            let edge = body.edges.insert(Edge {
                curve,
                start_parameter: FRAC_PI_2 * quarter as f64,
                end_parameter: FRAC_PI_2 * (quarter + 1) as f64,
                start,
                end,
                coedges: Vec::new(),
                provenance: Provenance::Synthesized,
            });
            edges.insert((start, end), edge);
        }
    }

    let (lump, shell) = add_shell(&mut body);
    for signs in 0..8 {
        let sign = |axis: usize| if signs >> axis & 1 == 0 { 1.0 } else { -1.0 };
        let middle = (axes[0] * sign(0) + axes[1] * sign(1) + axes[2] * sign(2)).normalize()?;
        // The pole axis lies square to the face's middle, leaning away from
        // all three of its corners, and the seam sits behind it.
        let pole = axes[0] * sign(0) - axes[1] * (0.8 * sign(1)) - axes[2] * (0.2 * sign(2));
        let frame = Plane::orthonormal(centre, (-middle).to_array(), pole.to_array())?;
        let surface = body
            .surfaces
            .insert(Surface::Sphere(Sphere { frame, radius }));
        // Round the corners anticlockwise seen from outside.
        let mut ring = [(0, sign(0) > 0.0), (1, sign(1) > 0.0), (2, sign(2) > 0.0)];
        if sign(0) * sign(1) * sign(2) < 0.0 {
            ring.swap(1, 2);
        }
        let mut senses = Vec::with_capacity(3);
        for index in 0..3 {
            let from = corners[&ring[index]];
            let to = corners[&ring[(index + 1) % 3]];
            senses.push(match (edges.get(&(from, to)), edges.get(&(to, from))) {
                (Some(edge), _) => (*edge, true),
                (None, Some(edge)) => (*edge, false),
                (None, None) => return None,
            });
        }
        close_shell(&mut body, lump, shell, surface, true, &senses)?;
    }
    body.validate().is_empty().then_some(body)
}

/// A circular or elliptical cylinder standing on `base`.
///
/// The two radii are the semi-axes of both planar end faces. Circular input
/// keeps the analytic cylinder representation; elliptical input stays
/// analytic too ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â the cone-surface whose section keeps the base's aspect,
/// with closed `ellipse-curve` rims ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â which is the form a native ACIS
/// author writes for an elliptical cylinder, rather than a spline
/// approximation of one.
pub fn elliptical_cylinder(
    base: [f64; 3],
    x_radius: f64,
    y_radius: f64,
    height: f64,
) -> Option<Body> {
    if [x_radius, y_radius, height]
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return None;
    }
    let scale = x_radius
        .abs()
        .max(y_radius.abs())
        .max(height.abs())
        .max(1.0);
    if (x_radius - y_radius).abs() <= 1e-12 * scale {
        return cylinder(base, x_radius, height);
    }
    let (major, ratio, direction) = canonical_ellipse_axes(x_radius, y_radius)?;
    elliptic_wall_body(base, direction, major, ratio, major, height)
}

/// The larger semi-axis, the smaller over it, and the world direction the
/// record every elliptical section is written from.
fn canonical_ellipse_axes(x_radius: f64, y_radius: f64) -> Option<(f64, f64, [f64; 3])> {
    let (major, minor, direction) = if x_radius >= y_radius {
        (x_radius, y_radius, [1.0, 0.0, 0.0])
    } else {
        (y_radius, x_radius, [0.0, 1.0, 0.0])
    };
    (minor > 0.0).then_some((major, minor / major, direction))
}

/// An elliptical cylinder, frustum or cone standing on `base` ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â the
/// analytic counterpart of the circular family, and the wire form native
/// authors give an elliptical solid: one wall face on a cone-surface whose
/// section keeps a constant aspect while its major shrinks along the axis.
///
/// The wall frame carries the major axis on `x`, so `ratio` is at most one.
/// `top_major` is the major semi-axis at the top; a frustum keeps both end
/// sections proportional, which is what makes one record cover it. A zero
/// `top_major` ends the solid in an apex point, mirroring the circular
/// cone's own construction.
fn elliptic_wall_body(
    base: [f64; 3],
    major_direction: [f64; 3],
    major: f64,
    ratio: f64,
    top_major: f64,
    height: f64,
) -> Option<Body> {
    if !(major.is_finite() && major > 0.0)
        || !(ratio.is_finite() && ratio > 0.0)
        || !(top_major.is_finite() && top_major >= 0.0)
        || !(height.is_finite() && height > 0.0)
    {
        return None;
    }
    let mut body = Body::new();
    let bottom_plane = Plane::orthonormal(base, major_direction, [0.0, 0.0, 1.0])?;
    let top_origin = (Vec3::from(base) + Vec3::Z * height).to_array();
    let top_plane = Plane::orthonormal(top_origin, major_direction, [0.0, 0.0, 1.0])?;

    // The seam sits where the frame's own major axis meets each rim, as on
    // the circular cylinder.
    let seam_low = bottom_plane.point_at([major, 0.0]);
    let seam_high = if top_major > 0.0 {
        top_plane.point_at([top_major, 0.0])
    } else {
        top_origin
    };
    let low = add_vertex(&mut body, seam_low);
    let high = add_vertex(&mut body, seam_high);

    let bottom_rim_curve = body.curves.insert(Curve3::Ellipse(Ellipse3 {
        plane: bottom_plane,
        major_radius: major,
        minor_radius: major * ratio,
    }));
    let bottom_rim = body.edges.insert(Edge {
        curve: bottom_rim_curve,
        start_parameter: 0.0,
        end_parameter: TAU,
        start: low,
        end: low,
        coedges: Vec::new(),
        provenance: Provenance::Synthesized,
    });
    let seam = add_line_edge(&mut body, low, high)?;

    let half_angle = if top_major > 0.0 {
        ((major - top_major) / height).atan()
    } else {
        (major / height).atan()
    };
    let wall = body.surfaces.insert(Surface::EllipticCone(EllipticCone {
        base: bottom_plane,
        radius: major,
        ratio,
        half_angle,
    }));

    if top_major > 0.0 {
        // A cylinder or frustum: two discs and a wall whose loop runs round
        // the bottom rim, up the seam, back round the top rim and down
        // again ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â the circular cylinder's own wiring, ellipses for circles.
        let top_rim_curve = body.curves.insert(Curve3::Ellipse(Ellipse3 {
            plane: top_plane,
            major_radius: top_major,
            minor_radius: top_major * ratio,
        }));
        let top_rim = body.edges.insert(Edge {
            curve: top_rim_curve,
            start_parameter: 0.0,
            end_parameter: TAU,
            start: high,
            end: high,
            coedges: Vec::new(),
            provenance: Provenance::Synthesized,
        });
        let (lump, shell) = add_shell(&mut body);
        disc(&mut body, shell, bottom_plane, bottom_rim, false)?;
        disc(&mut body, shell, top_plane, top_rim, true)?;

        let side = body.faces.insert(Face {
            surface: wall,
            forward: true,
            loops: Vec::new(),
            owner: shell,
            provenance: Provenance::Synthesized,
        });
        let ring = body.loops.insert(Loop {
            coedges: Vec::new(),
            owner: side,
            provenance: Provenance::Synthesized,
        });
        let mut coedges = Vec::with_capacity(4);
        for (edge, forward) in [
            (bottom_rim, true),
            (seam, true),
            (top_rim, false),
            (seam, false),
        ] {
            let coedge = body.coedges.insert(Coedge {
                edge,
                forward,
                pcurve: None,
                owner: ring,
                provenance: Provenance::Synthesized,
            });
            body.edges.get_mut(edge)?.coedges.push(coedge);
            coedges.push(coedge);
        }
        body.loops.get_mut(ring)?.coedges = coedges;
        body.faces.get_mut(side)?.loops = vec![ring];
        body.shells.get_mut(shell)?.faces.push(side);
        body.lumps.get_mut(lump)?.shells = vec![shell];
        body.roots = vec![lump];
    } else {
        // A cone: the wall's loop runs round the rim and up and back down
        // the seam, whose meeting point at the apex is the vertex the
        // surface is singular at ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â the circular cone's own construction.
        let (lump, shell) = add_shell(&mut body);
        disc(&mut body, shell, bottom_plane, bottom_rim, false)?;
        close_shell(
            &mut body,
            lump,
            shell,
            wall,
            true,
            &[(bottom_rim, true), (seam, true), (seam, false)],
        )?;
    }
    body.validate().is_empty().then_some(body)
}

/// A cone standing on `base`, `height` tall, with a base circle of `radius`
/// and a point at the top.
///
/// Two faces: the cone-surface wall and the disc it stands on. The wall's
/// loop runs round the rim, up the seam and back down it ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â the apex is where
/// the seam meets itself, a vertex the surface is singular at.
///
/// `None` for a size that is not positive. A truncated cone is a different
/// shape and is not this: pass a cylinder's radius twice, or build it from
/// the surface directly.
pub fn cone(base: [f64; 3], radius: f64, height: f64) -> Option<Body> {
    if radius.is_nan() || height.is_nan() || radius <= 0.0 || height <= 0.0 {
        return None;
    }
    let mut body = Body::new();
    let base_plane = Plane::orthonormal(base, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0])?;
    let rim_point = base_plane.point_at([radius, 0.0]);
    let apex_point = (Vec3::from(base) + Vec3::Z * height).to_array();
    let rim_vertex = add_vertex(&mut body, rim_point);
    let apex = add_vertex(&mut body, apex_point);

    let rim_curve = body.curves.insert(Curve3::Circle(Circle3 {
        plane: base_plane,
        radius,
    }));
    let rim = body.edges.insert(Edge {
        curve: rim_curve,
        start_parameter: 0.0,
        end_parameter: TAU,
        start: rim_vertex,
        end: rim_vertex,
        coedges: Vec::new(),
        provenance: Provenance::Synthesized,
    });
    let seam = add_line_edge(&mut body, rim_vertex, apex)?;

    let (lump, shell) = add_shell(&mut body);
    // The base looks down, so its own plane ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â which looks up ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â is the wrong
    // way round and the face says so.
    disc(&mut body, shell, base_plane, rim, false)?;
    let wall = body.surfaces.insert(Surface::Cone(Cone {
        base: base_plane,
        radius,
        // The radius falls to nothing over the height, which is what sets
        // the slope.
        half_angle: (radius / height).atan(),
    }));
    close_shell(
        &mut body,
        lump,
        shell,
        wall,
        true,
        &[(rim, true), (seam, true), (seam, false)],
    )?;
    body.validate().is_empty().then_some(body)
}

/// A circular or elliptical cone/frustum standing on `base`.
///
/// `base_x_radius` and `base_y_radius` are the base semi-axes. `top_radius`
/// is the top X semi-axis; its Y semi-axis keeps the base aspect ratio. A
/// zero top radius produces a pointed cone. Every tier of it stays analytic:
/// circular cones and frusta keep their own records, and elliptical ones
/// ride the same cone-surface record their circular special case does ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â a
/// native ACIS author's form for all of them ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â rather than a spline
/// approximation.
pub fn frustum(
    base: [f64; 3],
    base_x_radius: f64,
    base_y_radius: f64,
    top_radius: f64,
    height: f64,
) -> Option<Body> {
    if [base_x_radius, base_y_radius, height]
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
        || !top_radius.is_finite()
        || top_radius < 0.0
    {
        return None;
    }
    let scale = base_x_radius
        .abs()
        .max(base_y_radius.abs())
        .max(top_radius.abs())
        .max(height.abs())
        .max(1.0);
    let circular = (base_x_radius - base_y_radius).abs() <= 1e-12 * scale;
    if top_radius <= 1e-12 * scale {
        return if circular {
            cone(base, base_x_radius, height)
        } else {
            elliptical_cone(base, base_x_radius, base_y_radius, height)
        };
    }

    if base_x_radius == base_y_radius {
        let section = Plane::from_axes(base, Vec3::X.to_array(), Vec3::Z.to_array());
        let points = [
            [0.0, 0.0],
            [base_x_radius, 0.0],
            [top_radius, height],
            [0.0, height],
        ];
        let profile = (0..4)
            .map(|index| {
                Curve2::Line(Line2 {
                    start: points[index],
                    end: points[(index + 1) % 4],
                })
            })
            .collect::<Vec<_>>();
        return super::sweep::revolve(section, &profile, base, Vec3::Z.to_array(), TAU);
    }

    // An elliptical frustum whose end sections keep the base's aspect ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â one
    // analytic cone-surface, as every native-authored instance of the shape.
    let (major, ratio, direction) = canonical_ellipse_axes(base_x_radius, base_y_radius)?;
    let top_major = if base_x_radius >= base_y_radius {
        top_radius
    } else {
        top_radius * base_y_radius / base_x_radius
    };
    elliptic_wall_body(base, direction, major, ratio, top_major, height)
}


fn elliptical_cone(base: [f64; 3], x_radius: f64, y_radius: f64, height: f64) -> Option<Body> {
    let (major, ratio, direction) = canonical_ellipse_axes(x_radius, y_radius)?;
    elliptic_wall_body(base, direction, major, ratio, 0.0, height)
}

/// A torus about `centre`, its tube `minor_radius` thick at
/// `major_radius` out.
///
/// One face and two seams ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â one round the ring, one round the tube ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â which is
/// what a surface closed in both directions needs. Its Euler characteristic
/// is zero rather than two, because a torus has a hole through it and no
/// closed shell of genus one can have any other.
///
/// A tube smaller than the major radius produces the usual ring. At or past
/// the major radius the self-intersecting half of the generating circle is
/// trimmed at the axis, producing a closed horn or spindle solid.
pub fn torus(centre: [f64; 3], major_radius: f64, minor_radius: f64) -> Option<Body> {
    if centre.iter().any(|coordinate| !coordinate.is_finite())
        || !major_radius.is_finite()
        || !minor_radius.is_finite()
        || major_radius <= 0.0
        || minor_radius <= 0.0
    {
        return None;
    }
    if major_radius <= minor_radius {
        return closed_hole_torus(centre, major_radius, minor_radius);
    }
    let mut body = Body::new();
    let frame = Plane::orthonormal(centre, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0])?;
    // Both seams pass through the outermost point of the ring, so the whole
    // face has one vertex.
    let outer = frame.point_at([major_radius + minor_radius, 0.0]);
    let corner = add_vertex(&mut body, outer);

    // Round the ring, at the tube's outer edge.
    let around = body.curves.insert(Curve3::Circle(Circle3 {
        plane: frame,
        radius: major_radius + minor_radius,
    }));
    let ring = body.edges.insert(Edge {
        curve: around,
        start_parameter: 0.0,
        end_parameter: TAU,
        start: corner,
        end: corner,
        coedges: Vec::new(),
        provenance: Provenance::Synthesized,
    });
    // Round the tube, in the plane holding the axis and the frame's x.
    let tube_plane = Plane::orthonormal(
        frame.point_at([major_radius, 0.0]),
        [1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    )?;
    let across = body.curves.insert(Curve3::Circle(Circle3 {
        plane: tube_plane,
        radius: minor_radius,
    }));
    let tube = body.edges.insert(Edge {
        curve: across,
        start_parameter: 0.0,
        end_parameter: TAU,
        start: corner,
        end: corner,
        coedges: Vec::new(),
        provenance: Provenance::Synthesized,
    });

    let (lump, shell) = add_shell(&mut body);
    let surface = body.surfaces.insert(Surface::Torus(Torus {
        frame,
        major_radius,
        minor_radius,
    }));
    close_shell(
        &mut body,
        lump,
        shell,
        surface,
        true,
        &[(ring, true), (tube, true), (ring, false), (tube, false)],
    )?;
    body.validate().is_empty().then_some(body)
}

fn closed_hole_torus(centre: [f64; 3], major_radius: f64, minor_radius: f64) -> Option<Body> {
    let profile_plane = Plane::orthonormal(centre, [1.0, 0.0, 0.0], [0.0, -1.0, 0.0])?;
    let intersection_angle = (major_radius / minor_radius).clamp(0.0, 1.0).acos();
    let half_axis_length = (minor_radius * minor_radius - major_radius * major_radius)
        .max(0.0)
        .sqrt();
    let mut profile = vec![Curve2::Arc(Arc {
        centre: [major_radius, 0.0],
        radius: minor_radius,
        start_angle: PI + intersection_angle,
        end_angle: PI - intersection_angle,
    })];
    if half_axis_length > f64::EPSILON * minor_radius.max(1.0) {
        profile.push(Curve2::Line(Line2 {
            start: [0.0, half_axis_length],
            end: [0.0, 0.0],
        }));
        profile.push(Curve2::Line(Line2 {
            start: [0.0, 0.0],
            end: [0.0, -half_axis_length],
        }));
    }
    super::sweep::revolve(profile_plane, &profile, centre, [0.0, 0.0, 1.0], TAU)
}

/// A right triangular prism: a right triangle in the XZ plane at `origin`,
/// `width` deep along Y.
pub fn wedge(origin: [f64; 3], length: f64, width: f64, height: f64) -> Option<Body> {
    if [length, width, height]
        .iter()
        .any(|extent| extent.is_nan() || *extent <= 0.0)
    {
        return None;
    }
    // The profile lies in the XZ plane, so the sweep runs along Y.
    let plane = Plane::orthonormal(origin, [1.0, 0.0, 0.0], [0.0, -1.0, 0.0])?;
    let triangle = [[0.0, 0.0], [length, 0.0], [0.0, height]];
    let profile: Vec<crate::geom2d::Curve> = (0..3)
        .map(|index| {
            crate::geom2d::Curve::Line(crate::geom2d::Line {
                start: triangle[index],
                end: triangle[(index + 1) % 3],
            })
        })
        .collect();
    super::sweep::extrude(plane, &profile, [0.0, width, 0.0])
}

/// A pointed pyramid on a regular polygon of `sides` corners.
pub fn pyramid(base: [f64; 3], radius: f64, height: f64, sides: usize) -> Option<Body> {
    pyramid_frustum(base, radius, 0.0, height, sides)
}

/// A pyramid or polygonal frustum whose base and top are concentric regular
/// polygons with the same number of sides and orientation.
///
/// A zero `top_radius` produces a pointed pyramid. A positive value produces
/// a closed frustum with a planar top and one quadrilateral face per side.
///
/// `None` is returned for fewer than three sides, a non-positive base radius
/// or height, a negative top radius, or any non-finite dimension.
pub fn pyramid_frustum(
    base: [f64; 3],
    base_radius: f64,
    top_radius: f64,
    height: f64,
    sides: usize,
) -> Option<Body> {
    if sides < 3
        || [base_radius, top_radius, height]
            .iter()
            .any(|value| !value.is_finite())
        || base_radius <= 0.0
        || top_radius < 0.0
        || height <= 0.0
    {
        return None;
    }
    let mut body = Body::new();
    let ground = Plane::orthonormal(base, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0])?;
    let base_corners: Vec<VertexKey> = (0..sides)
        .map(|index| {
            let angle = TAU * index as f64 / sides as f64;
            add_vertex(
                &mut body,
                ground.point_at([base_radius * angle.cos(), base_radius * angle.sin()]),
            )
        })
        .collect();

    let base_rim: Vec<EdgeKey> = (0..sides)
        .map(|index| {
            add_line_edge(
                &mut body,
                base_corners[index],
                base_corners[(index + 1) % sides],
            )
        })
        .collect::<Option<Vec<_>>>()?;

    let (lump, shell) = add_shell(&mut body);

    // The base looks down, so its loop runs the rim backwards ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â which is also
    // what leaves each rim edge traversed one way by the base and the other
    // by the side standing on it.
    let ground_down = Plane::orthonormal(base, [1.0, 0.0, 0.0], [0.0, 0.0, -1.0])?;
    let surface = body.surfaces.insert(Surface::Plane(ground_down));
    let floor: Vec<(EdgeKey, bool)> = base_rim.iter().rev().map(|edge| (*edge, false)).collect();
    close_shell(&mut body, lump, shell, surface, true, &floor)?;

    if top_radius <= 1e-12 {
        let apex = add_vertex(&mut body, (Vec3::from(base) + Vec3::Z * height).to_array());
        let rails: Vec<EdgeKey> = base_corners
            .iter()
            .map(|corner| add_line_edge(&mut body, *corner, apex))
            .collect::<Option<Vec<_>>>()?;

        for index in 0..sides {
            let next = (index + 1) % sides;
            let here = Vec3::from(body.vertices.get(base_corners[index])?.point);
            let along = Vec3::from(body.vertices.get(base_corners[next])?.point) - here;
            let up = Vec3::from(body.vertices.get(apex)?.point) - here;
            let surface = body.surfaces.insert(Surface::Plane(Plane::orthonormal(
                here.to_array(),
                along.to_array(),
                along.cross(up).normalize()?.to_array(),
            )?));
            close_shell(
                &mut body,
                lump,
                shell,
                surface,
                true,
                &[
                    (base_rim[index], true),
                    (rails[next], true),
                    (rails[index], false),
                ],
            )?;
        }
        return body.validate().is_empty().then_some(body);
    }

    let top_center = Vec3::from(base) + Vec3::Z * height;
    let top_corners: Vec<VertexKey> = (0..sides)
        .map(|index| {
            let angle = TAU * index as f64 / sides as f64;
            add_vertex(
                &mut body,
                (top_center + Vec3::new(top_radius * angle.cos(), top_radius * angle.sin(), 0.0))
                    .to_array(),
            )
        })
        .collect();
    let top_rim: Vec<EdgeKey> = (0..sides)
        .map(|index| {
            add_line_edge(
                &mut body,
                top_corners[index],
                top_corners[(index + 1) % sides],
            )
        })
        .collect::<Option<Vec<_>>>()?;
    let rails: Vec<EdgeKey> = (0..sides)
        .map(|index| add_line_edge(&mut body, base_corners[index], top_corners[index]))
        .collect::<Option<Vec<_>>>()?;

    let top_plane = Plane::orthonormal(top_center.to_array(), [1.0, 0.0, 0.0], [0.0, 0.0, 1.0])?;
    let surface = body.surfaces.insert(Surface::Plane(top_plane));
    let ceiling: Vec<(EdgeKey, bool)> = top_rim.iter().map(|edge| (*edge, true)).collect();
    close_shell(&mut body, lump, shell, surface, true, &ceiling)?;

    for index in 0..sides {
        let next = (index + 1) % sides;
        let here = Vec3::from(body.vertices.get(base_corners[index])?.point);
        let along = Vec3::from(body.vertices.get(base_corners[next])?.point) - here;
        let up = Vec3::from(body.vertices.get(top_corners[index])?.point) - here;
        let surface = body.surfaces.insert(Surface::Plane(Plane::orthonormal(
            here.to_array(),
            along.to_array(),
            along.cross(up).normalize()?.to_array(),
        )?));
        close_shell(
            &mut body,
            lump,
            shell,
            surface,
            true,
            &[
                (base_rim[index], true),
                (rails[next], true),
                (top_rim[index], false),
                (rails[index], false),
            ],
        )?;
    }
    body.validate().is_empty().then_some(body)
}

/// A lump with one shell in it, ready for faces.
fn add_shell(body: &mut Body) -> (super::topology::LumpKey, super::topology::ShellKey) {
    let lump = body.lumps.insert(Lump {
        shells: Vec::new(),
        provenance: Provenance::Synthesized,
    });
    let shell = body.shells.insert(Shell {
        faces: Vec::new(),
        owner: lump,
        provenance: Provenance::Synthesized,
    });
    (lump, shell)
}

/// Adds one face from a list of edges and senses, then closes the body around
/// it.
fn close_shell(
    body: &mut Body,
    lump: super::topology::LumpKey,
    shell: super::topology::ShellKey,
    surface: super::topology::SurfaceKey,
    forward: bool,
    edges: &[(EdgeKey, bool)],
) -> Option<()> {
    let face = body.faces.insert(Face {
        surface,
        forward,
        loops: Vec::new(),
        owner: shell,
        provenance: Provenance::Synthesized,
    });
    let ring = body.loops.insert(Loop {
        coedges: Vec::new(),
        owner: face,
        provenance: Provenance::Synthesized,
    });
    let mut coedges = Vec::with_capacity(edges.len());
    for (edge, sense) in edges {
        let coedge = body.coedges.insert(Coedge {
            edge: *edge,
            forward: *sense,
            pcurve: None,
            owner: ring,
            provenance: Provenance::Synthesized,
        });
        body.edges.get_mut(*edge)?.coedges.push(coedge);
        coedges.push(coedge);
    }
    body.loops.get_mut(ring)?.coedges = coedges;
    body.faces.get_mut(face)?.loops = vec![ring];
    body.shells.get_mut(shell)?.faces.push(face);
    if body.lumps.get(lump)?.shells.is_empty() {
        body.lumps.get_mut(lump)?.shells = vec![shell];
    }
    if body.roots.is_empty() {
        body.roots = vec![lump];
    }
    Some(())
}


fn add_vertex(body: &mut Body, point: [f64; 3]) -> VertexKey {
    body.vertices.insert(Vertex {
        point,
        provenance: Provenance::Synthesized,
    })
}

fn add_line_edge(body: &mut Body, from: VertexKey, to: VertexKey) -> Option<EdgeKey> {
    let start = body.vertices.get(from)?.point;
    let end = body.vertices.get(to)?.point;
    let curve = body.curves.insert(Curve3::Line(Line3 {
        origin: start,
        direction: (Vec3::from(end) - Vec3::from(start)).to_array(),
    }));
    Some(body.edges.insert(Edge {
        curve,
        start_parameter: 0.0,
        end_parameter: 1.0,
        start: from,
        end: to,
        coedges: Vec::new(),
        provenance: Provenance::Synthesized,
    }))
}

/// Which of the twelve edges joins two corners, and whether the given order
/// runs along it or against it.
fn find_edge(from: usize, to: usize) -> Option<(usize, bool)> {
    const EDGES: [(usize, usize); 12] = [
        (0, 1),
        (2, 3),
        (4, 5),
        (6, 7),
        (0, 2),
        (1, 3),
        (4, 6),
        (5, 7),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];
    EDGES.iter().enumerate().find_map(|(index, (a, b))| {
        if (*a, *b) == (from, to) {
            Some((index, true))
        } else if (*a, *b) == (to, from) {
            Some((index, false))
        } else {
            None
        }
    })
}

/// The plane a face's four corners lie in, with its normal pointing out of
/// the box.
///
/// Built from the ring rather than from a table of axis directions, so the
/// two cannot disagree: the loop's own winding is what decides which way is
/// out.
fn face_plane(body: &Body, corners: &[Key<Vertex>], ring: [usize; 4]) -> Option<Plane> {
    let at = |index: usize| -> Option<Vec3> {
        Some(Vec3::from(body.vertices.get(corners[ring[index]])?.point))
    };
    let origin = at(0)?;
    let along = at(1)? - origin;
    let across = at(3)? - origin;
    // Counter-clockwise seen from outside, so along x across points out.
    // The frame runs along `across`, which keeps each face's parameterisation
    // what it was before the loops were turned the right way round.
    let normal = along.cross(across).normalize()?;
    Plane::orthonormal(origin.to_array(), across.to_array(), normal.to_array())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep::topology::Flaw;
    use std::collections::HashSet;

    fn unit_box() -> Body {
        cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).expect("a unit box")
    }

    #[test]
    fn a_closed_faceted_mesh_becomes_a_valid_solid() {
        let vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let faces = vec![vec![0, 2, 1], vec![0, 1, 3], vec![1, 2, 3], vec![2, 0, 3]];

        let body = faceted_solid(&vertices, &faces).expect("a tetrahedron");
        assert_eq!(body.faces.len(), 4);
        assert_eq!(body.euler_characteristic(), 2);
        assert!(body.validate().is_empty());
    }

    #[test]
    fn an_open_faceted_mesh_is_refused() {
        let vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let faces = vec![vec![0, 2, 1], vec![0, 1, 3], vec![1, 2, 3]];

        assert!(faceted_solid(&vertices, &faces).is_none());
    }

    #[test]
    fn a_box_has_the_parts_a_box_has() {
        let body = unit_box();
        assert_eq!(body.vertices.len(), 8);
        assert_eq!(body.edges.len(), 12);
        assert_eq!(body.faces.len(), 6);
        assert_eq!(body.coedges.len(), 24, "four per face");
        assert_eq!(body.loops.len(), 6);
        assert_eq!(body.shells.len(), 1);
        assert_eq!(body.lumps.len(), 1);
    }

    #[test]
    fn a_box_is_a_closed_surface() {
        // V ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€šÃ‚Â¹ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¾ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ E + F = 8 ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€šÃ‚Â¹ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¾ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ 12 + 6 = 2, which is what a shell of genus zero
        // must give. A face left out or an edge duplicated changes it.
        assert_eq!(unit_box().euler_characteristic(), 2);
    }

    #[test]
    fn a_box_has_nothing_wrong_with_it() {
        let flaws = unit_box().validate();
        assert!(flaws.is_empty(), "{flaws:?}");
    }

    #[test]
    fn every_edge_is_shared_by_two_faces_running_it_opposite_ways() {
        // The property that makes it a solid. Getting one face's winding
        // backwards leaves that face's edges traversed the same way on both
        // sides, and the box is inside out along that seam.
        let body = unit_box();
        for (key, edge) in body.edges.iter() {
            assert_eq!(edge.coedges.len(), 2, "edge {key:?}");
            let senses: Vec<bool> = edge
                .coedges
                .iter()
                .map(|c| body.coedges.get(*c).unwrap().forward)
                .collect();
            assert_ne!(senses[0], senses[1], "edge {key:?} is used twice one way");
        }
    }

    #[test]
    fn every_loop_closes() {
        let body = unit_box();
        for (key, ring) in body.loops.iter() {
            let count = ring.coedges.len();
            assert_eq!(count, 4, "loop {key:?}");
            for index in 0..count {
                let (_, ends) = body.coedge_vertices(ring.coedges[index]).unwrap();
                let (begins, _) = body
                    .coedge_vertices(ring.coedges[(index + 1) % count])
                    .unwrap();
                assert_eq!(ends, begins, "loop {key:?} breaks after {index}");
            }
        }
    }

    #[test]
    fn every_face_lies_on_its_own_surface() {
        let body = unit_box();
        for (key, face) in body.faces.iter() {
            let surface = body.surfaces.get(face.surface).unwrap();
            for coedge in body.face_coedges(key) {
                let (start, _) = body.coedge_vertices(coedge).unwrap();
                let point = body.vertices.get(start).unwrap().point;
                assert!(surface.contains(point, 1e-9), "{point:?} off face {key:?}");
            }
        }
    }

    #[test]
    fn every_face_normal_points_out_of_the_box() {
        // The centre is inside; a face's normal must lead away from it.
        let body = cuboid([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]).unwrap();
        let centre = Vec3::new(1.0, 2.0, 3.0);
        for (_, face) in body.faces.iter() {
            let Surface::Plane(plane) = body.surfaces.get(face.surface).unwrap() else {
                panic!("a box is planar");
            };
            let normal = Vec3::from(plane.normal().unwrap());
            let outward = Vec3::from(plane.origin) - centre;
            assert!(normal.dot(outward) > 0.0, "a face pointed inwards");
        }
    }

    #[test]
    fn every_vertex_sits_where_its_edges_end() {
        assert!(unit_box().worst_vertex_gap() < 1e-12);
    }

    #[test]
    fn the_corners_are_where_the_size_puts_them() {
        let body = cuboid([10.0, 20.0, 30.0], [1.0, 2.0, 3.0]).unwrap();
        let points: HashSet<[u64; 3]> = body
            .vertices
            .iter()
            .map(|(_, v)| {
                [
                    v.point[0].to_bits(),
                    v.point[1].to_bits(),
                    v.point[2].to_bits(),
                ]
            })
            .collect();
        assert_eq!(points.len(), 8, "no two corners coincide");
        for corner in [[10.0_f64, 20.0, 30.0], [11.0, 22.0, 33.0]] {
            let bits = [
                corner[0].to_bits(),
                corner[1].to_bits(),
                corner[2].to_bits(),
            ];
            assert!(points.contains(&bits), "{corner:?} missing");
        }
    }

    #[test]
    fn each_coedge_has_the_one_on_the_far_side() {
        let body = unit_box();
        for key in body.coedges.keys() {
            let partner = body.partner(key).expect("a closed box shares every edge");
            assert_ne!(partner, key);
            assert_eq!(body.partner(partner), Some(key));
        }
    }

    #[test]
    fn a_cylinder_has_the_parts_a_cylinder_has() {
        let solid = cylinder([0.0; 3], 5.0, 10.0).expect("a cylinder");
        assert_eq!(solid.faces.len(), 3, "two caps and a wall");
        assert_eq!(solid.edges.len(), 3, "two rims and a seam");
        assert_eq!(solid.vertices.len(), 2, "where the seam meets each rim");
        let flaws = solid.validate();
        assert!(flaws.is_empty(), "{flaws:?}");
        assert_eq!(solid.euler_characteristic(), 2);
    }

    #[test]
    fn a_cylinders_rims_are_each_shared_by_two_faces() {
        let solid = cylinder([0.0; 3], 3.0, 4.0).unwrap();
        for (key, edge) in solid.edges.iter() {
            assert_eq!(edge.coedges.len(), 2, "edge {key:?}");
            let senses: Vec<bool> = edge
                .coedges
                .iter()
                .map(|c| solid.coedges.get(*c).unwrap().forward)
                .collect();
            assert_ne!(senses[0], senses[1], "edge {key:?} runs one way twice");
        }
    }

    #[test]
    fn a_sphere_is_one_face_bounded_by_a_seam() {
        // ACIS's own shape for it: no edges of its own, so a meridian
        // traversed once each way is what bounds the face, and the poles are
        // vertices the surface is singular at.
        let solid = sphere([0.0; 3], 5.0).expect("a sphere");
        assert_eq!(solid.faces.len(), 1);
        assert_eq!(solid.edges.len(), 1, "one seam");
        assert_eq!(solid.vertices.len(), 2, "two poles");
        assert_eq!(solid.euler_characteristic(), 2);
        let flaws = solid.validate();
        assert!(flaws.is_empty(), "{flaws:?}");
    }

    #[test]
    fn a_spheres_poles_are_a_radius_apart_along_its_axis() {
        let solid = sphere([1.0, 2.0, 3.0], 4.0).unwrap();
        let mut heights: Vec<f64> = solid.vertices.iter().map(|(_, v)| v.point[2]).collect();
        heights.sort_by(f64::total_cmp);
        assert!((heights[0] - (-1.0)).abs() < 1e-9, "{heights:?}");
        assert!((heights[1] - 7.0).abs() < 1e-9, "{heights:?}");
    }

    #[test]
    fn a_spheres_seam_lies_on_the_sphere() {
        let solid = sphere([0.0; 3], 5.0).unwrap();
        let (_, face) = solid.faces.iter().next().unwrap();
        let surface = solid.surfaces.get(face.surface).unwrap();
        let (key, edge) = solid.edges.iter().next().unwrap();
        let curve = solid.curves.get(edge.curve).unwrap();
        for step in 0..=8 {
            let t = edge.start_parameter
                + (edge.end_parameter - edge.start_parameter) * step as f64 / 8.0;
            let point = curve.point_at(t);
            assert!(surface.contains(point, 1e-9), "{point:?} off the sphere");
        }
        let _ = key;
    }

    #[test]
    fn a_cone_is_a_wall_and_a_base_meeting_at_a_point() {
        let solid = cone([0.0; 3], 5.0, 12.0).expect("a cone");
        assert_eq!(solid.faces.len(), 2, "the wall and the disc");
        assert_eq!(solid.edges.len(), 2, "the rim and the seam");
        assert_eq!(solid.vertices.len(), 2, "the seam's foot and the apex");
        assert_eq!(solid.euler_characteristic(), 2);
        assert!(solid.validate().is_empty());
    }

    #[test]
    fn a_cones_surface_reaches_its_own_apex() {
        let solid = cone([0.0; 3], 5.0, 12.0).unwrap();
        let apex = solid
            .vertices
            .iter()
            .map(|(_, v)| v.point)
            .find(|p| p[2] > 1.0)
            .expect("an apex above the base");
        assert!((apex[2] - 12.0).abs() < 1e-9, "{apex:?}");
        // The wall's surface holds it: the radius really does fall to nothing
        // over the height, which the half-angle is what sets.
        let wall = solid
            .surfaces
            .iter()
            .find_map(|(_, s)| matches!(s, Surface::Cone(_)).then_some(s))
            .unwrap();
        assert!(wall.contains(apex, 1e-9), "the apex is off its own cone");
    }

    #[test]
    fn a_torus_has_a_hole_through_it() {
        // Euler zero rather than two: no closed shell of genus one can have
        // any other, so this is the check that the two seams are right.
        let solid = torus([0.0; 3], 10.0, 2.0).expect("a torus");
        assert_eq!(solid.faces.len(), 1);
        assert_eq!(solid.edges.len(), 2, "one seam each way");
        assert_eq!(solid.vertices.len(), 1, "where they cross");
        assert_eq!(solid.euler_characteristic(), 0);
        let flaws = solid.validate();
        assert!(flaws.is_empty(), "{flaws:?}");
    }

    #[test]
    fn horn_and_spindle_tori_are_trimmed_at_the_axis() {
        for (major, minor) in [(2.0, 2.0), (1.0, 2.0)] {
            let solid = torus([0.0; 3], major, minor).expect("a closed torus");
            let flaws = solid.validate();
            assert!(flaws.is_empty(), "{major}/{minor}: {flaws:?}");
        }
        assert!(torus([0.0; 3], 10.0, 0.0).is_none());
    }

    #[test]
    fn both_of_a_toruss_seams_lie_on_it() {
        let solid = torus([0.0; 3], 10.0, 2.0).unwrap();
        let (_, face) = solid.faces.iter().next().unwrap();
        let surface = solid.surfaces.get(face.surface).unwrap();
        for (_, edge) in solid.edges.iter() {
            let curve = solid.curves.get(edge.curve).unwrap();
            for step in 0..8 {
                let point = curve.point_at(TAU * step as f64 / 8.0);
                assert!(surface.contains(point, 1e-9), "{point:?} off the torus");
            }
        }
    }

    #[test]
    fn a_wedge_is_a_prism_on_a_right_triangle() {
        let solid = wedge([0.0; 3], 4.0, 3.0, 5.0).expect("a wedge");
        assert_eq!(solid.faces.len(), 5, "two triangles and three walls");
        assert_eq!(solid.edges.len(), 9);
        assert_eq!(solid.vertices.len(), 6);
        assert_eq!(solid.euler_characteristic(), 2);
        assert!(solid.validate().is_empty());
        let bounds = crate::brep::body_bounds(&solid).unwrap();
        assert_eq!(bounds.max, [4.0, 3.0, 5.0]);
    }

    /// How much a meshed solid encloses, by the divergence theorem.
    ///
    /// Negative when the faces are wound inwards, and short of the true
    /// figure by however much the flat facets cut off ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â never over, since
    /// every chord lies inside the surface it approximates.
    fn meshed_volume(solid: &Body) -> f64 {
        let mesh = crate::brep::mesh::body(solid, crate::tessellation::DEFAULT_ANGLE, 1e-9);
        assert!(!mesh.is_empty(), "nothing meshed");
        mesh.triangles
            .iter()
            .map(|triangle| {
                let at = |index: usize| Vec3::from(mesh.positions[triangle[index]]);
                at(0).cross(at(1)).dot(at(2)) / 6.0
            })
            .sum()
    }

    #[test]
    fn every_primitive_that_meshes_encloses_what_it_should() {
        // Measuring the whole mesh rather than looking at each triangle. A
        // face left out entirely passes any per-triangle check ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â the ones
        // that are there are still wound correctly ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â and for a long while a
        // cylinder's wall was missing from every mesh this way, with only
        // its two discs drawn.
        //
        // A volume notices both faults at once: a missing face reads far too
        // small, and an inside-out one reads negative.
        let cases: [(Body, f64); 7] = [
            (cuboid([0.0; 3], [4.0; 3]).unwrap(), 64.0),
            (
                wedge([0.0; 3], 4.0, 3.0, 5.0).unwrap(),
                0.5 * 4.0 * 5.0 * 3.0,
            ),
            (
                cylinder([0.0; 3], 3.0, 6.0).unwrap(),
                std::f64::consts::PI * 9.0 * 6.0,
            ),
            (
                cone([0.0; 3], 5.0, 12.0).unwrap(),
                std::f64::consts::PI * 25.0 * 12.0 / 3.0,
            ),
            (
                torus([0.0; 3], 10.0, 2.0).unwrap(),
                2.0 * std::f64::consts::PI.powi(2) * 10.0 * 4.0,
            ),
            (
                sphere([0.0; 3], 5.0).unwrap(),
                4.0 / 3.0 * std::f64::consts::PI * 125.0,
            ),
            (
                pyramid([0.0; 3], 4.0, 9.0, 6).unwrap(),
                // A regular hexagon of circumradius four: six equilateral
                // triangles, each with area rÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â²ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€šÃ‚Â¹ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€šÃ‚Â¦ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡3/4 at side r.
                6.0 * (16.0 * 3.0_f64.sqrt() / 4.0) * 9.0 / 3.0,
            ),
        ];
        for (solid, expected) in cases {
            let volume = meshed_volume(&solid);
            assert!(volume > 0.0, "wound inwards: {volume}");
            assert!(volume <= expected * 1.000_001, "reads over: {volume}");
            assert!(volume > 0.98 * expected, "{volume} vs {expected}");
        }
    }

    #[test]
    fn the_primitives_that_mesh_do_so_with_their_normals_out() {
        // The check that each one's senses hang together, for the shapes
        // convex enough that every normal must lead away from one point
        // inside: a face wound the wrong way lights the solid inside out.
        let solids = [
            (
                cuboid([0.0; 3], [4.0; 3]).unwrap(),
                Vec3::new(2.0, 2.0, 2.0),
            ),
            (
                cylinder([0.0; 3], 3.0, 6.0).unwrap(),
                Vec3::new(0.0, 0.0, 3.0),
            ),
            (cone([0.0; 3], 5.0, 12.0).unwrap(), Vec3::new(0.0, 0.0, 1.0)),
            (
                wedge([0.0; 3], 4.0, 3.0, 5.0).unwrap(),
                Vec3::new(1.0, 1.5, 1.0),
            ),
        ];
        for (solid, inside) in solids {
            let mesh = crate::brep::mesh::body(&solid, crate::tessellation::DEFAULT_ANGLE, 1e-9);
            assert!(!mesh.is_empty(), "nothing meshed");
            for triangle in &mesh.triangles {
                let corner = Vec3::from(mesh.positions[triangle[0]]);
                let normal = Vec3::from(mesh.normals[triangle[0]]);
                assert!(
                    normal.dot(corner - inside) > 0.0,
                    "a face pointed inwards at {corner:?}"
                );
            }
        }
    }

    #[test]
    fn a_pyramid_is_a_base_and_a_triangle_per_side() {
        for sides in [3, 4, 5, 8] {
            let solid = pyramid([0.0; 3], 4.0, 9.0, sides).expect("a pyramid");
            assert_eq!(solid.faces.len(), sides + 1, "{sides}");
            assert_eq!(solid.edges.len(), 2 * sides, "the rim and one rail each");
            assert_eq!(solid.vertices.len(), sides + 1);
            assert_eq!(solid.euler_characteristic(), 2, "{sides}");
            assert!(solid.validate().is_empty(), "{sides}");
        }
        assert!(
            pyramid([0.0; 3], 4.0, 9.0, 2).is_none(),
            "a sliver is not a solid"
        );
        assert!(pyramid([0.0; 3], 0.0, 9.0, 4).is_none());
        assert!(pyramid([0.0; 3], 4.0, -1.0, 4).is_none());
    }

    #[test]
    fn a_pyramid_encloses_a_third_of_the_prism_around_it() {
        // Base area times height over three, whatever the base is. A square
        // one is the case worth naming, since its base area is exact.
        let solid = pyramid([0.0; 3], 4.0, 9.0, 4).unwrap();
        // A square inscribed in a circle of radius four has diagonal eight.
        let expected = (0.5 * 8.0 * 8.0) * 9.0 / 3.0;
        let volume = meshed_volume(&solid);
        assert!(
            (volume - expected).abs() < 1e-9 * expected,
            "{volume} vs {expected}"
        );
    }

    #[test]
    fn a_sphere_meshes_from_what_bounds_it_rather_than_where() {
        // Its seam ends at the poles, where longitude has no single value ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â
        // every meridian passes through them ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â so where that seam sits in
        // (u, v) cannot be worked out from the geometry at all. What can be
        // worked out is that the seam is the *only* thing bounding the face,
        // and a face bounded by nothing but its own seams covers the whole
        // of its surface. That is enough, and it needs no pcurve.
        let solid = sphere([0.0; 3], 5.0).unwrap();
        let volume = meshed_volume(&solid);
        let expected = 4.0 / 3.0 * std::f64::consts::PI * 125.0;
        assert!(volume > 0.0, "wound inwards: {volume}");
        assert!(
            volume < expected,
            "a tessellation cannot read over: {volume}"
        );
        assert!(volume > 0.98 * expected, "{volume} vs {expected}");
    }

    #[test]
    fn a_sphere_with_no_size_is_refused() {
        assert!(sphere([0.0; 3], 0.0).is_none());
        assert!(sphere([0.0; 3], -1.0).is_none());
        assert!(cone([0.0; 3], 1.0, 0.0).is_none());
        assert!(wedge([0.0; 3], 1.0, 0.0, 1.0).is_none());
    }

    #[test]
    fn a_cylinder_with_no_size_is_refused() {
        assert!(cylinder([0.0; 3], 0.0, 1.0).is_none());
        assert!(cylinder([0.0; 3], 1.0, -1.0).is_none());
        assert!(cylinder([0.0; 3], f64::NAN, 1.0).is_none());
    }

    #[test]
    fn a_box_with_no_thickness_is_refused() {
        assert!(cuboid([0.0; 3], [1.0, 0.0, 1.0]).is_none());
        assert!(cuboid([0.0; 3], [1.0, -2.0, 1.0]).is_none());
        assert!(cuboid([0.0; 3], [1.0, f64::NAN, 1.0]).is_none());
    }

    #[test]
    fn a_box_at_survey_coordinates_is_still_a_box() {
        let body = cuboid([512_345.678, 4_512_345.678, 91.5], [0.5, 0.5, 0.5]).unwrap();
        assert!(body.validate().is_empty());
        assert_eq!(body.euler_characteristic(), 2);
        assert!(body.worst_vertex_gap() < 1e-9);
    }

    #[test]
    fn validation_notices_a_face_taken_away() {
        let mut body = unit_box();
        let victim = body.faces.keys().next().unwrap();
        body.faces.remove(victim);
        let flaws = body.validate();
        assert!(
            flaws
                .iter()
                .any(|flaw| matches!(flaw, Flaw::DanglingKey(_))),
            "{flaws:?}"
        );
        assert_ne!(body.euler_characteristic(), 2);
    }

    #[test]
    fn validation_notices_a_loop_that_no_longer_closes() {
        let mut body = unit_box();
        let ring = body.loops.keys().next().unwrap();
        body.loops.get_mut(ring).unwrap().coedges.swap(0, 1);
        let flaws = body.validate();
        assert!(
            flaws.iter().any(|flaw| matches!(flaw, Flaw::OpenLoop(_))),
            "{flaws:?}"
        );
    }

    #[test]
    fn validation_notices_a_face_wound_the_wrong_way() {
        let mut body = unit_box();
        let ring = body.loops.keys().next().unwrap();
        let coedges = body.loops.get(ring).unwrap().coedges.clone();
        for key in coedges {
            let coedge = body.coedges.get_mut(key).unwrap();
            coedge.forward = !coedge.forward;
        }
        let flaws = body.validate();
        assert!(
            flaws
                .iter()
                .any(|flaw| matches!(flaw, Flaw::SameSidedEdge(_))),
            "{flaws:?}"
        );
    }

    #[test]
    fn moving_a_vertex_dirties_what_names_it_and_nothing_else() {
        let mut body = unit_box();
        // Pretend it came from a file, so there is something to dirty.
        for node in body.vertices.values_mut() {
            node.provenance = Provenance::Clean(crate::brep::SourceRef::new(0));
        }
        for node in body.edges.values_mut() {
            node.provenance = Provenance::Clean(crate::brep::SourceRef::new(0));
        }
        for node in body.faces.values_mut() {
            node.provenance = Provenance::Clean(crate::brep::SourceRef::new(0));
        }
        let corner = body.vertices.keys().next().unwrap();
        body.soil_vertex(corner);

        // Three edges meet at a corner of a box, and three faces.
        let dirty_edges = body
            .edges
            .iter()
            .filter(|(_, e)| !e.provenance.is_reusable())
            .count();
        assert_eq!(dirty_edges, 3, "only the edges that end on it");
        let dirty_faces = body
            .faces
            .iter()
            .filter(|(_, f)| !f.provenance.is_reusable())
            .count();
        assert_eq!(dirty_faces, 3, "only the faces those bound");
    }

    /// The wall of `body`'s shell, which primitives keep as the one face
    /// that is not a plane.
    fn wall_surface(body: &Body) -> (&crate::brep::SurfaceKey, &Surface) {
        let face = body
            .faces
            .iter()
            .find(|(_, face)| !matches!(body.surfaces.get(face.surface), Some(Surface::Plane(_))))
            .expect("a non-planar wall face");
        (&face.1.surface, body.surfaces.get(face.1.surface).unwrap())
    }

    #[test]
    fn an_elliptical_cylinder_is_the_circular_census_with_a_ratio() {
        let circular = cylinder([0.0, 0.0, 0.0], 5.0, 10.0).expect("the circular twin");
        let body =
            elliptical_cylinder([0.0, 0.0, 0.0], 5.0, 3.0, 10.0).expect("an elliptical cylinder");
        assert!(body.validate().is_empty());
        assert_eq!(body.euler_characteristic(), 2);
        // The same census the circular cylinder carries, counted against it
        // rather than restated: two closed rims, one seam, the wiring both
        // walls take.
        assert_eq!(body.vertices.len(), circular.vertices.len());
        assert_eq!(body.edges.len(), circular.edges.len());
        assert_eq!(body.faces.len(), circular.faces.len());
        assert_eq!(body.loops.len(), circular.loops.len());
        assert_eq!(body.coedges.len(), circular.coedges.len());

        let (_, wall) = wall_surface(&body);
        match wall {
            Surface::EllipticCone(cone) => {
                assert_eq!(cone.radius, 5.0);
                assert!((cone.ratio - 0.6).abs() < 1e-12);
                assert_eq!(cone.half_angle, 0.0, "a cylinder does not taper");
                assert!((cone.base.point_at([5.0, 0.0])[0] - 5.0).abs() < 1e-12);
            }
            other => panic!("the wall left the analytic family: {other:?}"),
        }

        let ellipse_count = body
            .curves
            .iter()
            .filter(|(_, curve)| {
                matches!(
                    curve,
                    Curve3::Ellipse(ellipse)
                        if (ellipse.major_radius - 5.0).abs() < 1e-12
                            && (ellipse.minor_radius - 3.0).abs() < 1e-12
                )
            })
            .count();
        assert_eq!(ellipse_count, 2, "both rims stay analytic ellipses");
    }

    #[test]
    fn an_elliptical_cylinder_keeps_the_major_on_its_frame() {
        // Major on Y: the frame swaps rather than the record stretching,
        // which is the direction the SAT basis vector carries.
        let body = elliptical_cylinder([0.0, 0.0, 0.0], 3.0, 5.0, 10.0).unwrap();
        let (_, wall) = wall_surface(&body);
        match wall {
            Surface::EllipticCone(cone) => {
                assert_eq!(cone.radius, 5.0);
                assert!((cone.ratio - 0.6).abs() < 1e-12);
                // The frame's own x axis points along world Y.
                let [x, y, _] = cone.base.x_axis;
                assert!(y > x);
            }
            other => panic!("the wall left the analytic family: {other:?}"),
        }
    }

    #[test]
    fn an_elliptical_cone_tapers_like_its_circular_twin() {
        let circular = cone([0.0, 0.0, 0.0], 5.0, 10.0).unwrap();
        let elliptical = frustum([0.0, 0.0, 0.0], 5.0, 3.0, 0.0, 10.0).unwrap();
        assert!(elliptical.validate().is_empty());
        assert_eq!(
            (elliptical.vertices.len(), elliptical.edges.len(), elliptical.faces.len()),
            (circular.vertices.len(), circular.edges.len(), circular.faces.len()),
            "the elliptical cone carries the circular cone's census"
        );
        let (_, wall) = wall_surface(&elliptical);
        match wall {
            Surface::EllipticCone(cone) => {
                assert_eq!(cone.radius, 5.0);
                assert!((cone.ratio - 0.6).abs() < 1e-12);
                assert!(
                    (cone.half_angle - (5.0f64 / 10.0).atan()).abs() < 1e-12,
                    "the major falls to nothing over the height, as on a circular cone"
                );
            }
            other => panic!("the wall left the analytic family: {other:?}"),
        }
    }

    #[test]
    fn an_elliptical_frustum_keeps_its_end_sections_proportional() {
        let body = frustum([0.0, 0.0, 0.0], 5.0, 3.0, 2.0, 10.0).unwrap();
        assert!(body.validate().is_empty());
        assert_eq!(body.euler_characteristic(), 2);
        let rims: Vec<f64> = body
            .curves
            .iter()
            .filter_map(|(_, curve)| match curve {
                Curve3::Ellipse(ellipse) => Some(ellipse.major_radius),
                _ => None,
            })
            .collect();
        assert_eq!(rims.len(), 2, "two closed rims");
        let top = rims.iter().cloned().fold(f64::INFINITY, f64::min);
        let bottom = rims.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!((top - 2.0).abs() < 1e-12, "top major {top}");
        assert!((bottom - 5.0).abs() < 1e-12, "bottom major {bottom}");
        let (_, wall) = wall_surface(&body);
        match wall {
            Surface::EllipticCone(cone) => {
                assert!(
                    (cone.half_angle - (3.0f64 / 10.0).atan()).abs() < 1e-12,
                    "the major runs 5 to 2 over the height"
                );
                assert!((cone.ratio - 0.6).abs() < 1e-12);
            }
            other => panic!("the wall left the analytic family: {other:?}"),
        }
    }

    fn of_elliptic(built: Option<Body>) -> (Body, usize, usize, (f64, f64)) {
        let body = built.expect("a body");
        let expected = match wall_surface(&body).1 {
            Surface::EllipticCone(cone) => (cone.radius, cone.ratio),
            other => panic!("the wall left the analytic family: {other:?}"),
        };
        let faces = body.faces.len();
        let edges = body.edges.len();
        (body, faces, edges, expected)
    }

    #[test]
    fn elliptical_solids_round_trip_through_sat() {
        for (body, faces, edges, expected) in [
            of_elliptic(elliptical_cylinder([0.0, 0.0, 0.0], 5.0, 3.0, 10.0)),
            of_elliptic(frustum([0.0, 0.0, 0.0], 5.0, 3.0, 2.0, 10.0)),
            of_elliptic(frustum([0.0, 0.0, 0.0], 5.0, 3.0, 0.0, 10.0)),
        ] {
            let mut document = cadcodec::entities::acis::SatDocument::new();
            assert!(crate::acis::append(&body, &mut document).is_ok());
            let parsed =
                cadcodec::entities::acis::SatDocument::parse(&document.to_sat_string()).unwrap();
            let (restored, loss) = crate::acis::lift(&parsed);
            assert!(loss.is_empty(), "{loss:?}");
            assert_eq!(restored.len(), 1);
            assert!(restored[0].validate().is_empty());
            assert_eq!(restored[0].faces.len(), faces, "face census changed");
            assert_eq!(restored[0].edges.len(), edges, "edge census changed");
            let restored_pair = match wall_surface(&restored[0]).1 {
                Surface::EllipticCone(cone) => (cone.radius, cone.ratio),
                other => panic!("the lifted wall left the analytic family: {other:?}"),
            };
            assert!((restored_pair.0 - expected.0).abs() < 1e-12);
            assert!((restored_pair.1 - expected.1).abs() < 1e-12);
        }
    }

    #[test]
    fn the_emitted_record_is_the_authored_elliptical_wire_form() {
        // A native ACIS author writes an elliptical cylinder as one
        // cone-surface whose major basis carries the semi-axis, whose
        // fourth token is the section ratio, and whose half-angle tokens
        // read sine 0, cosine 1 ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â measured from the authored gold
        // fixtures this form was built against.
        let body = elliptical_cylinder([0.0, 0.0, 0.0], 5.0, 3.0, 10.0).unwrap();
        let mut document = cadcodec::entities::acis::SatDocument::new();
        crate::acis::append(&body, &mut document).unwrap();
        let parsed = cadcodec::entities::acis::SatDocument::parse(&document.to_sat_string()).unwrap();
        for record in &parsed.records {
            if let Some(cone) = cadcodec::entities::acis::SatConeSurface::from_record(record) {
                let (mx, my, mz) = cone.major_axis();
                assert!(
                    ((mx.hypot(my).hypot(mz)) - 5.0).abs() < 1e-12,
                    "the basis is the major semi-axis, not a unit vector"
                );
                assert!((cone.ratio() - 0.6).abs() < 1e-12);
                assert_eq!(cone.sin_half_angle().abs(), 0.0);
                assert_eq!(cone.cos_half_angle(), 1.0);
                assert_eq!(cone.radius(), 5.0);
                return;
            }
        }
        panic!("no cone-surface record was emitted for an elliptical cylinder");
    }
}
