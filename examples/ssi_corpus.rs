//! Runner for jyviko's SSI corpus (github.com/jyviko/ssi-corpus) through the
//! kernel's own surfaces and `brep::ssi::meet`.
//!
//! ```sh
//! cargo build --release --features brep --example ssi_corpus
//! SSI_ENGINE=opencadkernel SSI_RUNNER=target/release/examples/ssi_corpus \
//!     python3 -I <ssi-corpus>/scoring/score.py fs 1e-6
//! ```
//!
//! Each surface is built as the kernel's own `Surface` over the case's box,
//! and every curve point is written with each surface evaluated at the
//! parameters `meet` gave back, so the scorer checks those too.

#[path = "../crates/opencadkernel-ssi/src/bin/ssi_corpus_runner/json.rs"]
mod json;

use std::fmt::Write as _;
use std::time::Instant;

use json::Json;
use opencadkernel::brep::ssi::{meet, Bounded, Refusal, Shared};
use opencadkernel::brep::{Cone, Cylinder, Sphere, Surface, Torus};
use opencadkernel::space::{NurbsSurface3, Plane};

fn nums(j: Option<&Json>, what: &str) -> Result<Vec<f64>, String> {
    j.and_then(Json::nums).ok_or_else(|| format!("{what} is not a number array"))
}
fn v3(j: Option<&Json>, what: &str) -> Result<[f64; 3], String> {
    nums(j, what)?.try_into().map_err(|_| format!("{what} is not a 3-vector"))
}
fn pair(j: Option<&Json>, what: &str) -> Result<[f64; 2], String> {
    nums(j, what)?.try_into().map_err(|_| format!("{what} is not a pair"))
}
fn number(j: Option<&Json>, what: &str) -> Result<f64, String> {
    j.and_then(Json::num).ok_or_else(|| format!("{what} is not a number"))
}
fn unit(a: [f64; 3]) -> [f64; 3] {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    [a[0] / l, a[1] / l, a[2] / l]
}
fn radians(r: [f64; 2]) -> [f64; 2] {
    [r[0].to_radians(), r[1].to_radians()]
}

/// A case's surface and its box, in the kernel's parameters.
fn surface(s: &Json) -> Result<(Surface, [f64; 2], [f64; 2]), String> {
    let kind = s.get("type").and_then(Json::str).ok_or("a surface has no type")?;
    if kind == "nurbs" {
        let degree = pair(s.get("degree"), "degree")?;
        let (nu, nv) = (number(s.get("nu"), "nu")? as usize, number(s.get("nv"), "nv")? as usize);
        let points: Vec<[f64; 3]> = s
            .get("points")
            .and_then(Json::arr)
            .ok_or("points")?
            .iter()
            .map(|p| v3(Some(p), "a control point"))
            .collect::<Result<_, _>>()?;
        let weights = nums(s.get("weights"), "weights")?;
        let net = (0..nu).map(|i| points[i * nv..(i + 1) * nv].to_vec()).collect();
        let w = (0..nu).map(|i| weights[i * nv..(i + 1) * nv].to_vec()).collect();
        let spline = NurbsSurface3::new(
            degree[0] as usize,
            degree[1] as usize,
            net,
            nums(s.get("knots_u"), "knots_u")?,
            nums(s.get("knots_v"), "knots_v")?,
            Some(w),
        )
        .ok_or("a spline the kernel cannot take")?;
        let ((u0, u1), (v0, v1)) = spline.domain();
        return Ok((Surface::Nurbs(spline), [u0, u1], [v0, v1]));
    }
    let (x, y) = (unit(v3(s.get("xaxis"), "xaxis")?), unit(v3(s.get("yaxis"), "yaxis")?));
    let frame = |origin| Plane { origin, x_axis: x, y_axis: y };
    Ok(match kind {
        "plane" => {
            let plane = Plane {
                origin: v3(s.get("origin"), "origin")?,
                x_axis: v3(s.get("xaxis"), "xaxis")?,
                y_axis: v3(s.get("yaxis"), "yaxis")?,
            };
            (Surface::Plane(plane), pair(s.get("u"), "u")?, pair(s.get("v"), "v")?)
        }
        "cylinder" => (
            Surface::Cylinder(Cylinder {
                base: frame(v3(s.get("origin"), "origin")?),
                radius: number(s.get("radius"), "radius")?,
            }),
            radians(pair(s.get("angle"), "angle")?),
            pair(s.get("height"), "height")?,
        ),
        "cone" => (
            // Radius 0 at the apex, widening along the axis.
            Surface::Cone(Cone {
                base: frame(v3(s.get("apex"), "apex")?),
                radius: 0.0,
                half_angle: -number(s.get("half_angle_deg"), "half_angle_deg")?.to_radians(),
            }),
            [0.0, std::f64::consts::TAU],
            pair(s.get("height"), "height")?,
        ),
        "sphere" => (
            Surface::Sphere(Sphere {
                frame: frame(v3(s.get("center"), "center")?),
                radius: number(s.get("radius"), "radius")?,
            }),
            radians(pair(s.get("lon"), "lon")?),
            radians(pair(s.get("lat"), "lat")?),
        ),
        "torus" => (
            Surface::Torus(Torus {
                frame: frame(v3(s.get("center"), "center")?),
                major_radius: number(s.get("R"), "R")?,
                minor_radius: number(s.get("r"), "r")?,
            }),
            radians(pair(s.get("major"), "major")?),
            radians(pair(s.get("minor"), "minor")?),
        ),
        other => return Err(format!("unknown surface type: {other}")),
    })
}

fn run(args: &[String]) -> Result<(), String> {
    let [case_path, _mode, tol, out_path] = args else {
        return Err("usage: ssi_corpus case.json mode approx_tol out.json".into());
    };
    let text = std::fs::read_to_string(case_path).map_err(|e| format!("{case_path}: {e}"))?;
    let case = json::parse(&text)?;
    let [first, second] = case.get("surfaces").and_then(Json::arr).ok_or("no surfaces")? else {
        return Err("a case has two surfaces".into());
    };
    let (a, b) = (surface(first)?, surface(second)?);
    let feature: f64 = tol.parse().unwrap_or(1e-6);
    let (one, other) = (
        Bounded { surface: &a.0, u: a.1, v: a.2 },
        Bounded { surface: &b.0, u: b.1, v: b.2 },
    );
    let started = Instant::now();
    let result = meet(&one, &other, feature);
    let ms = started.elapsed().as_secs_f64() * 1e3;

    let mut out = String::from("{\"case\":");
    json::string(&mut out, case.get("id").and_then(Json::str).unwrap_or(""));
    let (status, name) = match &result {
        Ok(_) => (1000, "SUCCESS".to_string()),
        Err(Refusal::Tangent(why)) => (2000, format!("TANGENT: {why}")),
        Err(Refusal::Unsupported(why)) => (3000, format!("UNSUPPORTED: {why}")),
        Err(Refusal::Failed(why)) => (3000, format!("FAILED: {why}")),
    };
    let _ = write!(out, ",\"mode\":\"fs\",\"status\":{status},\"status_name\":");
    json::string(&mut out, &name);
    out.push_str(",\"time_ms\":");
    json::num(&mut out, ms);
    out.push_str(",\"curves\":[");
    if let Ok(found) = &result {
        for (index, branch) in found.branches.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            let mut points = branch.points.clone();
            if branch.closed {
                points.extend(branch.points.first().copied());
            }
            let _ = write!(out, "{{\"type\":\"CURVE\",\"closed\":{}", branch.closed);
            let lists: [(&str, &dyn Fn(&Shared) -> Vec<f64>); 5] = [
                ("xyz", &|q| q.point.to_vec()),
                ("uv1", &|q| q.first.to_vec()),
                ("s1_at_uv", &|q| a.0.point_at(q.first[0], q.first[1]).to_vec()),
                ("uv2", &|q| q.second.to_vec()),
                ("s2_at_uv", &|q| b.0.point_at(q.second[0], q.second[1]).to_vec()),
            ];
            for (key, value) in lists {
                let _ = write!(out, ",\"{key}\":[");
                for (k, q) in points.iter().enumerate() {
                    if k > 0 {
                        out.push(',');
                    }
                    json::point(&mut out, &value(q));
                }
                out.push(']');
            }
            out.push('}');
        }
    }
    out.push_str("],\"points\":[");
    if let Ok(found) = &result {
        let points = found.touches.iter().map(|q| q.point).chain(found.crossings.iter().copied());
        for (index, p) in points.enumerate() {
            if index > 0 {
                out.push(',');
            }
            json::point(&mut out, &p);
        }
    }
    out.push_str("]}\n");
    std::fs::write(out_path, out).map_err(|e| format!("{out_path}: {e}"))
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if let Err(message) = run(&args) {
        eprintln!("ssi_corpus: {message}");
        std::process::exit(2);
    }
}
