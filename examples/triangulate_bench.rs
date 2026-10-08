// Synthetic triangulation workloads, including mostly-reflex and many-hole rings.
// Run on both revisions: cargo run --release --example triangulate_bench
// Optional substring selects a workload: -- three-convex
use std::{
    hash::{DefaultHasher, Hash, Hasher},
    hint::black_box,
    time::Instant,
};
fn ellipse(n: usize, y: f64) -> Vec<[f64; 2]> {
    (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            [10.0 * a.cos(), y * a.sin()]
        })
        .collect()
}
fn main() {
    let mut cases = Vec::new();
    for n in [4, 8, 16, 17, 32, 128, 512, 1024, 2048] {
        cases.push((format!("ellipse-{n}"), vec![ellipse(n, 3.0)]));
    }
    for n in [16, 128, 512] {
        let mut star = ellipse(n, 10.0);
        for (i, p) in star.iter_mut().enumerate() {
            if i % 2 == 0 {
                p[0] *= 0.1;
                p[1] *= 0.1;
            }
        }
        cases.push((format!("star-{n}"), vec![star]));
        cases.push((format!("thin-{n}"), vec![ellipse(n, 0.00001)]));
        let mut comb = vec![[0.0, 0.0], [n as f64, 0.0], [n as f64, 10.0]];
        for i in (0..n).rev() {
            comb.extend([
                [i as f64 + 0.8, 10.0],
                [i as f64 + 0.8, 1.0],
                [i as f64 + 0.2, 1.0],
                [i as f64 + 0.2, 10.0],
            ]);
        }
        comb.push([0.0, 10.0]);
        cases.push((format!("comb-{n}"), vec![comb]));
    }
    for n in [1, 10, 50] {
        let mut rings = vec![[[0., 0.], [100., 0.], [100., 100.], [0., 100.]].to_vec()];
        for i in 0..n {
            let x = 3.0 + (i % 10) as f64 * 9.0;
            let y = 3.0 + (i / 10) as f64 * 9.0;
            rings.push(
                ellipse(12, 0.1)
                    .iter()
                    .map(|p| [x + p[0] * 0.01, y + p[1]])
                    .collect(),
            );
        }
        cases.push((format!("holes-{n}"), rings));
    }
    cases.push((
        "self-crossing".into(),
        vec![vec![[0., 0.], [1., 1.], [0., 1.], [1., 0.]]],
    ));
    for steps in [4, 16, 64, 256] {
        let tips = [[-10.0, -5.0], [10.0, -5.0], [0.0, 10.0]];
        let mut ring = Vec::new();
        for side in 0..3 {
            let a = tips[side];
            let b = tips[(side + 1) % 3];
            let control = [(a[0] + b[0]) * 0.1, (a[1] + b[1]) * 0.1];
            for i in 0..steps {
                let t = i as f64 / steps as f64;
                let s = 1.0 - t;
                ring.push([
                    s * s * a[0] + 2.0 * s * t * control[0] + t * t * b[0],
                    s * s * a[1] + 2.0 * s * t * control[1] + t * t * b[1],
                ]);
            }
        }
        cases.push((format!("three-convex-{steps}"), vec![ring]));
    }
    cases.push((
        "collinear".into(),
        vec![vec![[0., 0.], [1., 0.], [2., 0.], [3., 0.]]],
    ));
    for (name, rings) in cases {
        if std::env::args()
            .nth(1)
            .is_some_and(|filter| !name.contains(&filter))
        {
            continue;
        }
        let (vertices, triangles) = opencadkernel::geom2d::triangulate_rings(&rings);
        // Compare output hashes across revisions using the same Rust toolchain.
        let mut mesh_hash = DefaultHasher::new();
        for point in &vertices {
            point.map(f64::to_bits).hash(&mut mesh_hash);
        }
        triangles.hash(&mut mesh_hash);
        let n: usize = rings.iter().map(Vec::len).sum();
        let repeats = if n <= 32 {
            5000
        } else if n <= 150 {
            100
        } else if n <= 700 {
            5
        } else {
            1
        };
        let mut times = Vec::new();
        for _ in 0..5 {
            let start = Instant::now();
            for _ in 0..repeats {
                black_box(opencadkernel::geom2d::triangulate_rings(black_box(&rings)));
            }
            times.push(start.elapsed().as_nanos() as f64 / repeats as f64);
        }
        times.sort_by(f64::total_cmp);
        println!(
            "{name:20} vertices={n:5} triangles={:5} mesh={:016x} median_us={:.3}",
            triangles.len(),
            mesh_hash.finish(),
            times[2] / 1000.0
        );
    }
}
