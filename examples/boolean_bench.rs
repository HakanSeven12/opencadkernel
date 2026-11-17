// Boolean timings: a 100-hole plate and a few curved pairs.
// cargo run --release --features brep --example boolean_bench [-- <filter>]
use opencadkernel::brep::{
    combine, combine_all, make, mass_properties, transform, Body, Operation, Placement,
};
use std::f64::consts::PI;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

/// Each repetition moves every input a little, so caches keyed by geometry
/// do not carry one run's work into the next.
static SHIFT: AtomicU64 = AtomicU64::new(0);

fn shifted(body: Body) -> Body {
    let step = SHIFT.load(Ordering::Relaxed) as f64 * 1.37e-3;
    transform(&body, &Placement::at([step, step * 0.7, step * 0.3])).expect("moved")
}

fn cb(a: Body, b: Body, how: Operation) -> Option<Body> {
    combine(shifted(a), shifted(b), how, 1e-6).ok()
}

fn turned(body: Body, x_axis: [f64; 3], y_axis: [f64; 3], z_axis: [f64; 3]) -> Body {
    let place = Placement {
        x_axis,
        y_axis,
        z_axis,
        origin: [0.0; 3],
    };
    transform(&body, &place).expect("placed")
}

fn plate(holes: usize) -> (Body, Vec<Body>, f64) {
    let side = 20.0 * (holes as f64).sqrt().ceil();
    let plate = shifted(make::cuboid([0.0, 0.0, 0.0], [side, side, 10.0]).unwrap());
    let per_row = (holes as f64).sqrt().ceil() as usize;
    let cutters: Vec<Body> = (0..holes)
        .map(|n| {
            let (i, j) = (n % per_row, n / per_row);
            let centre = [10.0 + 20.0 * i as f64, 10.0 + 20.0 * j as f64, -1.0];
            shifted(make::cylinder(centre, 4.0, 12.0).unwrap())
        })
        .collect();
    let volume = side * side * 10.0 - holes as f64 * PI * 16.0 * 10.0;
    (plate, cutters, volume)
}

/// Volume two perpendicular cylinders share, radii `big` >= `small`.
fn crossing(big: f64, small: f64) -> f64 {
    let steps = 20000;
    let width = 2.0 * small / steps as f64;
    (0..steps)
        .map(|i| {
            let y = -small + (i as f64 + 0.5) * width;
            4.0 * (small * small - y * y).sqrt() * (big * big - y * y).sqrt() * width
        })
        .sum()
}

fn run(name: &str, truth: f64, job: impl FnOnce() -> Option<Body>) {
    let start = Instant::now();
    let result = job();
    let took = start.elapsed().as_secs_f64() * 1e3;
    match result.as_ref().and_then(mass_properties) {
        Some(mass) => {
            let error = (mass.volume - truth).abs() / truth.abs().max(1e-12);
            println!("{name:<28} {took:>9.1} ms  volume {:.6} (rel err {error:.1e})", mass.volume);
        }
        None => println!("{name:<28} {took:>9.1} ms  FAILED ({})", result.is_some()),
    }
}

#[cfg(target_os = "linux")]
extern "C" {
    fn prctl(option: i32, a2: u64, a3: u64, a4: u64, a5: u64) -> i32;
}

fn main() {
    // Lets a sampling profiler (eu-stack) attach: PR_SET_PTRACER, any.
    #[cfg(target_os = "linux")]
    if std::env::var_os("BENCH_PROFILE").is_some() {
        unsafe { prctl(0x59616d61, u64::MAX, 0, 0, 0) };
    }
    let repeat: usize = std::env::var("BENCH_REPEAT").ok().and_then(|n| n.parse().ok()).unwrap_or(1);
    let filter = std::env::args().nth(1).unwrap_or_default();
    let tolerance = 1e-6;
    let cases = || -> Vec<(&str, f64, Box<dyn FnOnce() -> Option<Body>>)> { vec![
        ("box-cylinder", 8.0 - PI * 0.25 * 2.0, Box::new(move || {
            let a = make::cuboid([-1.0, -1.0, -1.0], [2.0, 2.0, 2.0])?;
            let b = make::cylinder([0.0, 0.0, -2.0], 0.5, 4.0)?;
            cb(a, b, Operation::Difference)
        })),
        ("sphere+box", {
            // sphere r=1 at origin united with the box over z >= 0.5 of it.
            let cap = PI * 0.5 * 0.5 * (3.0 - 0.5) / 3.0;
            4.0 / 3.0 * PI - cap + 4.0 * 1.5
        }, Box::new(move || {
            let a = make::sphere([0.0, 0.0, 0.0], 1.0)?;
            let b = make::cuboid([-1.0, -1.0, 0.5], [2.0, 2.0, 1.5])?;
            cb(a, b, Operation::Union)
        })),
        ("steinmetz", 16.0 / 3.0, Box::new(move || {
            let a = make::cylinder([0.0, 0.0, -2.0], 1.0, 4.0)?;
            let b = turned(make::cylinder([0.0, 0.0, -2.0], 1.0, 4.0)?, [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
            cb(a, b, Operation::Intersection)
        })),
        ("tpipe", PI * 0.25 * 6.0 + PI * 0.16 * 4.0 - crossing(0.5, 0.4), Box::new(move || {
            // r=0.5 along z (length 6) united with r=0.4 along x (length 4).
            let a = make::cylinder([0.0, 0.0, -3.0], 0.5, 6.0)?;
            let b = turned(make::cylinder([0.0, 0.0, -2.0], 0.4, 4.0)?, [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
            cb(a, b, Operation::Union)
        })),
        ("torus-box", 2.0 * PI * PI * 2.0 * 0.25 * 0.5, Box::new(move || {
            // half the torus (R=2, r=0.5) removed by a box over y >= 0.
            let a = make::torus([0.0, 0.0, 0.0], 2.0, 0.5)?;
            let b = make::cuboid([-3.0, 0.0, -1.0], [6.0, 3.0, 2.0])?;
            cb(a, b, Operation::Difference)
        })),
        ("sphere-sphere", {
            // two unit spheres a unit apart, intersected: a lens.
            let h = 0.5;
            2.0 * PI * h * h * (3.0 - h) / 3.0
        }, Box::new(move || {
            let a = make::sphere([0.0, 0.0, 0.0], 1.0)?;
            let b = make::sphere([1.0, 0.0, 0.0], 1.0)?;
            cb(a, b, Operation::Intersection)
        })),
        ("cone-cylinder", PI * 4.0 * 2.0 / 3.0 - PI * 0.25 * (1.5 + 1.0 / 6.0), Box::new(move || {
            let a = make::cone([0.0, 0.0, 0.0], 2.0, 2.0)?;
            let b = make::cylinder([0.0, 0.0, -1.0], 0.5, 4.0)?;
            cb(a, b, Operation::Difference)
        })),
        ("viviani", {
            // unit sphere minus a cylinder r=0.5 through its centre, off axis.
            4.0 / 3.0 * PI - (2.0 * PI / 3.0 - 8.0 / 9.0)
        }, Box::new(move || {
            let a = make::sphere([0.0, 0.0, 0.0], 1.0)?;
            let b = make::cylinder([0.5, 0.0, -2.0], 0.5, 4.0)?;
            cb(a, b, Operation::Difference)
        })),
        ("touch-box-edge", 16.0, Box::new(move || {
            let a = make::cuboid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0])?;
            let b = make::cuboid([2.0, 2.0, 0.0], [2.0, 2.0, 2.0])?;
            cb(a, b, Operation::Union)
        })),
        ("touch-sphere-box", 8.0 + 4.0 / 3.0 * PI, Box::new(move || {
            let a = make::cuboid([-1.0, -1.0, -2.0], [2.0, 2.0, 2.0])?;
            let b = make::sphere([0.0, 0.0, 1.0], 1.0)?;
            cb(a, b, Operation::Union)
        })),
        ("sunk-sphere-box", 8.0 + 4.0 / 3.0 * PI, Box::new(move || {
            let a = make::cuboid([-1.0, -1.0, -2.0], [2.0, 2.0, 2.0])?;
            let b = make::sphere([0.0, 0.0, 1.0 - 2e-7], 1.0)?;
            cb(a, b, Operation::Union)
        })),
        ("sunk-sphere-box-diff", 8.0, Box::new(move || {
            let a = make::cuboid([-1.0, -1.0, -2.0], [2.0, 2.0, 2.0])?;
            let b = make::sphere([0.0, 0.0, 1.0 - 2e-7], 1.0)?;
            cb(a, b, Operation::Difference)
        })),
        ("touch-cyl-box", 8.0 + PI * 2.0, Box::new(move || {
            let a = make::cuboid([-1.0, -1.0, -2.0], [2.0, 2.0, 2.0])?;
            let b = turned(make::cylinder([0.0, 1.0, -1.0], 1.0, 2.0)?, [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
            cb(a, b, Operation::Union)
        })),
        ("touch-cyl-cyl", 4.0 * PI, Box::new(move || {
            let a = make::cylinder([0.0, 0.0, 0.0], 1.0, 2.0)?;
            let b = make::cylinder([2.0, 0.0, 0.0], 1.0, 2.0)?;
            cb(a, b, Operation::Union)
        })),
        ("touch-sphere-sphere", 8.0 / 3.0 * PI, Box::new(move || {
            let a = make::sphere([0.0, 0.0, 0.0], 1.0)?;
            let b = make::sphere([2.0, 0.0, 0.0], 1.0)?;
            cb(a, b, Operation::Union)
        })),
        ("plate-25-seq", plate(25).2, Box::new(move || {
            let (mut body, cutters, _) = plate(25);
            for cutter in cutters {
                body = combine(body, cutter, Operation::Difference, tolerance).ok()?;
            }
            Some(body)
        })),
        ("plate-100-seq", plate(100).2, Box::new(move || {
            let (mut body, cutters, _) = plate(100);
            for cutter in cutters {
                body = combine(body, cutter, Operation::Difference, tolerance).ok()?;
            }
            Some(body)
        })),
        ("plate-100-batch", plate(100).2, Box::new(move || {
            let (body, cutters, _) = plate(100);
            combine_all(body, cutters, Operation::Difference, tolerance).ok()
        })),
    ] };
    for round in 0..repeat {
        SHIFT.store(round as u64, Ordering::Relaxed);
        for (name, truth, job) in cases() {
            if name.contains(&filter) {
                run(name, truth, job);
            }
        }
    }
}
