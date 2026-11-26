//! Moving an imported solid's surfaces, as little as possible, where its faces miss each other.
//!
//! Where two faces touch rather than cross, or a B-spline face stops short of its corner, no
//! movement of vertices and edges alone brings them within [`tol::LINEAR`] of every face they
//! bound. Each chosen face's surface gets a few degrees of freedom: a rigid motion and its radii,
//! or its control points. The vertices and sample points of every edge on those faces are free
//! too, and weighted least-change Gauss–Newton steps move them all until every point lies on each
//! of its faces. Moving a surface costs more than moving a point, so surfaces move only where the
//! points cannot close a gap on their own.

use std::collections::HashMap;

use super::topology::Raw;
use crate::{
    geom::{NurbsSurface, Surface},
    math::{Point3, Transform, Vec2, Vec3},
    numeric::solve_dense,
    tol,
};

/// What moving a surface by a length costs, against moving a point by the same length (1).
const SURFACE_COST: f64 = 100.0;
const ROUNDS: usize = 12;
/// The steps tried in turn until one shrinks the worst gap: the damping added to the step's
/// normal equations, and how much of the step is taken. Touching faces make those equations
/// nearly singular, and an undamped step then runs far along a direction they barely fix.
const TRIALS: [(f64, f64); 8] = [
    (1e-12, 1.0),
    (1e-12, 0.25),
    (1e-6, 1.0),
    (1e-6, 0.25),
    (1e-3, 1.0),
    (1e-3, 0.25),
    (1e-12, 0.0625),
    (1e-12, 0.015625),
];
/// The conjugate-gradient iterations without the residual halving after which a solve stops. A
/// stall on these nearly singular systems is lasting: past this many, more iterations only cost
/// time.
const CG_STALL: usize = 25;
/// The fewest knot spans a moved B-spline surface has along each parameter.
const MIN_SPANS: usize = 16;

/// A point the step moves, and the faces it must lie on.
pub(super) struct Target {
    pub at: Point3,
    pub faces: Vec<usize>,
    /// The vertex it is, if it is one rather than a sample of an edge.
    pub vertex: Option<usize>,
    /// The edge it samples.
    pub edge: Option<usize>,
    /// The directions its three unknowns move it in, and what moving it each way costs.
    frame: [Vec3; 3],
    cost: [f64; 3],
}

/// Below this sine of the angle between two faces' normals, the faces touch rather than cross.
const TOUCHING: f64 = 0.1;

impl Target {
    /// A vertex. Where all its faces touch (a corner of rounds that blend into each other), only
    /// how far along their common normal it lies is fixed to first order; sliding across that
    /// normal costs as much as moving a surface.
    pub(super) fn vertex(raw: &Raw, at: Point3, faces: Vec<usize>, vertex: usize) -> Target {
        let normals = normals_at(raw, &faces, at);
        let n = common_normal(&normals);
        let touching = normals
            .iter()
            .all(|m| n.is_some_and(|n| n.cross(*m).norm() < TOUCHING));
        let (frame, cost) = match n.filter(|_| touching) {
            Some(n) => {
                let a = n.any_perpendicular().unit_or_zero();
                ([n, a, n.cross(a)], [1.0, SURFACE_COST, SURFACE_COST])
            }
            None => ([Vec3::X, Vec3::Y, Vec3::Z], [1.0; 3]),
        };
        Target {
            at,
            faces,
            vertex: Some(vertex),
            edge: None,
            frame,
            cost,
        }
    }

    /// A point of edge `edge`, which runs along `along` there. Where along its edge a sample lies
    /// is arbitrary, so moving it that way costs as much as moving a surface: free, the step's
    /// linear model slides it far along an edge between faces that touch.
    pub(super) fn sample(
        raw: &Raw,
        at: Point3,
        faces: Vec<usize>,
        edge: usize,
        along: Vec3,
    ) -> Target {
        let (frame, cost) = match along.unit() {
            Some(a) => {
                let n = common_normal(&normals_at(raw, &faces, at))
                    .map(|n| n - a * a.dot(n))
                    .and_then(Vec3::unit)
                    .unwrap_or(a.any_perpendicular().unit_or_zero());
                ([a, a.cross(n), n], [SURFACE_COST, 1.0, 1.0])
            }
            None => ([Vec3::X, Vec3::Y, Vec3::Z], [1.0; 3]),
        };
        Target {
            at,
            faces,
            vertex: None,
            edge: Some(edge),
            frame,
            cost,
        }
    }

    fn moved_by(&self, x: &[f64]) -> Vec3 {
        let f = &self.frame;
        f[0] * x[0] + f[1] * x[1] + f[2] * x[2]
    }
}

/// The normal of each of `faces`' surfaces at the point nearest `at`.
fn normals_at(raw: &Raw, faces: &[usize], at: Point3) -> Vec<Vec3> {
    faces
        .iter()
        .map(|&f| {
            let s = &raw.faces[f].surface;
            let uv = s.project(at);
            s.normal(uv.x, uv.y)
        })
        .collect()
}

/// The mean of `normals`, each turned to agree with the first.
fn common_normal(normals: &[Vec3]) -> Option<Vec3> {
    normals
        .iter()
        .fold(Vec3::ZERO, |acc, &n| {
            acc + if acc.dot(n) < 0.0 { -n } else { n }
        })
        .unit()
}

/// One face's unknowns, from column `first`.
enum Dofs {
    /// Translation, rotation about `center` (times `lever`, so every unknown is a length), and
    /// one unknown per radius.
    Rigid {
        center: Point3,
        lever: f64,
        sizes: usize,
    },
    /// Three per group of control points that coincide (a pole, or a closed surface's seam), so
    /// moving them keeps the surface closed.
    Control { group_of: Vec<usize>, groups: usize },
}

struct Vars {
    first: usize,
    dofs: Dofs,
}

impl Vars {
    fn count(&self) -> usize {
        match &self.dofs {
            Dofs::Rigid { sizes, .. } => 6 + sizes,
            Dofs::Control { groups, .. } => 3 * groups,
        }
    }
}

/// Directions in which the point of `s` at `uv`, `q`, moves when each radius grows.
fn size_dirs(s: &Surface, uv: Vec2, q: Point3) -> Vec<Vec3> {
    match s {
        Surface::Cylinder(c) => vec![c.frame.radial(uv.x)],
        Surface::Cone(c) => vec![c.frame.radial(uv.x)],
        Surface::Sphere(c) => vec![(q - c.frame.origin).unit_or_zero()],
        Surface::Torus(t) => {
            let (sv, cv) = uv.y.sin_cos();
            vec![
                t.frame.radial(uv.x),
                t.frame.radial(uv.x) * cv + t.frame.z * sv,
            ]
        }
        _ => Vec::new(),
    }
}

fn with_sizes(s: &Surface, d: &[f64]) -> Option<Surface> {
    let positive = |r: f64| (r > tol::LINEAR).then_some(r);
    Some(match s {
        Surface::Cylinder(c) => {
            let mut c = c.clone();
            c.radius = positive(c.radius + d[0])?;
            Surface::Cylinder(c)
        }
        Surface::Cone(c) => {
            let mut c = c.clone();
            c.radius += d[0];
            Surface::Cone(c)
        }
        Surface::Sphere(c) => {
            let mut c = c.clone();
            c.radius = positive(c.radius + d[0])?;
            Surface::Sphere(c)
        }
        Surface::Torus(t) => {
            let mut t = t.clone();
            t.major += d[0];
            t.minor = positive(t.minor + d[1])?;
            Surface::Torus(t)
        }
        other => other.clone(),
    })
}

/// The parameter ranges a surface ends at, where it is bounded and not periodic.
fn bounds(s: &Surface) -> [Option<(f64, f64)>; 2] {
    let curve_bound = |c: &crate::geom::Curve| {
        let d = c.domain();
        (c.period().is_none() && d.0.is_finite() && d.1.is_finite()).then_some(d)
    };
    match s {
        Surface::Nurbs(n) => {
            let (du, dv) = n.domain();
            [
                n.u_period().is_none().then_some(du),
                n.v_period().is_none().then_some(dv),
            ]
        }
        Surface::Extrusion(e) => [curve_bound(&e.curve), None],
        Surface::Revolution(r) => [None, curve_bound(&r.curve)],
        _ => [None, None],
    }
}

/// The directions in which a point's distance from `s` at `uv` counts: the normal, and along any
/// parameter clamped at the end of the surface the tangent across that end.
fn directions(s: &Surface, uv: Vec2) -> Vec<Vec3> {
    let (_, su, sv) = s.d1(uv.x, uv.y);
    let n = s.normal(uv.x, uv.y);
    let b = bounds(s);
    let at_end = |b: Option<(f64, f64)>, t: f64| {
        b.is_some_and(|(lo, hi)| {
            let eps = 1e-9 * (hi - lo);
            t <= lo + eps || t >= hi - eps
        })
    };
    let (cu, cv) = (at_end(b[0], uv.x), at_end(b[1], uv.y));
    let flat = |v: Vec3| (v - n * n.dot(v)).unit();
    match (cu, cv, n.unit()) {
        (false, false, Some(n)) => vec![n],
        (true, false, Some(n)) => match (flat(su), flat(sv)) {
            (Some(a), Some(b)) => match (a - b * b.dot(a)).unit() {
                Some(t) => vec![n, t],
                None => vec![Vec3::X, Vec3::Y, Vec3::Z],
            },
            _ => vec![Vec3::X, Vec3::Y, Vec3::Z],
        },
        (false, true, Some(n)) => match (flat(sv), flat(su)) {
            (Some(a), Some(b)) => match (a - b * b.dot(a)).unit() {
                Some(t) => vec![n, t],
                None => vec![Vec3::X, Vec3::Y, Vec3::Z],
            },
            _ => vec![Vec3::X, Vec3::Y, Vec3::Z],
        },
        _ => vec![Vec3::X, Vec3::Y, Vec3::Z],
    }
}

/// Knot `t` inserted once into the curve of homogeneous points `pts` over `knots` (Boehm).
fn insert_knot(
    knots: &[f64],
    degree: usize,
    pts: &[[f64; 4]],
    t: f64,
) -> (Vec<f64>, Vec<[f64; 4]>) {
    let k = (degree..pts.len())
        .rev()
        .find(|&i| knots[i] <= t)
        .unwrap_or(degree);
    let mut q = Vec::with_capacity(pts.len() + 1);
    for i in 0..=pts.len() {
        if i + degree <= k {
            q.push(pts[i]);
        } else if i > k {
            q.push(pts[i - 1]);
        } else {
            let den = knots[i + degree] - knots[i];
            let a = if den > 0.0 { (t - knots[i]) / den } else { 0.0 };
            let (p0, p1) = (pts[i - 1], pts[i]);
            q.push(std::array::from_fn(|c| (1.0 - a) * p0[c] + a * p1[c]));
        }
    }
    let mut kn = knots.to_vec();
    kn.insert(k + 1, t);
    (kn, q)
}

/// The knots that split the spans of `knots` over `domain` until there are at least `spans`, each
/// at the middle of the longest span left.
fn splits(knots: &[f64], domain: (f64, f64), spans: usize) -> Vec<f64> {
    let mut cuts: Vec<f64> = knots
        .iter()
        .copied()
        .filter(|&k| k >= domain.0 && k <= domain.1)
        .collect();
    cuts.dedup();
    let mut out = Vec::new();
    while cuts.len() - 1 < spans {
        let i = (0..cuts.len() - 1)
            .max_by(|&a, &b| (cuts[a + 1] - cuts[a]).total_cmp(&(cuts[b + 1] - cuts[b])))
            .unwrap_or(0);
        let m = 0.5 * (cuts[i] + cuts[i + 1]);
        out.push(m);
        cuts.insert(i + 1, m);
    }
    out
}

/// The same surface with knots inserted until each parameter has at least `spans` spans, so the
/// step can bend it along an edge as closely as the file's own curves bend.
fn refined(n: &NurbsSurface, spans: usize) -> Option<NurbsSurface> {
    let (du, dv) = n.degrees();
    let (mut ku, mut kv) = (n.knots().0.to_vec(), n.knots().1.to_vec());
    let (mut nu, mut nv) = n.counts();
    let w = n.weights();
    let mut net: Vec<[f64; 4]> = n
        .control_points()
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let wi = w.map_or(1.0, |w| w[i]);
            [p.x * wi, p.y * wi, p.z * wi, wi]
        })
        .collect();
    let (dom_u, dom_v) = n.domain();
    for t in splits(&ku, dom_u, spans) {
        let mut next = vec![[0.0; 4]; (nu + 1) * nv];
        let mut knots = Vec::new();
        for j in 0..nv {
            let col: Vec<[f64; 4]> = (0..nu).map(|i| net[i * nv + j]).collect();
            let (k, q) = insert_knot(&ku, du, &col, t);
            for (i, p) in q.into_iter().enumerate() {
                next[i * nv + j] = p;
            }
            knots = k;
        }
        ku = knots;
        nu += 1;
        net = next;
    }
    for t in splits(&kv, dom_v, spans) {
        let mut next = Vec::with_capacity(nu * (nv + 1));
        let mut knots = Vec::new();
        for i in 0..nu {
            let (k, q) = insert_knot(&kv, dv, &net[i * nv..(i + 1) * nv], t);
            next.extend(q);
            knots = k;
        }
        kv = knots;
        nv += 1;
        net = next;
    }
    let ctrl = net
        .iter()
        .map(|p| Point3::new(p[0] / p[3], p[1] / p[3], p[2] / p[3]))
        .collect();
    let weights = w.map(|_| net.iter().map(|p| p[3]).collect());
    NurbsSurface::new((du, dv), (ku, kv), (nu, nv), ctrl, weights).ok()
}

/// Groups of control points that coincide, as each point's group and the number of groups.
fn control_groups(ctrl: &[Point3]) -> (Vec<usize>, usize) {
    let mut group_of = vec![usize::MAX; ctrl.len()];
    let mut groups = 0;
    for i in 0..ctrl.len() {
        if group_of[i] != usize::MAX {
            continue;
        }
        group_of[i] = groups;
        for j in i + 1..ctrl.len() {
            if group_of[j] == usize::MAX && ctrl[i].dist(ctrl[j]) <= tol::LINEAR {
                group_of[j] = groups;
            }
        }
        groups += 1;
    }
    (group_of, groups)
}

/// The rows of a sparse matrix, stored one after another: row `i` is
/// `entries[start[i]..start[i + 1]]`, each entry a column and its value. Conjugate gradients read
/// every row several times an iteration, and rows laid end to end are read in one sweep of memory.
struct Rows {
    start: Vec<usize>,
    entries: Vec<(usize, f64)>,
}

impl Rows {
    fn of(rows: &[Vec<(usize, f64)>]) -> Rows {
        let mut start = Vec::with_capacity(rows.len() + 1);
        let mut entries = Vec::with_capacity(rows.iter().map(Vec::len).sum());
        start.push(0);
        for row in rows {
            entries.extend_from_slice(row);
            start.push(entries.len());
        }
        Rows { start, entries }
    }

    fn len(&self) -> usize {
        self.start.len() - 1
    }

    fn row(&self, i: usize) -> &[(usize, f64)] {
        &self.entries[self.start[i]..self.start[i + 1]]
    }
}

/// A sparse system `A x = b`, solved for the `x` of least weighted norm `Σ wⱼ xⱼ²`.
struct System {
    rows: Rows,
    rhs: Vec<f64>,
    inv_w: Vec<f64>,
    /// Consecutive rows that share a point's unknowns, as `start..end`.
    blocks: Vec<(usize, usize)>,
}

impl System {
    /// The inverse of each block of `A W⁻¹ Aᵀ + λ`, row-major. A point on two faces that touch
    /// has two nearly equal rows, which only their surfaces' columns tell apart; inverting the
    /// block keeps that from stalling conjugate gradients.
    fn block_inverses(&self, lambda: f64) -> Vec<Vec<f64>> {
        let mut dense = vec![0.0; self.inv_w.len()];
        self.blocks
            .iter()
            .map(|&(a, b)| {
                let k = b - a;
                let mut m = vec![0.0; k * k];
                for i in 0..k {
                    for &(j, x) in self.rows.row(a + i) {
                        dense[j] += x;
                    }
                    for l in 0..k {
                        m[i * k + l] = (self.rows.row(a + l).iter())
                            .map(|&(j, x)| x * dense[j] * self.inv_w[j])
                            .sum::<f64>();
                    }
                    m[i * k + i] += lambda;
                    for &(j, _) in self.rows.row(a + i) {
                        dense[j] = 0.0;
                    }
                }
                let mut eye = vec![0.0; k * k];
                for i in 0..k {
                    eye[i * k + i] = 1.0;
                }
                solve_dense(m.clone(), eye, k, k).unwrap_or_else(|_| {
                    let mut d = vec![0.0; k * k];
                    for i in 0..k {
                        d[i * k + i] = 1.0 / m[i * k + i].max(lambda);
                    }
                    d
                })
            })
            .collect()
    }

    /// `z = W⁻¹ Aᵀ y`.
    fn at(&self, y: &[f64], z: &mut [f64]) {
        z.fill(0.0);
        for (i, &yi) in y.iter().enumerate() {
            for &(j, a) in self.rows.row(i) {
                z[j] += a * yi;
            }
        }
        for (zj, w) in z.iter_mut().zip(&self.inv_w) {
            *zj *= w;
        }
    }

    /// Block-preconditioned conjugate gradients on `(A W⁻¹ Aᵀ + λ) y = b`, then
    /// `x = W⁻¹ Aᵀ y`.
    fn solve(&self, lambda: f64) -> Vec<f64> {
        let m = self.rows.len();
        let n = self.inv_w.len();
        // One buffer for `W⁻¹ Aᵀ p` and one for the product, reused by every iteration.
        let mut z_n = vec![0.0; n];
        let apply = |y: &[f64], z: &mut [f64], out: &mut [f64]| {
            self.at(y, z);
            for (i, (o, &yi)) in out.iter_mut().zip(y).enumerate() {
                *o = self.rows.row(i).iter().map(|&(j, a)| a * z[j]).sum::<f64>() + lambda * yi;
            }
        };
        let inverses = self.block_inverses(lambda);
        let precondition = |r: &[f64], z: &mut [f64]| {
            z.fill(0.0);
            for (&(a, b), inv) in self.blocks.iter().zip(&inverses) {
                let k = b - a;
                for i in 0..k {
                    z[a + i] = (0..k).map(|l| inv[i * k + l] * r[a + l]).sum();
                }
            }
        };
        let mut y = vec![0.0; m];
        let mut r = self.rhs.clone();
        let mut z = vec![0.0; m];
        precondition(&r, &mut z);
        let mut p = z.clone();
        let mut ap = vec![0.0; m];
        let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        let goal = 1e-4 * tol::LINEAR;
        // On nearly singular systems (faces that touch) the residual stalls far above the goal and
        // stays there: once it has not halved for a while, more iterations only cost time.
        let (mut mark, mut since) = (f64::INFINITY, 0usize);
        for _ in 0..(4 * m).clamp(50, 4000) {
            let worst = r.iter().fold(0.0, |a: f64, x| a.max(x.abs()));
            if worst <= goal {
                break;
            }
            if worst <= 0.5 * mark {
                (mark, since) = (worst, 0);
            } else {
                since += 1;
                if since >= CG_STALL {
                    break;
                }
            }
            apply(&p, &mut z_n, &mut ap);
            let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
            if pap <= 0.0 || !pap.is_finite() {
                break;
            }
            let alpha = rz / pap;
            for i in 0..m {
                y[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            precondition(&r, &mut z);
            let rz_next: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
            let beta = rz_next / rz;
            rz = rz_next;
            for i in 0..m {
                p[i] = z[i] + beta * p[i];
            }
        }
        self.at(&y, &mut z_n);
        z_n
    }
}

/// Moves the surfaces of `faces` and every target until each target lies on each of its faces,
/// or until the rounds run out. Returns where each target vertex went, and whether the targets
/// all lie on their faces within a small part of the tolerance.
pub(super) fn adjust(
    raw: &mut Raw,
    faces: &[usize],
    targets: &mut [Target],
) -> (HashMap<usize, Point3>, bool) {
    let surfaces: Vec<Surface> = faces
        .iter()
        .map(|&f| raw.faces[f].surface.clone())
        .collect();
    let points: Vec<Point3> = targets.iter().map(|t| t.at).collect();
    let mut met = gauss_newton(raw, faces, targets);
    // A B-spline face's own control points often cannot follow a curve of another surface
    // closely enough (a cubic along a circle), so where they fall short the surfaces get more.
    // Only then: more control points also let a surface bend locally between the targets.
    if !met {
        for (&f, s) in faces.iter().zip(&surfaces) {
            raw.faces[f].surface = s.clone();
        }
        for (t, &p) in targets.iter_mut().zip(&points) {
            t.at = p;
        }
        let mut more = false;
        for &f in faces {
            if let Surface::Nurbs(n) = &raw.faces[f].surface
                && let Some(r) = refined(n, MIN_SPANS)
            {
                raw.faces[f].surface = Surface::Nurbs(Box::new(r));
                more = true;
            }
        }
        if more {
            met = gauss_newton(raw, faces, targets);
        }
    }
    let at = targets
        .iter()
        .filter_map(|t| t.vertex.map(|v| (v, t.at)))
        .collect();
    (at, met)
}

/// Damped Gauss–Newton steps on the surfaces of `faces` and on `targets`; whether every target
/// ends on its faces.
fn gauss_newton(raw: &mut Raw, faces: &[usize], targets: &mut [Target]) -> bool {
    let mut vars: HashMap<usize, Vars> = HashMap::new();
    let mut next = 3 * targets.len();
    for &f in faces {
        let pts: Vec<Point3> = targets
            .iter()
            .filter(|t| t.faces.contains(&f))
            .map(|t| t.at)
            .collect();
        if pts.is_empty() {
            continue;
        }
        let center = pts.iter().fold(Vec3::ZERO, |a, &p| a + p) / pts.len() as f64;
        let lever = pts
            .iter()
            .map(|p| p.dist(center))
            .fold(0.0, f64::max)
            .max(1e-3);
        let dofs = match &raw.faces[f].surface {
            Surface::Nurbs(n) => {
                let (group_of, groups) = control_groups(n.control_points());
                Dofs::Control { group_of, groups }
            }
            s => Dofs::Rigid {
                center,
                lever,
                sizes: size_dirs(s, Vec2::new(0.0, 0.0), center).len(),
            },
        };
        let v = Vars { first: next, dofs };
        next += v.count();
        vars.insert(f, v);
    }
    let n_unknowns = next;
    let mut inv_w = vec![1.0; n_unknowns];
    for w in inv_w.iter_mut().skip(3 * targets.len()) {
        *w = 1.0 / SURFACE_COST;
    }
    let mut hints: HashMap<(usize, usize), Vec2> = HashMap::new();
    let Some(mut state) = linearize(raw, targets, &vars, &mut hints, None, 0) else {
        unreachable!("rows without a bound are always built")
    };
    let mut converged = false;
    for _round in 0..ROUNDS {
        if state.worst <= 0.05 * tol::LINEAR {
            converged = true;
            break;
        }
        let mut sys = System {
            rows: Rows::of(&std::mem::take(&mut state.rows)),
            rhs: std::mem::take(&mut state.rhs),
            inv_w: inv_w.clone(),
            blocks: std::mem::take(&mut state.blocks),
        };
        let mut steps: Vec<(f64, Vec<f64>)> = Vec::new();
        let mut taken = None;
        for (k, t) in targets.iter().enumerate() {
            for (i, c) in t.cost.iter().enumerate() {
                sys.inv_w[3 * k + i] = 1.0 / c;
            }
        }
        for (damp, scale) in TRIALS {
            if !steps.iter().any(|s| s.0 == damp) {
                steps.push((damp, sys.solve(damp)));
            }
            let Some((_, x)) = steps.iter().find(|s| s.0 == damp) else {
                break;
            };
            let mut faces_at: Vec<(usize, Surface)> = Vec::with_capacity(vars.len());
            for (&f, v) in &vars {
                let xs: Vec<f64> = x[v.first..v.first + v.count()]
                    .iter()
                    .map(|a| a * scale)
                    .collect();
                match moved(&raw.faces[f].surface, &v.dofs, &xs) {
                    Some(m) => faces_at.push((f, m)),
                    None => return false,
                }
            }
            let mut trial_raw_faces: Vec<(usize, Surface)> = Vec::with_capacity(faces_at.len());
            for (f, m) in faces_at {
                trial_raw_faces.push((f, std::mem::replace(&mut raw.faces[f].surface, m)));
            }
            let before: Vec<Point3> = targets.iter().map(|t| t.at).collect();
            for (k, t) in targets.iter_mut().enumerate() {
                t.at += t.moved_by(&x[3 * k..3 * k + 3]) * scale;
            }
            let mut trial_hints = hints.clone();
            let next = linearize(
                raw,
                targets,
                &vars,
                &mut trial_hints,
                Some(state.worst),
                state.worst_at,
            );
            if let Some(next) = next {
                hints = trial_hints;
                taken = Some(next);
                break;
            }
            for (f, old) in trial_raw_faces {
                raw.faces[f].surface = old;
            }
            for (t, p) in targets.iter_mut().zip(before) {
                t.at = p;
            }
        }
        match taken {
            Some(next) => state = next,
            None => break,
        }
    }
    converged
}

/// The rows of one Gauss–Newton step and the worst distance of any target from any of its faces.
struct Linear {
    rows: Vec<Vec<(usize, f64)>>,
    rhs: Vec<f64>,
    blocks: Vec<(usize, usize)>,
    worst: f64,
    /// The target that lies `worst` from one of its faces.
    worst_at: usize,
}

/// The rows of one Gauss–Newton step at `targets`, or `None` as soon as a target is found as far
/// as `below` (where given) from one of its faces or farther: a step is kept only where it brings the worst
/// gap below the last one, so the rows of a step that does not are never needed. The targets are
/// measured first and the rows built only for a step that is kept, starting with target `first`,
/// the worst before, which a step that overshoots usually leaves worse still.
fn linearize(
    raw: &Raw,
    targets: &[Target],
    vars: &HashMap<usize, Vars>,
    hints: &mut HashMap<(usize, usize), Vec2>,
    below: Option<f64>,
    first: usize,
) -> Option<Linear> {
    let mut worst: f64 = 0.0;
    let mut worst_at = 0;
    // Each target's nearest point on each of its faces, in the order of its faces.
    let mut feet: Vec<Vec<(Vec2, Point3)>> = vec![Vec::new(); targets.len()];
    let order = (first < targets.len())
        .then_some(first)
        .into_iter()
        .chain((0..targets.len()).filter(|&k| k != first));
    for k in order {
        let t = &targets[k];
        for &f in &t.faces {
            let s = &raw.faces[f].surface;
            let uv = match hints.get(&(k, f)) {
                Some(&h) => s.project_near(t.at, h),
                None => s.project(t.at),
            };
            hints.insert((k, f), uv);
            let q = s.eval(uv.x, uv.y);
            let gap = (t.at - q).norm();
            if gap > worst {
                (worst, worst_at) = (gap, k);
            }
            if below.is_some_and(|b| worst >= b) {
                return None;
            }
            feet[k].push((uv, q));
        }
    }
    let mut rows = Vec::new();
    let mut rhs = Vec::new();
    let mut blocks = Vec::with_capacity(targets.len());
    for (k, t) in targets.iter().enumerate() {
        let start = rows.len();
        let fr = t.frame;
        for (&f, &(uv, q)) in t.faces.iter().zip(&feet[k]) {
            let s = &raw.faces[f].surface;
            let r = t.at - q;
            for d in directions(s, uv) {
                let mut row = vec![
                    (3 * k, d.dot(fr[0])),
                    (3 * k + 1, d.dot(fr[1])),
                    (3 * k + 2, d.dot(fr[2])),
                ];
                if let Some(v) = vars.get(&f) {
                    columns(s, v, uv, q, d, &mut row);
                }
                rows.push(row);
                rhs.push(-d.dot(r));
            }
        }
        blocks.push((start, rows.len()));
    }
    Some(Linear {
        rows,
        rhs,
        blocks,
        worst,
        worst_at,
    })
}

/// The columns of face unknowns `v` in the row for the distance of a point from `s` along `d`,
/// where its nearest point of `s` is `q`, at `uv`.
fn columns(s: &Surface, v: &Vars, uv: Vec2, q: Point3, d: Vec3, row: &mut Vec<(usize, f64)>) {
    match &v.dofs {
        Dofs::Rigid {
            center,
            lever,
            sizes,
        } => {
            let arm = (q - *center).cross(d);
            for (i, x) in d.to_array().into_iter().enumerate() {
                row.push((v.first + i, -x));
            }
            for (i, x) in arm.to_array().into_iter().enumerate() {
                row.push((v.first + 3 + i, -x / lever));
            }
            for (i, dir) in size_dirs(s, uv, q).into_iter().take(*sizes).enumerate() {
                row.push((v.first + 6 + i, -d.dot(dir)));
            }
        }
        Dofs::Control { group_of, .. } => {
            if let Surface::Nurbs(n) = s {
                for (idx, b) in n.basis(uv.x, uv.y) {
                    let g = v.first + 3 * group_of[idx];
                    for (i, x) in d.to_array().into_iter().enumerate() {
                        row.push((g + i, -b * x));
                    }
                }
            }
        }
    }
}

/// Surface `s` moved by its unknowns `xs`; `None` where that would change what it is (a radius
/// through zero, a closed B-spline opened).
fn moved(s: &Surface, dofs: &Dofs, xs: &[f64]) -> Option<Surface> {
    match dofs {
        Dofs::Rigid { center, lever, .. } => {
            let sized = with_sizes(s, &xs[6..])?;
            let w = Vec3::new(xs[3], xs[4], xs[5]) / *lever;
            let rot = match w.unit() {
                Some(axis) => Transform::rotation(*center, axis, w.norm()),
                None => Transform::IDENTITY,
            };
            let t = rot.then(&Transform::translation(Vec3::new(xs[0], xs[1], xs[2])));
            Some(sized.transformed(&t))
        }
        Dofs::Control { group_of, .. } => {
            let Surface::Nurbs(n) = s else { return None };
            let ctrl: Vec<Point3> = n
                .control_points()
                .iter()
                .zip(group_of)
                .map(|(&p, &g)| p + Vec3::new(xs[3 * g], xs[3 * g + 1], xs[3 * g + 2]))
                .collect();
            let m = n.with_control_points(ctrl).ok()?;
            let same = m.u_period().is_some() == n.u_period().is_some()
                && m.v_period().is_some() == n.v_period().is_some();
            same.then(|| Surface::Nurbs(Box::new(m)))
        }
    }
}
