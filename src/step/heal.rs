//! Making a file's geometry agree with itself within the kernel tolerance.
//!
//! Exporters write each vertex, edge curve and surface separately, and many hold them together
//! only to their own, looser tolerance (often 1e-6 to 1e-2 mm). The kernel's own tolerance is
//! [`tol::LINEAR`], so where a vertex or edge misses its faces by more, the vertex
//! is moved onto every surface that meets there and the edge is refitted as the intersection of
//! its two faces' surfaces, through its vertices. Where faces miss each other, so that no vertex
//! or edge can lie on both (two tangent faces the file holds apart, a B-spline face that stops
//! short of its corner), their surfaces are moved too ([`adjust`]). Every solid records how far
//! its geometry moved, and one that would move more than [`MAX_HEAL`] is an error. Where healing
//! fails, or moves geometry further than the file says it meets, [`tolerate`] keeps the file's
//! geometry and gives its edges and vertices the local tolerances it does meet within.

use std::collections::{BTreeSet, HashMap};

use super::adjust::{Target, adjust};
use super::topology::{Raw, edge_range, param_of};
use crate::{
    Error, Result,
    geom::{BSpline, Curve, IntersectionCurve, Surface, refine_onto_both},
    math::{Point3, Vec2, Vec3},
    numeric::solve_dense,
    threads, tol,
};

/// How far healing may move a vertex, in mm, and in multiples of how far it was from meeting its
/// constraints: a vertex that has to go further has found some other meeting of its surfaces.
const MAX_MOVE: f64 = 0.5;
const MOVE_PER_GAP: f64 = 10.0;

fn fail<T>(id: u64, what: &str) -> Result<T> {
    Err(Error::InvalidInput(format!("STEP #{id}: {what}")))
}

/// The distance from `p` to `s`, projecting from `hint` when it is near.
fn surface_gap(s: &Surface, p: Point3) -> f64 {
    let uv = s.project(p);
    s.eval(uv.x, uv.y).dist(p)
}

/// Distances to a surface from points taken in turn along a curve, each found from the foot of
/// the one before. Where the surface's nearest point is found by search (a B-spline, a swept or
/// ruled surface), [`Surface::project_near`] keeps the foot Newton's method reaches from that hint
/// only where it lies within the tolerance of the point (on a B-spline, within a share of its seed
/// spacing, which global projection also stops at), and otherwise projects globally, as
/// [`surface_gap`] does. A point that near one sheet of the surface can be nearer another only
/// where two sheets come within the tolerance of each other, so the distance is the one global
/// projection finds, without a global search at every point of a walk that stays on the
/// surface. The elementary surfaces project in closed form, so they are projected afresh.
struct Walk<'a> {
    s: &'a Surface,
    hint: Option<Vec2>,
}

impl<'a> Walk<'a> {
    fn on(s: &'a Surface) -> Walk<'a> {
        Walk { s, hint: None }
    }

    /// The distance from `p` to the surface.
    fn gap(&mut self, p: Point3) -> f64 {
        let s = self.s;
        let uv = match (s, self.hint) {
            (
                Surface::Plane(_)
                | Surface::Cylinder(_)
                | Surface::Cone(_)
                | Surface::Sphere(_)
                | Surface::Torus(_),
                _,
            )
            | (_, None) => s.project(p),
            (_, Some(h)) => s.project_near(p, h),
        };
        self.hint = Some(uv);
        s.eval(uv.x, uv.y).dist(p)
    }
}

/// The even steps [`crate::check`] samples an edge at.
const CHECK_SAMPLES: usize = 16;

fn sample_params(c: &Curve, r: (f64, f64), per_span: usize) -> Vec<f64> {
    let mut cuts = vec![r.0];
    cuts.extend(c.breakpoints(r.0, r.1));
    cuts.push(r.1);
    let per = if matches!(c, Curve::Line(_)) {
        2
    } else {
        per_span.max(16 / (cuts.len() - 1))
    };
    let mut ts = Vec::new();
    for w in cuts.windows(2) {
        for k in 0..per {
            ts.push(w[0] + (w[1] - w[0]) * k as f64 / per as f64);
        }
    }
    // The checker samples every edge at even steps; an edge that fits must fit there too.
    if per > 2 {
        ts.extend((1..CHECK_SAMPLES).map(|k| r.0 + (r.1 - r.0) * k as f64 / CHECK_SAMPLES as f64));
        ts.sort_by(f64::total_cmp);
    }
    ts.push(r.1);
    ts
}

/// Whether edge `e` over `range` ends at its vertices and lies on each of its faces.
fn edge_fits(raw: &Raw, e: usize, range: Option<(f64, f64)>) -> bool {
    let Some(r) = range else { return false };
    let edge = &raw.edges[e];
    let c = &edge.curve;
    if c.eval(r.0).dist(raw.points[edge.start]) > tol::LINEAR
        || c.eval(r.1).dist(raw.points[edge.end]) > tol::LINEAR
    {
        return false;
    }
    let ts = sample_params(c, r, 4);
    let mut faces = edge.faces.clone();
    faces.dedup();
    faces.iter().all(|&f| {
        let s = &raw.faces[f].surface;
        let mut hint = s.project(c.eval(ts[0]));
        ts.iter().all(|&t| {
            let p = c.eval(t);
            let uv = s.project_near(p, hint);
            hint = uv;
            s.eval(uv.x, uv.y).dist(p) <= tol::LINEAR
        })
    })
}

/// A constraint on a vertex: lie on a surface, on a curve, or on the plane through a point
/// square to a direction.
enum On<'a> {
    Surface(&'a Surface),
    Curve(&'a Curve),
    Plane(Point3, Vec3),
}

/// How far from the file's curve a point of a refitted edge may be found: edges whose faces meet
/// away from it are moved onto their meeting.
const MEET_REACH: f64 = 0.05;

/// The point nearest `p0` that lies on every constraint within a small part of the tolerance, by
/// damped Gauss–Newton on the distances; `None` if they do not meet within `reach` of `p0`
/// (by default, a few times the distance `p0` misses them by).
fn solve_vertex(p0: Point3, on: &[On], reach: Option<f64>) -> Option<Point3> {
    let mut x = p0;
    let mut reach = reach;
    let mut best = (f64::INFINITY, p0);
    for _ in 0..100 {
        let mut rows: Vec<(Vec3, f64)> = Vec::new();
        for c in on {
            match c {
                On::Surface(s) => {
                    let uv = s.project(x);
                    let q = s.eval(uv.x, uv.y);
                    let d = x - q;
                    let n = d
                        .unit()
                        .filter(|_| d.norm() > 1e-14)
                        .unwrap_or(s.normal(uv.x, uv.y));
                    rows.push((n, n.dot(d)));
                }
                On::Curve(c) => {
                    let t = param_of(c, x);
                    let (q, tan) = c.d1(t);
                    let d = x - q;
                    let tan = tan.unit_or_zero();
                    let e1 = (d - tan * d.dot(tan))
                        .unit()
                        .filter(|_| d.norm() > 1e-14)
                        .unwrap_or(tan.any_perpendicular().unit_or_zero());
                    let e2 = tan.cross(e1);
                    rows.push((e1, e1.dot(d)));
                    rows.push((e2, e2.dot(d)));
                }
                On::Plane(o, n) => rows.push((*n, n.dot(x - *o))),
            }
        }
        let worst = rows.iter().map(|r| r.1.abs()).fold(0.0, f64::max);
        let reach = *reach.get_or_insert((MOVE_PER_GAP * (worst + tol::LINEAR)).min(MAX_MOVE));
        if worst < best.0 {
            best = (worst, x);
        }
        if worst <= 1e-3 * tol::LINEAR {
            break;
        }
        let mut a = vec![0.0; 9];
        let mut b = vec![0.0; 3];
        for (n, r) in &rows {
            let n = n.to_array();
            for i in 0..3 {
                for j in 0..3 {
                    a[i * 3 + j] += n[i] * n[j];
                }
                b[i] -= n[i] * r;
            }
        }
        // Where the constraints nearly agree in direction (faces that touch), the plain step runs
        // far along their common tangent; more damping shortens it until it stays within reach.
        let trace = a[0] + a[4] + a[8];
        let mut damp = 1e-10 * trace + 1e-300;
        let step = loop {
            let mut ad = a.clone();
            for i in 0..3 {
                ad[i * 3 + i] += damp;
            }
            let st = solve_dense(ad, b.clone(), 3, 1).ok()?;
            let st = Vec3::new(st[0], st[1], st[2]);
            if (x + st).dist(p0) <= reach || damp > trace {
                break st;
            }
            damp *= 100.0;
        };
        if (x + step).dist(p0) > reach {
            break;
        }
        x += step;
        if step.norm() <= 1e-16 * (1.0 + x.norm()) {
            break;
        }
    }
    let x = best.1;
    let worst = on
        .iter()
        .map(|c| match c {
            On::Surface(s) => surface_gap(s, x),
            On::Curve(c) => c.eval(param_of(c, x)).dist(x),
            On::Plane(o, n) => n.dot(x - *o).abs(),
        })
        .fold(0.0, f64::max);
    (worst <= 0.25 * tol::LINEAR && x.dist(p0) <= reach.unwrap_or(0.0)).then_some(x)
}

/// The most points a refit puts through an edge before it gives up on the faces meeting along one
/// curve: halving further only finds the same band where they touch.
const MAX_REFIT_POINTS: usize = 4096;

/// The point on both surfaces nearest `p`: by their tangent planes where they cross, by
/// alternating projection where they touch, or onto the one surface of a seam. Where neither
/// reaches them, the point on both in the plane through `p` square to `along`, the edge's
/// direction there.
fn onto_both(a: &Surface, b: &Surface, same: bool, p: Point3, along: Vec3) -> Option<Point3> {
    let ok = |x: Point3| {
        surface_gap(a, x) <= 0.25 * tol::LINEAR && surface_gap(b, x) <= 0.25 * tol::LINEAR
    };
    if same {
        let uv = a.project(p);
        return Some(a.eval(uv.x, uv.y));
    }
    if ok(p) {
        return Some(p);
    }
    if let Some(x) = refine_onto_both(a, b, p).filter(|&x| ok(x)) {
        return Some(x);
    }
    let mut x = p;
    let (mut ha, mut hb) = (a.project(x), b.project(x));
    for _ in 0..400 {
        ha = a.project_near(x, ha);
        x = a.eval(ha.x, ha.y);
        hb = b.project_near(x, hb);
        let y = b.eval(hb.x, hb.y);
        let moved = y.dist(x);
        x = y;
        if moved <= 1e-3 * tol::LINEAR {
            break;
        }
    }
    if ok(x) {
        return Some(x);
    }
    let along = along.unit()?;
    let on = [On::Surface(a), On::Surface(b), On::Plane(p, along)];
    solve_vertex(p, &on, Some(MEET_REACH))
}

fn chord_params(points: &[Point3], closed: bool) -> Vec<f64> {
    let n = points.len();
    let count = if closed { n } else { n - 1 };
    let mut params = vec![0.0];
    for i in 0..count {
        params.push(params[i] + points[i].dist(points[(i + 1) % n]));
    }
    params
}

/// Edge `e` refitted onto both its faces' surfaces through its vertices, following its file
/// curve; why not, where the surfaces do not meet along it within tolerance.
fn refit(
    raw: &Raw,
    e: usize,
    range: Option<(f64, f64)>,
    samples: Option<&Vec<Point3>>,
) -> std::result::Result<(Curve, (f64, f64)), String> {
    let edge = &raw.edges[e];
    let mut faces = edge.faces.clone();
    faces.sort_unstable();
    if faces.len() != 2 {
        return Err(format!("the edge is used by {} faces, not 2", faces.len()));
    }
    let (a, b) = (&raw.faces[faces[0]].surface, &raw.faces[faces[1]].surface);
    let same = faces[0] == faces[1];
    let closed = edge.start == edge.end;
    let (ps, pe) = (raw.points[edge.start], raw.points[edge.end]);
    let c = &edge.curve;
    let r = match range {
        Some(r) => r,
        None => {
            let (t0, t1) = (param_of(c, ps), param_of(c, pe));
            match c.period() {
                Some(p) if closed => (t0, t0 + p),
                Some(p) => (t0, t0 + (t1 - t0).rem_euclid(p)),
                None => (t0.min(t1), t0.max(t1)),
            }
        }
    };
    let length: f64 = {
        let ts = sample_params(c, r, 8);
        ts.windows(2).map(|w| c.eval(w[0]).dist(c.eval(w[1]))).sum()
    };
    let n = ((length / (100.0 * tol::LINEAR)) as usize).clamp(1, 8);
    // Where the surface step placed points along the edge, they lie on both faces already: two
    // faces that touch meet anywhere across a band, and refining points onto both picks a line.
    let placed = samples.filter(|s| !s.is_empty());
    let mut pts = match placed {
        Some(s) => {
            let mut p = vec![ps];
            p.extend(s.iter().copied());
            if !closed {
                p.push(pe);
            }
            p
        }
        None => {
            let mut p = bent(c, r, n, (ps, pe));
            if closed {
                p.pop();
            }
            p
        }
    };
    if closed && pts.len() < 3 {
        return Err("the closed edge is too short to refit".into());
    }
    let kinds = format!("{} and {}", a.kind(), b.kind());
    let refine = |p: Point3, along: Vec3| onto_both(a, b, same, p, along);
    let last = pts.len() - 1;
    let count = pts.len();
    let before = pts.clone();
    for (i, p) in pts.iter_mut().enumerate() {
        if closed || (i != 0 && i != last) {
            let t = r.0 + (r.1 - r.0) * i as f64 / (if closed { count } else { last }) as f64;
            let along = match placed {
                Some(_) => before[(i + 1) % count] - before[(i + count - 1) % count],
                None => c.d1(t).1,
            };
            *p = refine(*p, along)
                .ok_or_else(|| format!("its {kinds} faces do not meet near it"))?;
        }
    }
    for _ in 0..16 {
        let approx = if closed {
            BSpline::interpolate_periodic(&pts).ok()
        } else {
            BSpline::interpolate(&pts).ok()
        };
        let Some(approx) = approx else {
            return Err("its points on both faces could not be interpolated".into());
        };
        let params = chord_params(&pts, closed);
        let n_seg = params.len() - 1;
        // The samples are taken in order along the spline, so each is projected from the last.
        let (mut wa, mut wb) = (Walk::on(a), Walk::on(b));
        let mut gap = |m: Point3| wa.gap(m).max(wb.gap(m));
        let err: Vec<f64> = (0..n_seg)
            .map(|i| {
                [0.25, 0.5, 0.75]
                    .iter()
                    .map(|f| gap(approx.eval(params[i] + (params[i + 1] - params[i]) * f)))
                    .fold(0.0, f64::max)
            })
            .collect();
        if err.iter().all(|&x| x <= 0.5 * tol::LINEAR) {
            let dom = approx.domain();
            // A spline that stays on both faces everywhere is the edge as it is: evaluating it
            // costs nothing, where an intersection curve refines every point onto both faces.
            let (mut wa, mut wb) = (Walk::on(a), Walk::on(b));
            let mut gap = |m: Point3| wa.gap(m).max(wb.gap(m));
            let on_both = (0..n_seg).all(|i| {
                (1..8).all(|k| {
                    let m = approx.eval(params[i] + (params[i + 1] - params[i]) * k as f64 / 8.0);
                    gap(m) <= 0.5 * tol::LINEAR
                })
            }) && {
                let (mut wa, mut wb) = (Walk::on(a), Walk::on(b));
                (1..CHECK_SAMPLES).all(|k| {
                    let m = approx.eval(dom.0 + (dom.1 - dom.0) * k as f64 / CHECK_SAMPLES as f64);
                    wa.gap(m).max(wb.gap(m)) <= 0.5 * tol::LINEAR
                })
            };
            let curve = if on_both {
                Curve::BSpline(approx)
            } else {
                Curve::Intersection(Box::new(IntersectionCurve {
                    approx,
                    surfaces: [a.clone(), b.clone()],
                }))
            };
            return Ok((curve, dom));
        }
        let mut next = Vec::with_capacity(pts.len() * 2);
        for i in 0..n_seg {
            next.push(pts[i]);
            if err[i] > 0.125 * tol::LINEAR {
                let m = approx.eval(0.5 * (params[i] + params[i + 1]));
                let along = pts[(i + 1) % pts.len()] - pts[i];
                let q =
                    refine(m, along).ok_or_else(|| format!("its {kinds} faces part along it"))?;
                // A midpoint that lands on its neighbours means the faces meet along no single
                // curve here: they touch, or meet the file's curve somewhere else.
                if q.dist(pts[i]) <= 4.0 * tol::LINEAR
                    || q.dist(pts[(i + 1) % pts.len()]) <= 4.0 * tol::LINEAR
                {
                    return Err(format!(
                        "its {kinds} faces meet {:.1e} mm off the curve through its ends",
                        err[i]
                    ));
                }
                next.push(q);
            }
        }
        if !closed {
            next.push(pts[pts.len() - 1]);
        }
        pts = next;
        if pts.len() > MAX_REFIT_POINTS {
            break;
        }
    }
    Err(format!(
        "no curve through {} points stays within tolerance of its {kinds} faces, which touch \
         rather than cross",
        pts.len()
    ))
}

/// The range of edge `e` between its (moved) vertices nearest its range `old` before they moved:
/// on a periodic curve a vertex moved back past its end must not turn a short arc into a long one.
fn range_near(raw: &Raw, e: usize, old: (f64, f64)) -> (f64, f64) {
    let edge = &raw.edges[e];
    let c = &edge.curve;
    let near = |p: Point3, to: f64| {
        let t = param_of(c, p);
        match c.period() {
            Some(per) => t + ((to - t) / per).round() * per,
            None => t,
        }
    };
    let t0 = near(raw.points[edge.start], old.0);
    match c.period() {
        Some(per) if edge.start == edge.end => (t0, t0 + per),
        _ => (t0, near(raw.points[edge.end], old.1)),
    }
}

/// How far healing moved a solid's geometry: the largest movement, where, and what moved.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct Healed {
    pub moved: f64,
    pub at: Point3,
    pub what: String,
    /// The largest distance any surface moved, which bounds how far the solid's faces moved.
    pub surfaces: f64,
}

/// Each edge's parameter range, and how far healing moved the solid's geometry.
pub(super) type Ranges = (Vec<(f64, f64)>, Option<Healed>);

/// The most healing may move any of a solid's geometry, in mm.
pub const MAX_HEAL: f64 = 0.05;

/// How many times surfaces are moved before a solid that still does not close is rejected.
const ADJUST_ROUNDS: usize = 3;

/// A vertex or edge that would not fit, and why.
enum Failure {
    Vertex(usize, String),
    Edge(usize, String),
}

fn range_of(raw: &Raw, e: usize) -> Option<(f64, f64)> {
    let edge = &raw.edges[e];
    edge_range(
        &edge.curve,
        raw.points[edge.start],
        raw.points[edge.end],
        edge.start == edge.end,
    )
}

/// The fewest refits worth a thread of their own: a refit takes milliseconds, a thread far less.
const REFITS_PER_RUN: usize = 8;

/// What [`refit`] made of an edge, from what: the version of each face's surface, the edge's ends
/// and range, and the points the surface step placed along it.
type RefitKey = (
    [usize; 2],
    [Point3; 2],
    Option<(f64, f64)>,
    Option<Vec<Point3>>,
);
type Refitted = std::result::Result<(Curve, (f64, f64)), String>;

/// The edges [`heal`] has refitted in earlier rounds. A round moves only the surfaces of the faces
/// it adjusts, so an edge whose two faces kept their surfaces, and whose ends and range came out
/// the same, refits as it did before: [`refit`] reads nothing else.
#[derive(Default)]
struct Refits {
    /// How many times each face's surface has moved.
    versions: Vec<usize>,
    done: HashMap<usize, (RefitKey, Refitted)>,
}

/// Moves each vertex that misses its faces onto them, starting from `start` where given, and
/// refits each edge that still misses. Returns each edge's range, or every vertex and edge that
/// would not fit. `loose` refits the edges at a vertex that would not move too, from where it is.
fn attempt(
    raw: &mut Raw,
    start: &HashMap<usize, Point3>,
    samples: &HashMap<usize, Vec<Point3>>,
    loose: bool,
    mut refits: Option<&mut Refits>,
) -> std::result::Result<Vec<(f64, f64)>, Vec<Failure>> {
    let ranges: Vec<Option<(f64, f64)>> = (0..raw.edges.len()).map(|e| range_of(raw, e)).collect();
    let fits: Vec<bool> = (0..raw.edges.len())
        .map(|e| edge_fits(raw, e, ranges[e]))
        .collect();
    let mut at_vertex: Vec<Vec<usize>> = vec![Vec::new(); raw.points.len()];
    for (e, edge) in raw.edges.iter().enumerate() {
        at_vertex[edge.start].push(e);
        if edge.end != edge.start {
            at_vertex[edge.end].push(e);
        }
    }
    let mut failures = Vec::new();
    let mut moved = raw.points.clone();
    let mut stuck = vec![false; raw.points.len()];
    for (v, edges) in at_vertex.iter().enumerate() {
        let p0 = start.get(&v).copied().unwrap_or(raw.points[v]);
        if edges.iter().all(|&e| fits[e]) && p0 == raw.points[v] {
            continue;
        }
        let mut faces: Vec<usize> = edges
            .iter()
            .flat_map(|&e| raw.edges[e].faces.iter().copied())
            .collect();
        faces.sort_unstable();
        faces.dedup();
        let mut on: Vec<On> = faces
            .iter()
            .map(|&f| On::Surface(&raw.faces[f].surface))
            .collect();
        let surfaces_only = on.len();
        on.extend(
            edges
                .iter()
                .filter(|&&e| fits[e])
                .map(|&e| On::Curve(&raw.edges[e].curve)),
        );
        let p =
            solve_vertex(p0, &on, None).or_else(|| solve_vertex(p0, &on[..surfaces_only], None));
        match p {
            Some(p) => moved[v] = p,
            None => {
                let gap = faces
                    .iter()
                    .map(|&f| surface_gap(&raw.faces[f].surface, p0))
                    .fold(0.0, f64::max);
                stuck[v] = true;
                failures.push(Failure::Vertex(
                    v,
                    format!(
                        "the {} faces at an end of this edge have no common point near it (the \
                         file's vertex is {gap:.1e} mm off one)",
                        faces.len()
                    ),
                ));
            }
        }
    }
    let before = std::mem::replace(&mut raw.points, moved);
    // A vertex moved without its edge's curve keeps the edge only where the curve still ends
    // close by: a face's pcurves meet there only as closely as the curves' ends, measured in
    // parameters, which on a thin cylinder is stricter than the tolerance in space.
    let ends_close = |raw: &Raw, e: usize, r: Option<(f64, f64)>| {
        let edge = &raw.edges[e];
        r.is_some_and(|r| {
            [(r.0, edge.start), (r.1, edge.end)].iter().all(|&(t, v)| {
                raw.points[v] == before[v]
                    || edge.curve.eval(t).dist(raw.points[v]) <= 0.25 * tol::LINEAR
            })
        })
    };

    // Each edge is kept, or refitted from its own curve, its faces and its ends alone, so the
    // refits are independent of each other: they are made on as many threads as there are, and
    // taken in the order of the edges.
    let mut out = vec![(0.0, 0.0); raw.edges.len()];
    let mut todo: Vec<(usize, Option<(f64, f64)>)> = Vec::new();
    for (e, &old) in ranges.iter().enumerate() {
        let edge = &raw.edges[e];
        if !loose && (stuck[edge.start] || stuck[edge.end]) {
            continue;
        }
        let r = match old {
            Some(old) => Some(range_near(raw, e, old)),
            None => range_of(raw, e),
        };
        if edge_fits(raw, e, r) && ends_close(raw, e, r) {
            out[e] = r.unwrap_or((0.0, 0.0));
            continue;
        }
        todo.push((e, r.or(old)));
    }
    let key = |memo: &Refits, e: usize, r: Option<(f64, f64)>| -> RefitKey {
        let edge = &raw.edges[e];
        let (fa, fb) = match edge.faces[..] {
            [fa, fb] => (fa.min(fb), fa.max(fb)),
            _ => (usize::MAX, usize::MAX),
        };
        let version = |f: usize| memo.versions.get(f).copied().unwrap_or(usize::MAX);
        (
            [version(fa), version(fb)],
            [raw.points[edge.start], raw.points[edge.end]],
            r,
            samples.get(&e).cloned(),
        )
    };
    let keys: Vec<Option<RefitKey>> = todo
        .iter()
        .map(|&(e, r)| refits.as_deref().map(|memo| key(memo, e, r)))
        .collect();
    let kept: Vec<Option<Refitted>> = todo
        .iter()
        .zip(&keys)
        .map(|(&(e, _), k)| match (refits.as_deref(), k) {
            (Some(memo), Some(k)) => memo
                .done
                .get(&e)
                .filter(|(old, _)| old == k)
                .map(|(_, rf)| rf.clone()),
            _ => None,
        })
        .collect();
    let fresh: Vec<usize> = (0..todo.len()).filter(|&i| kept[i].is_none()).collect();
    let made = {
        let raw = &*raw;
        threads::map(fresh.len(), REFITS_PER_RUN, |i| {
            let (e, r) = todo[fresh[i]];
            match samples.get(&e) {
                Some(s) => refit(raw, e, r, Some(s)).or_else(|_| refit(raw, e, r, None)),
                None => refit(raw, e, r, None),
            }
        })
    };
    let mut made = made.into_iter();
    for ((&(e, _), k), rf) in todo.iter().zip(keys).zip(kept) {
        let rf = match rf {
            Some(rf) => rf,
            None => {
                let rf = made.next().expect("one refit made for each edge not kept");
                if let (Some(memo), Some(k)) = (refits.as_deref_mut(), k) {
                    memo.done.insert(e, (k, rf.clone()));
                }
                rf
            }
        };
        match rf {
            Ok((curve, r)) => {
                raw.edges[e].curve = curve;
                out[e] = r;
            }
            Err(why) => failures.push(Failure::Edge(e, why)),
        }
    }
    if failures.is_empty() {
        Ok(out)
    } else {
        Err(failures)
    }
}

/// The parameter range of edge `e` between its vertices, however far they are from its curve.
fn curve_range(raw: &Raw, e: usize) -> (f64, f64) {
    if let Some(r) = range_of(raw, e) {
        return r;
    }
    let edge = &raw.edges[e];
    let c = &edge.curve;
    let (t0, t1) = (
        param_of(c, raw.points[edge.start]),
        param_of(c, raw.points[edge.end]),
    );
    match c.period() {
        Some(p) if edge.start == edge.end => (t0, t0 + p),
        Some(p) => (t0, t0 + (t1 - t0).rem_euclid(p)),
        None => (t0.min(t1), t0.max(t1)),
    }
}

/// Points along edge `e` strictly between its ends.
fn inner_samples(raw: &Raw, e: usize, count: usize) -> Vec<Point3> {
    let c = &raw.edges[e].curve;
    let r = curve_range(raw, e);
    (1..=count)
        .map(|i| c.eval(r.0 + (r.1 - r.0) * i as f64 / (count + 1) as f64))
        .collect()
}

/// `n + 1` points along curve `c` over `r`, its two ends moved to `ends` and the points between
/// moved by the blend of those two moves. A file's vertex can sit well off its edge's curve, to
/// the side where two faces touch; following the curve unbent would lead a refit, or the surface
/// step, to a second line where the faces meet beside the vertex's.
fn bent(c: &Curve, r: (f64, f64), n: usize, ends: (Point3, Point3)) -> Vec<Point3> {
    let (d0, d1) = (ends.0 - c.eval(r.0), ends.1 - c.eval(r.1));
    let mut pts: Vec<Point3> = (0..=n)
        .map(|i| {
            let s = i as f64 / n as f64;
            c.eval(r.0 + (r.1 - r.0) * s) + d0 * (1.0 - s) + d1 * s
        })
        .collect();
    pts[0] = ends.0;
    pts[n] = ends.1;
    pts
}

/// The points the surface step moves: every vertex and edge on `faces`, more densely along the
/// edges that would not fit.
fn targets(
    raw: &Raw,
    faces: &BTreeSet<usize>,
    failing: &BTreeSet<usize>,
    start: &HashMap<usize, Point3>,
) -> Vec<Target> {
    let mut out = Vec::new();
    let mut vertices = BTreeSet::new();
    let mut at_vertex: HashMap<usize, Vec<usize>> = HashMap::new();
    for (e, edge) in raw.edges.iter().enumerate() {
        at_vertex.entry(edge.start).or_default().push(e);
        at_vertex.entry(edge.end).or_default().push(e);
        if !edge.faces.iter().any(|f| faces.contains(f)) {
            continue;
        }
        vertices.insert(edge.start);
        vertices.insert(edge.end);
        let mut fs = edge.faces.clone();
        fs.sort_unstable();
        fs.dedup();
        // Faces that touch meet only where the step puts a point, so a failing edge gets several
        // points per span of its curve: between them the moved surfaces could part again.
        let count = if failing.contains(&e) {
            let r = curve_range(raw, e);
            (4 * (edge.curve.breakpoints(r.0, r.1).len() + 1)).clamp(33, 96)
        } else {
            3
        };
        let at = |v: usize| start.get(&v).copied().unwrap_or(raw.points[v]);
        let pts = bent(
            &edge.curve,
            curve_range(raw, e),
            count + 1,
            (at(edge.start), at(edge.end)),
        );
        for i in 1..=count {
            out.push(Target::sample(
                raw,
                pts[i],
                fs.clone(),
                e,
                pts[i + 1] - pts[i - 1],
            ));
        }
    }
    for v in vertices {
        let mut fs: Vec<usize> = at_vertex[&v]
            .iter()
            .flat_map(|&e| raw.edges[e].faces.iter().copied())
            .collect();
        fs.sort_unstable();
        fs.dedup();
        out.push(Target::vertex(
            raw,
            start.get(&v).copied().unwrap_or(raw.points[v]),
            fs,
            v,
        ));
    }
    out
}

/// The largest distance between `file`'s geometry and `healed`'s: vertices, the surfaces of
/// `faces` along their edges, and edges given new curves.
fn movement(file: &Raw, healed: &Raw, faces: &BTreeSet<usize>) -> Option<Healed> {
    let mut best: Option<Healed> = None;
    let mut surfaces: f64 = 0.0;
    let mut note = |moved: f64, at: Point3, what: String| {
        if moved > best.as_ref().map_or(0.0, |b| b.moved) {
            best = Some(Healed {
                moved,
                at,
                what,
                surfaces: 0.0,
            });
        }
    };
    for (v, (&a, &b)) in file.points.iter().zip(&healed.points).enumerate() {
        note(a.dist(b), b, format!("vertex #{}", file.vertex_ids[v]));
    }
    for &f in faces {
        let (s0, s1) = (&file.faces[f].surface, &healed.faces[f].surface);
        for (e, edge) in file.edges.iter().enumerate() {
            if !edge.faces.contains(&f) {
                continue;
            }
            for p in inner_samples(file, e, 5) {
                let uv = s0.project(p);
                let a = s0.eval(uv.x, uv.y);
                let uv = s1.project_near(a, uv);
                let b = s1.eval(uv.x, uv.y);
                surfaces = surfaces.max(a.dist(b));
                note(
                    a.dist(b),
                    b,
                    format!("the surface of face #{}", file.faces[f].id),
                );
            }
        }
    }
    for (e, (a, b)) in file.edges.iter().zip(&healed.edges).enumerate() {
        if a.curve == b.curve {
            continue;
        }
        for p in inner_samples(healed, e, 7) {
            let q = a.curve.eval(param_of(&a.curve, p));
            note(p.dist(q), p, format!("edge #{}", a.id));
        }
    }
    best.map(|b| Healed { surfaces, ..b })
}

/// Makes `raw` agree with itself within the kernel tolerance, and returns each edge's parameter
/// range and how far its geometry moved. A file that already agrees is left exactly as it is.
///
/// Vertices move onto their faces and edges are refitted onto theirs; where faces miss each
/// other so that cannot work, [`adjust`] moves the surfaces of the faces involved. A solid whose
/// geometry would move more than [`MAX_HEAL`] is an error naming the distance.
pub(super) fn heal(raw: &mut Raw) -> Result<Ranges> {
    let ranges: Vec<Option<(f64, f64)>> = (0..raw.edges.len()).map(|e| range_of(raw, e)).collect();
    if (0..raw.edges.len()).all(|e| edge_fits(raw, e, ranges[e])) {
        let ranges = ranges
            .into_iter()
            .map(|r| r.unwrap_or((0.0, 0.0)))
            .collect();
        return Ok((ranges, None));
    }
    let file = raw.clone();
    let mut surfaces: Vec<Surface> = file.faces.iter().map(|f| f.surface.clone()).collect();
    let mut chosen = BTreeSet::new();
    let mut failing = BTreeSet::new();
    let mut start = HashMap::new();
    let mut samples: HashMap<usize, Vec<Point3>> = HashMap::new();
    let mut settled = true;
    let mut failed_before = usize::MAX;
    let mut refits = Refits {
        versions: vec![0; file.faces.len()],
        ..Refits::default()
    };
    let with_surfaces = |surfaces: &[Surface]| {
        let mut r = file.clone();
        for (f, s) in r.faces.iter_mut().zip(surfaces) {
            f.surface = s.clone();
        }
        r
    };
    for round in 0..=ADJUST_ROUNDS {
        let mut work = with_surfaces(&surfaces);
        let failures = match attempt(&mut work, &start, &samples, false, Some(&mut refits)) {
            Ok(ranges) => {
                let healed = movement(&file, &work, &chosen);
                if let Some(h) = healed.as_ref().filter(|h| h.moved > MAX_HEAL) {
                    return fail(
                        file.edges.first().map_or(0, |e| e.id),
                        &format!(
                            "healing would move {} {:.1e} mm, more than the {MAX_HEAL} mm allowed",
                            h.what, h.moved
                        ),
                    );
                }
                *raw = work;
                return Ok((ranges, healed));
            }
            Err(f) => f,
        };
        let before = chosen.len();
        for f in &failures {
            match f {
                Failure::Vertex(v, _) => {
                    for edge in file.edges.iter().filter(|e| e.start == *v || e.end == *v) {
                        chosen.extend(edge.faces.iter().copied());
                    }
                }
                Failure::Edge(e, _) => {
                    failing.insert(*e);
                    chosen.extend(file.edges[*e].faces.iter().copied());
                }
            }
        }
        // Moving the surfaces again helps only where the last move closed every gap or the
        // failures reached new faces, and only from a start no worse than the last one: a move
        // that missed its targets and left more failing than before only spreads to more faces.
        let worse = round > 0 && !settled && failures.len() > failed_before;
        failed_before = failures.len();
        let stuck = round > 0 && chosen.len() == before && !settled;
        if round == ADJUST_ROUNDS || stuck || worse {
            let (id, why) = match &failures[0] {
                Failure::Vertex(v, why) => (
                    file.edges
                        .iter()
                        .find(|e| e.start == *v || e.end == *v)
                        .map_or(0, |e| e.id),
                    why.clone(),
                ),
                Failure::Edge(e, why) => (file.edges[*e].id, why.clone()),
            };
            return fail(id, &why);
        }
        let mut shaped = with_surfaces(&surfaces);
        let mut pts = targets(&file, &chosen, &failing, &start);
        let faces: Vec<usize> = chosen.iter().copied().collect();
        let (at, met) = adjust(&mut shaped, &faces, &mut pts);
        settled = met;
        for (f, face) in shaped.faces.iter().enumerate() {
            if face.surface != surfaces[f] {
                refits.versions[f] += 1;
            }
        }
        surfaces = shaped.faces.into_iter().map(|f| f.surface).collect();
        start = at;
        samples.clear();
        for t in &pts {
            if let Some(e) = t.edge.filter(|e| failing.contains(e)) {
                samples.entry(e).or_default().push(t.at);
            }
        }
    }
    unreachable!("the last round returns")
}

/// Moves each vertex of `raw` onto its faces and refits each edge onto its faces wherever that
/// works, leaving the rest as the file has it and every surface where it is, for [`tolerate`] to
/// hold what is left: where faces meet nowhere exactly, as four faces at a corner can miss each
/// other by far less than the file's uncertainty, healing fails, but an edge elsewhere may still
/// be further off than the file promises until it is refitted. Returns how far it moved.
pub(super) fn heal_loosely(raw: &mut Raw) -> Option<Healed> {
    let file = raw.clone();
    let _ = attempt(raw, &HashMap::new(), &HashMap::new(), true, None);
    movement(&file, raw, &BTreeSet::new())
}

/// What [`tolerate`] gives a solid it keeps as the file has it: each edge's range, each vertex's
/// and edge's tolerance (at least [`tol::LINEAR`]), and how far the curves it put [`midway`]
/// between their faces moved those edges.
pub(super) struct Tolerated {
    pub ranges: Vec<(f64, f64)>,
    pub vertices: Vec<f64>,
    pub edges: Vec<f64>,
    pub moved: Option<Healed>,
}

/// How much wider than the largest gap sampled along an edge its tolerance is: the gap between
/// samples can be larger than at them.
const GAP_MARGIN: f64 = 1.5;

/// The range of edge `e` between its vertices as the file places them, which may lie off its
/// curve; `None` where they give no stretch of it.
fn file_range(raw: &Raw, e: usize) -> Option<(f64, f64)> {
    let edge = &raw.edges[e];
    let c = &edge.curve;
    if let Some(r) = range_of(raw, e) {
        return Some(r);
    }
    let (t0, t1) = (
        param_of(c, raw.points[edge.start]),
        param_of(c, raw.points[edge.end]),
    );
    let r = match c.period() {
        Some(p) if edge.start == edge.end => (t0, t0 + p),
        Some(p) => (t0, t0 + (t1 - t0).rem_euclid(p)),
        None if edge.start == edge.end => c.domain(),
        None => (t0.min(t1), t0.max(t1)),
    };
    (r.0.is_finite() && r.1.is_finite() && r.1 > r.0).then_some(r)
}

/// The largest distance from curve `c` over `r` to surface `s`, sampled densely.
fn curve_gap(c: &Curve, r: (f64, f64), s: &Surface) -> f64 {
    let mut ts = sample_params(c, r, 8);
    ts.extend((1..64).map(|k| r.0 + (r.1 - r.0) * k as f64 / 64.0));
    ts.sort_by(f64::total_cmp);
    ts.dedup();
    // The walk starts in the middle: an edge's end may sit at a pole of the surface, where the
    // point nearest it is found only by arriving from along the edge.
    let mid = ts.len() / 2;
    let start = s.project(c.eval(ts[mid]));
    let mut worst: f64 = 0.0;
    for half in [&ts[mid..], &ts[..=mid]] {
        let mut hint = start;
        let walk: Vec<f64> = if std::ptr::eq(half.as_ptr(), ts[mid..].as_ptr()) {
            half.to_vec()
        } else {
            half.iter().rev().copied().collect()
        };
        for t in walk {
            let p = c.eval(t);
            hint = s.project_near(p, hint);
            let d = s.eval(hint.x, hint.y).dist(p);
            // A projection that stalls short of the nearest point is tried again from scratch.
            let d = if d > tol::LINEAR {
                d.min(surface_gap(s, p))
            } else {
                d
            };
            worst = worst.max(d);
        }
    }
    worst
}

/// A curve for edge `e` over its file range `r` halfway between its two faces, through its
/// vertices as they are: where the faces lie apart along it rather than crossing, the file's curve
/// on one of them misses the other by the whole gap, and this one by half.
fn midway(raw: &Raw, e: usize, r: (f64, f64)) -> Option<Curve> {
    let edge = &raw.edges[e];
    let mut faces = edge.faces.clone();
    faces.sort_unstable();
    faces.dedup();
    let [fa, fb] = faces.as_slice() else {
        return None;
    };
    let (a, b) = (&raw.faces[*fa].surface, &raw.faces[*fb].surface);
    let closed = edge.start == edge.end;
    let mut ts = sample_params(&edge.curve, r, 8);
    ts.sort_by(f64::total_cmp);
    ts.dedup_by(|x, y| (*x - *y).abs() <= 1e-9 * (r.1 - r.0));
    if ts.len() < 2 {
        return None;
    }
    let last = ts.len() - 1;
    let mut pts: Vec<Point3> = ts
        .iter()
        .enumerate()
        .map(|(i, &t)| match i {
            0 => raw.points[edge.start],
            i if i == last => raw.points[edge.end],
            _ => {
                let mut p = edge.curve.eval(t);
                for _ in 0..50 {
                    let (ua, ub) = (a.project(p), b.project(p));
                    let (qa, qb) = (a.eval(ua.x, ua.y), b.eval(ub.x, ub.y));
                    let next = qa + (qb - qa) * 0.5;
                    let done = next.dist(p) <= 1e-12;
                    p = next;
                    if done {
                        break;
                    }
                }
                p
            }
        })
        .collect();
    let spline = if closed {
        pts.pop();
        BSpline::interpolate_periodic(&pts)
    } else {
        BSpline::interpolate(&pts)
    };
    spline.ok().map(Curve::BSpline)
}

/// The share of an edge's length that its tolerance, or that of a vertex it ends at, may reach:
/// the spheres of its two ends then cover at most half of it. Past that the ends cannot be told
/// apart, and the solid has no definite shape.
const MAX_SHARE: f64 = 0.25;

/// An edge's length over `r`, near enough to judge a tolerance against.
fn edge_length(c: &Curve, r: (f64, f64)) -> f64 {
    let at = |k: usize| c.eval(r.0 + (r.1 - r.0) * k as f64 / 32.0);
    (0..32).map(|k| at(k).dist(at(k + 1))).sum()
}

/// Keeps a solid's geometry as the file has it where healing could not make it meet within
/// [`tol::LINEAR`]: each edge and vertex that misses its faces by more gets a local tolerance that
/// holds the gap it measures, as the exporting system's own does. What the file states as its
/// uncertainty does not bound that: exporters write gaps wider than they state. An edge that
/// misses by more than `stated` takes a curve [`midway`] between its faces where that one misses
/// by less. A tolerance past [`MAX_SHARE`] of an edge it holds is an error.
pub(super) fn tolerate(
    raw: &mut Raw,
    stated: Option<f64>,
) -> std::result::Result<Tolerated, String> {
    let file = raw.clone();
    let mut ranges = Vec::with_capacity(raw.edges.len());
    let mut edges = Vec::with_capacity(raw.edges.len());
    let mut vertices = vec![tol::LINEAR; raw.points.len()];
    // The shortest edge ending at each vertex.
    let mut shortest = vec![f64::INFINITY; raw.points.len()];
    let widened = |gap: f64| {
        if gap <= tol::LINEAR {
            tol::LINEAR
        } else {
            GAP_MARGIN * gap
        }
    };
    for e in 0..raw.edges.len() {
        let Some(mut r) = file_range(raw, e) else {
            return Err(format!(
                "edge #{} has no stretch between its vertices",
                raw.edges[e].id
            ));
        };
        let mut faces = raw.edges[e].faces.clone();
        faces.sort_unstable();
        faces.dedup();
        let gap_of = |c: &Curve, r: (f64, f64)| {
            faces
                .iter()
                .map(|&f| curve_gap(c, r, &raw.faces[f].surface))
                .fold(0.0, f64::max)
        };
        let mut gap = gap_of(&raw.edges[e].curve, r);
        if gap > stated.unwrap_or(tol::LINEAR)
            && let Some(c) = midway(raw, e, r)
        {
            let (rm, gm) = (c.domain(), gap_of(&c, c.domain()));
            if gm < gap {
                (r, gap) = (rm, gm);
                raw.edges[e].curve = c;
            }
        }
        let edge = &raw.edges[e];
        let length = edge_length(&edge.curve, r);
        let t = widened(gap);
        if t > tol::LINEAR && t > MAX_SHARE * length {
            return Err(format!(
                "edge #{} misses its faces by {gap:.1e} mm, too far for an edge {length:.1e} mm long",
                edge.id
            ));
        }
        for (v, at) in [(edge.start, r.0), (edge.end, r.1)] {
            let p = raw.points[v];
            let far = faces
                .iter()
                .map(|&f| surface_gap(&raw.faces[f].surface, p))
                .fold(edge.curve.eval(at).dist(p), f64::max);
            vertices[v] = vertices[v].max(t).max(widened(far));
            shortest[v] = shortest[v].min(length);
        }
        ranges.push(r);
        edges.push(t);
    }
    for (v, &t) in vertices.iter().enumerate() {
        if t > tol::LINEAR && t > MAX_SHARE * shortest[v] {
            return Err(format!(
                "vertex #{} lies {:.1e} mm off its edges and faces, too far for an edge {:.1e} mm \
                 long that ends at it",
                raw.vertex_ids[v],
                t / GAP_MARGIN,
                shortest[v]
            ));
        }
    }
    Ok(Tolerated {
        ranges,
        vertices,
        edges,
        moved: movement(&file, raw, &BTreeSet::new()),
    })
}

/// Where `raw` first disagrees with itself by more than [`tol::LINEAR`]: an edge that does not end
/// at its vertices or does not lie on one of its faces, sampled as [`heal`] samples it to decide
/// whether a file needs healing at all. `None` where the file agrees with itself.
pub(super) fn misfit(raw: &Raw) -> Option<String> {
    for e in 0..raw.edges.len() {
        let edge = &raw.edges[e];
        let Some(r) = range_of(raw, e) else {
            return Some(format!(
                "edge #{} does not run between its vertices",
                edge.id
            ));
        };
        if edge_fits(raw, e, Some(r)) {
            continue;
        }
        let c = &edge.curve;
        let ends = [(r.0, edge.start), (r.1, edge.end)]
            .iter()
            .map(|&(t, v)| (c.eval(t).dist(raw.points[v]), raw.vertex_ids[v]))
            .fold((0.0, 0), |a, b| if b.0 > a.0 { b } else { a });
        if ends.0 > tol::LINEAR {
            return Some(format!(
                "edge #{} ends {:.1e} mm from vertex #{}",
                edge.id, ends.0, ends.1
            ));
        }
        let (gap, f) = edge
            .faces
            .iter()
            .map(|&f| (curve_gap(c, r, &raw.faces[f].surface), raw.faces[f].id))
            .fold((0.0, 0), |a, b| if b.0 > a.0 { b } else { a });
        return Some(format!("edge #{} lies {gap:.1e} mm off face #{f}", edge.id));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{NurbsSurface, Revolution};
    use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, TAU};

    /// A quarter of the cylinder of radius 1 about the z axis, 0 ≤ z ≤ 2, as a rational quadratic
    /// by linear patch.
    fn quarter_cylinder() -> Surface {
        let arc = [(1.0, 0.0, 1.0), (1.0, 1.0, FRAC_1_SQRT_2), (0.0, 1.0, 1.0)];
        let mut ctrl = Vec::new();
        let mut weights = Vec::new();
        for &(x, y, w) in &arc {
            for z in [0.0, 2.0] {
                ctrl.push(Vec3::new(x, y, z));
                weights.push(w);
            }
        }
        let knots = (vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], vec![0.0, 0.0, 1.0, 1.0]);
        let s = NurbsSurface::new((2, 1), knots, (3, 2), ctrl, Some(weights)).unwrap();
        Surface::Nurbs(Box::new(s))
    }

    /// The cylinder of radius 2 about the z axis, 0 ≤ z ≤ 2, as a straight meridian turned about
    /// the axis.
    fn turned_cylinder() -> Surface {
        let meridian = [0.0, 1.0, 2.0].map(|z| Vec3::new(2.0, 0.0, z));
        Surface::Revolution(Box::new(Revolution {
            curve: Curve::BSpline(BSpline::interpolate(&meridian).unwrap()),
            origin: Vec3::ZERO,
            axis: Vec3::Z,
        }))
    }

    // A point at radius r + d from a cylinder's axis, between its ends and within the angle the
    // surface spans, lies |d| from it. A walk projects each point from the last one's foot, so it
    // must measure the same distance there whether the point is on the surface, a hair off it, or
    // too far off for the foot from a hint to be taken on trust.
    #[test]
    fn a_walk_measures_each_point_as_projecting_afresh_does() {
        for (s, r, turn) in [
            (quarter_cylinder(), 1.0, FRAC_PI_2),
            (turned_cylinder(), 2.0, TAU),
        ] {
            let mut walk = Walk::on(&s);
            for k in 0..=64 {
                let t = k as f64 / 64.0;
                let d = [0.0, 3e-8, -2e-7, 1e-3, -0.25][k % 5];
                let a = 0.02 + (turn - 0.04) * t;
                let p = Point3::new((r + d) * a.cos(), (r + d) * a.sin(), 0.1 + 1.8 * t);
                let gap = walk.gap(p);
                assert!((gap - d.abs()).abs() <= 1e-12, "{}: {k}: {gap:e}", s.kind());
                assert!(
                    (gap - surface_gap(&s, p)).abs() <= 1e-12,
                    "{}: {k}",
                    s.kind()
                );
            }
        }
    }
}
