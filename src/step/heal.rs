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
    tol,
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
    let (x, worst, reach, solved) = nearest_meeting(p0, on, reach);
    (solved && worst <= 0.25 * tol::LINEAR && x.dist(p0) <= reach).then_some(x)
}

/// The point [`solve_vertex`] finds, whether or not it is near enough: the point of its search
/// within `reach` of `p0` nearest all of `on`, how far it lies from the farthest of them, the
/// reach, and false where a step could not be solved for, which ends the search.
fn nearest_meeting(p0: Point3, on: &[On], reach: Option<f64>) -> (Point3, f64, f64, bool) {
    let mut x = p0;
    let mut reach = reach;
    let mut best = (f64::INFINITY, p0);
    let mut solved = true;
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
            let Ok(st) = solve_dense(ad, b.clone(), 3, 1) else {
                break None;
            };
            let st = Vec3::new(st[0], st[1], st[2]);
            if (x + st).dist(p0) <= reach || damp > trace {
                break Some(st);
            }
            damp *= 100.0;
        };
        let Some(step) = step else {
            solved = false;
            break;
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
    (x, worst, reach.unwrap_or(0.0), solved)
}

/// Below this sine of the angle between two faces' normals at a vertex they touch there rather
/// than cross, as [`adjust`] has it.
const TOUCHING: f64 = 0.1;

/// The point the surfaces `on` fix near `p0`, and how far it lies from the farthest of them, where
/// some of them touch there. Surfaces whose normals agree within [`TOUCHING`] are a stack of nearly
/// parallel sheets: to first order they fix the point only along their common normal, where the
/// point nearest the farthest of them lies midway between the outermost two (the centre of an
/// interval, which minimises the largest distance). Stacks that cross fix it along each of their
/// normals; the point is the least move from `p0` that meets every stack's middle, the
/// minimum-norm solution of at most three equations in three unknowns. More than three stacks
/// overdetermine the point, which is then no such point: those faces cross and miss.
fn touching_centre(p0: Point3, on: &[&Surface]) -> Option<(Point3, f64)> {
    let mut x = p0;
    let mut hints: Vec<Vec2> = on.iter().map(|s| s.project(x)).collect();
    let mut touching = false;
    for _ in 0..16 {
        // Each stack's unit normal and the least and greatest signed offset of `x` from its sheets.
        let mut stacks: Vec<(Vec3, f64, f64)> = Vec::new();
        for (s, h) in on.iter().zip(hints.iter_mut()) {
            *h = s.project_near(x, *h);
            let q = s.eval(h.x, h.y);
            let n = s.normal(h.x, h.y);
            if n.norm() < 0.5 {
                return None;
            }
            match stacks.iter_mut().find(|g| g.0.cross(n).norm() < TOUCHING) {
                Some(g) => {
                    let n = if g.0.dot(n) < 0.0 { -n } else { n };
                    let o = n.dot(x - q);
                    g.1 = g.1.min(o);
                    g.2 = g.2.max(o);
                    touching = true;
                }
                None => {
                    let o = n.dot(x - q);
                    stacks.push((n, o, o));
                }
            }
        }
        let m = stacks.len();
        if m > 3 {
            return None;
        }
        // δ = Nᵀ (N Nᵀ)⁻¹ c, with N the stacks' normals and c the moves onto their middles.
        let mut a = vec![0.0; m * m];
        let mut c = vec![0.0; m];
        for i in 0..m {
            for j in 0..m {
                a[i * m + j] = stacks[i].0.dot(stacks[j].0);
            }
            c[i] = -0.5 * (stacks[i].1 + stacks[i].2);
        }
        let y = solve_dense(a, c, m, 1).ok()?;
        let step = (0..m).fold(Vec3::ZERO, |d, i| d + stacks[i].0 * y[i]);
        x += step;
        if step.norm() <= 1e-3 * tol::LINEAR {
            break;
        }
    }
    if !touching {
        return None;
    }
    let worst = on.iter().map(|s| surface_gap(s, x)).fold(0.0, f64::max);
    Some((x, worst))
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
    (held, reach): (&[bool], f64),
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
    let (mut ps, mut pe) = (raw.points[edge.start], raw.points[edge.end]);
    let c = &edge.curve;
    // At a held vertex the edge ends on its own two faces, as near the vertex as they meet, within
    // `reach` of it: the vertex's tolerance, found after, must hold that end.
    let end_on = |p: Point3| {
        if same {
            return p;
        }
        touching_centre(p, &[a, b])
            .filter(|&(x, g)| g <= 0.25 * tol::LINEAR && x.dist(p) <= reach)
            .map(|(x, _)| x)
            .or_else(|| {
                let along = c.d1(param_of(c, p)).1;
                onto_both(a, b, false, p, along).filter(|x| x.dist(p) <= reach)
            })
            .unwrap_or(p)
    };
    if held[edge.start] {
        ps = end_on(ps);
    }
    if held[edge.end] {
        pe = if closed { ps } else { end_on(pe) };
    }
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
        let err: Vec<f64> = (0..n_seg)
            .map(|i| {
                [0.25, 0.5, 0.75]
                    .iter()
                    .map(|f| {
                        let m = approx.eval(params[i] + (params[i + 1] - params[i]) * f);
                        surface_gap(a, m).max(surface_gap(b, m))
                    })
                    .fold(0.0, f64::max)
            })
            .collect();
        if err.iter().all(|&x| x <= 0.5 * tol::LINEAR) {
            let dom = approx.domain();
            // A spline that stays on both faces everywhere is the edge as it is: evaluating it
            // costs nothing, where an intersection curve refines every point onto both faces.
            let on_both = (0..n_seg).all(|i| {
                (1..8).all(|k| {
                    let m = approx.eval(params[i] + (params[i + 1] - params[i]) * k as f64 / 8.0);
                    surface_gap(a, m).max(surface_gap(b, m)) <= 0.5 * tol::LINEAR
                })
            }) && (1..CHECK_SAMPLES).all(|k| {
                let m = approx.eval(dom.0 + (dom.1 - dom.0) * k as f64 / CHECK_SAMPLES as f64);
                surface_gap(a, m).max(surface_gap(b, m)) <= 0.5 * tol::LINEAR
            });
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

/// Moves each vertex that misses its faces onto them, starting from `start` where given, and
/// refits each edge that still misses. Returns each edge's range, or every vertex and edge that
/// would not fit. `loose` refits the edges at a vertex that would not move too, from where it is.
/// `hold`, the file's stated uncertainty, is how far apart a vertex or an edge whose faces touch
/// there may hold its edges' ends or its faces ([`touching_centre`], [`held_edge`]), and no
/// further than healing may move geometry ([`MAX_HEAL`]); without it, as for a file that states
/// none, nothing is held.
fn attempt(
    raw: &mut Raw,
    start: &HashMap<usize, Point3>,
    samples: &HashMap<usize, Vec<Point3>>,
    loose: bool,
    hold: Option<f64>,
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
    // Vertices whose faces touch there and have no common point: each edge ends on its own faces,
    // and the vertex's tolerance holds those ends.
    let mut held = vec![false; raw.points.len()];
    let hold = hold.map(|u| u.min(MAX_HEAL));
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
        // Where faces touch, the point midway between their sheets, within the reach a vertex
        // may move.
        let centre = || {
            let surfaces: Vec<&Surface> = faces.iter().map(|&f| &raw.faces[f].surface).collect();
            let gap = surfaces
                .iter()
                .map(|s| surface_gap(s, p0))
                .fold(0.0, f64::max);
            let reach = (MOVE_PER_GAP * (gap + tol::LINEAR)).min(MAX_MOVE);
            touching_centre(p0, &surfaces).filter(|&(x, _)| x.dist(p0) <= reach)
        };
        match p {
            Some(p) => moved[v] = p,
            None => match hold.and_then(|cap| Some((cap, centre()?))) {
                Some((cap, (x, r))) if r <= cap => {
                    moved[v] = x;
                    held[v] = true;
                }
                c => {
                    stuck[v] = true;
                    // The reason is the best point's distance from the faces, not the file's
                    // vertex's: where faces touch, the point midway between them; else the
                    // point nearest them all that the search found.
                    let why = match c {
                        Some((_, (_, r))) => format!(
                            "the {} faces at an end of this edge touch there and come no closer \
                             than {r:.1e} mm to one point",
                            faces.len()
                        ),
                        None => format!(
                            "the {} faces at an end of this edge have no common point near it (the \
                             nearest point found is {:.1e} mm off one)",
                            faces.len(),
                            nearest_meeting(p0, &on[..surfaces_only], None).1
                        ),
                    };
                    failures.push(Failure::Vertex(v, why));
                }
            },
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
                    || held[v]
                    || edge.curve.eval(t).dist(raw.points[v]) <= 0.25 * tol::LINEAR
            })
        })
    };

    let mut out = Vec::with_capacity(raw.edges.len());
    let mut held_edges: HashMap<usize, f64> = HashMap::new();
    let reach = hold.unwrap_or(tol::LINEAR);
    for (e, &old) in ranges.iter().enumerate() {
        let edge = &raw.edges[e];
        if !loose && (stuck[edge.start] || stuck[edge.end]) {
            out.push((0.0, 0.0));
            continue;
        }
        let r = match old {
            Some(old) => Some(range_near(raw, e, old)),
            None => range_of(raw, e),
        };
        if edge_fits(raw, e, r) && ends_close(raw, e, r) {
            out.push(r.unwrap_or((0.0, 0.0)));
            continue;
        }
        // Where the faces touch along the edge, the curve midway between them is the edge; a
        // refit could only find the same band, and fails where they lie apart.
        let touching = match hold {
            Some(cap) => held_edge(raw, e, r.or(old), cap),
            None => Err(None),
        };
        if let Ok((curve, r, t)) = touching {
            raw.edges[e].curve = curve;
            out.push(r);
            if t > tol::LINEAR {
                held_edges.insert(e, t);
            }
            continue;
        }
        // Touching sheets further apart than half the tolerance have no point within a quarter
        // of it of both, which is what a refit's points must be.
        if let Err(Some(gap)) = touching
            && gap > 0.25 * tol::LINEAR
        {
            failures.push(Failure::Edge(
                e,
                format!(
                    "its faces touch along it and lie {:.1e} mm apart",
                    2.0 * gap
                ),
            ));
            out.push((0.0, 0.0));
            continue;
        }
        let rf = match samples.get(&e) {
            Some(s) => refit(raw, e, r.or(old), Some(s), (&held, reach))
                .or_else(|_| refit(raw, e, r.or(old), None, (&held, reach))),
            None => refit(raw, e, r.or(old), None, (&held, reach)),
        };
        match rf {
            Ok((curve, r)) => {
                raw.edges[e].curve = curve;
                out.push(r);
            }
            Err(why) => {
                let held = match touching {
                    Err(Some(gap)) => {
                        format!(", and a curve midway between them misses them by {gap:.1e} mm")
                    }
                    _ => String::new(),
                };
                failures.push(Failure::Edge(e, format!("{why}{held}")));
                out.push((0.0, 0.0));
            }
        }
    }
    // A held vertex's tolerance is the widest gap between it and its edges' ends, widened as
    // [`tolerate`] widens a gap, and no more than the file's stated uncertainty or a quarter of an
    // edge it holds. Only a file that states one holds any.
    let cap = hold.unwrap_or(tol::LINEAR);
    raw.held.clear();
    for (v, edges) in at_vertex.iter().enumerate() {
        let widest = edges
            .iter()
            .filter_map(|e| held_edges.get(e))
            .fold(0.0, |a: f64, &b| a.max(b));
        if !held[v] && widest == 0.0 {
            continue;
        }
        let x = raw.points[v];
        let mut far: f64 = 0.0;
        let mut shortest = f64::INFINITY;
        for &e in edges {
            let edge = &raw.edges[e];
            let r = out[e];
            if r == (0.0, 0.0) {
                continue;
            }
            for (t, w) in [(r.0, edge.start), (r.1, edge.end)] {
                if w == v {
                    far = far.max(edge.curve.eval(t).dist(x));
                }
            }
            shortest = shortest.min(edge_length(&edge.curve, r));
        }
        let t = if far <= tol::LINEAR {
            tol::LINEAR
        } else {
            GAP_MARGIN * far
        }
        .max(widest);
        if t > cap || (t > tol::LINEAR && t > MAX_SHARE * shortest) {
            failures.push(Failure::Vertex(
                v,
                format!(
                    "the faces at an end of this edge touch there, and its edges end {far:.1e} mm \
                     apart on them"
                ),
            ));
        } else if t > tol::LINEAR {
            raw.held.insert(v, t);
        }
    }
    raw.held_edges = held_edges;
    if failures.is_empty() {
        Ok(out)
    } else {
        Err(failures)
    }
}

/// Edge `e` over `range` held with a tolerance where its two faces touch along it and no curve
/// lies on both within the kernel tolerance: the curve [`midway`] between them through its
/// vertices, and the widened largest gap from it to either face, if that is within `cap` and
/// a quarter of the edge; else that gap, where the faces touch. Faces that cross are left to the
/// surface step.
fn held_edge(
    raw: &Raw,
    e: usize,
    range: Option<(f64, f64)>,
    cap: f64,
) -> std::result::Result<(Curve, (f64, f64), f64), Option<f64>> {
    let edge = &raw.edges[e];
    let mut faces = edge.faces.clone();
    faces.sort_unstable();
    faces.dedup();
    let [fa, fb] = faces.as_slice() else {
        return Err(None);
    };
    let (a, b) = (&raw.faces[*fa].surface, &raw.faces[*fb].surface);
    let r = range.or_else(|| file_range(raw, e)).ok_or(None)?;
    let m = edge.curve.eval(0.5 * (r.0 + r.1));
    let (ua, ub) = (a.project(m), b.project(m));
    if a.normal(ua.x, ua.y).cross(b.normal(ub.x, ub.y)).norm() >= TOUCHING {
        return Err(None);
    }
    let c = midway_by(raw, e, r, true).ok_or(None)?;
    let dom = c.domain();
    let gap = curve_gap(&c, dom, a).max(curve_gap(&c, dom, b));
    let t = if gap <= tol::LINEAR {
        tol::LINEAR
    } else {
        GAP_MARGIN * gap
    };
    let length = edge_length(&c, dom);
    if t <= cap && (t <= tol::LINEAR || t <= MAX_SHARE * length) {
        Ok((c, dom, t))
    } else {
        Err(Some(gap))
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
    let with_surfaces = |surfaces: &[Surface]| {
        let mut r = file.clone();
        for (f, s) in r.faces.iter_mut().zip(surfaces) {
            f.surface = s.clone();
        }
        r
    };
    // Faces that touch are held within the file's stated uncertainty rather than moved: moving
    // their surfaces cannot make nearly parallel sheets meet without moving them by their gap. A
    // file that states none is healed as before.
    for round in 0..=ADJUST_ROUNDS {
        let mut work = with_surfaces(&surfaces);
        let failures = match attempt(&mut work, &start, &samples, false, file.stated) {
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
            let (id, why) = first_failure(&file, &failures);
            return fail(id, &why);
        }
        let mut shaped = with_surfaces(&surfaces);
        let mut pts = targets(&file, &chosen, &failing, &start);
        let faces: Vec<usize> = chosen.iter().copied().collect();
        let (at, met) = adjust(&mut shaped, &faces, &mut pts);
        settled = met;
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

/// The STEP id an error names for the first of `failures`, and why it failed.
fn first_failure(file: &Raw, failures: &[Failure]) -> (u64, String) {
    match &failures[0] {
        Failure::Vertex(v, why) => (
            file.edges
                .iter()
                .find(|e| e.start == *v || e.end == *v)
                .map_or(0, |e| e.id),
            why.clone(),
        ),
        Failure::Edge(e, why) => (file.edges[*e].id, why.clone()),
    }
}

/// Moves each vertex of `raw` onto its faces and refits each edge onto its faces wherever that
/// works, leaving the rest as the file has it and every surface where it is, for [`tolerate`] to
/// hold what is left: where faces meet nowhere exactly, as four faces at a corner can miss each
/// other by far less than the file's uncertainty, healing fails, but an edge elsewhere may still
/// be further off than the file promises until it is refitted. Returns how far it moved.
pub(super) fn heal_loosely(raw: &mut Raw) -> Option<Healed> {
    let file = raw.clone();
    let _ = attempt(raw, &HashMap::new(), &HashMap::new(), true, raw.stated);
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
    midway_by(raw, e, r, false)
}

/// [`midway`], with each sample's feet after the first found by walking from the last ones where
/// `walk` is set: the point moves by less each time, and a projection onto a spline from scratch
/// searches the whole surface.
fn midway_by(raw: &Raw, e: usize, r: (f64, f64), walk: bool) -> Option<Curve> {
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
                let (mut ua, mut ub) = (a.project(p), b.project(p));
                for k in 0..50 {
                    if k > 0 {
                        (ua, ub) = if walk {
                            (a.project_near(p, ua), b.project_near(p, ub))
                        } else {
                            (a.project(p), b.project(p))
                        };
                    }
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
    use crate::geom::Plane;
    use crate::math::Frame;

    fn plane(o: Point3, n: Vec3) -> Surface {
        Surface::Plane(Plane {
            frame: Frame::from_axis(o, n).unwrap(),
        })
    }

    #[test]
    fn touching_sheets_fix_a_point_midway_between_them_and_crossing_faces_fix_it_exactly() {
        // Planes z = 0 and z = d (one facing down) touch; x = 0.5 crosses them. The point is
        // midway between the sheets on the crossing plane, the least move from p0, d/2 from both.
        let d = 1e-7;
        let a = plane(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
        let b = plane(Point3::new(0.0, 0.0, d), Vec3::new(0.0, 0.0, -1.0));
        let c = plane(Point3::new(0.5, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let p0 = Point3::new(0.5 + 3e-8, 0.2, 0.9 * d);
        let (x, r) = touching_centre(p0, &[&a, &b, &c]).unwrap();
        assert!(x.dist(Point3::new(0.5, 0.2, 0.5 * d)) <= 1e-15, "{x:?}");
        assert!((r - 0.5 * d).abs() <= 1e-15, "{r}");
        // Three crossing planes are not touching faces: the solver of crossing faces finds them.
        let e = plane(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
        assert!(touching_centre(p0, &[&a, &c, &e]).is_none());
        // Four stacks overdetermine the point.
        let f = plane(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
        assert!(touching_centre(p0, &[&a, &b, &c, &e, &f]).is_none());
    }
}
