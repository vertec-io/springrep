//! STEP B-rep topology into a kernel solid: vertices and edges shared across faces, each face's
//! loops given pcurves that close in its parameter plane, with the seams and poles the kernel's
//! representation needs (doc/design.md, "Representation").
//!
//! A solid is read in three passes: the file's topology into [`Raw`], [`super::heal`] where the
//! file's geometry disagrees with itself by more than the kernel tolerance, then the arena.

use std::collections::{HashMap, HashSet};

use super::WRITER;
use super::geometry::Geometry;
use super::heal;
use super::part21::Param;
use crate::{
    Error, Result,
    boolean::classify,
    check,
    geom::{Circle, Cone, Curve, Line, Surface},
    math::{Frame, Point3, Transform, Vec2},
    tol,
    topo::{Brep, CoedgeSpec, EdgeId, FaceId, PCurve, ShellId, SolidId, VertexId, project_on},
};

fn bad<T>(id: u64, what: &str) -> Result<T> {
    Err(Error::InvalidInput(format!("STEP #{id}: {what}")))
}

fn unsupported<T>(id: u64, what: &str) -> Result<T> {
    Err(Error::Unsupported(format!("STEP #{id}: {what}")))
}

/// A face's use of an edge, in the kernel edge's direction or against it.
#[derive(Clone, Copy, Debug)]
pub(super) struct Use {
    pub edge: usize,
    pub reversed: bool,
}

/// An edge as the file gives it, turned to run along its curve: from `start` to `end`.
#[derive(Clone, Debug)]
pub(super) struct RawEdge {
    pub id: u64,
    pub curve: Curve,
    pub start: usize,
    pub end: usize,
    /// The faces that use it, once per use.
    pub faces: Vec<usize>,
}

#[derive(Clone, Debug)]
pub(super) struct RawFace {
    pub id: u64,
    pub name: String,
    pub surface: Surface,
    pub reversed: bool,
    /// Each bound: whether the file calls it the outer one, and its uses in loop order.
    pub bounds: Vec<(bool, Vec<Use>)>,
    /// The points its `VERTEX_LOOP`s hold inside it, where something touches it there.
    pub points: Vec<usize>,
}

/// One solid's topology and geometry as read, before any pcurve.
#[derive(Clone, Debug, Default)]
pub(super) struct Raw {
    pub points: Vec<Point3>,
    /// Each point's `VERTEX_POINT` instance.
    pub vertex_ids: Vec<u64>,
    pub edges: Vec<RawEdge>,
    pub faces: Vec<RawFace>,
    pub shells: Vec<Vec<usize>>,
}

/// The parameter of the point of `c` nearest `p`; a periodic curve's within its first period.
pub(super) fn param_of(c: &Curve, p: Point3) -> f64 {
    let (d0, d1) = c.domain();
    match c.period() {
        // A search over one period can stall at its ends, so a periodic curve is searched over
        // two periods half a period apart: the nearest point is inside at least one of them.
        Some(per) => {
            let t = [d0, d0 + 0.5 * per]
                .map(|a| c.project(p, a, a + per))
                .into_iter()
                .min_by(|&a, &b| c.eval(a).dist(p).total_cmp(&c.eval(b).dist(p)))
                .unwrap_or(d0);
            d0 + (t - d0).rem_euclid(per)
        }
        None if d0.is_finite() => c.project(p, d0, d1),
        None => c.project(p, -1e12, 1e12),
    }
}

/// The parameter range of an edge of `c` from `ps` to `pe`; `closed` when both are one vertex.
pub(super) fn edge_range(c: &Curve, ps: Point3, pe: Point3, closed: bool) -> Option<(f64, f64)> {
    let (d0, d1) = c.domain();
    let t0 = param_of(c, ps);
    let mut t1 = param_of(c, pe);
    match c.period() {
        Some(p) if closed => t1 = t0 + p,
        Some(p) => {
            t1 = t0 + (t1 - t0).rem_euclid(p);
            // A vertex at the join of the period projects to either end of it.
            if (t1 - t0).abs() <= tol::param(c.d1(t0).1.norm()) {
                t1 = t0 + p;
            }
        }
        None if closed => return None,
        None => {}
    }
    // A parameter whose point lies within `LINEAR` of an end of the curve is that end. The test is
    // on the points themselves, not on the parameter distance over the speed: where the curve
    // stands still at its end, as where its last two control points coincide, the speed there is
    // zero and the point moves only with the square of the parameter (Piegl & Tiller, §3.3), so a
    // first-order estimate of the parameter tolerance is unbounded and would take any parameter
    // for the other end.
    let snap = |t: f64| {
        if c.period().is_some() || !d0.is_finite() {
            return t;
        }
        let at = c.eval(t);
        let end = if t - d0 <= d1 - t { d0 } else { d1 };
        let mid = c.eval(0.5 * (t + end));
        if c.eval(end).dist(at) <= tol::LINEAR && c.eval(end).dist(mid) <= tol::LINEAR {
            end
        } else {
            t
        }
    };
    let (t0, t1) = (snap(t0), snap(t1));
    (t1 > t0).then_some((t0, t1))
}

/// `raw` without the edges split in two after they were read, which no face uses any more.
fn without_split(mut raw: Raw, splits: &HashMap<usize, [usize; 2]>) -> Raw {
    if splits.is_empty() {
        return raw;
    }
    let mut index = vec![usize::MAX; raw.edges.len()];
    let mut kept = Vec::with_capacity(raw.edges.len());
    for (i, e) in raw.edges.drain(..).enumerate() {
        if !splits.contains_key(&i) {
            index[i] = kept.len();
            kept.push(e);
        }
    }
    raw.edges = kept;
    for f in &mut raw.faces {
        for (_, uses) in &mut f.bounds {
            for u in uses {
                u.edge = index[u.edge];
            }
        }
    }
    raw
}

/// Whether the stretch of `c` between the points nearest `a` and `b` is no longer than
/// [`tol::LINEAR`], where `c` does not close: on a closed curve two coincident ends bound either
/// none of it or all of it, and the edge is read as all of it.
fn no_length(c: &Curve, a: Point3, b: Point3) -> bool {
    let (d0, d1) = c.domain();
    let closes = d0.is_finite() && d1.is_finite() && c.eval(d0).dist(c.eval(d1)) <= tol::LINEAR;
    if c.period().is_some() || closes {
        return false;
    }
    let (t0, t1) = (param_of(c, a), param_of(c, b));
    let at = |k: usize| c.eval(t0 + (t1 - t0) * k as f64 / 16.0);
    (0..16).map(|k| at(k).dist(at(k + 1))).sum::<f64>() <= tol::LINEAR
}

/// `raw` without the edges `null` marks, and their uses. An edge is kept where it is all of a
/// loop, which would be left with nothing to bound.
fn without_null(mut raw: Raw, null: &[bool]) -> Raw {
    let mut null = null.to_vec();
    for f in &raw.faces {
        for (_, uses) in &f.bounds {
            if uses.iter().all(|u| null[u.edge]) {
                for u in uses {
                    null[u.edge] = false;
                }
            }
        }
    }
    if !null.contains(&true) {
        return raw;
    }
    let mut index = vec![usize::MAX; raw.edges.len()];
    let mut kept = Vec::with_capacity(raw.edges.len());
    for (i, e) in raw.edges.drain(..).enumerate() {
        if !null[i] {
            index[i] = kept.len();
            kept.push(e);
        }
    }
    raw.edges = kept;
    for f in &mut raw.faces {
        for (_, uses) in &mut f.bounds {
            uses.retain(|u| !null[u.edge]);
            for u in uses {
                u.edge = index[u.edge];
            }
        }
    }
    raw
}

/// `raw` with its coincident vertices made one, and then its coincident edges. A solid that
/// touches itself is written, as STEP's manifold solids require, with each side of the contact
/// given its own vertices and edges; two vertices of one valid solid are never that close, so
/// coincident ones are such a contact, and joining them gives it back.
fn join_contacts(mut raw: Raw) -> Raw {
    let n = raw.points.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn root(p: &mut [usize], mut i: usize) -> usize {
        while p[i] != i {
            p[i] = p[p[i]];
            i = p[i];
        }
        i
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| raw.points[a].x.total_cmp(&raw.points[b].x));
    let mut joined = false;
    for (k, &i) in order.iter().enumerate() {
        for &j in &order[k + 1..] {
            if raw.points[j].x - raw.points[i].x > tol::LINEAR {
                break;
            }
            if raw.points[i].dist(raw.points[j]) <= tol::LINEAR {
                let (a, b) = (root(&mut parent, i), root(&mut parent, j));
                if a != b {
                    parent[a.max(b)] = a.min(b);
                    joined = true;
                }
            }
        }
    }
    // An edge whose two ends are now one vertex, on a curve that does not close, and no longer
    // than the kernel tolerance between them, is no edge: its ends were written apart only to be
    // joined. Each loop that uses it goes straight on from the edge before it to the edge after.
    let null: Vec<bool> = raw
        .edges
        .iter()
        .map(|e| {
            root(&mut parent, e.start) == root(&mut parent, e.end)
                && no_length(&e.curve, raw.points[e.start], raw.points[e.end])
        })
        .collect();
    if joined {
        // The points kept, renumbered in order.
        let mut index = vec![usize::MAX; n];
        let (mut points, mut ids) = (Vec::new(), Vec::new());
        for (i, slot) in index.iter_mut().enumerate() {
            if root(&mut parent, i) == i {
                *slot = points.len();
                points.push(raw.points[i]);
                ids.push(raw.vertex_ids[i]);
            }
        }
        let at: Vec<usize> = (0..n).map(|i| index[root(&mut parent, i)]).collect();
        raw.points = points;
        raw.vertex_ids = ids;
        for e in &mut raw.edges {
            e.start = at[e.start];
            e.end = at[e.end];
        }
        for f in &mut raw.faces {
            for v in &mut f.points {
                *v = at[*v];
            }
        }
    }
    let mut raw = without_null(raw, &null);
    // Edges between the same vertices that lie along each other: each later one is the first,
    // its uses turned where it runs the other way.
    let mut ends: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (i, e) in raw.edges.iter().enumerate() {
        ends.entry((e.start.min(e.end), e.start.max(e.end)))
            .or_default()
            .push(i);
    }
    let m = raw.edges.len();
    let mut into: Vec<(usize, bool)> = (0..m).map(|i| (i, false)).collect();
    let mut any = false;
    for group in ends.values().filter(|g| g.len() > 1) {
        for (k, &i) in group.iter().enumerate() {
            if into[i].0 != i {
                continue;
            }
            let a = &raw.edges[i];
            let (pa, pz) = (raw.points[a.start], raw.points[a.end]);
            let Some((t0, t1)) = edge_range(&a.curve, pa, pz, a.start == a.end) else {
                continue;
            };
            let (q, d) = a.curve.d1(0.5 * (t0 + t1));
            for &j in &group[k + 1..] {
                let b = &raw.edges[j];
                if into[j].0 != j {
                    continue;
                }
                // On `b`'s curve within its range: two arcs of one circle share their ends too.
                let (pb, pe) = (raw.points[b.start], raw.points[b.end]);
                let Some((u0, u1)) = edge_range(&b.curve, pb, pe, b.start == b.end) else {
                    continue;
                };
                let mut t = param_of(&b.curve, q);
                if let Some(p) = b.curve.period() {
                    t = u0 + (t - u0).rem_euclid(p);
                }
                let (r, dr) = b.curve.d1(t);
                if (u0..=u1).contains(&t) && r.dist(q) <= tol::LINEAR {
                    into[j] = (i, d.dot(dr) < 0.0);
                    any = true;
                }
            }
        }
    }
    if !any {
        return raw;
    }
    let mut kept = vec![usize::MAX; m];
    let mut edges = Vec::new();
    for i in 0..m {
        if into[i].0 == i {
            kept[i] = edges.len();
            edges.push(raw.edges[i].clone());
        }
    }
    for (i, &(to, _)) in into.iter().enumerate() {
        if to != i {
            let faces = raw.edges[i].faces.clone();
            edges[kept[to]].faces.extend(faces);
        }
    }
    raw.edges = edges;
    for f in &mut raw.faces {
        for (_, uses) in &mut f.bounds {
            for u in uses {
                let (to, flip) = into[u.edge];
                u.edge = kept[to];
                u.reversed ^= flip;
            }
        }
    }
    raw
}

/// The angles a run of edges covers about an axis, from `lo` to `hi` in turn, and `pad`, the
/// most angle that a length of [`tol::LINEAR`] round the axis subtends at any point of it: a line
/// within `pad` of the run passes within the linear tolerance of it, as near as two points are
/// taken to lie on one line.
struct Arc {
    lo: f64,
    hi: f64,
    pad: f64,
}

impl Arc {
    /// Whether the line at angle `a` passes beside the run, farther from it than the linear
    /// tolerance, neither through it nor along it.
    fn clear_of(&self, a: f64) -> bool {
        let from = (a - self.lo).rem_euclid(std::f64::consts::TAU);
        from > self.hi - self.lo + self.pad && from < std::f64::consts::TAU - self.pad
    }
}

/// The middle of the widest gap round an axis between the angles `points` and the arcs `arcs`:
/// the line there is farthest round from all of them. `None` where they leave no gap.
fn widest_gap(points: &[f64], arcs: &[Arc]) -> Option<f64> {
    const TAU: f64 = std::f64::consts::TAU;
    let mut spans: Vec<(f64, f64)> = points
        .iter()
        .map(|&a| a.rem_euclid(TAU))
        .map(|a| (a, a))
        .collect();
    spans.extend(arcs.iter().map(|h| {
        let lo = (h.lo - h.pad).rem_euclid(TAU);
        (lo, lo + h.hi - h.lo + 2.0 * h.pad)
    }));
    spans.sort_by(|x, y| x.0.total_cmp(&y.0));
    // Twice round, from the farthest any span reaches past a turn, so that every gap, the one
    // across angle 0 too, is one between neighbours.
    let laps = spans
        .iter()
        .copied()
        .chain(spans.iter().map(|&(lo, hi)| (lo + TAU, hi + TAU)));
    let mut end = spans.iter().map(|s| s.1).fold(f64::NEG_INFINITY, f64::max) - TAU;
    let mut best: Option<(f64, f64)> = None;
    for (lo, hi) in laps {
        let gap = lo - end;
        if gap > tol::ANGULAR && best.is_none_or(|b| gap > b.0) {
            best = Some((gap, end + 0.5 * gap));
        }
        end = end.max(hi);
    }
    best.map(|b| b.1)
}

/// The vertex a use starts at.
/// A point's angle about `frame`'s axis.
fn about(frame: &Frame) -> impl Fn(Point3) -> f64 + '_ {
    move |p| {
        let l = frame.to_local(p);
        l.y.atan2(l.x)
    }
}

fn use_start(raw: &Raw, u: &Use) -> usize {
    let e = &raw.edges[u.edge];
    if u.reversed { e.end } else { e.start }
}

/// Reads one solid's topology, sharing vertices and edges between its faces.
struct Reader<'g, 'a> {
    geo: &'g mut Geometry<'a>,
    raw: Raw,
    vertices: HashMap<u64, usize>,
    /// Each edge's index and whether it runs against the STEP edge's direction.
    edges: HashMap<u64, (usize, bool)>,
    /// A `POLY_LOOP`'s corners, by their coordinates, and its sides, by their two corners.
    corners: HashMap<[u64; 3], usize>,
    sides: HashMap<(usize, usize), usize>,
    /// Edges split in two after they were read, to give a seam a vertex: each one's two halves.
    splits: HashMap<usize, [usize; 2]>,
}

/// One bound of a face as the file gives it.
enum Bound {
    /// Whether the file calls it the outer bound, and its uses in loop order.
    Edges(bool, Vec<Use>),
    /// A `VERTEX_LOOP`: the face closes on this one vertex, the apex of a cone or the pole of a
    /// sphere, or, anywhere else, holds it as a point where something touches the face.
    Vertex(usize),
}

impl<'g, 'a> Reader<'g, 'a> {
    fn params(&self, id: u64, kinds: &[&str]) -> Result<(&'a str, &'a [Param])> {
        let Some(e) = self.geo.file.get(id) else {
            return bad(id, "missing instance");
        };
        for (k, p) in &e.parts {
            if kinds.contains(&k.as_str()) {
                return Ok((k.as_str(), p.as_slice()));
            }
        }
        unsupported(
            id,
            &format!("{} where {} was expected", e.kind(), kinds.join(" or ")),
        )
    }

    fn vertex(&mut self, id: u64) -> Result<usize> {
        if let Some(&v) = self.vertices.get(&id) {
            return Ok(v);
        }
        let (_, a) = self.params(id, &["VERTEX_POINT"])?;
        let p = self.geo.point(self.geo.reference(id, a.get(1))?)?;
        self.raw.points.push(p);
        self.raw.vertex_ids.push(id);
        let v = self.raw.points.len() - 1;
        self.vertices.insert(id, v);
        Ok(v)
    }

    fn edge(&mut self, id: u64) -> Result<(usize, bool)> {
        if let Some(&e) = self.edges.get(&id) {
            return Ok(e);
        }
        let (_, a) = self.params(id, &["EDGE_CURVE"])?;
        let vs = self.vertex(self.geo.reference(id, a.get(1))?)?;
        let ve = self.vertex(self.geo.reference(id, a.get(2))?)?;
        let flip = a.get(4).and_then(Param::as_bool) == Some(false);
        let (start, end) = if flip { (ve, vs) } else { (vs, ve) };
        let curve_id = self.geo.reference(id, a.get(3))?;
        let curve =
            match self
                .geo
                .conic_arc(curve_id, self.raw.points[start], self.raw.points[end])?
            {
                Some(c) => c,
                None => self.geo.curve(curve_id)?,
            };
        self.raw.edges.push(RawEdge {
            id,
            curve,
            start,
            end,
            faces: Vec::new(),
        });
        let e = (self.raw.edges.len() - 1, flip);
        self.edges.insert(id, e);
        Ok(e)
    }

    /// The vertex at a `POLY_LOOP` corner; corners at one place are one vertex.
    fn corner(&mut self, id: u64) -> Result<usize> {
        let p = self.geo.point(id)?;
        let key = [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
        if let Some(&v) = self.corners.get(&key) {
            return Ok(v);
        }
        self.raw.points.push(p);
        self.raw.vertex_ids.push(id);
        let v = self.raw.points.len() - 1;
        self.corners.insert(key, v);
        Ok(v)
    }

    /// The uses of a `POLY_LOOP`'s straight sides, each side shared with the face across it.
    fn poly_loop(&mut self, lp: u64, la: &[Param], forward: bool) -> Result<Vec<Use>> {
        let mut corners = Vec::new();
        for p in la.get(1).and_then(Param::as_list).unwrap_or(&[]) {
            corners.push(self.corner(self.geo.reference(lp, Some(p))?)?);
        }
        corners.dedup();
        while corners.len() > 1 && corners.first() == corners.last() {
            corners.pop();
        }
        if corners.len() < 3 {
            return bad(lp, "a polygon needs three distinct corners");
        }
        let mut uses = Vec::new();
        for k in 0..corners.len() {
            let (a, b) = (corners[k], corners[(k + 1) % corners.len()]);
            let edge = match self.sides.get(&(a.min(b), a.max(b))) {
                Some(&e) => e,
                None => {
                    let (pa, pb) = (self.raw.points[a], self.raw.points[b]);
                    let Some(dir) = (pb - pa).unit() else {
                        return bad(lp, "a polygon side has no length");
                    };
                    self.raw.edges.push(RawEdge {
                        id: lp,
                        curve: Curve::Line(Line { origin: pa, dir }),
                        start: a,
                        end: b,
                        faces: Vec::new(),
                    });
                    let e = self.raw.edges.len() - 1;
                    self.sides.insert((a.min(b), a.max(b)), e);
                    e
                }
            };
            let along = self.raw.edges[edge].start == a;
            uses.push(Use {
                edge,
                reversed: along != forward,
            });
        }
        if !forward {
            uses.reverse();
        }
        Ok(uses)
    }

    /// The edge uses of one `FACE_BOUND`, in the kernel's direction for the face.
    fn bound(&mut self, id: u64, flip_face: bool) -> Result<Bound> {
        let (kind, a) = self.params(id, &["FACE_OUTER_BOUND", "FACE_BOUND"])?;
        let outer = kind == "FACE_OUTER_BOUND";
        let lp = self.geo.reference(id, a.get(1))?;
        let orient = a.get(2).and_then(Param::as_bool) != Some(false);
        let (lk, la) = self.params(lp, &["EDGE_LOOP", "VERTEX_LOOP", "POLY_LOOP"])?;
        match lk {
            "VERTEX_LOOP" => {
                return Ok(Bound::Vertex(
                    self.vertex(self.geo.reference(lp, la.get(1))?)?,
                ));
            }
            "POLY_LOOP" => {
                return Ok(Bound::Edges(
                    outer,
                    self.poly_loop(lp, la, orient != flip_face)?,
                ));
            }
            _ => {}
        }
        let mut uses = Vec::new();
        for oe in la.get(1).and_then(Param::as_list).unwrap_or(&[]) {
            let oe = self.geo.reference(lp, Some(oe))?;
            let (_, p) = self.params(oe, &["ORIENTED_EDGE"])?;
            let (edge, flip) = self.edge(self.geo.reference(oe, p.get(3))?)?;
            let along = p.get(4).and_then(Param::as_bool) != Some(false);
            uses.push(Use {
                edge,
                reversed: ((along == orient) ^ flip_face) == flip,
            });
        }
        if !orient ^ flip_face {
            uses.reverse();
        }
        if uses.is_empty() {
            return bad(lp, "empty loop");
        }
        Ok(Bound::Edges(outer, self.expand(&uses)))
    }

    /// A seam from `from` to `apex` along `surface`, for a face that closes on its apex: a
    /// generator of a cone, a meridian of a sphere.
    fn seam(&mut self, id: u64, surface: &Surface, from: usize, apex: usize) -> Result<usize> {
        let (p, q) = (self.raw.points[from], self.raw.points[apex]);
        let curve = match surface {
            Surface::Cone(_) => match (q - p).unit() {
                Some(dir) => Curve::Line(Line { origin: p, dir }),
                None => return bad(id, "the vertex loop sits on its face's other bound"),
            },
            Surface::Sphere(s) => {
                let c = s.frame.origin;
                let x = p - c;
                let y = (q - c) - x.unit_or_zero() * (q - c).dot(x.unit_or_zero());
                match Frame::new(c, x.cross(y), x) {
                    Ok(frame) => Curve::Circle(Circle {
                        frame,
                        radius: s.radius,
                    }),
                    Err(_) => return bad(id, "the vertex loop is not at a pole of its sphere"),
                }
            }
            other => {
                return unsupported(id, &format!("VERTEX_LOOP bound on a {} face", other.kind()));
            }
        };
        self.raw.edges.push(RawEdge {
            id,
            curve,
            start: from,
            end: apex,
            faces: Vec::new(),
        });
        Ok(self.raw.edges.len() - 1)
    }

    /// `uses` with each use of a split edge as the uses of its halves.
    fn expand(&self, uses: &[Use]) -> Vec<Use> {
        let mut out = Vec::with_capacity(uses.len());
        for &u in uses {
            match self.splits.get(&u.edge) {
                Some(&[a, b]) => {
                    let halves = [Use { edge: a, ..u }, Use { edge: b, ..u }];
                    let ordered = if u.reversed {
                        [halves[1], halves[0]]
                    } else {
                        halves
                    };
                    out.extend(self.expand(&ordered));
                }
                None => out.push(u),
            }
        }
        out
    }

    /// Splits edge `e` at its parameter `t` into two edges meeting at a new vertex there, and
    /// gives every face read so far, and `bounds`, the halves in its place. Returns the vertex.
    fn split_edge(&mut self, e: usize, t: f64, bounds: &mut [(bool, Vec<Use>)]) -> usize {
        let edge = self.raw.edges[e].clone();
        let p = edge.curve.eval(t);
        self.raw.points.push(p);
        self.raw.vertex_ids.push(edge.id);
        let v = self.raw.points.len() - 1;
        self.raw.edges.push(RawEdge {
            end: v,
            ..edge.clone()
        });
        self.raw.edges.push(RawEdge { start: v, ..edge });
        let n = self.raw.edges.len();
        self.splits.insert(e, [n - 2, n - 1]);
        for f in 0..self.raw.faces.len() {
            for b in 0..self.raw.faces[f].bounds.len() {
                let uses = self.expand(&self.raw.faces[f].bounds[b].1);
                self.raw.faces[f].bounds[b].1 = uses;
            }
        }
        for (_, uses) in bounds.iter_mut() {
            *uses = self.expand(uses);
        }
        v
    }

    /// Where `uses` cross the angle `at`, as an edge and its parameter, `angle` giving a point's.
    fn crossing(
        &self,
        angle: &dyn Fn(Point3) -> f64,
        uses: &[Use],
        at: f64,
    ) -> Option<(usize, f64)> {
        let off = |p: Point3| {
            (angle(p) - at + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
                - std::f64::consts::PI
        };
        for u in uses {
            let e = &self.raw.edges[u.edge];
            let (ps, pe) = (self.raw.points[e.start], self.raw.points[e.end]);
            let (t0, t1) = edge_range(&e.curve, ps, pe, e.start == e.end)?;
            let n = 64;
            for k in 0..n {
                let (mut lo, mut hi) = (
                    t0 + (t1 - t0) * k as f64 / n as f64,
                    t0 + (t1 - t0) * (k + 1) as f64 / n as f64,
                );
                let (dl, dh) = (off(e.curve.eval(lo)), off(e.curve.eval(hi)));
                // A sign change away from the wrap at half a turn.
                if dl.signum() == dh.signum() || dl.abs().max(dh.abs()) > 1.0 {
                    continue;
                }
                for _ in 0..100 {
                    let mid = 0.5 * (lo + hi);
                    if off(e.curve.eval(mid)).signum() == dl.signum() {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                return Some((u.edge, 0.5 * (lo + hi)));
            }
        }
        None
    }

    /// The angle `uses` sweep about `frame`'s axis, from samples along each edge in its use's
    /// direction; `None` where an edge's span between its vertices cannot be found, or where the
    /// loop passes over the axis, as a seam through a pole does, and its angle jumps.
    fn swept(&self, frame: &Frame, uses: &[Use]) -> Option<f64> {
        self.swept_by(&about(frame), uses)
    }

    /// The angle `uses` sweep, `angle` giving a point's: [`Reader::swept`] about any measure of it.
    fn swept_by(&self, angle: &dyn Fn(Point3) -> f64, uses: &[Use]) -> Option<f64> {
        let mut total = 0.0;
        for u in uses {
            let e = &self.raw.edges[u.edge];
            let (ps, pe) = (self.raw.points[e.start], self.raw.points[e.end]);
            let (t0, t1) = edge_range(&e.curve, ps, pe, e.start == e.end)?;
            let n = 64;
            let mut ts: Vec<f64> = (0..=n)
                .map(|k| t0 + (t1 - t0) * k as f64 / n as f64)
                .collect();
            if u.reversed {
                ts.reverse();
            }
            for w in ts.windows(2) {
                let d = angle(e.curve.eval(w[1])) - angle(e.curve.eval(w[0]));
                let d = (d + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
                    - std::f64::consts::PI;
                if d.abs() > 1.0 {
                    return None;
                }
                total += d;
            }
        }
        Some(total)
    }

    /// The angles about an axis that `uses`, in order, cover, `angle` giving a point's and `pad`
    /// the angle the linear tolerance subtends there: from samples along each edge, as
    /// [`Reader::swept_by`] takes them. `None` where an edge's span cannot be found.
    fn arc(
        &self,
        angle: &dyn Fn(Point3) -> f64,
        pad: &dyn Fn(Point3) -> f64,
        uses: &[Use],
    ) -> Option<Arc> {
        let start = self.raw.points[use_start(&self.raw, uses.first()?)];
        let first = angle(start);
        let (mut at, mut lo, mut hi, mut most) = (first, first, first, pad(start));
        for u in uses {
            let e = &self.raw.edges[u.edge];
            let (ps, pe) = (self.raw.points[e.start], self.raw.points[e.end]);
            let (t0, t1) = edge_range(&e.curve, ps, pe, e.start == e.end)?;
            let n = 64;
            for i in 0..=n {
                let f = if u.reversed { n - i } else { i } as f64 / n as f64;
                let p = e.curve.eval(t0 + (t1 - t0) * f);
                let d = angle(p) - at;
                at += (d + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
                    - std::f64::consts::PI;
                (lo, hi, most) = (lo.min(at), hi.max(at), most.max(pad(p)));
            }
        }
        // Near the axis the tolerance subtends any angle: a line round it passes everything
        // there.
        let pad = most.min(std::f64::consts::PI);
        Some(Arc { lo, hi, pad })
    }

    /// What a seam of a face that closes round an axis must keep clear of, as the angles about
    /// the axis each covers, `angle` giving a point's: every bound but the two `rims` the seam
    /// joins, holes through the face, and every edge of those two that lies along one line of
    /// the surface, as a slot's side does, since a seam along it would be that edge a second
    /// time, to within the linear tolerance, `pad` giving the angle it subtends at a point.
    /// `None` where an edge's span cannot be found, or a hole goes round the axis.
    fn obstacles(
        &self,
        angle: &dyn Fn(Point3) -> f64,
        pad: &dyn Fn(Point3) -> f64,
        bounds: &[(bool, Vec<Use>)],
        rims: [usize; 2],
    ) -> Option<Vec<Arc>> {
        let mut out = Vec::new();
        for (k, (_, uses)) in bounds.iter().enumerate() {
            if !rims.contains(&k) {
                let hole = self.arc(angle, pad, uses)?;
                if hole.hi - hole.lo >= std::f64::consts::TAU {
                    return None;
                }
                out.push(hole);
                continue;
            }
            for u in uses {
                let edge = self.arc(angle, pad, std::slice::from_ref(u))?;
                if edge.hi - edge.lo <= edge.pad {
                    out.push(edge);
                }
            }
        }
        Some(out)
    }

    /// A face of revolution (cylinder, cone, sphere, torus) with two loops that each go once round
    /// its axis, the ways opposite, and no seam: those two become one loop, the first, a seam
    /// along the surface to a vertex of the second at the same angle, the second, and the seam
    /// back. Its pcurves then close in the parameter plane, as the kernel's faces do. Other loops,
    /// holes through the band, stay, and the seam runs clear of them: where the line through the
    /// two loops' vertices meets a hole, or runs along an edge of a hole or of the two loops, the
    /// seam would cross the hole or lie on that edge a second time, and the face's boundary would
    /// overlap itself. It then runs where it is farthest round from every hole and every vertex
    /// of the two loops, each loop's edge split there. `reversed` is the face's orientation
    /// against its surface's. `None` when the face is not such a band, or no line round it is
    /// clear of [`Reader::obstacles`].
    fn seamed_band(
        &mut self,
        id: u64,
        surface: &Surface,
        reversed: bool,
        bounds: &mut Vec<(bool, Vec<Use>)>,
    ) -> Option<()> {
        let frame = match surface {
            Surface::Cylinder(c) => c.frame,
            Surface::Cone(c) => c.frame,
            Surface::Sphere(s) => s.frame,
            Surface::Torus(t) => t.frame,
            _ => return None,
        };
        let full = |w: f64| (w.abs() - std::f64::consts::TAU).abs() < 0.1;
        let mut winding = Vec::new();
        for (k, (_, uses)) in bounds.iter().enumerate() {
            let w = self.swept(&frame, uses)?;
            if full(w) {
                winding.push((k, w));
            }
        }
        let &[(ka, wa), (kb, wb)] = winding.as_slice() else {
            return None;
        };
        if wa * wb >= 0.0 {
            return None;
        }
        let (a, b) = (&bounds[ka].1, &bounds[kb].1);
        // A vertex of each loop on one line of the surface: a cylinder's lines run along its
        // axis, a cone's through its apex.
        let on_line = |p: Point3, q: Point3| {
            let (lp, lq) = (frame.to_local(p), frame.to_local(q));
            let radial = |l: Point3| (l.x * l.x + l.y * l.y).sqrt();
            let (rp, rq) = (radial(lp), radial(lq));
            if rp <= tol::LINEAR || rq <= tol::LINEAR {
                return false;
            }
            let (ux, uy) = (lp.x / rp, lp.y / rp);
            (lq.x - ux * rq).hypot(lq.y - uy * rq) <= tol::LINEAR
        };
        let angle = about(&frame);
        // The angle a length of the linear tolerance round the axis subtends at a point, as
        // `on_line` measures how far a point lies off a line.
        let pad = |p: Point3| {
            let l = frame.to_local(p);
            tol::LINEAR / l.x.hypot(l.y)
        };
        let obstacles = self.obstacles(&angle, &pad, bounds, [ka, kb])?;
        let clear = |p: Point3| obstacles.iter().all(|h| h.clear_of(angle(p)));
        let at = |u: &Use| self.raw.points[use_start(&self.raw, u)];
        let aligned = (0..a.len())
            .flat_map(|i| (0..b.len()).map(move |j| (i, j)))
            .find(|&(i, j)| on_line(at(&a[i]), at(&b[j])) && clear(at(&a[i])));
        let (a, b, ia, ib) = match aligned {
            Some((i, j)) => (a.clone(), b.clone(), i, j),
            None if !clear(at(&a[0])) => {
                // The line through the first loop's vertex meets an obstacle: the seam runs in
                // the widest gap between the obstacles and the two loops' vertices, each loop
                // split there where it has no vertex on that line.
                let vertices: Vec<f64> = [ka, kb]
                    .iter()
                    .flat_map(|&k| bounds[k].1.iter())
                    .map(|u| angle(self.raw.points[use_start(&self.raw, u)]))
                    .collect();
                let mid = widest_gap(&vertices, &obstacles)?;
                let mut ends = [0usize; 2];
                for (slot, k) in [ka, kb].into_iter().enumerate() {
                    let on_mid = |u: &Use| {
                        let p = self.raw.points[use_start(&self.raw, u)];
                        let d = angle(p) - mid;
                        let d = (d + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
                            - std::f64::consts::PI;
                        d.abs() <= pad(p)
                    };
                    ends[slot] = match bounds[k].1.iter().find(|u| on_mid(u)) {
                        Some(u) => use_start(&self.raw, u),
                        None => {
                            let uses = bounds[k].1.clone();
                            let (e, t) = self.crossing(&angle, &uses, mid)?;
                            self.split_edge(e, t, bounds)
                        }
                    };
                }
                let (a, b) = (bounds[ka].1.clone(), bounds[kb].1.clone());
                let ia = a.iter().position(|u| use_start(&self.raw, u) == ends[0])?;
                let ib = b.iter().position(|u| use_start(&self.raw, u) == ends[1])?;
                (a, b, ia, ib)
            }
            None => {
                // No vertex of the second loop at the first's angle: split the second's edge
                // there, in every face that uses it.
                let va = use_start(&self.raw, &a[0]);
                let la = frame.to_local(self.raw.points[va]);
                let (e, t) = self.crossing(&about(&frame), b, la.y.atan2(la.x))?;
                let v = self.split_edge(e, t, bounds);
                let (a, b) = (bounds[ka].1.clone(), bounds[kb].1.clone());
                let ia = a.iter().position(|u| use_start(&self.raw, u) == va)?;
                let ib = b.iter().position(|u| use_start(&self.raw, u) == v)?;
                (a, b, ia, ib)
            }
        };
        let (va, vb) = (use_start(&self.raw, &a[ia]), use_start(&self.raw, &b[ib]));
        let (pa, pb) = (self.raw.points[va], self.raw.points[vb]);
        // The seam's line at the vertices' angle about the axis.
        let la = frame.to_local(pa);
        let radial = (frame.x * la.x + frame.y * la.y).unit()?;
        let meridian = |centre: Point3, radius: f64, toward_v: bool| {
            // A circle whose parameter is the surface's v there, or runs against it.
            let normal = radial.cross(frame.z) * if toward_v { 1.0 } else { -1.0 };
            Frame::new(centre, normal, radial)
                .ok()
                .map(|frame| Curve::Circle(Circle { frame, radius }))
        };
        let curve = match surface {
            Surface::Cylinder(_) | Surface::Cone(_) => Curve::Line(Line {
                origin: pa,
                dir: (pb - pa).unit()?,
            }),
            // A latitude band: the meridian between the two loops, which v crosses one way.
            Surface::Sphere(s) => {
                let v = |p: Point3| {
                    frame
                        .to_local(p)
                        .z
                        .atan2(frame.to_local(p).x.hypot(frame.to_local(p).y))
                };
                meridian(s.frame.origin, s.radius, v(pb) > v(pa))?
            }
            // A ring of the tube: the first loop's material lies on its left, seen from outside,
            // which is the way v grows from it when it runs the way u grows on a face along its
            // surface's normal.
            Surface::Torus(t) => {
                let toward_v = (wa > 0.0) != reversed;
                meridian(frame.origin + radial * t.major, t.minor, toward_v)?
            }
            _ => return None,
        };
        self.raw.edges.push(RawEdge {
            id,
            curve,
            start: va,
            end: vb,
            faces: Vec::new(),
        });
        let seam = self.raw.edges.len() - 1;
        let mut out: Vec<Use> = a[ia..].iter().chain(&a[..ia]).copied().collect();
        out.push(Use {
            edge: seam,
            reversed: false,
        });
        out.extend(b[ib..].iter().chain(&b[..ib]).copied());
        out.push(Use {
            edge: seam,
            reversed: true,
        });
        bounds.remove(ka.max(kb));
        bounds.remove(ka.min(kb));
        bounds.insert(0, (true, out));
        Some(())
    }

    /// A face closed round one of its surface's parameters, and not round the other, with two
    /// loops that each go once round it, the ways opposite, and no seam: [`Reader::seamed_band`]
    /// with that parameter for the angle, and the surface's own line of the other parameter for the
    /// seam: a spline's iso curve, a revolution's profile turned to that angle, an extrusion's
    /// line along its direction. `None` when the face is not such a band.
    ///
    /// A loop lifts to a closed path in the parameter plane exactly when it winds round the period
    /// zero times (Hatcher, *Algebraic Topology*, §1.3, Propositions 1.30 and 1.31). Each rim winds
    /// once, and the two the ways opposite, joined by the seam used once each way, wind zero times.
    fn seamed_parameter_band(
        &mut self,
        id: u64,
        surface: &Surface,
        bounds: &mut Vec<(bool, Vec<Use>)>,
    ) -> Option<()> {
        let (in_v, period) = match (surface.u_period(), surface.v_period()) {
            (None, Some(p)) => (true, p),
            (Some(p), None) => (false, p),
            _ => return None,
        };
        // The line of the other parameter where the closed one is `at`, as a curve whose own
        // parameter is the surface's along it.
        let line_at = |at: f64| -> Option<Curve> {
            match (surface, in_v) {
                (Surface::Nurbs(n), _) => Some(Curve::Nurbs(Box::new(n.iso_curve(!in_v, at)))),
                (Surface::Revolution(r), false) => Some(
                    r.curve
                        .transformed(&Transform::rotation(r.origin, r.axis, at)),
                ),
                (Surface::Extrusion(x), false) => Some(Curve::Line(Line {
                    origin: x.curve.eval(at),
                    dir: x.dir.unit()?,
                })),
                _ => None,
            }
        };
        // Only these surfaces have such a line; ask before any loop is looked at.
        line_at(0.0)?;
        const TAU: f64 = std::f64::consts::TAU;
        let angle = |p: Point3| {
            let uv = surface.project(p);
            (if in_v { uv.y } else { uv.x }) / period * TAU
        };
        let full = |w: f64| (w.abs() - TAU).abs() < 0.1;
        let mut winding = Vec::new();
        for (k, (_, uses)) in bounds.iter().enumerate() {
            let w = self.swept_by(&angle, uses)?;
            if full(w) {
                winding.push((k, w));
            }
        }
        let &[(ka, wa), (kb, wb)] = winding.as_slice() else {
            return None;
        };
        if wa * wb >= 0.0 {
            return None;
        }
        // The seam runs where it is farthest round from every vertex of the two loops and every
        // hole, the middle of the widest gap between their angles, so that splitting a loop for
        // it there leaves no sliver beside a vertex a file placed a little off the line, and the
        // seam neither crosses a hole nor lies along its edge.
        let angles: Vec<f64> = [ka, kb]
            .iter()
            .flat_map(|&k| bounds[k].1.iter())
            .map(|u| angle(self.raw.points[use_start(&self.raw, u)]))
            .collect();
        // The angle a length of the linear tolerance along the closed parameter subtends at a
        // point.
        let pad = |p: Point3| {
            let uv = surface.project(p);
            let (_, du, dv) = surface.d1(uv.x, uv.y);
            tol::LINEAR / ((if in_v { dv } else { du }).norm() * period / TAU)
        };
        let obstacles = self.obstacles(&angle, &pad, bounds, [ka, kb])?;
        let at = widest_gap(&angles, &obstacles)?;
        let seam = line_at(at / TAU * period)?;
        let on_seam = |raw: &Raw, u: &Use| {
            let q = raw.points[use_start(raw, u)];
            seam.eval(param_of(&seam, q)).dist(q) <= tol::LINEAR
        };
        // A vertex of each loop on the seam, the loop's edge split there where it has none.
        let mut ends = [0usize; 2];
        for (slot, k) in [ka, kb].into_iter().enumerate() {
            ends[slot] = match bounds[k].1.iter().find(|u| on_seam(&self.raw, u)) {
                Some(u) => use_start(&self.raw, u),
                None => {
                    let uses = bounds[k].1.clone();
                    let (e, t) = self.crossing(&angle, &uses, at)?;
                    self.split_edge(e, t, bounds)
                }
            };
        }
        let [va, vb] = ends;
        let (a, b) = (bounds[ka].1.clone(), bounds[kb].1.clone());
        let ia = a.iter().position(|u| use_start(&self.raw, u) == va)?;
        let ib = b.iter().position(|u| use_start(&self.raw, u) == vb)?;
        // An edge runs the way its curve does: from whichever end comes first along the iso.
        let back = param_of(&seam, self.raw.points[va]) > param_of(&seam, self.raw.points[vb]);
        let (start, end) = if back { (vb, va) } else { (va, vb) };
        self.raw.edges.push(RawEdge {
            id,
            curve: seam,
            start,
            end,
            faces: Vec::new(),
        });
        let edge = self.raw.edges.len() - 1;
        let mut out: Vec<Use> = a[ia..].iter().chain(&a[..ia]).copied().collect();
        out.push(Use {
            edge,
            reversed: back,
        });
        out.extend(b[ib..].iter().chain(&b[..ib]).copied());
        out.push(Use {
            edge,
            reversed: !back,
        });
        bounds.remove(ka.max(kb));
        bounds.remove(ka.min(kb));
        bounds.insert(0, (true, out));
        Some(())
    }

    /// A sphere or cone face with one loop that goes once round its axis and no seam: a cap over
    /// the pole on the loop's left seen from outside, or over the cone's apex, the only side that
    /// closes. As a vertex loop there would, it gives the loop a seam up to the pole and back. A
    /// loop through a pole of a sphere's frame, such as a hemisphere's great circle, is a latitude
    /// of the frame about its circle's normal, which then becomes the surface's frame. Other loops,
    /// holes through the cap, stay. False when the face is not such a cap.
    fn capped(
        &mut self,
        id: u64,
        surface: &mut Surface,
        reversed: bool,
        bounds: &mut [(bool, Vec<Use>)],
    ) -> Result<bool> {
        let (base, sphere) = match &*surface {
            Surface::Sphere(s) => (s.frame, true),
            Surface::Cone(c) if c.half_angle.abs() > tol::ANGULAR => (c.frame, false),
            _ => return Ok(false),
        };
        let mut found = None;
        for (k, (_, uses)) in bounds.iter().enumerate() {
            // A loop through the apex, or a pole, already closes on it by its own seam.
            let at_pole = uses.iter().any(|u| {
                let l = base.to_local(self.raw.points[use_start(&self.raw, u)]);
                l.x.hypot(l.y) <= tol::LINEAR
            });
            let frame = match uses.as_slice() {
                _ if !at_pole => base,
                [u] if sphere => match &self.raw.edges[u.edge].curve {
                    Curve::Circle(c) => match Frame::new(base.origin, c.frame.z, c.frame.x) {
                        Ok(f) => f,
                        Err(_) => continue,
                    },
                    _ => continue,
                },
                _ => continue,
            };
            let Some(w) = self.swept(&frame, uses) else {
                continue;
            };
            if (w.abs() - std::f64::consts::TAU).abs() < 0.1 {
                if found.is_some() {
                    return Ok(false);
                }
                found = Some((k, frame, w));
            }
        }
        let Some((k, frame, w)) = found else {
            return Ok(false);
        };
        let uses = &bounds[k].1;
        let pole = match surface {
            Surface::Sphere(s) => {
                s.frame = frame;
                // Along the outward normal, a loop running the way u grows has +z on its left.
                let up = (w > 0.0) != reversed;
                frame.origin + frame.z * if up { s.radius } else { -s.radius }
            }
            Surface::Cone(c) => c.frame.origin - c.frame.z * (c.radius / c.half_angle.tan()),
            _ => return Ok(false),
        };
        self.raw.points.push(pole);
        self.raw.vertex_ids.push(id);
        let apex = self.raw.points.len() - 1;
        let from = use_start(&self.raw, &uses[0]);
        let seam = self.seam(id, surface, from, apex)?;
        for reversed in [false, true] {
            bounds[k].1.push(Use {
                edge: seam,
                reversed,
            });
        }
        Ok(true)
    }

    /// Splits each edge of `bounds` that passes through a pole of `surface` between its vertices,
    /// at that pole. A pole is a whole line of the parameter plane, so a loop through one leaves
    /// it at another value of the free parameter than it arrived with (a meridian over a sphere's
    /// pole turns by half a turn there). Only at a vertex can the loop's parameters jump, along
    /// the pole's own line, as a degenerate edge does, and only a vertex there gives it one. An
    /// edge is split at the parameter nearest the pole, so one that passes through the same pole
    /// twice is split at one of the two only.
    fn split_at_poles(&mut self, surface: &Surface, bounds: &mut [(bool, Vec<Use>)]) {
        for pole in poles(surface) {
            loop {
                let mut found = None;
                'search: for (_, uses) in bounds.iter() {
                    for u in uses {
                        let e = &self.raw.edges[u.edge];
                        let (ps, pe) = (self.raw.points[e.start], self.raw.points[e.end]);
                        if ps.dist(pole) <= tol::LINEAR || pe.dist(pole) <= tol::LINEAR {
                            continue;
                        }
                        let Some((t0, t1)) = edge_range(&e.curve, ps, pe, e.start == e.end) else {
                            continue;
                        };
                        let t = param_of(&e.curve, pole);
                        let shifts = match e.curve.period() {
                            Some(p) => [0.0, p, -p],
                            None => [0.0; 3],
                        };
                        for d in shifts {
                            let t = t + d;
                            if t > t0 && t < t1 && e.curve.eval(t).dist(pole) <= tol::LINEAR {
                                found = Some((u.edge, t));
                                break 'search;
                            }
                        }
                    }
                }
                let Some((e, t)) = found else { break };
                self.split_edge(e, t, bounds);
            }
        }
    }

    fn face(&mut self, id: u64, flip_face: bool) -> Result<usize> {
        let (_, a) = self.params(id, &["ADVANCED_FACE", "FACE_SURFACE"])?;
        let name = a.first().and_then(Param::as_str).unwrap_or("").to_string();
        let mut surface = self.geo.surface(self.geo.reference(id, a.get(2))?)?;
        let same_sense = a.get(3).and_then(Param::as_bool) != Some(false);
        let mut bounds = Vec::new();
        let mut apexes = Vec::new();
        for bd in a.get(1).and_then(Param::as_list).unwrap_or(&[]) {
            match self.bound(self.geo.reference(id, Some(bd))?, flip_face)? {
                Bound::Edges(outer, uses) => bounds.push((outer, uses)),
                Bound::Vertex(v) => apexes.push(v),
            }
        }
        let (mut apexes, points): (Vec<usize>, Vec<usize>) = apexes
            .into_iter()
            .partition(|&v| closes_at(&surface, self.raw.points[v]));
        for &v in &points {
            let p = self.raw.points[v];
            let uv = surface.project(p);
            if surface.eval(uv.x, uv.y).dist(p) > tol::LINEAR {
                return bad(id, "a VERTEX_LOOP off its face's surface");
            }
        }
        // A whole sphere, bounded only by a vertex loop at one pole: the kernel's seam, a
        // meridian from the other pole to it, used both ways.
        if let ([apex], 0, Surface::Sphere(sp)) = (apexes.as_slice(), bounds.len(), &surface) {
            let (c, p) = (sp.frame.origin, self.raw.points[*apex]);
            let other = c * 2.0 - p;
            self.raw.points.push(other);
            self.raw.vertex_ids.push(id);
            let from = self.raw.points.len() - 1;
            let x = (other - c).unit_or_zero();
            let Ok(frame) = Frame::new(c, x.cross(sp.frame.x), x) else {
                return bad(id, "the vertex loop is not at a pole of its sphere");
            };
            self.raw.edges.push(RawEdge {
                id,
                curve: Curve::Circle(Circle {
                    frame,
                    radius: sp.radius,
                }),
                start: from,
                end: *apex,
                faces: Vec::new(),
            });
            let seam = self.raw.edges.len() - 1;
            bounds.push((
                true,
                vec![
                    Use {
                        edge: seam,
                        reversed: false,
                    },
                    Use {
                        edge: seam,
                        reversed: true,
                    },
                ],
            ));
            apexes.clear();
        }
        // A face that closes on a vertex is, in the kernel, its other bound followed by a seam
        // up to the vertex and back, the pole between the seam's two uses.
        match (apexes.as_slice(), bounds.len()) {
            ([], _) => {}
            (&[apex], 1) => {
                let uses = &bounds[0].1;
                let first = uses[0];
                let e = &self.raw.edges[first.edge];
                let from = if first.reversed { e.end } else { e.start };
                let seam = self.seam(id, &surface, from, apex)?;
                bounds[0].1.push(Use {
                    edge: seam,
                    reversed: false,
                });
                bounds[0].1.push(Use {
                    edge: seam,
                    reversed: true,
                });
            }
            _ => {
                return unsupported(
                    id,
                    "a vertex loop other than one face's apex beside one other bound",
                );
            }
        }
        if bounds.is_empty() {
            return unsupported(id, "face without bounds");
        }
        let reversed = !same_sense ^ flip_face;
        if self
            .seamed_band(id, &surface, reversed, &mut bounds)
            .is_none()
            && self
                .seamed_parameter_band(id, &surface, &mut bounds)
                .is_none()
            && apexes.is_empty()
        {
            self.capped(id, &mut surface, reversed, &mut bounds)?;
        }
        self.split_at_poles(&surface, &mut bounds);
        let f = self.raw.faces.len();
        for (_, uses) in &bounds {
            for u in uses {
                self.raw.edges[u.edge].faces.push(f);
            }
        }
        self.raw.faces.push(RawFace {
            id,
            name,
            surface,
            reversed: !same_sense ^ flip_face,
            bounds,
            points,
        });
        Ok(f)
    }

    /// Reads the shells of solid `id` into `raw`, and returns their instances as the solid lists
    /// them, outer first.
    fn solid(&mut self, id: u64) -> Result<Vec<u64>> {
        let Some(e) = self.geo.file.get(id) else {
            return bad(id, "missing instance");
        };
        let (kind, a) = (&e.parts[0].0, &e.parts[0].1);
        let mut shells: Vec<u64> = Vec::new();
        match kind.as_str() {
            "MANIFOLD_SOLID_BREP" | "FACETED_BREP" => {
                shells.push(self.geo.reference(id, a.get(1))?)
            }
            "BREP_WITH_VOIDS" => {
                shells.push(self.geo.reference(id, a.get(1))?);
                for v in a.get(2).and_then(Param::as_list).unwrap_or(&[]) {
                    shells.push(self.geo.reference(id, Some(v))?);
                }
            }
            other => return unsupported(id, &format!("solid {other}")),
        }
        for &listed in &shells {
            let (mut sh, mut flip) = (listed, false);
            if let Some(p) = self
                .geo
                .file
                .get(sh)
                .and_then(|e| e.part("ORIENTED_CLOSED_SHELL"))
            {
                flip = p.get(3).and_then(Param::as_bool) == Some(false);
                sh = self.geo.reference(sh, p.get(2))?;
            }
            let (_, p) = self.params(sh, &["CLOSED_SHELL"])?;
            let mut faces = Vec::new();
            for f in p.get(1).and_then(Param::as_list).unwrap_or(&[]) {
                faces.push(self.face(self.geo.reference(sh, Some(f))?, flip)?);
            }
            self.raw.shells.push(faces);
        }
        Ok(shells)
    }
}

/// One coedge's samples in the parameter plane, in its own direction.
#[derive(Clone, Debug)]
struct Track {
    ts: Vec<f64>,
    uvs: Vec<Vec2>,
}

/// Builds the arena of a read (and, where needed, healed) solid.
struct Builder {
    b: Brep,
    vertices: Vec<VertexId>,
    edges: Vec<EdgeId>,
}

impl Builder {
    fn new(raw: &Raw, ranges: &[(f64, f64)]) -> Builder {
        let mut b = Brep::new();
        let vertices: Vec<VertexId> = raw.points.iter().map(|&p| b.add_vertex(p)).collect();
        let edges = raw
            .edges
            .iter()
            .zip(ranges)
            .map(|(e, &r)| b.add_edge(e.curve.clone(), r, vertices[e.start], vertices[e.end]))
            .collect();
        Builder { b, vertices, edges }
    }

    fn face(&mut self, f: &RawFace) -> Result<()> {
        let surface = &f.surface;
        let want = if f.reversed { -1.0 } else { 1.0 };
        let single = f.bounds.len() == 1;
        let mut tracks: Vec<(bool, &[Use], Vec<Track>)> = Vec::new();
        for (outer, uses) in &f.bounds {
            let sign = (*outer || single).then_some(want);
            let t = self
                .track_loop(surface, uses, sign)
                .or_else(|| self.track_loop(surface, uses, None));
            let Some(t) = t else {
                return Err(Error::InvalidInput(format!(
                    "STEP #{}: a loop of this {} face does not close in its parameters",
                    f.id,
                    surface.kind()
                )));
            };
            tracks.push((*outer, uses, t));
        }
        let areas: Vec<f64> = tracks.iter().map(|(_, _, t)| area2(t)).collect();
        let outer = (0..tracks.len())
            .filter(|&i| areas[i] * want > 0.0)
            .max_by(|&i, &j| {
                let key = |k: usize| (tracks[k].0, areas[k].abs());
                key(i)
                    .partial_cmp(&key(j))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(0);
        tracks.swap(0, outer);
        let (lo, hi) = uv_box(&tracks[0].2);
        for (_, _, t) in tracks.iter_mut().skip(1) {
            shift_into(surface, t, lo, hi);
        }
        let face = self.b.add_face(surface.clone(), f.reversed, []);
        for (_, uses, ts) in &tracks {
            let specs = self.specs(surface, uses, ts);
            self.b.add_loop(face, specs);
        }
        for &p in &f.points {
            let v = self.vertices[p];
            let mut at = [Track {
                ts: vec![0.0],
                uvs: vec![surface.project(self.b.vertex(v).point)],
            }];
            shift_into(surface, &mut at, lo, hi);
            let edge = self.b.add_degenerate_edge(v, (0.0, 1.0));
            let spec = CoedgeSpec {
                edge,
                reversed: false,
                pcurve: PCurve::Affine {
                    origin: at[0].uvs[0],
                    dir: Vec2::ZERO,
                },
            };
            self.b.add_loop(face, vec![spec]);
        }
        Ok(())
    }

    fn specs(&mut self, s: &Surface, uses: &[Use], tracks: &[Track]) -> Vec<CoedgeSpec> {
        let mut out = Vec::new();
        let n = uses.len();
        for k in 0..n {
            let (edge, reversed) = (self.edges[uses[k].edge], uses[k].reversed);
            let tr = &tracks[k];
            let pcurve = match (s, &self.b.edge(edge).curve) {
                (Surface::Plane(_), Some(Curve::Line(_))) => {
                    let (t0, t1) = (tr.ts[0], tr.ts[tr.ts.len() - 1]);
                    let (a, z) = (tr.uvs[0], tr.uvs[tr.uvs.len() - 1]);
                    let dir = (z - a) / (t1 - t0);
                    PCurve::Affine {
                        origin: a - dir * t0,
                        dir,
                    }
                }
                _ => PCurve::Projected { seed: tr.uvs[0] },
            };
            out.push(CoedgeSpec {
                edge,
                reversed,
                pcurve,
            });
            let end = tr.uvs[tr.uvs.len() - 1];
            let next = tracks[(k + 1) % n].uvs[0];
            if gap_is_pole(s, end, next) {
                let e = self.b.edge(edge);
                let v = if reversed { e.start } else { e.end };
                let d = self.b.add_degenerate_edge(v, (0.0, 1.0));
                out.push(CoedgeSpec {
                    edge: d,
                    reversed: false,
                    pcurve: PCurve::Affine {
                        origin: end,
                        dir: next - end,
                    },
                });
            }
        }
        out
    }

    /// The 3-D samples of a use of an edge, in its direction.
    fn samples(&self, u: Use) -> (Vec<f64>, Vec<Point3>) {
        let e = self.b.edge(self.edges[u.edge]);
        let c = e.curve.as_ref().expect("STEP edges have curves");
        let (t0, t1) = e.range;
        let mut cuts = vec![t0];
        cuts.extend(c.breakpoints(t0, t1));
        cuts.push(t1);
        let per = match c {
            Curve::Line(_) => 2,
            _ => (16 / (cuts.len() - 1)).max(4),
        };
        let mut ts = Vec::new();
        for w in cuts.windows(2) {
            for k in 0..per {
                ts.push(w[0] + (w[1] - w[0]) * k as f64 / per as f64);
            }
        }
        ts.push(t1);
        if u.reversed {
            ts.reverse();
        }
        let pts = ts.iter().map(|&t| c.eval(t)).collect();
        (ts, pts)
    }

    /// Each use's parameter-plane samples, continuous along the loop and closing back on its
    /// start; `None` if no choice of periods closes it. Where a loop leaves a pole the period it
    /// leaves in is free, and `sign`, when given, picks the choice that encloses area that way.
    fn track_loop(&self, s: &Surface, uses: &[Use], sign: Option<f64>) -> Option<Vec<Track>> {
        let samples: Vec<(Vec<f64>, Vec<Point3>)> = uses.iter().map(|&u| self.samples(u)).collect();
        let start = s.project(samples[0].1[0]);
        // A search over each use's choices in turn, backing up where a later use finds none, kept
        // on a stack of its own: a loop of hundreds of edges, as a gear's outline is, would
        // otherwise nest a call per edge, past the 1 MiB a WebAssembly program's stack holds.
        let mut acc: Vec<Track> = Vec::with_capacity(uses.len());
        let mut choices = vec![self.choices(s, uses, &samples, &acc, start).into_iter()];
        while let Some(next) = choices.last_mut() {
            let Some(track) = next.next() else {
                // Every choice for this use failed: back up to the one before.
                choices.pop();
                acc.pop();
                continue;
            };
            let end = track.uvs[track.uvs.len() - 1];
            acc.push(track);
            let k = acc.len();
            if k < uses.len() {
                choices.push(self.choices(s, uses, &samples, &acc, end).into_iter());
                continue;
            }
            let (last, first) = (acc[k - 1].uvs[acc[k - 1].uvs.len() - 1], acc[0].uvs[0]);
            let closed = closes(s, last, first) || gap_is_pole(s, last, first);
            if closed && sign.is_none_or(|w| area2(&acc) * w > 0.0) {
                return Some(acc);
            }
            acc.pop();
        }
        None
    }

    /// The samples use `acc.len()` may take, starting near `hint` where the uses before it end, in
    /// the order they are tried: those that are not a slit, then those that are.
    fn choices(
        &self,
        s: &Surface,
        uses: &[Use],
        samples: &[(Vec<f64>, Vec<Point3>)],
        acc: &[Track],
        hint: Vec2,
    ) -> Vec<Track> {
        let k = acc.len();
        let (ts, pts) = &samples[k];
        let first = s.project_near(pts[0], hint);
        let free = pole_free(s, first);
        let mut starts = vec![first];
        if free != (false, false) && pts.len() > 1 {
            // A use leaving a pole may leave at any angle: take the one its next point implies,
            // in each period near the hint.
            let inner = s.project_near(pts[1], hint);
            starts.clear();
            let shifts = |on: bool, p: Option<f64>| match p {
                Some(p) if on => vec![0.0, -p, p],
                _ => vec![0.0],
            };
            for du in shifts(free.0, s.u_period()) {
                for dv in shifts(free.1, s.v_period()) {
                    let mut c = first;
                    if free.0 {
                        c.x = inner.x + du;
                    }
                    if free.1 {
                        c.y = inner.y + dv;
                    }
                    starts.push(c);
                }
            }
            starts.sort_by(|a, b| a.dist(hint).total_cmp(&b.dist(hint)));
        }
        let curve = self.b.edge(self.edges[uses[k].edge]).curve.as_ref();
        let (mut plain, mut slits) = (Vec::new(), Vec::new());
        for st in starts {
            let mut uvs = Vec::with_capacity(pts.len());
            uvs.push(st);
            for p in &pts[1..] {
                let prev = uvs[uvs.len() - 1];
                uvs.push(s.project_near(*p, prev));
            }
            let n = uvs.len();
            // At a pole the free parameter is the one the edge leaves it by, found as the
            // checker and every later reader of the pcurve finds it.
            if let Some(c) = curve.filter(|_| n > 1) {
                if pole_free(s, uvs[0]) != (false, false) {
                    uvs[0] = project_on(s, c, ts[0], ts[1], st);
                }
                if pole_free(s, uvs[n - 1]) != (false, false) {
                    uvs[n - 1] = project_on(s, c, ts[n - 1], ts[n - 2], uvs[n - 2]);
                }
            }
            // A use running back over an earlier use of its edge is a slit, which a face may hold
            // where something touched it; a seam's second use runs a period away, so that is tried
            // first.
            let slit = coincides_with_twin(uses, acc, k, &uvs);
            let track = Track {
                ts: ts.clone(),
                uvs,
            };
            if slit {
                slits.push(track);
            } else {
                plain.push(track);
            }
        }
        plain.extend(slits);
        plain
    }
}

/// Whether an earlier use of the same edge in this loop runs back over the same parameters: a slit,
/// where a seam's would run a period away.
fn coincides_with_twin(uses: &[Use], acc: &[Track], k: usize, uvs: &[Vec2]) -> bool {
    let n = uvs.len();
    (0..k).any(|j| {
        uses[j].edge == uses[k].edge && {
            let other = &acc[j].uvs;
            other.len() == n && other[n / 2].dist(uvs[n - 1 - n / 2]) < 1e-6
        }
    })
}

/// Whether `uv` is on a pole of the surface: `(u free, v free)` where the whole line of that
/// parameter through `uv` is one point. The derivative is the filter, [`Surface::line_is_point`]
/// the test: along a line where the surface only stands still the derivative vanishes too, and
/// the line is a curve that the loop runs along as along any other.
fn pole_free(s: &Surface, uv: Vec2) -> (bool, bool) {
    let (_, su, sv) = s.d1(uv.x, uv.y);
    (
        su.norm() <= tol::LINEAR && s.line_is_point(uv, false),
        sv.norm() <= tol::LINEAR && s.line_is_point(uv, true),
    )
}

/// The poles of `s` through which an edge of a face on it may pass between its vertices: a
/// sphere's two, `origin ± radius·z`. The other surfaces with poles are not split for: an edge
/// through a cone's apex passes from one nappe to the other, and the reader keeps a cone face on
/// one; loops through the pole of a surface of revolution, a spindle torus or a spline are not yet
/// tracked there.
fn poles(s: &Surface) -> Vec<Point3> {
    match s {
        Surface::Sphere(sp) => {
            let (o, z) = (sp.frame.origin, sp.frame.z);
            vec![o + z * sp.radius, o - z * sp.radius]
        }
        _ => Vec::new(),
    }
}

/// Whether `a` and `b` are the same point of the parameter plane, not a period apart.
fn closes(s: &Surface, a: Vec2, b: Vec2) -> bool {
    let near = |d: f64, p: Option<f64>| p.is_none_or(|p| d.abs() < 0.25 * p);
    near(a.x - b.x, s.u_period()) && near(a.y - b.y, s.v_period())
}

/// Whether the step from `a` to `b` is a pole's span: one point in space, apart in parameters.
fn gap_is_pole(s: &Surface, a: Vec2, b: Vec2) -> bool {
    let free = pole_free(s, a);
    if free == (false, false) {
        return false;
    }
    let (du, dv) = ((a.x - b.x).abs(), (a.y - b.y).abs());
    let (moved, held) = (1e-9, 1e-6);
    (free.0 && du > moved && dv <= held) || (free.1 && dv > moved && du <= held)
}

fn area2(tracks: &[Track]) -> f64 {
    let pts: Vec<Vec2> = tracks
        .iter()
        .flat_map(|t| t.uvs[..t.uvs.len() - 1].iter().copied())
        .collect();
    (0..pts.len())
        .map(|i| pts[i].cross(pts[(i + 1) % pts.len()]))
        .sum()
}

/// Whether a face on `s` bounded by a `VERTEX_LOOP` at `p` closes there: at a cone's apex or a
/// sphere's pole, where every parameter line meets.
fn closes_at(s: &Surface, p: Point3) -> bool {
    let near = |q: Point3| q.dist(p) <= tol::LINEAR;
    match s {
        Surface::Cone(c) => near(c.frame.origin - c.frame.z * (c.radius / c.half_angle.tan())),
        Surface::Sphere(sp) => {
            let (o, z) = (sp.frame.origin, sp.frame.z);
            near(o + z * sp.radius) || near(o - z * sp.radius)
        }
        _ => false,
    }
}

fn uv_box(tracks: &[Track]) -> (Vec2, Vec2) {
    let mut lo = Vec2::new(f64::INFINITY, f64::INFINITY);
    let mut hi = Vec2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in tracks.iter().flat_map(|t| t.uvs.iter()) {
        lo = Vec2::new(lo.x.min(p.x), lo.y.min(p.y));
        hi = Vec2::new(hi.x.max(p.x), hi.y.max(p.y));
    }
    (lo, hi)
}

/// Moves a hole's samples by whole periods so its centre lies inside the outer loop's box.
fn shift_into(s: &Surface, tracks: &mut [Track], lo: Vec2, hi: Vec2) {
    let (a, b) = uv_box(tracks);
    let c = (a + b) * 0.5;
    let shift = |x: f64, lo: f64, hi: f64, p: Option<f64>| -> f64 {
        match p {
            Some(p) if x < lo || x > hi => ((0.5 * (lo + hi) - x) / p).round() * p,
            _ => 0.0,
        }
    };
    let d = Vec2::new(
        shift(c.x, lo.x, hi.x, s.u_period()),
        shift(c.y, lo.y, hi.y, s.v_period()),
    );
    if d.x != 0.0 || d.y != 0.0 {
        for t in tracks {
            for uv in &mut t.uvs {
                *uv = *uv + d;
            }
        }
    }
}

/// A solid read, its faces' STEP ids and names in order, how far healing moved the geometry of
/// each of its solids it moved, by the solid's position, the instances of the void shells that
/// lie outside its outer shell and were read as solids of their own, in the order the solids follow
/// the first, and the faces whose bounds it read against the file's flags ([`agree_senses`]): each
/// face's position among its faces and how many of its bounds.
pub(super) type ReadSolid = (
    Brep,
    Vec<(u64, String)>,
    Vec<(usize, heal::Healed)>,
    Vec<u64>,
    Vec<(usize, usize)>,
);

/// A solid read from some of a solid's shells, its faces' STEP ids and names in order, how far
/// healing moved its geometry, and the faces whose bounds it read against the file's flags.
type Shells = (
    Brep,
    Vec<(u64, String)>,
    Option<heal::Healed>,
    Vec<(usize, usize)>,
);

/// The larger of two movements.
fn farther(a: Option<heal::Healed>, b: Option<heal::Healed>) -> Option<heal::Healed> {
    match (a, b) {
        (Some(a), Some(b)) => Some(if b.moved > a.moved {
            heal::Healed {
                surfaces: a.surfaces.max(b.surfaces),
                ..b
            }
        } else {
            heal::Healed {
                surfaces: a.surfaces.max(b.surfaces),
                ..a
            }
        }),
        (a, b) => a.or(b),
    }
}

/// Each cone face on its cone's own nappe. A file may place a face on the nappe past the apex from
/// its cone's placement, which no cone the kernel makes does and which its intersections and
/// projections onto a face's cone do not follow. The same cone about its axis turned over, from the
/// point as far past the apex as the placement's origin is before it, holds that face on its own
/// nappe; the face turns over where the natural normal there turns with it. A face with points on
/// both nappes is an error: it would pass through the apex, which no face does.
fn own_nappes(mut raw: Raw) -> Result<Raw> {
    for f in 0..raw.faces.len() {
        let Surface::Cone(c) = &raw.faces[f].surface else {
            continue;
        };
        let (sa, ca) = c.half_angle.sin_cos();
        if sa.abs() <= tol::ANGULAR || ca.abs() <= tol::ANGULAR {
            continue;
        }
        // The face's point furthest from the apex says which nappe it lies on: a face reaches the
        // apex at most at a point.
        let mut pts = Vec::new();
        for (_, uses) in &raw.faces[f].bounds {
            for u in uses {
                let e = &raw.edges[u.edge];
                let (ps, pe) = (raw.points[e.start], raw.points[e.end]);
                pts.extend([ps, pe]);
                if let Some((t0, t1)) = edge_range(&e.curve, ps, pe, e.start == e.end) {
                    pts.push(e.curve.eval(0.5 * (t0 + t1)));
                }
            }
        }
        let radius_at = |q: Point3| c.radius + c.frame.to_local(q).z / ca * sa;
        let sides = |below: bool| {
            pts.iter().any(|&q| {
                if below {
                    radius_at(q) < -tol::LINEAR
                } else {
                    radius_at(q) > tol::LINEAR
                }
            })
        };
        if sides(true) && sides(false) {
            return bad(raw.faces[f].id, "a face lies on both nappes of its cone");
        }
        let Some(far) = pts
            .iter()
            .copied()
            .max_by(|&a, &b| radius_at(a).abs().total_cmp(&radius_at(b).abs()))
        else {
            continue;
        };
        if radius_at(far) >= -tol::LINEAR {
            continue;
        }
        let apex = c.frame.origin + c.frame.z * (-c.radius / sa * ca);
        let Ok(frame) = Frame::new(apex * 2.0 - c.frame.origin, -c.frame.z, c.frame.x) else {
            continue;
        };
        let own = Surface::Cone(Cone {
            frame,
            radius: c.radius,
            half_angle: c.half_angle,
        });
        let old = &raw.faces[f].surface;
        let (a, b) = (old.project(far), own.project(far));
        let turns = old.normal(a.x, a.y).dot(own.normal(b.x, b.y)) < 0.0;
        let face = &mut raw.faces[f];
        face.surface = own;
        face.reversed ^= turns;
    }
    Ok(raw)
}

/// Whether either sense of a loop on `s` bounds a region of it, as on a surface closed in both
/// directions: a torus, or a sphere, closed round its axis and at its poles.
fn closed_both_ways(s: &Surface) -> bool {
    matches!(s, Surface::Sphere(_)) || (s.u_period().is_some() && s.v_period().is_some())
}

/// Moves the cap seam of bound `k` of face `f`, just turned, to the other pole. A sphere face
/// with one loop round its axis and no seam is read as the cap on the loop's left, given a
/// meridian up to that cap's pole and back ([`Reader::capped`]); turned, the loop has the other
/// cap on its left, so the meridian must run to the antipode instead. The seam is told from a
/// `VERTEX_LOOP`'s, whose pole is the file's vertex and stays, by its pole: a point the read made
/// for the face, which no other edge uses. The bound is rotated so the seam's two uses close it,
/// as the read lays them. Nothing else changes: other bounds and surfaces have no seam of the read's.
fn recap(raw: &mut Raw, f: usize, k: usize) {
    let Surface::Sphere(s) = &raw.faces[f].surface else {
        return;
    };
    let (c, z) = (s.frame.origin, s.frame.z);
    let bound = &raw.faces[f].bounds[k].1;
    let Some(i) = (0..bound.len()).find(|&i| {
        let e = bound[i].edge;
        let pole = raw.edges[e].end;
        bound.iter().filter(|u| u.edge == e).count() == 2
            && raw.vertex_ids[pole] == raw.faces[f].id
            && (raw.points[pole] - c).cross(z).norm() <= tol::LINEAR
            && raw
                .edges
                .iter()
                .enumerate()
                .all(|(j, x)| j == e || (x.start != pole && x.end != pole))
    }) else {
        return;
    };
    let e = bound[i].edge;
    let (from, pole) = (raw.edges[e].start, raw.edges[e].end);
    let q = c * 2.0 - raw.points[pole];
    // A meridian from the loop's vertex to the new pole, as `Reader::seam` draws one.
    let x = raw.points[from] - c;
    let y = (q - c) - x.unit_or_zero() * (q - c).dot(x.unit_or_zero());
    let Ok(frame) = Frame::new(c, x.cross(y), x) else {
        return;
    };
    let radius = s.radius;
    raw.points[pole] = q;
    raw.edges[e].curve = Curve::Circle(Circle { frame, radius });
    // The seam's uses go last: out along it from `from`, and back.
    let bound = &mut raw.faces[f].bounds[k].1;
    if let Some(j) = bound.iter().position(|u| u.edge == e && !u.reversed)
        && bound[(j + 1) % bound.len()].edge == e
    {
        let n = bound.len();
        bound.rotate_left((j + 2) % n);
    }
}

/// Takes each bound of a face on a surface closed both ways in the sense the faces across its
/// edges agree with. On such a surface both senses of a loop bound a region of it, the face or the
/// rest of the surface, so the `FACE_BOUND` flag alone says which the face is, and some files flag
/// a loop against the way they list it. In a closed shell the two uses of an edge run opposite
/// ways; where every use of a bound with a face across it runs the same way as that face's, the
/// bound is turned. A bound settles once it has a settled face across one of its edges, every
/// bound of any other face being settled as written, so the sense passes from face to face.
/// Returns each face with bounds turned and how many.
fn agree_senses(raw: &mut Raw) -> Vec<(usize, usize)> {
    let mut settled: Vec<Vec<bool>> = raw
        .faces
        .iter()
        .map(|f| vec![!closed_both_ways(&f.surface); f.bounds.len()])
        .collect();
    if settled.iter().flatten().all(|&s| s) {
        return Vec::new();
    }
    // Each edge's uses: the face, its bound and whether the use runs against the edge.
    let mut uses: Vec<Vec<(usize, usize, bool)>> = vec![Vec::new(); raw.edges.len()];
    for (f, face) in raw.faces.iter().enumerate() {
        for (k, (_, bound)) in face.bounds.iter().enumerate() {
            for u in bound {
                uses[u.edge].push((f, k, u.reversed));
            }
        }
    }
    let mut turned = vec![0; raw.faces.len()];
    let mut changed = true;
    while changed {
        changed = false;
        for f in 0..raw.faces.len() {
            for k in 0..raw.faces[f].bounds.len() {
                if settled[f][k] {
                    continue;
                }
                let (mut with, mut against) = (0, 0);
                for u in &raw.faces[f].bounds[k].1 {
                    // Only an edge two faces use has a face across it to agree with.
                    let &[a, b] = &uses[u.edge][..] else {
                        continue;
                    };
                    let (mine, other) = if (a.0, a.1) == (f, k) { (a, b) } else { (b, a) };
                    if (other.0, other.1) == (f, k) || !settled[other.0][other.1] {
                        continue;
                    }
                    if mine.2 == other.2 {
                        against += 1;
                    } else {
                        with += 1;
                    }
                }
                if with + against == 0 {
                    continue;
                }
                settled[f][k] = true;
                changed = true;
                if with > 0 {
                    continue;
                }
                let bound = &mut raw.faces[f].bounds[k].1;
                bound.reverse();
                for u in bound.iter_mut() {
                    u.reversed = !u.reversed;
                }
                let mut edges: Vec<usize> = bound.iter().map(|u| u.edge).collect();
                edges.sort_unstable();
                edges.dedup();
                for e in edges {
                    for x in uses[e].iter_mut().filter(|x| (x.0, x.1) == (f, k)) {
                        x.2 = !x.2;
                    }
                }
                recap(raw, f, k);
                turned[f] += 1;
            }
        }
    }
    turned
        .into_iter()
        .enumerate()
        .filter(|&(_, n)| n > 0)
        .collect()
}

/// Builds the solid of a `MANIFOLD_SOLID_BREP` or `BREP_WITH_VOIDS`, healed where its geometry
/// disagrees with itself, and checked. Returns the solid, its faces' STEP ids and names in order,
/// how far healing moved the geometry of each solid read, and the void shells read as solids of
/// their own. `strict`
/// builds it only where the file's geometry already agrees with itself, and is an error naming
/// where it does not.
///
/// A void is a hole in the material its outer shell bounds (ISO 10303-42, `brep_with_voids`). Some
/// files write a compound of bodies as one solid whose "voids" are the other bodies, outside the
/// outer shell and often touching it, or passing through it. Each void with a point outside the
/// outer shell is read as a solid of its own instead, and returned to be reported; a strict read
/// rejects the solid, saying why.
pub(super) fn solid(geo: &mut Geometry, id: u64, strict: bool) -> Result<ReadSolid> {
    if geo
        .file
        .get(id)
        .is_some_and(|e| e.part("TESSELLATED_SOLID").is_some())
    {
        return tessellated(geo, id).map(|b| {
            let names = vec![(0, String::new()); b.face_count()];
            (b, names, Vec::new(), Vec::new(), Vec::new())
        });
    }
    let mut rd = Reader {
        geo,
        raw: Raw::default(),
        vertices: HashMap::new(),
        edges: HashMap::new(),
        corners: HashMap::new(),
        sides: HashMap::new(),
        splits: HashMap::new(),
    };
    let listed = rd.solid(id)?;
    let uncertainty = rd.geo.uncertainty();
    let ours = rd.geo.file.originating_system() == Some(WRITER);
    let raw = without_split(rd.raw, &rd.splits);
    let read = |raw: Raw| built(raw, uncertainty, ours, strict);
    let first = |healed: Option<heal::Healed>| healed.map(|h| (0, h)).into_iter().collect();
    if raw.shells.len() == 1 {
        return read(raw)
            .map(|(b, names, healed, turned)| (b, names, first(healed), Vec::new(), turned));
    }
    // A void lies inside the outer shell, so one with a point outside it, beyond the tolerance, is
    // no void: it lies beside the outer shell, touching it or not, or passes through it, as one
    // body of a compound does through another. Read as voids, such shells make edges four faces
    // use where they touch, or a solid that passes through itself where they cross. Each void is
    // placed against the outer shell of the solid read whole, or, where that read is rejected, of
    // the outer shell read alone; a solid whose voids all lie inside is read once.
    let whole = read(raw.clone());
    let mut outer = None;
    let against = |b: &Brep, faces: &[FaceId]| -> Vec<usize> {
        (1..raw.shells.len())
            .filter(|&k| classify::any_outside(b, faces, &shell_points(&raw, k)))
            .collect()
    };
    let outside = match &whole {
        Ok((b, names, _, _)) => {
            let ids: HashSet<u64> = raw.shells[0].iter().map(|&f| raw.faces[f].id).collect();
            let faces: Vec<FaceId> = b
                .face_ids()
                .filter(|f| ids.contains(&names[f.index()].0))
                .collect();
            against(b, &faces)
        }
        Err(_) => {
            outer = read(shells_of(&raw, &[0])).ok();
            outer.as_ref().map_or_else(Vec::new, |(b, _, _, _)| {
                against(b, &b.face_ids().collect::<Vec<_>>())
            })
        }
    };
    if outside.is_empty() {
        return whole
            .map(|(b, names, healed, turned)| (b, names, first(healed), Vec::new(), turned));
    }
    let named: Vec<String> = outside.iter().map(|&k| format!("#{}", listed[k])).collect();
    if strict {
        return Err(Error::InvalidInput(format!(
            "void shell{} {} lie{} outside the outer shell, which a void of a solid does not; a \
             strict read does not read them as solids of their own",
            if outside.len() == 1 { "" } else { "s" },
            named.join(", "),
            if outside.len() == 1 { "s" } else { "" },
        )));
    }
    let kept: Vec<usize> = (0..raw.shells.len())
        .filter(|k| !outside.contains(k))
        .collect();
    let (mut b, mut names, healed, mut turned) = match outer {
        Some(outer) if kept.len() == 1 => outer,
        _ => read(shells_of(&raw, &kept))?,
    };
    let mut healed: Vec<(usize, heal::Healed)> = first(healed);
    for (&k, name) in outside.iter().zip(&named) {
        let (part, more, moved, also) = read(shells_of(&raw, &[k])).map_err(|e| {
            Error::InvalidInput(format!(
                "void shell {name} lies outside the outer shell, and read as a solid of its own: \
                 {e}"
            ))
        })?;
        healed.extend(moved.map(|h| (b.solid_ids().count(), h)));
        // The part's faces follow those read before it.
        turned.extend(also.into_iter().map(|(f, n)| (names.len() + f, n)));
        b.append(&part);
        names.extend(more);
    }
    Ok((
        b,
        names,
        healed,
        outside.iter().map(|&k| listed[k]).collect(),
        turned,
    ))
}

/// The part of `raw` bounded by its shells `keep`, in that order, with its own copies of the
/// points and edges they use.
fn shells_of(raw: &Raw, keep: &[usize]) -> Raw {
    let mut out = Raw::default();
    let mut face_at = vec![usize::MAX; raw.faces.len()];
    let mut edge_at = vec![usize::MAX; raw.edges.len()];
    let mut point_at = vec![usize::MAX; raw.points.len()];
    let mut point = |out: &mut Raw, p: usize| {
        if point_at[p] == usize::MAX {
            point_at[p] = out.points.len();
            out.points.push(raw.points[p]);
            out.vertex_ids.push(raw.vertex_ids[p]);
        }
        point_at[p]
    };
    let mut faces = Vec::new();
    for &k in keep {
        for &f in &raw.shells[k] {
            face_at[f] = faces.len();
            faces.push(raw.faces[f].clone());
        }
    }
    for f in &mut faces {
        for (_, uses) in &mut f.bounds {
            for u in uses {
                if edge_at[u.edge] == usize::MAX {
                    let e = &raw.edges[u.edge];
                    let (start, end) = (point(&mut out, e.start), point(&mut out, e.end));
                    edge_at[u.edge] = out.edges.len();
                    out.edges.push(RawEdge {
                        id: e.id,
                        curve: e.curve.clone(),
                        start,
                        end,
                        faces: e
                            .faces
                            .iter()
                            .map(|&g| face_at[g])
                            .filter(|&g| g != usize::MAX)
                            .collect(),
                    });
                }
                u.edge = edge_at[u.edge];
            }
        }
        for v in &mut f.points {
            *v = point(&mut out, *v);
        }
    }
    out.faces = faces;
    let mut next = 0;
    for &k in keep {
        let n = raw.shells[k].len();
        out.shells.push((next..next + n).collect());
        next += n;
    }
    out
}

/// The vertices of shell `k` of `raw` and the middles of its edges.
fn shell_points(raw: &Raw, k: usize) -> Vec<Point3> {
    let mut seen = vec![false; raw.edges.len()];
    let mut out = Vec::new();
    for &f in &raw.shells[k] {
        for (_, uses) in &raw.faces[f].bounds {
            for u in uses {
                if std::mem::replace(&mut seen[u.edge], true) {
                    continue;
                }
                let e = &raw.edges[u.edge];
                let (ps, pe) = (raw.points[e.start], raw.points[e.end]);
                out.extend([ps, pe]);
                if let Some((t0, t1)) = edge_range(&e.curve, ps, pe, e.start == e.end) {
                    out.push(e.curve.eval(0.5 * (t0 + t1)));
                }
            }
        }
        out.extend(raw.faces[f].points.iter().map(|&v| raw.points[v]));
    }
    out
}

/// The checked solid of the shells of `raw`, as [`solid`] reads them, with its faces' STEP ids and
/// names and how far healing moved it; `uncertainty` is the one its file states, and `ours`
/// whether this kernel wrote the file.
fn built(raw: Raw, uncertainty: Option<f64>, ours: bool, strict: bool) -> Result<Shells> {
    let mut file = own_nappes(join_contacts(raw))?;
    let turned = agree_senses(&mut file);
    let names = || file.faces.iter().map(|f| (f.id, f.name.clone())).collect();
    if strict {
        if let Some(why) = heal::misfit(&file) {
            return Err(Error::InvalidInput(format!(
                "{why}, and a strict read moves nothing to close it"
            )));
        }
        // A file that agrees with itself is left exactly as it is.
        let mut raw = file.clone();
        let (ranges, _) = heal::heal(&mut raw)?;
        return assemble(&raw, &ranges, None).map(|b| (b, names(), None, turned));
    }
    let exact = || {
        let mut raw = file.clone();
        heal::heal(&mut raw)
            .and_then(|(ranges, healed)| assemble(&raw, &ranges, None).map(|b| (b, healed)))
    };
    let within = uncertainty.filter(|&u| u > tol::LINEAR);
    // This kernel states an uncertainty above its own only for geometry it kept with local
    // tolerances, which a read of its file keeps too rather than healing within them; it heals
    // only if that fails.
    let kept = within.is_some() && ours;
    // Healing that moves geometry further than the file says its geometry meets has changed the
    // shape, not closed it: the file's geometry is kept instead where tolerances can hold it.
    let why = if kept {
        None
    } else {
        match exact() {
            Ok((b, healed))
                if within.is_none_or(|w| healed.as_ref().is_none_or(|h| h.moved <= w)) =>
            {
                return Ok((b, names(), healed, turned));
            }
            why => Some(why),
        }
    };
    // Geometry that cannot be made to meet within the kernel tolerance is kept as the file has
    // it, its edges and vertices given the tolerances it does meet within, whatever the file
    // states.
    let tolerant = |raw: &Raw| {
        let mut raw = raw.clone();
        heal::tolerate(&mut raw, within)
            .map_err(|t| format!("kept as the file has it, {t}"))
            .and_then(|t| {
                assemble(&raw, &t.ranges, Some(&t))
                    .map(|b| (b, t.moved))
                    .map_err(|e| format!("kept as the file has it, with local tolerances: {e}"))
            })
    };
    let t = match tolerant(&file) {
        Ok((b, moved)) => return Ok((b, names(), moved, turned)),
        Err(t) => t,
    };
    match why.unwrap_or_else(exact) {
        Ok((b, healed)) => Ok((b, names(), healed, turned)),
        Err(e) => {
            // Healed where it can be, and held by tolerances where it cannot.
            let mut loose = file.clone();
            let healed = heal::heal_loosely(&mut loose);
            let l = match (tolerant(&loose), healed.as_ref(), within) {
                (Ok(_), Some(h), Some(within)) if h.moved > within => format!(
                    "healing would move {} {:.1e} mm, more than the file's {within:.1e} mm",
                    h.what, h.moved
                ),
                (Ok((b, moved)), _, _) => {
                    return Ok((b, names(), farther(healed, moved), turned));
                }
                (Err(l), _, _) => l,
            };
            Err(Error::InvalidInput(format!(
                "{e}; {t}; healed where it can be and {l}"
            )))
        }
    }
}

/// The checked arena of a solid read as `raw`, its edges over `ranges`, and its vertices and edges
/// given the local tolerances of `tolerated` where it is given.
fn assemble(raw: &Raw, ranges: &[(f64, f64)], tolerated: Option<&heal::Tolerated>) -> Result<Brep> {
    let mut b = build(raw, ranges, tolerated)?;
    // Some systems write a void's faces facing into it and turn its shell as well, which reads
    // inside out: a shell is the way the volume it encloses says, outer positive, void negative.
    let turned = check::turned_shells(&b);
    if !turned.is_empty() {
        b = build(&turned_round(raw, &turned), ranges, tolerated)?;
    }
    let issues = check::check_structure(&b);
    if !issues.is_empty() {
        return Err(Error::InvalidResult(
            issues.into_iter().map(|i| i.to_string()).collect(),
        ));
    }
    Ok(b)
}

/// `raw` with the faces of each of `turned`'s shells facing the other way.
fn turned_round(raw: &Raw, turned: &[(SolidId, ShellId, bool)]) -> Raw {
    let mut raw = raw.clone();
    for &(_, sh, _) in turned {
        for &f in &raw.shells[sh.index()].clone() {
            let face = &mut raw.faces[f];
            face.reversed = !face.reversed;
            for (_, uses) in &mut face.bounds {
                uses.reverse();
                for u in uses.iter_mut() {
                    u.reversed = !u.reversed;
                }
            }
        }
    }
    raw
}

/// The arena of a solid read as `raw`, as the file orients it, unchecked.
fn build(raw: &Raw, ranges: &[(f64, f64)], tolerated: Option<&heal::Tolerated>) -> Result<Brep> {
    let mut bld = Builder::new(raw, ranges);
    if let Some(t) = tolerated {
        for (v, &x) in t.vertices.iter().enumerate() {
            bld.b.set_vertex_tolerance(bld.vertices[v], x);
        }
        for (e, &x) in t.edges.iter().enumerate() {
            bld.b.set_edge_tolerance(bld.edges[e], x);
        }
    }
    for f in &raw.faces {
        bld.face(f)?;
    }
    let mut b = bld.b;
    let mut base = 0;
    let mut shells = Vec::new();
    for sh in &raw.shells {
        let ids: Vec<_> = b.face_ids().skip(base).take(sh.len()).collect();
        base += sh.len();
        shells.push(b.add_shell(ids));
    }
    b.add_solid(shells);
    // A solid that touches itself along an edge uses it four times; the file does not say which
    // faces pair across the material there, the faces round it do.
    b.pair_radial();
    Ok(b)
}

/// An AP242 `TESSELLATED_SOLID`: each triangle of its `TRIANGULATED_FACE`s and
/// `COMPLEX_TRIANGULATED_FACE`s a plane face, the triangles sharing the edges between their
/// corners. The triangles run counter-clockwise seen from outside, as the schema has them.
fn tessellated(geo: &mut Geometry, id: u64) -> Result<Brep> {
    let unit = geo.length_unit();
    let Some(solid) = geo
        .file
        .get(id)
        .and_then(|e| e.part("TESSELLATED_SOLID"))
        .map(<[Param]>::to_vec)
    else {
        return bad(id, "expected a tessellated solid");
    };
    // Corners index their own face's list; one list's point is one vertex wherever it is used.
    let mut points: Vec<Point3> = Vec::new();
    let mut index: HashMap<(u64, usize), usize> = HashMap::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    for item in solid.get(1).and_then(Param::as_list).unwrap_or(&[]) {
        let fid = geo.reference(id, Some(item))?;
        let entity = geo.file.get(fid);
        let (face, complex) = match (
            entity.and_then(|e| e.part("TRIANGULATED_FACE")),
            entity.and_then(|e| e.part("COMPLEX_TRIANGULATED_FACE")),
        ) {
            (Some(f), _) => (f.to_vec(), false),
            (None, Some(f)) => (f.to_vec(), true),
            _ => return unsupported(fid, "a tessellated face other than a triangulated one"),
        };
        let list = geo.reference(fid, face.get(1))?;
        let Some(coords) = geo
            .file
            .get(list)
            .and_then(|e| e.part("COORDINATES_LIST"))
            .map(<[Param]>::to_vec)
        else {
            return bad(list, "expected a coordinates list");
        };
        let mut xyz = Vec::new();
        for c in coords.get(2).and_then(Param::as_list).unwrap_or(&[]) {
            let v = geo.numbers(list, Some(c))?;
            if v.len() != 3 || v.iter().any(|x| !x.is_finite()) {
                return bad(list, "expected three finite coordinates");
            }
            xyz.push(Point3::new(v[0], v[1], v[2]) * unit);
        }
        // Triangle corners index the face's own points, which index the coordinates list.
        let own: Vec<usize> = face
            .get(5)
            .and_then(Param::as_list)
            .unwrap_or(&[])
            .iter()
            .map(|p| p.as_i64().map_or(0, |i| i as usize))
            .collect();
        let corner = |k: i64| -> Option<usize> {
            let k = usize::try_from(k).ok()?.checked_sub(1)?;
            let i = if own.is_empty() {
                k
            } else {
                own.get(k)?.checked_sub(1)?
            };
            (i < xyz.len()).then_some(i)
        };
        let lists = |at: usize| -> Vec<Vec<i64>> {
            face.get(at)
                .and_then(Param::as_list)
                .unwrap_or(&[])
                .iter()
                .map(|l| {
                    l.as_list()
                        .unwrap_or(&[])
                        .iter()
                        .filter_map(Param::as_i64)
                        .collect()
                })
                .collect()
        };
        let tris = if complex {
            strips_and_fans(&lists(6), &lists(7))
        } else {
            lists(6)
        };
        for ks in tris {
            let idx: Option<Vec<usize>> = ks.iter().map(|&k| corner(k)).collect();
            let Some(idx) = idx.filter(|v| v.len() == 3) else {
                return bad(fid, "a triangle needs three corners in its list");
            };
            let mut corners = [0; 3];
            for (corner, &i) in corners.iter_mut().zip(&idx) {
                *corner = *index.entry((list, i)).or_insert_with(|| {
                    points.push(xyz[i]);
                    points.len() - 1
                });
            }
            triangles.push(corners);
        }
    }
    if triangles.is_empty() {
        return bad(id, "a tessellated solid with no triangles");
    }
    crate::mesh::from_indexed(&points, &triangles)
}

/// The triangles of a `COMPLEX_TRIANGULATED_FACE`'s strips and fans, each keeping the face's
/// turn: in a strip each point after the first two makes a triangle with the two before it, every
/// second one taken the other way round; in a fan each point after the second makes one with the
/// point before it and the first. A triangle with a corner twice, which joins two strips, has no
/// area and is left out.
fn strips_and_fans(strips: &[Vec<i64>], fans: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let strip = strips.iter().flat_map(|s| {
        s.windows(3).enumerate().map(|(i, w)| match i % 2 {
            0 => vec![w[0], w[1], w[2]],
            _ => vec![w[1], w[0], w[2]],
        })
    });
    let fan = fans
        .iter()
        .flat_map(|f| f.windows(2).skip(1).map(|w| vec![f[0], w[0], w[1]]));
    strip
        .chain(fan)
        .filter(|t| t[0] != t[1] && t[1] != t[2] && t[0] != t[2])
        .collect()
}

/// Every face of the shells in `file`, each with its own vertices and edges, sewn into one solid:
/// vertices within `tolerance` of each other become one, edges that then share their ends and lie
/// together become one, faces are turned so each shared edge runs opposite ways in its two faces,
/// and the solid is healed onto its faces' surfaces as a read one is.
pub(super) fn sew(geo: &mut Geometry, tolerance: f64) -> Result<Brep> {
    let shells: Vec<u64> = geo
        .file
        .ids
        .iter()
        .zip(&geo.file.entities)
        .filter(|(_, e)| matches!(e.kind(), "OPEN_SHELL" | "CLOSED_SHELL"))
        .map(|(&id, _)| id)
        .collect();
    let mut rd = Reader {
        geo,
        raw: Raw::default(),
        vertices: HashMap::new(),
        edges: HashMap::new(),
        corners: HashMap::new(),
        sides: HashMap::new(),
        splits: HashMap::new(),
    };
    let mut faces = Vec::new();
    for sh in shells {
        let (_, p) = rd.params(sh, &["OPEN_SHELL", "CLOSED_SHELL"])?;
        for f in p.get(1).and_then(Param::as_list).unwrap_or(&[]) {
            faces.push(rd.face(rd.geo.reference(sh, Some(f))?, false)?);
        }
    }
    if faces.is_empty() {
        return Err(Error::InvalidInput(
            "STEP: the file holds no faces to sew".into(),
        ));
    }
    sew_raw(own_nappes(without_split(rd.raw, &rd.splits))?, tolerance)
}

/// The faces of `raw`, each with its own vertices and edges, sewn into one solid as [`sew`] does.
fn sew_raw(mut raw: Raw, tolerance: f64) -> Result<Brep> {
    // Points within the tolerance of each other are one, at their mean.
    let n = raw.points.len();
    let mut rep: Vec<usize> = (0..n).collect();
    fn root(rep: &mut [usize], i: usize) -> usize {
        let mut r = i;
        while rep[r] != r {
            r = rep[r];
        }
        rep[i] = r;
        r
    }
    for i in 0..n {
        for j in i + 1..n {
            if raw.points[i].dist(raw.points[j]) <= tolerance {
                let (a, b) = (root(&mut rep, i), root(&mut rep, j));
                if a != b {
                    rep[b] = a;
                }
            }
        }
    }
    let mut sum: HashMap<usize, (Point3, f64)> = HashMap::new();
    for i in 0..n {
        let r = root(&mut rep, i);
        let e = sum.entry(r).or_insert((Point3::ZERO, 0.0));
        e.0 += raw.points[i];
        e.1 += 1.0;
    }
    for i in 0..n {
        let r = root(&mut rep, i);
        let (s, k) = sum[&r];
        raw.points[i] = s * (1.0 / k);
    }
    // Edges that now run between one pair of points, and lie together along it, are one: the
    // first stands for the rest, and their uses are turned to its direction.
    let ends = |e: &RawEdge, rep: &mut [usize]| (root(rep, e.start), root(rep, e.end));
    let mut keep: Vec<usize> = (0..raw.edges.len()).collect();
    let mut turned = vec![false; raw.edges.len()];
    for i in 0..raw.edges.len() {
        if keep[i] != i {
            continue;
        }
        let (a, b) = ends(&raw.edges[i], &mut rep);
        // Points a quarter, half and three quarters along an edge: two curves between one pair of
        // points can both pass nearest the middle of their chord at one of its ends.
        let inside = |e: &RawEdge, closed: bool| -> Vec<Point3> {
            let (s, z) = (raw.points[e.start], raw.points[e.end]);
            match edge_range(&e.curve, s, z, closed) {
                Some((t0, t1)) => [0.25, 0.5, 0.75]
                    .map(|f| e.curve.eval(t0 + (t1 - t0) * f))
                    .to_vec(),
                None => vec![e.curve.eval(param_of(&e.curve, s.lerp(z, 0.5)))],
            }
        };
        let on = |ps: &[Point3], e: &RawEdge| {
            ps.iter()
                .all(|&p| e.curve.eval(param_of(&e.curve, p)).dist(p) <= tolerance)
        };
        let mi = inside(&raw.edges[i], a == b);
        for j in i + 1..raw.edges.len() {
            if keep[j] != j {
                continue;
            }
            let (c, d) = ends(&raw.edges[j], &mut rep);
            let same = (a, b) == (c, d);
            if !(same || (a, b) == (d, c))
                || !on(&mi, &raw.edges[j])
                || !on(&inside(&raw.edges[j], c == d), &raw.edges[i])
            {
                continue;
            }
            keep[j] = i;
            turned[j] = !same;
        }
    }
    for f in &mut raw.faces {
        for (_, uses) in &mut f.bounds {
            for u in uses.iter_mut() {
                if keep[u.edge] != u.edge {
                    u.reversed ^= turned[u.edge];
                    u.edge = keep[u.edge];
                }
            }
        }
    }
    for e in &mut raw.edges {
        e.start = root(&mut rep, e.start);
        e.end = root(&mut rep, e.end);
    }
    // Each face's neighbours across its edges, and whether they agree: a shared edge runs one way
    // in one face and the other way in the other.
    let mut uses: HashMap<usize, Vec<(usize, bool)>> = HashMap::new();
    for (fi, f) in raw.faces.iter().enumerate() {
        for (_, us) in &f.bounds {
            for u in us {
                uses.entry(u.edge).or_default().push((fi, u.reversed));
            }
        }
    }
    let mut flip: Vec<Option<bool>> = vec![None; raw.faces.len()];
    for start in 0..raw.faces.len() {
        if flip[start].is_some() {
            continue;
        }
        flip[start] = Some(false);
        let mut stack = vec![start];
        while let Some(f) = stack.pop() {
            let mine = flip[f].unwrap_or(false);
            for (_, us) in &raw.faces[f].bounds {
                for u in us {
                    for &(g, rev) in uses.get(&u.edge).into_iter().flatten() {
                        if g == f || flip[g].is_some() {
                            continue;
                        }
                        // Agreeing faces use the edge opposite ways once both are turned as set.
                        let want = (u.reversed ^ mine) == rev;
                        flip[g] = Some(want);
                        stack.push(g);
                    }
                }
            }
        }
    }
    for (f, t) in raw.faces.iter_mut().zip(&flip) {
        if *t == Some(true) {
            f.reversed = !f.reversed;
            for (_, us) in &mut f.bounds {
                us.reverse();
                for u in us.iter_mut() {
                    u.reversed = !u.reversed;
                }
            }
        }
    }
    // Only the edges and points still used are kept, renumbered.
    let mut edge_at = vec![usize::MAX; raw.edges.len()];
    let mut edges = Vec::new();
    for (i, e) in raw.edges.drain(..).enumerate() {
        if keep[i] == i {
            edge_at[i] = edges.len();
            edges.push(e);
        }
    }
    raw.edges = edges;
    for f in &mut raw.faces {
        for (_, us) in &mut f.bounds {
            for u in us.iter_mut() {
                u.edge = edge_at[u.edge];
            }
        }
    }
    let mut point_at = vec![usize::MAX; raw.points.len()];
    let (mut points, mut ids) = (Vec::new(), Vec::new());
    for e in &mut raw.edges {
        for v in [&mut e.start, &mut e.end] {
            if point_at[*v] == usize::MAX {
                point_at[*v] = points.len();
                points.push(raw.points[*v]);
                ids.push(raw.vertex_ids[*v]);
            }
            *v = point_at[*v];
        }
    }
    for f in &mut raw.faces {
        for p in &mut f.points {
            *p = root(&mut rep, *p);
            if point_at[*p] == usize::MAX {
                point_at[*p] = points.len();
                points.push(raw.points[*p]);
                ids.push(raw.vertex_ids[*p]);
            }
            *p = point_at[*p];
        }
    }
    raw.points = points;
    raw.vertex_ids = ids;
    raw.edges.iter_mut().for_each(|e| e.faces.clear());
    for (fi, f) in raw.faces.iter().enumerate() {
        for (_, us) in &f.bounds {
            for u in us {
                raw.edges[u.edge].faces.push(fi);
            }
        }
    }
    raw.shells = vec![(0..raw.faces.len()).collect()];
    let (ranges, _) = heal::heal(&mut raw)?;
    assemble(&raw, &ranges, None)
}

/// `b` with each vertex at `points[v]`, each edge on `curves[e]` and each face on `surfaces[f]`,
/// its pcurves and parameter ranges found again as a file's are; `turned` runs every loop the
/// other way, for geometry that a mirror has turned inside out.
pub(crate) fn rebuilt(
    b: &Brep,
    points: Vec<Point3>,
    curves: Vec<Option<Curve>>,
    surfaces: Vec<Surface>,
    turned: bool,
) -> Result<Brep> {
    let mut raw = Raw {
        vertex_ids: (0..points.len() as u64).collect(),
        points,
        ..Raw::default()
    };
    let mut edge_at = vec![usize::MAX; curves.len()];
    for ((e, id), c) in b.edge_ids().enumerate().zip(curves) {
        if let Some(curve) = c {
            let x = b.edge(id);
            edge_at[e] = raw.edges.len();
            raw.edges.push(RawEdge {
                id: e as u64,
                curve,
                start: x.start.index(),
                end: x.end.index(),
                faces: Vec::new(),
            });
        }
    }
    let faces: Vec<FaceId> = b
        .face_ids()
        .filter(|&f| b.face(f).shell.is_some())
        .collect();
    for (&f, surface) in faces.iter().zip(surfaces) {
        let face = b.face(f);
        let mut bounds = Vec::new();
        for (k, &l) in face.loops.iter().enumerate() {
            let mut uses: Vec<Use> = b
                .lp(l)
                .coedges
                .iter()
                .map(|&c| b.coedge(c))
                .filter(|c| edge_at[c.edge.index()] != usize::MAX)
                .map(|c| Use {
                    edge: edge_at[c.edge.index()],
                    reversed: c.reversed ^ turned,
                })
                .collect();
            if turned {
                uses.reverse();
            }
            for u in &uses {
                raw.edges[u.edge].faces.push(raw.faces.len());
            }
            bounds.push((k == 0, uses));
        }
        raw.faces.push(RawFace {
            id: f.index() as u64,
            name: String::new(),
            surface,
            reversed: face.reversed ^ turned,
            bounds,
            points: Vec::new(),
        });
    }
    for sh in b.shell_ids() {
        let ids = &b.shell(sh).faces;
        raw.shells.push(
            ids.iter()
                .filter_map(|f| faces.iter().position(|g| g == f))
                .collect(),
        );
    }
    let (ranges, _) = heal::heal(&mut raw)?;
    let mut out = assemble(&raw, &ranges, None)?;
    for (g, &f) in out.face_ids().collect::<Vec<_>>().into_iter().zip(&faces) {
        out.set_face_names(g, b.face(f).names.clone());
    }
    Ok(out)
}

/// A solid as loose geometry and who uses what: the form a file gives one in, and the way to
/// build one whose faces are known by their surfaces and bounding edges alone.
#[derive(Clone, Debug, Default)]
pub(crate) struct Parts {
    pub points: Vec<Point3>,
    /// Each edge's curve and its ends, as indices into `points`.
    pub edges: Vec<(Curve, usize, usize)>,
    pub faces: Vec<PartFace>,
}

#[derive(Clone, Debug)]
pub(crate) struct PartFace {
    pub surface: Surface,
    pub reversed: bool,
    /// Each loop's uses of edges in order, the outer loop first: the edge, and whether the loop
    /// runs against it.
    pub loops: Vec<Vec<(usize, bool)>>,
    pub names: std::collections::BTreeSet<crate::topo::FaceName>,
    /// Which of the solid's shells it lies in, the outer one first.
    pub shell: usize,
}

/// The parts of `b`, one solid bounded by one outer shell and any voids. Degenerate edges, the
/// kernel's own, are left out.
pub(crate) fn parts_of(b: &Brep) -> Result<Parts> {
    let mut solids = b.solid_ids();
    let (Some(solid), None) = (solids.next(), solids.next()) else {
        return Err(Error::Unsupported("the parts of several solids".into()));
    };
    let solid = b.solid(solid);
    if solid.outer != 1 {
        return Err(Error::Unsupported(
            "the parts of a solid that touches itself".into(),
        ));
    }
    let mut parts = Parts {
        points: b.vertex_ids().map(|v| b.vertex(v).point).collect(),
        ..Parts::default()
    };
    let mut edge_at = vec![usize::MAX; b.edge_ids().count()];
    for e in b.edge_ids() {
        let x = b.edge(e);
        if let Some(curve) = &x.curve {
            edge_at[e.index()] = parts.edges.len();
            parts
                .edges
                .push((curve.clone(), x.start.index(), x.end.index()));
        }
    }
    for (shell, &sh) in solid.shells.iter().enumerate() {
        for &f in &b.shell(sh).faces {
            let face = b.face(f);
            let loops = face
                .loops
                .iter()
                .map(|&l| {
                    b.lp(l)
                        .coedges
                        .iter()
                        .map(|&c| b.coedge(c))
                        .filter(|c| edge_at[c.edge.index()] != usize::MAX)
                        .map(|c| (edge_at[c.edge.index()], c.reversed))
                        .collect()
                })
                .collect();
            parts.faces.push(PartFace {
                surface: face.surface.clone(),
                reversed: face.reversed,
                loops,
                names: face.names.clone(),
                shell,
            });
        }
    }
    Ok(parts)
}

/// The solid `parts` describe, its pcurves and parameter ranges found as a file's are. Points and
/// edges no face uses are left out.
pub(crate) fn solid_of(parts: &Parts) -> Result<Brep> {
    let (mut raw, order) = raw_of(parts);
    let (ranges, _) = heal::heal(&mut raw)?;
    let mut out = assemble(&raw, &ranges, None)?;
    named(&mut out, parts, &order);
    Ok(out)
}

/// The faces `parts` describe, each with its own points and edges, sewn into one solid as a STEP
/// file's open shells are ([`sew`]): points within `tolerance` become one, and so do the edges
/// that then lie together.
pub(crate) fn sewn(parts: &Parts, tolerance: f64) -> Result<Brep> {
    let (raw, order) = raw_of(parts);
    let mut out = sew_raw(raw, tolerance)?;
    named(&mut out, parts, &order);
    Ok(out)
}

/// Gives each face of `out`, built in `order` from `parts`, its part's names.
fn named(out: &mut Brep, parts: &Parts, order: &[usize]) {
    for (g, &at) in out.face_ids().collect::<Vec<_>>().into_iter().zip(order) {
        out.set_face_names(g, parts.faces[at].names.clone());
    }
}

/// `parts` as a read solid's topology, its faces in shell order, and the part each face came from.
fn raw_of(parts: &Parts) -> (Raw, Vec<usize>) {
    let mut raw = Raw::default();
    let mut edge_at = vec![usize::MAX; parts.edges.len()];
    let mut point_at = vec![usize::MAX; parts.points.len()];
    let mut order: Vec<usize> = (0..parts.faces.len()).collect();
    order.sort_by_key(|&f| parts.faces[f].shell);
    for (fi, &at) in order.iter().enumerate() {
        let f = &parts.faces[at];
        let mut bounds = Vec::new();
        for (k, l) in f.loops.iter().enumerate() {
            let mut uses = Vec::new();
            for &(e, reversed) in l {
                if edge_at[e] == usize::MAX {
                    let (curve, start, end) = &parts.edges[e];
                    let mut ends = [*start, *end];
                    for v in &mut ends {
                        if point_at[*v] == usize::MAX {
                            point_at[*v] = raw.points.len();
                            raw.vertex_ids.push(raw.points.len() as u64);
                            raw.points.push(parts.points[*v]);
                        }
                        *v = point_at[*v];
                    }
                    edge_at[e] = raw.edges.len();
                    raw.edges.push(RawEdge {
                        id: e as u64,
                        curve: curve.clone(),
                        start: ends[0],
                        end: ends[1],
                        faces: Vec::new(),
                    });
                }
                raw.edges[edge_at[e]].faces.push(fi);
                uses.push(Use {
                    edge: edge_at[e],
                    reversed,
                });
            }
            bounds.push((k == 0, uses));
        }
        raw.faces.push(RawFace {
            id: fi as u64,
            name: String::new(),
            surface: f.surface.clone(),
            reversed: f.reversed,
            bounds,
            points: Vec::new(),
        });
        if raw.shells.len() <= f.shell {
            raw.shells.resize(f.shell + 1, Vec::new());
        }
        raw.shells[f.shell].push(fi);
    }
    (raw, order)
}
