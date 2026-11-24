//! Rational B-spline curves and surfaces of any degree up to [`MAX_DEGREE`], with explicit knot
//! vectors, as STEP files carry them. A curve or surface whose two ends meet is periodic in that
//! parameter: evaluation wraps into the knot domain, so an edge or loop may run across the join.

use std::{
    cell::RefCell,
    sync::{
        Arc, OnceLock,
        atomic::{AtomicU32, Ordering},
    },
};

use super::steady::Swept;
use crate::{
    Result,
    error::invalid,
    math::{Aabb, Point3, Transform, Vec2, Vec3},
    mesh::bvh::Bvh,
    numeric::{self, solve2},
    tol,
};

pub(super) const MAX_DEGREE: usize = 15;
/// How near the surface, as a share of the seed spacing, a local projection is taken as nearest.
const NEAR_SHARE: f64 = 0.05;
const M: usize = MAX_DEGREE + 1;
/// The highest degree the span polynomials evaluate (see [`NurbsSurface::by_polys`]).
const MAX_POLY_DEGREE: usize = 6;

/// A knot vector and the degree of the basis it defines.
#[derive(Clone, Debug, PartialEq)]
struct Basis {
    degree: usize,
    knots: Vec<f64>,
    /// The number of basis functions (control points along this parameter).
    n: usize,
    /// `1 / (uₖ₊d − uₖ)` for each `k`, zero where the knots coincide, for `d` the degree and
    /// the degree less one: the derivative recurrence divides by these at every evaluation.
    inverse: [Vec<f64>; 2],
}

impl Basis {
    fn new(degree: usize, knots: Vec<f64>, n: usize) -> Result<Basis> {
        if degree == 0 || degree > MAX_DEGREE {
            return invalid(format!("B-spline degree {degree} is not 1 to {MAX_DEGREE}"));
        }
        if n <= degree {
            return invalid("a B-spline needs more control points than its degree");
        }
        if knots.len() != n + degree + 1 {
            return invalid("B-spline knot count must be control points + degree + 1");
        }
        if knots.iter().any(|k| !k.is_finite()) || knots.windows(2).any(|w| w[1] < w[0]) {
            return invalid("B-spline knots must be finite and not decrease");
        }
        let inverse = [degree, degree - 1].map(|d| {
            (0..knots.len() - d)
                .map(|k| {
                    let den = knots[k + d] - knots[k];
                    if den == 0.0 { 0.0 } else { 1.0 / den }
                })
                .collect()
        });
        let b = Basis {
            degree,
            knots,
            n,
            inverse,
        };
        let (a, z) = b.domain();
        if z - a <= 0.0 {
            return invalid("B-spline has an empty domain");
        }
        Ok(b)
    }

    fn domain(&self) -> (f64, f64) {
        (self.knots[self.degree], self.knots[self.n])
    }

    fn wrap(&self, t: f64, period: Option<f64>) -> f64 {
        let (a, b) = self.domain();
        match period {
            Some(p) => {
                let w = a + (t - a).rem_euclid(p);
                if w >= b { a } else { w }
            }
            None => t,
        }
    }

    fn span(&self, t: f64) -> usize {
        let (p, m) = (self.degree, self.n);
        let u = &self.knots;
        if t >= u[m] {
            let mut i = m - 1;
            while i > p && u[i] >= u[m] {
                i -= 1;
            }
            return i;
        }
        if t <= u[p] {
            let mut i = p;
            while i + 1 < m && u[i + 1] <= u[p] {
                i += 1;
            }
            return i;
        }
        let (mut lo, mut hi) = (p, m);
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if t < u[mid] { hi = mid } else { lo = mid }
        }
        lo
    }

    /// Basis values and derivatives up to `nd` (at most 2) at `t` in span `i`: `out[k][j]` is the
    /// `k`-th derivative of function `i − degree + j`. Values are raised a degree at a time,
    /// keeping the two lower degrees for the derivative recurrence.
    /// `N` is any length past the degree: the arrays are only that long, so a low degree does
    /// not clear and copy the room a degree-15 basis needs.
    fn ders<const N: usize>(&self, i: usize, t: f64, nd: usize) -> [[f64; N]; 3] {
        let p = self.degree;
        // The values at the degree reached, then at one and two degrees below it.
        let mut n = [0.0; N];
        let mut lower = [[0.0; N]; 2];
        n[0] = 1.0;
        for q in 1..=p {
            lower[1] = lower[0];
            lower[0] = n;
            n = raised(&self.knots, i, t, q, &lower[0]);
        }
        let mut out = [n, [0.0; N], [0.0; N]];
        if nd >= 1 {
            out[1] = self.deriv(&lower[0], p, i);
        }
        if nd >= 2 && p >= 2 {
            let below = self.deriv(&lower[1], p - 1, i);
            out[2] = self.deriv(&below, p, i);
        }
        out
    }

    /// `f` of [`Basis::ders`] in arrays sized for the degree.
    fn ders_then<R>(&self, i: usize, t: f64, nd: usize, f: impl FnOnce([&[f64]; 3]) -> R) -> R {
        match self.degree {
            0..4 => {
                let d = self.ders::<4>(i, t, nd);
                f([&d[0], &d[1], &d[2]])
            }
            4..8 => {
                let d = self.ders::<8>(i, t, nd);
                f([&d[0], &d[1], &d[2]])
            }
            _ => {
                let d = self.ders::<M>(i, t, nd);
                f([&d[0], &d[1], &d[2]])
            }
        }
    }

    /// The derivatives of the degree-`deg` functions on span `i` from the values (or derivatives)
    /// `a` of the degree-`deg − 1` ones:
    /// `N'ₖ,d = d·(Nₖ,d−1 / (uₖ₊d − uₖ) − Nₖ₊₁,d−1 / (uₖ₊d₊₁ − uₖ₊₁))`.
    fn deriv<const N: usize>(&self, a: &[f64; N], deg: usize, i: usize) -> [f64; N] {
        let inv = &self.inverse[self.degree - deg];
        let mut out = [0.0; N];
        let d = deg as f64;
        for (j, o) in out.iter_mut().enumerate().take(deg + 1) {
            let k = i + j - deg;
            let mut v = 0.0;
            if j >= 1 {
                v += a[j - 1] * inv[k];
            }
            if j < deg {
                v -= a[j] * inv[k + 1];
            }
            *o = d * v;
        }
        out
    }

    fn span_poly(&self, i: usize, end: bool) -> Vec<Vec<f64>> {
        span_poly(self.degree, &self.knots, i, end)
    }

    /// The distinct knots strictly inside (t0, t1), repeated each period when periodic.
    fn breakpoints(&self, t0: f64, t1: f64, period: Option<f64>) -> Vec<f64> {
        let (a, b) = self.domain();
        let inner: Vec<f64> = self.knots[self.degree..=self.n].to_vec();
        let mut out: Vec<f64> = Vec::new();
        let mut push = |k: f64| {
            if k > t0 && k < t1 && out.last().is_none_or(|&l| k > l) {
                out.push(k);
            }
        };
        match period {
            Some(p) => {
                let start = ((t0 - a) / p).floor() as i64;
                let end = ((t1 - a) / p).ceil() as i64;
                for w in start..=end {
                    for &k in &inner {
                        push(k + w as f64 * p);
                    }
                }
            }
            None => {
                for &k in &inner {
                    if k > a && k < b {
                        push(k);
                    }
                }
            }
        }
        out
    }

    /// The distinct knots of the domain, ends included.
    fn spans(&self) -> Vec<f64> {
        let mut out: Vec<f64> = Vec::new();
        for &k in &self.knots[self.degree..=self.n] {
            if out.last().is_none_or(|&l| k > l) {
                out.push(k);
            }
        }
        out
    }
}

fn check_weights(weights: &Option<Vec<f64>>, n: usize) -> Result<()> {
    if let Some(w) = weights {
        if w.len() != n {
            return invalid("B-spline weight count must match its control points");
        }
        if w.iter().any(|&x| !(x > 0.0 && x.is_finite())) {
            return invalid("B-spline weights must be positive");
        }
    }
    Ok(())
}

/// The values at `t` of the degree-`q` B-spline functions of knots `u` nonzero on span `i`, from
/// those of degree `q − 1` in `below`: entry `j` is function `i − q + j`. A function of degree q
/// is the blend of its two neighbours of degree q − 1, each weighted by where `t` lies in its
/// support; a neighbour that is zero on this span drops out, and a quotient with a zero
/// denominator, at a repeated knot, counts as zero.
fn raised<const N: usize>(u: &[f64], i: usize, t: f64, q: usize, below: &[f64; N]) -> [f64; N] {
    let ratio = |a: f64, b: f64| if b == 0.0 { 0.0 } else { a / b };
    let mut out = [0.0; N];
    for (j, o) in out.iter_mut().enumerate().take(q + 1) {
        let f = i + j - q;
        if j > 0 {
            *o += ratio(t - u[f], u[f + q] - u[f]) * below[j - 1];
        }
        if j < q {
            *o += ratio(u[f + q + 1] - t, u[f + q + 1] - u[f + 1]) * below[j];
        }
    }
    out
}

/// The values and first `nd` derivatives at `t` of the degree-`p` B-spline functions of knots `u`
/// that are nonzero on span `i`: `out[k][j]` is the `k`-th derivative of function `i − p + j`.
pub(crate) fn basis_derivatives<const N: usize>(
    p: usize,
    u: &[f64],
    i: usize,
    t: f64,
    nd: usize,
) -> [[f64; N]; N] {
    let ratio = |a: f64, b: f64| if b == 0.0 { 0.0 } else { a / b };
    // by_degree[q][j] is the value of the degree-q function `i − q + j`.
    let mut by_degree = [[0.0; N]; N];
    by_degree[0][0] = 1.0;
    for q in 1..=p {
        by_degree[q] = raised(u, i, t, q, &by_degree[q - 1]);
    }
    // A derivative of a degree-q function is q times the difference of its two neighbours of
    // degree q − 1, each over the width of its support. So the k-th derivative at degree p is
    // the values at degree p − k, differenced k times on the way back up.
    let mut out = [[0.0; N]; N];
    for (k, slot) in out.iter_mut().enumerate().take(nd.min(p) + 1) {
        let mut row = by_degree[p - k];
        for q in p - k + 1..=p {
            let mut raised = [0.0; N];
            for (j, r) in raised.iter_mut().enumerate().take(q + 1) {
                let f = i + j - q;
                let mut d = 0.0;
                if j > 0 {
                    d += ratio(row[j - 1], u[f + q] - u[f]);
                }
                if j < q {
                    d -= ratio(row[j], u[f + q + 1] - u[f + 1]);
                }
                *r = q as f64 * d;
            }
            row = raised;
        }
        *slot = row;
    }
    out
}

/// The degree-`p` B-spline functions of knots `u` nonzero on span `i`, as polynomials in
/// `s = (t − e)/(tᵢ₊₁ − tᵢ)`, `e` the span's start, or its end when `end`: `out[j][k]` is the
/// coefficient of `sᵏ` in function `i − p + j`. Each is the Taylor series at `e`, from every
/// derivative there, scaled by the span's width so the coefficients stay of the functions' own
/// size.
pub(crate) fn span_poly(p: usize, u: &[f64], i: usize, end: bool) -> Vec<Vec<f64>> {
    let e = if end { u[i + 1] } else { u[i] };
    let width = u[i + 1] - u[i];
    let ders = basis_derivatives::<M>(p, u, i, e, p);
    // The series term of order k is the k-th derivative over k!, times width^k.
    let mut out = vec![vec![0.0; p + 1]; p + 1];
    let mut scale = 1.0;
    for k in 0..=p {
        for j in 0..=p {
            out[j][k] = ders[k][j] * scale;
        }
        scale *= width / (k + 1) as f64;
    }
    out
}

/// The point and its derivatives on a rational surface from the homogeneous ones `h`, the first
/// `n` of `(0,0), (1,0), (0,1), (2,0), (1,1), (0,2)` filled; the rest come back zero.
fn rational(h: &[[f64; 4]; 6], n: usize) -> [Vec3; 6] {
    let a = |i: usize| Vec3::new(h[i][0], h[i][1], h[i][2]);
    let w = |i: usize| h[i][3];
    let s = a(0) / w(0);
    let mut out = [
        s,
        Vec3::ZERO,
        Vec3::ZERO,
        Vec3::ZERO,
        Vec3::ZERO,
        Vec3::ZERO,
    ];
    if n >= 3 {
        out[1] = (a(1) - s * w(1)) / w(0);
        out[2] = (a(2) - s * w(2)) / w(0);
    }
    if n >= 6 {
        let (su, sv) = (out[1], out[2]);
        out[3] = (a(3) - su * (2.0 * w(1)) - s * w(3)) / w(0);
        out[4] = (a(4) - su * w(2) - sv * w(1) - s * w(4)) / w(0);
        out[5] = (a(5) - sv * (2.0 * w(2)) - s * w(5)) / w(0);
    }
    out
}

/// A piece of a surface's parameter domain and a box that holds the surface over it.
#[derive(Clone, Copy, Debug)]
pub struct Patch {
    pub bbox: Aabb,
    pub lo: Vec2,
    pub hi: Vec2,
}

/// A piece of a surface over a rectangle of its parameters as a Bézier patch: its control points
/// in homogeneous coordinates, `rows` along u of `cols` along v. The patch lies in their convex
/// hull, and so does each half a split leaves, in its own.
#[derive(Clone, Debug)]
pub(crate) struct Bezier {
    pub lo: Vec2,
    pub hi: Vec2,
    rows: usize,
    cols: usize,
    net: Vec<[f64; 4]>,
}

impl Bezier {
    /// The box of the control points, which holds the patch.
    pub(crate) fn bbox(&self) -> Aabb {
        let mut bb = Aabb::EMPTY;
        for q in &self.net {
            bb.add(Vec3::new(q[0] / q[3], q[1] / q[3], q[2] / q[3]));
        }
        bb
    }

    /// The patch's two halves, cut across u or across v at the middle (de Casteljau).
    pub(crate) fn halves(&self, across_u: bool) -> [Bezier; 2] {
        let (lines, len) = if across_u {
            (self.cols, self.rows)
        } else {
            (self.rows, self.cols)
        };
        let at = |line: usize, k: usize| {
            if across_u {
                k * self.cols + line
            } else {
                line * self.cols + k
            }
        };
        let (mut first, mut second) = (self.net.clone(), self.net.clone());
        for line in 0..lines {
            let mut level: Vec<[f64; 4]> = (0..len).map(|k| self.net[at(line, k)]).collect();
            for step in 0..len {
                first[at(line, step)] = level[0];
                second[at(line, len - 1 - step)] = level[len - 1 - step];
                for k in 0..len - 1 - step {
                    level[k] = std::array::from_fn(|c| 0.5 * (level[k][c] + level[k + 1][c]));
                }
            }
        }
        let mid = (self.lo + self.hi) * 0.5;
        let (first_hi, second_lo) = if across_u {
            (Vec2::new(mid.x, self.hi.y), Vec2::new(mid.x, self.lo.y))
        } else {
            (Vec2::new(self.hi.x, mid.y), Vec2::new(self.lo.x, mid.y))
        };
        let half = |lo: Vec2, hi: Vec2, net: Vec<[f64; 4]>| Bezier {
            lo,
            hi,
            rows: self.rows,
            cols: self.cols,
            net,
        };
        [
            half(self.lo, first_hi, first),
            half(second_lo, self.hi, second),
        ]
    }

    /// The patch with each control point, in homogeneous coordinates, put through a linear map:
    /// the same piece of surface in other coordinates.
    pub(crate) fn mapped(&self, f: impl Fn([f64; 4]) -> [f64; 4]) -> Bezier {
        Bezier {
            lo: self.lo,
            hi: self.hi,
            rows: self.rows,
            cols: self.cols,
            net: self.net.iter().map(|&q| f(q)).collect(),
        }
    }

    /// The control points.
    pub(crate) fn points(&self) -> &[[f64; 4]] {
        &self.net
    }

    /// The least and greatest slope of homogeneous coordinate `c` across u and across v, each
    /// over the patch's own unit square: its derivative is a Bézier patch whose control values
    /// are the steps between neighbouring control points, so it lies between their extremes.
    pub(crate) fn slopes(&self, c: usize) -> [(f64, f64); 2] {
        let at = |i: usize, j: usize| self.net[i * self.cols + j][c];
        let span = |values: &mut dyn Iterator<Item = f64>| {
            values.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                (lo.min(x), hi.max(x))
            })
        };
        let (pu, pv) = ((self.rows - 1) as f64, (self.cols - 1) as f64);
        let across_u = span(
            &mut (0..self.rows - 1)
                .flat_map(|i| (0..self.cols).map(move |j| (i, j)))
                .map(|(i, j)| pu * (at(i + 1, j) - at(i, j))),
        );
        let across_v = span(
            &mut (0..self.rows)
                .flat_map(|i| (0..self.cols - 1).map(move |j| (i, j)))
                .map(|(i, j)| pv * (at(i, j + 1) - at(i, j))),
        );
        [across_u, across_v]
    }

    /// The patch at `(s, t)` of its own unit square, and its slopes across u and across v there,
    /// all in homogeneous coordinates (de Casteljau).
    pub(crate) fn at(&self, s: f64, t: f64) -> [[f64; 4]; 3] {
        let mix = |a: [f64; 4], b: [f64; 4], x: f64| -> [f64; 4] {
            std::array::from_fn(|c| (1.0 - x) * a[c] + x * b[c])
        };
        // A line of control points down to its last two, which give the value and the slope.
        let fold = |mut level: Vec<[f64; 4]>, x: f64| -> ([f64; 4], [f64; 4]) {
            let n = (level.len() - 1) as f64;
            while level.len() > 2 {
                level = level.windows(2).map(|w| mix(w[0], w[1], x)).collect();
            }
            if level.len() == 1 {
                return (level[0], [0.0; 4]);
            }
            (
                mix(level[0], level[1], x),
                std::array::from_fn(|c| n * (level[1][c] - level[0][c])),
            )
        };
        let along: Vec<([f64; 4], [f64; 4])> = (0..self.rows)
            .map(|i| fold(self.net[i * self.cols..(i + 1) * self.cols].to_vec(), t))
            .collect();
        let (value, across_u) = fold(along.iter().map(|r| r.0).collect(), s);
        let (across_v, _) = fold(along.iter().map(|r| r.1).collect(), s);
        [value, across_u, across_v]
    }

    /// Whether the patch runs further across u than across v, by its control points' spread.
    pub(crate) fn longer_in_u(&self) -> bool {
        let point = |i: usize, j: usize| {
            let q = self.net[i * self.cols + j];
            Vec3::new(q[0] / q[3], q[1] / q[3], q[2] / q[3])
        };
        let across_u = point(0, 0).dist(point(self.rows - 1, 0))
            + point(0, self.cols - 1).dist(point(self.rows - 1, self.cols - 1));
        let across_v = point(0, 0).dist(point(0, self.cols - 1))
            + point(self.rows - 1, 0).dist(point(self.rows - 1, self.cols - 1));
        across_u >= across_v
    }
}

/// The matrix that turns a polynomial's values at `0, 1/n, …, 1` into its Bézier coefficients.
fn bernstein_from_values(n: usize) -> Vec<Vec<f64>> {
    let binomial = |k: usize| (0..k).fold(1.0, |b, i| b * (n - i) as f64 / (i + 1) as f64);
    let size = n + 1;
    // Each row: the Bernstein functions at one node, then that row of the identity.
    let mut m: Vec<Vec<f64>> = (0..size)
        .map(|j| {
            let s = j as f64 / n.max(1) as f64;
            let mut row: Vec<f64> = (0..size)
                .map(|k| binomial(k) * s.powi(k as i32) * (1.0 - s).powi((n - k) as i32))
                .collect();
            row.extend((0..size).map(|k| if k == j { 1.0 } else { 0.0 }));
            row
        })
        .collect();
    for col in 0..size {
        let pivot = (col..size)
            .max_by(|&a, &b| m[a][col].abs().total_cmp(&m[b][col].abs()))
            .expect("a row is left");
        m.swap(col, pivot);
        let by = m[col].clone();
        for (row, values) in m.iter_mut().enumerate() {
            let f = if row == col {
                continue;
            } else {
                values[col] / by[col]
            };
            for (x, p) in values.iter_mut().zip(&by) {
                *x -= f * p;
            }
        }
    }
    (0..size)
        .map(|i| (0..size).map(|k| m[i][size + k] / m[i][i]).collect())
        .collect()
}

/// Inserts knot `t` once into a curve of homogeneous control points (Boehm's algorithm).
fn insert_knot(knots: &mut Vec<f64>, p: usize, pts: &mut Vec<[f64; 4]>, t: f64) {
    let n = pts.len();
    let mut k = p;
    while k + 1 < n && knots[k + 1] <= t {
        k += 1;
    }
    let mut out = Vec::with_capacity(n + 1);
    out.extend_from_slice(&pts[..=k - p]);
    for i in k - p + 1..=k {
        let den = knots[i + p] - knots[i];
        let a = if den > 0.0 { (t - knots[i]) / den } else { 0.0 };
        let (q, r) = (pts[i], pts[i - 1]);
        out.push(std::array::from_fn(|c| a * q[c] + (1.0 - a) * r[c]));
    }
    out.extend_from_slice(&pts[k..]);
    *pts = out;
    knots.insert(k + 1, t);
}

/// The stretch of a curve of homogeneous control points from `t0` to `t1`, both inside its domain,
/// as a curve of its own that starts and ends on its first and last control points.
pub(crate) fn piece(
    p: usize,
    knots: &[f64],
    pts: &[[f64; 4]],
    (t0, t1): (f64, f64),
) -> (Vec<f64>, Vec<[f64; 4]>) {
    debug_assert!(
        t0 < t1 && t0 >= knots[p] && t1 <= knots[knots.len() - 1 - p],
        "{t0}..{t1} is not a stretch of the curve's domain"
    );
    let (mut knots, mut pts) = (knots.to_vec(), pts.to_vec());
    for t in [t1, t0] {
        let have = knots.iter().filter(|&&k| k == t).count();
        for _ in have..p {
            insert_knot(&mut knots, p, &mut pts, t);
        }
    }
    let k = knots.iter().rposition(|&x| x == t0).expect("inserted");
    let j = knots.iter().position(|&x| x == t1).expect("inserted");
    let mut out = Vec::with_capacity(j - k + 2 * p + 1);
    out.push(t0);
    out.extend_from_slice(&knots[k + 1 - p..j + p]);
    out.push(t1);
    (out, pts[k - p..j].to_vec())
}

/// A rational (or, without weights, polynomial) B-spline curve.
#[derive(Clone, Debug, PartialEq)]
pub struct NurbsCurve {
    basis: Basis,
    ctrl: Vec<Point3>,
    weights: Option<Vec<f64>>,
    period: Option<f64>,
    stationary: StationaryCache,
}

/// Along each of x, y and z, the parameters where that coordinate of a curve is stationary, each
/// with its point.
type Stationary = [Vec<(f64, Point3)>; 3];

/// Each knot span's [`Stationary`] over the whole span, by span number, worked out when a box
/// first reaches the span and shared with the curve's copies. It depends only on the curve, so it
/// never takes part in equality.
#[derive(Clone, Debug, Default)]
struct StationaryCache(Arc<OnceLock<Vec<OnceLock<Stationary>>>>);

impl PartialEq for StationaryCache {
    fn eq(&self, _: &StationaryCache) -> bool {
        true
    }
}

impl NurbsCurve {
    /// The curve from its parts; weights all 1 when `weights` is `None`. The curve is periodic when
    /// its two ends meet within the linear tolerance.
    pub fn new(
        degree: usize,
        knots: Vec<f64>,
        ctrl: Vec<Point3>,
        weights: Option<Vec<f64>>,
    ) -> Result<NurbsCurve> {
        let basis = Basis::new(degree, knots, ctrl.len())?;
        check_weights(&weights, ctrl.len())?;
        let weights = weights.filter(|w| w.iter().any(|&x| x != w[0]));
        let mut c = NurbsCurve {
            basis,
            ctrl,
            weights,
            period: None,
            stationary: StationaryCache::default(),
        };
        let (a, b) = c.domain();
        if c.eval(a).dist(c.eval(b)) <= tol::LINEAR {
            c.period = Some(b - a);
        }
        Ok(c)
    }

    /// The curve exactly as [`NurbsCurve::new`] made it, from its parts and its period.
    pub(crate) fn restore(
        degree: usize,
        knots: Vec<f64>,
        ctrl: Vec<Point3>,
        weights: Option<Vec<f64>>,
        period: Option<f64>,
    ) -> Result<NurbsCurve> {
        let basis = Basis::new(degree, knots, ctrl.len())?;
        check_weights(&weights, ctrl.len())?;
        Ok(NurbsCurve {
            basis,
            ctrl,
            weights,
            period,
            stationary: StationaryCache::default(),
        })
    }

    pub fn degree(&self) -> usize {
        self.basis.degree
    }

    pub fn knots(&self) -> &[f64] {
        &self.basis.knots
    }

    pub fn control_points(&self) -> &[Point3] {
        &self.ctrl
    }

    /// The weights, or `None` when every weight is equal (the curve is polynomial).
    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    pub fn domain(&self) -> (f64, f64) {
        self.basis.domain()
    }

    pub fn period(&self) -> Option<f64> {
        self.period
    }

    pub fn breakpoints(&self, t0: f64, t1: f64) -> Vec<f64> {
        self.basis.breakpoints(t0, t1, self.period)
    }

    /// The box of the control points that shape the curve at `t`. Its weights being positive, the
    /// curve over the knot span there lies in their convex hull, and so in the box.
    pub(crate) fn span_box(&self, t: f64) -> Aabb {
        let t = self.basis.wrap(t, self.period);
        let (p, i) = (self.basis.degree, self.basis.span(t));
        let mut bb = Aabb::EMPTY;
        for &q in &self.ctrl[i - p..=i] {
            bb.add(q);
        }
        bb
    }

    /// The point and its first two derivatives at `t`.
    pub fn d2(&self, t: f64) -> (Point3, Vec3, Vec3) {
        self.ders(t, 2)
    }

    /// The point and its derivatives up to `nd` (at most 2) at `t`, the rest zero: each is what
    /// [`NurbsCurve::d2`] gives, without the sums the higher ones take.
    fn ders(&self, t: f64, nd: usize) -> (Point3, Vec3, Vec3) {
        let t = self.basis.wrap(t, self.period);
        let p = self.basis.degree;
        let i = self.basis.span(t);
        self.basis
            .ders_then(i, t, nd, |n| self.ders_from(n, p, i, nd))
    }

    #[allow(clippy::needless_range_loop)]
    fn ders_from(&self, n: [&[f64]; 3], p: usize, i: usize, nd: usize) -> (Point3, Vec3, Vec3) {
        let mut a = [Vec3::ZERO; 3];
        let mut w = [0.0; 3];
        for j in 0..=p {
            let idx = i - p + j;
            let wt = self.weights.as_ref().map_or(1.0, |w| w[idx]);
            let pw = self.ctrl[idx] * wt;
            for k in 0..=nd {
                a[k] += pw * n[k][j];
                w[k] += wt * n[k][j];
            }
        }
        let c = a[0] / w[0];
        if nd == 0 {
            return (c, Vec3::ZERO, Vec3::ZERO);
        }
        let d1 = (a[1] - c * w[1]) / w[0];
        if nd == 1 {
            return (c, d1, Vec3::ZERO);
        }
        let d2 = (a[2] - d1 * (2.0 * w[1]) - c * w[2]) / w[0];
        (c, d1, d2)
    }

    pub fn eval(&self, t: f64) -> Point3 {
        self.ders(t, 0).0
    }

    /// The point and its first derivative at `t`.
    pub(crate) fn d1(&self, t: f64) -> (Point3, Vec3) {
        let (c, d1, _) = self.ders(t, 1);
        (c, d1)
    }

    /// The curve's denominator at `t`, the blend of its weights there: 1 on a polynomial curve.
    /// The curve times it is a polynomial over each piece, where the curve itself is not.
    pub(crate) fn weight(&self, t: f64) -> f64 {
        let Some(weights) = &self.weights else {
            return 1.0;
        };
        let t = self.basis.wrap(t, self.period);
        let p = self.basis.degree;
        let i = self.basis.span(t);
        self.basis.ders_then(i, t, 0, |n| {
            (0..=p).map(|j| weights[i - p + j] * n[0][j]).sum()
        })
    }

    pub fn transformed(&self, t: &Transform) -> NurbsCurve {
        NurbsCurve {
            ctrl: self.ctrl.iter().map(|&p| t.point(p)).collect(),
            stationary: StationaryCache::default(),
            ..self.clone()
        }
    }

    /// Where each of x, y and z is stationary on [a, b], a stretch of one knot span: the roots of
    /// the derivative's coordinates, the three found from one pass of samples.
    fn stationary_on(&self, a: f64, b: f64) -> Stationary {
        let mut rates = |t: f64| {
            let d = self.d1(t).1;
            [d.dot(Vec3::X), d.dot(Vec3::Y), d.dot(Vec3::Z)]
        };
        numeric::roots_each(&mut rates, a, b, 12)
            .map(|ts| ts.into_iter().map(|t| (t, self.eval(t))).collect())
    }

    /// [`NurbsCurve::stationary_on`] over [a, b] when that is one whole knot span, bit for bit,
    /// kept from the first call: a periodic curve's span a period on is not one.
    fn stationary_on_span(&self, a: f64, b: f64) -> Option<&Stationary> {
        let i = self.basis.span(a);
        let u = &self.basis.knots;
        if a.to_bits() != u[i].to_bits() || b.to_bits() != u[i + 1].to_bits() {
            return None;
        }
        let spans = self
            .stationary
            .0
            .get_or_init(|| (0..self.basis.n).map(|_| OnceLock::new()).collect());
        Some(spans[i].get_or_init(|| self.stationary_on(a, b)))
    }

    /// Parameters strictly inside (t0, t1) where the curve's coordinate along `dir` is stationary.
    pub(crate) fn extremes(&self, dir: Vec3, t0: f64, t1: f64) -> Vec<f64> {
        let axis = [Vec3::X, Vec3::Y, Vec3::Z].iter().position(|&a| a == dir);
        let mut cuts = vec![t0];
        cuts.extend(self.breakpoints(t0, t1));
        cuts.push(t1);
        let mut out = Vec::new();
        let inside = |r: &f64| *r > t0 && *r < t1;
        for w in cuts.windows(2) {
            match axis.and_then(|k| Some(&self.stationary_on_span(w[0], w[1])?[k])) {
                Some(kept) => out.extend(kept.iter().map(|&(t, _)| t).filter(inside)),
                None => {
                    let mut f = |t: f64| self.d1(t).1.dot(dir);
                    out.extend(
                        numeric::roots(&mut f, w[0], w[1], 12)
                            .into_iter()
                            .filter(inside),
                    );
                }
            }
        }
        out
    }

    /// The tight box of the curve over [t0, t1]: its ends and every coordinate's extremes.
    pub(crate) fn bbox(&self, t0: f64, t1: f64) -> Aabb {
        let mut bb = Aabb::EMPTY;
        bb.add(self.eval(t0));
        bb.add(self.eval(t1));
        let mut cuts = vec![t0];
        cuts.extend(self.breakpoints(t0, t1));
        cuts.push(t1);
        for w in cuts.windows(2) {
            let fresh;
            let found = match self.stationary_on_span(w[0], w[1]) {
                Some(kept) => kept,
                None => {
                    fresh = self.stationary_on(w[0], w[1]);
                    &fresh
                }
            };
            for &(t, p) in found.iter().flatten() {
                if t > t0 && t < t1 {
                    bb.add(p);
                }
            }
        }
        bb
    }
}

/// A NURBS surface along one parameter, from [`NurbsSurface::iso`]: the span polynomials with the
/// held parameter put in where they evaluate the surface, and otherwise, for each control point
/// along the free parameter, the held parameter's weighted sums of the control net and of its
/// derivative.
pub struct IsoLine<'a> {
    surface: &'a NurbsSurface,
    along_v: bool,
    rows: Vec<([Vec3; 2], [f64; 2])>,
    /// Where the span polynomials evaluate the surface: the held parameter's span, its local
    /// parameter about the span's nearer end, that end, and the span's width.
    held: Option<(usize, f64, usize, f64)>,
    /// For each free span and end used so far, the span polynomial with the held parameter put
    /// in: the homogeneous point's series and its held derivative's, by power of the free one.
    series: RefCell<Vec<(usize, usize, IsoSeries)>>,
}

/// An iso line's homogeneous point and held derivative over one free span, by power of the free
/// parameter.
type IsoSeries = Vec<[[f64; 4]; 2]>;

impl IsoLine<'_> {
    /// The point and the partial derivatives along `u` and `v` at free parameter `t`.
    pub fn d1(&self, t: f64) -> (Point3, Vec3, Vec3) {
        let s = self.surface;
        let (free, period) = if self.along_v {
            (&s.bv, s.v_period)
        } else {
            (&s.bu, s.u_period)
        };
        let t = free.wrap(t, period);
        let (p, i) = (free.degree, free.span(t));
        if self.held.is_some() {
            return self.d1_poly(free, i, t);
        }
        free.ders_then(i, t, 1, |nf| self.d1_from(nf, p, i))
    }

    /// [`IsoLine::d1`] by a Horner pass over the span polynomial in the free parameter.
    fn d1_poly(&self, free: &Basis, i: usize, t: f64) -> (Point3, Vec3, Vec3) {
        let h = free.knots[i + 1] - free.knots[i];
        let x = (t - free.knots[i]) / h;
        let (x, end) = if x > 0.5 { (x - 1.0, 1) } else { (x, 0) };
        let mut cache = self.series.borrow_mut();
        let k = match cache.iter().position(|c| c.0 == i && c.1 == end) {
            Some(k) => k,
            None => {
                cache.push((i, end, self.collapse(i, end)));
                cache.len() - 1
            }
        };
        let c = &cache[k].2;
        let p = c.len() - 1;
        let (mut a0, mut a1, mut a2) = (c[p][0], c[p][1], [0.0; 4]);
        for k in (0..p).rev() {
            for j in 0..4 {
                a2[j] = a2[j] * x + a0[j];
                a0[j] = a0[j] * x + c[k][0][j];
                a1[j] = a1[j] * x + c[k][1][j];
            }
        }
        let v = |a: [f64; 4]| Vec3::new(a[0], a[1], a[2]);
        let pt = v(a0) / a0[3];
        let held = (v(a1) - pt * a1[3]) / a0[3];
        let along = (v(a2) - pt * a2[3]) / (a0[3] * h);
        if self.along_v {
            (pt, held, along)
        } else {
            (pt, along, held)
        }
    }

    /// The span polynomial over free span `i` about its end `end`, the held parameter put in.
    fn collapse(&self, i: usize, end: usize) -> IsoSeries {
        let s = self.surface;
        let (hi, hx, he, hw) = self.held.expect("a polynomial iso line");
        let (pu, pv) = (s.bu.degree, s.bv.degree);
        let c = if self.along_v {
            s.cell(hi, i, he, end)
        } else {
            s.cell(i, hi, end, he)
        };
        let (pf, ph) = if self.along_v { (pv, pu) } else { (pu, pv) };
        let at = |f: usize, g: usize| {
            if self.along_v {
                c[g * (pv + 1) + f]
            } else {
                c[f * (pv + 1) + g]
            }
        };
        (0..=pf)
            .map(|f| {
                let (mut val, mut der) = (at(f, ph), [0.0; 4]);
                for g in (0..ph).rev() {
                    let next = at(f, g);
                    for j in 0..4 {
                        der[j] = der[j] * hx + val[j];
                        val[j] = val[j] * hx + next[j];
                    }
                }
                [val, der.map(|d| d / hw)]
            })
            .collect()
    }

    #[allow(clippy::needless_range_loop)]
    fn d1_from(&self, nf: [&[f64]; 3], p: usize, i: usize) -> (Point3, Vec3, Vec3) {
        // a[0], w[0]: the point; [1]: along the held parameter; [2]: along the free one.
        let mut a = [Vec3::ZERO; 3];
        let mut w = [0.0; 3];
        for l in 0..=p {
            let (r, rw) = &self.rows[i - p + l];
            a[0] += r[0] * nf[0][l];
            w[0] += rw[0] * nf[0][l];
            a[1] += r[1] * nf[0][l];
            w[1] += rw[1] * nf[0][l];
            a[2] += r[0] * nf[1][l];
            w[2] += rw[0] * nf[1][l];
        }
        let pt = a[0] / w[0];
        let held = (a[1] - pt * w[1]) / w[0];
        let along = (a[2] - pt * w[2]) / w[0];
        if self.along_v {
            (pt, held, along)
        } else {
            (pt, along, held)
        }
    }
}

/// A rational (or polynomial) B-spline surface. Control points are stored row by row: index
/// `i·nv + j` is the `i`-th along `u` and the `j`-th along `v`.
#[derive(Clone, Debug)]
pub struct NurbsSurface {
    bu: Basis,
    bv: Basis,
    ctrl: Vec<Point3>,
    weights: Option<Vec<f64>>,
    u_period: Option<f64>,
    v_period: Option<f64>,
    /// The curve and steady motion that swept the surface, where it was made so and not changed.
    swept: Option<Arc<Swept>>,
    /// What is worked out from the fields above on first use, shared with the surface's copies.
    kept: Arc<Kept>,
}

/// What a surface works out of itself on first use. It depends only on the surface, so copies of
/// the surface share it, whichever of them first asks.
#[derive(Debug, Default)]
struct Kept {
    /// Sample points of the surface that seed global projection.
    seeds: OnceLock<SeedGrid>,
    /// Boxes over runs of the seeds, built once global projections have scanned them a few times.
    seed_blocks: OnceLock<SeedBlocks>,
    seed_scans: AtomicU32,
    /// Small patches whose boxes together hold the surface.
    patches: OnceLock<Vec<Patch>>,
    /// A hierarchy over the patches' boxes.
    patch_tree: OnceLock<Bvh>,
    /// Each patch as a Bézier patch.
    beziers: OnceLock<Vec<Bezier>>,
    /// The surface over each pair of spans as a polynomial.
    polys: OnceLock<SpanPolys>,
}

/// The seeds of global projection: points of the surface with their parameters, row by row along
/// `u`, `row` to a row, and the spacing they were laid out at along the surface.
#[derive(Clone, Debug)]
struct SeedGrid {
    seeds: Vec<(Vec2, Point3)>,
    spacing: f64,
    row: usize,
}

/// How many seeds one box of [`SeedBlocks`] holds.
const SEED_RUN: usize = 8;
/// How many scans of every seed a surface takes before its [`SeedBlocks`] are built: building
/// them costs about as much.
const SCANS_BEFORE_BLOCKS: u32 = 16;

/// A box round each row of the seeds and round each run of [`SEED_RUN`] seeds of a row, which
/// say where no seed can be among the nearest, and the seeds that have another within the
/// distance global projection takes as one point.
#[derive(Clone, Debug)]
struct SeedBlocks {
    rows: Vec<Aabb>,
    runs: Vec<Aabb>,
    runs_per_row: usize,
    crowded: Vec<bool>,
}

/// The square of the distance from `p` to box `b`, never more than the squared distance, as
/// `(q − p).norm2()` works it out, to any point `q` in the box: each coordinate's gap is no
/// longer than that of `q − p`, and squaring and summing keep the order.
fn box_gap2(b: &Aabb, p: Point3) -> f64 {
    Vec3::new(
        (b.min.x - p.x).max(p.x - b.max.x).max(0.0),
        (b.min.y - p.y).max(p.y - b.max.y).max(0.0),
        (b.min.z - p.z).max(p.z - b.max.z).max(0.0),
    )
    .norm2()
}

/// The surface over each pair of nonempty spans as a polynomial in the spans' local parameters,
/// in homogeneous coordinates: evaluating it is a Horner pass with no basis recurrence.
#[derive(Clone, Debug)]
struct SpanPolys {
    /// Each span's slot among the nonempty ones, by span number; `usize::MAX` for an empty span.
    u_slot: Vec<usize>,
    v_slot: Vec<usize>,
    /// Each nonempty span's functions as series about its start and about its end.
    u_basis: Vec<[Vec<Vec<f64>>; 2]>,
    v_basis: Vec<[Vec<Vec<f64>>; 2]>,
    /// Row by row along u, `(degree_u + 1)·(degree_v + 1)` coefficients each, by slot pair and
    /// the end each series is about.
    cells: Vec<OnceLock<Box<[[f64; 4]]>>>,
}

impl PartialEq for NurbsSurface {
    fn eq(&self, o: &NurbsSurface) -> bool {
        self.bu == o.bu && self.bv == o.bv && self.ctrl == o.ctrl && self.weights == o.weights
    }
}

impl NurbsSurface {
    pub fn new(
        degrees: (usize, usize),
        knots: (Vec<f64>, Vec<f64>),
        counts: (usize, usize),
        ctrl: Vec<Point3>,
        weights: Option<Vec<f64>>,
    ) -> Result<NurbsSurface> {
        if ctrl.len() != counts.0 * counts.1 {
            return invalid("B-spline surface control net has the wrong size");
        }
        let bu = Basis::new(degrees.0, knots.0, counts.0)?;
        let bv = Basis::new(degrees.1, knots.1, counts.1)?;
        check_weights(&weights, ctrl.len())?;
        let weights = weights.filter(|w| w.iter().any(|&x| x != w[0]));
        let mut s = NurbsSurface {
            bu,
            bv,
            ctrl,
            weights,
            u_period: None,
            v_period: None,
            swept: None,
            kept: Arc::default(),
        };
        let ((u0, u1), (v0, v1)) = (s.bu.domain(), s.bv.domain());
        let n = 12;
        let closed_u = (0..=n).all(|k| {
            let v = v0 + (v1 - v0) * k as f64 / n as f64;
            s.eval(u0, v).dist(s.eval(u1, v)) <= tol::LINEAR
        });
        let closed_v = (0..=n).all(|k| {
            let u = u0 + (u1 - u0) * k as f64 / n as f64;
            s.eval(u, v0).dist(s.eval(u, v1)) <= tol::LINEAR
        });
        s.u_period = closed_u.then_some(u1 - u0);
        s.v_period = closed_v.then_some(v1 - v0);
        Ok(s)
    }

    /// The surface exactly as [`NurbsSurface::new`] made it, from its parts and its periods.
    pub(crate) fn restore(
        degrees: (usize, usize),
        knots: (Vec<f64>, Vec<f64>),
        counts: (usize, usize),
        ctrl: Vec<Point3>,
        weights: Option<Vec<f64>>,
        periods: (Option<f64>, Option<f64>),
    ) -> Result<NurbsSurface> {
        if ctrl.len() != counts.0 * counts.1 {
            return invalid("B-spline surface control net has the wrong size");
        }
        check_weights(&weights, ctrl.len())?;
        Ok(NurbsSurface {
            bu: Basis::new(degrees.0, knots.0, counts.0)?,
            bv: Basis::new(degrees.1, knots.1, counts.1)?,
            ctrl,
            weights,
            u_period: periods.0,
            v_period: periods.1,
            swept: None,
            kept: Arc::default(),
        })
    }

    pub fn degrees(&self) -> (usize, usize) {
        (self.bu.degree, self.bv.degree)
    }

    pub fn knots(&self) -> (&[f64], &[f64]) {
        (&self.bu.knots, &self.bv.knots)
    }

    /// The number of control points along `u` and along `v`.
    pub fn counts(&self) -> (usize, usize) {
        (self.bu.n, self.bv.n)
    }

    pub fn control_points(&self) -> &[Point3] {
        &self.ctrl
    }

    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    pub fn domain(&self) -> ((f64, f64), (f64, f64)) {
        (self.bu.domain(), self.bv.domain())
    }

    pub fn u_period(&self) -> Option<f64> {
        self.u_period
    }

    pub fn v_period(&self) -> Option<f64> {
        self.v_period
    }

    /// The number of knot spans along `u` and along `v`.
    pub fn span_counts(&self) -> (usize, usize) {
        (self.bu.spans().len() - 1, self.bv.spans().len() - 1)
    }

    pub fn u_breakpoints(&self, u0: f64, u1: f64) -> Vec<f64> {
        self.bu.breakpoints(u0, u1, self.u_period)
    }

    pub fn v_breakpoints(&self, v0: f64, v1: f64) -> Vec<f64> {
        self.bv.breakpoints(v0, v1, self.v_period)
    }

    /// The point and its two first partial derivatives.
    pub fn d1(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        match self.bu.degree.max(self.bv.degree) {
            0..4 => self.d1_n::<4>(u, v),
            4..8 => self.d1_n::<8>(u, v),
            _ => self.d1_n::<M>(u, v),
        }
    }

    #[allow(clippy::needless_range_loop)]
    fn d1_n<const N: usize>(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        let u = self.bu.wrap(u, self.u_period);
        let v = self.bv.wrap(v, self.v_period);
        let (pu, pv) = (self.bu.degree, self.bv.degree);
        let (iu, iv) = (self.bu.span(u), self.bv.span(v));
        let nu = self.bu.ders::<N>(iu, u, 1);
        let nv = self.bv.ders::<N>(iv, v, 1);
        let mut a = [Vec3::ZERO; 3];
        let mut w = [0.0; 3];
        let stride = self.bv.n;
        for k in 0..=pu {
            let row = (iu - pu + k) * stride + iv - pv;
            // The row summed along v for the value and its v derivative.
            let (mut r0, mut r1) = (Vec3::ZERO, Vec3::ZERO);
            let (mut w0, mut w1) = (0.0, 0.0);
            for l in 0..=pv {
                let idx = row + l;
                let wt = self.weights.as_ref().map_or(1.0, |w| w[idx]);
                let pw = self.ctrl[idx] * wt;
                r0 += pw * nv[0][l];
                w0 += wt * nv[0][l];
                r1 += pw * nv[1][l];
                w1 += wt * nv[1][l];
            }
            a[0] += r0 * nu[0][k];
            w[0] += w0 * nu[0][k];
            a[1] += r0 * nu[1][k];
            w[1] += w0 * nu[1][k];
            a[2] += r1 * nu[0][k];
            w[2] += w1 * nu[0][k];
        }
        let s = a[0] / w[0];
        (s, (a[1] - s * w[1]) / w[0], (a[2] - s * w[2]) / w[0])
    }

    /// The point alone: what [`NurbsSurface::d1`] gives, without the derivatives' sums.
    pub fn eval(&self, u: f64, v: f64) -> Point3 {
        let u = self.bu.wrap(u, self.u_period);
        let v = self.bv.wrap(v, self.v_period);
        let (iu, iv) = (self.bu.span(u), self.bv.span(v));
        self.bu.ders_then(iu, u, 0, |nu| {
            self.bv
                .ders_then(iv, v, 0, |nv| self.point_from(iu, iv, nu[0], nv[0]))
        })
    }

    /// The point over spans `iu` and `iv` from the basis values there, summed as
    /// [`NurbsSurface::d1`] sums its point.
    fn point_from(&self, iu: usize, iv: usize, nu: &[f64], nv: &[f64]) -> Point3 {
        let (pu, pv) = (self.bu.degree, self.bv.degree);
        let (mut a, mut w) = (Vec3::ZERO, 0.0);
        for (k, fu) in nu.iter().enumerate().take(pu + 1) {
            let row = (iu - pu + k) * self.bv.n + iv - pv;
            let (mut r0, mut w0) = (Vec3::ZERO, 0.0);
            for (l, fv) in nv.iter().enumerate().take(pv + 1) {
                let wt = self.weights.as_ref().map_or(1.0, |w| w[row + l]);
                r0 += self.ctrl[row + l] * wt * *fv;
                w0 += wt * *fv;
            }
            a += r0 * *fu;
            w += w0 * *fu;
        }
        a / w
    }

    /// [`NurbsSurface::eval`] at every pair of `us` and `vs`, row by row along `u`, each
    /// parameter's basis worked out once.
    fn grid(&self, us: &[f64], vs: &[f64]) -> Vec<Point3> {
        let basis = |b: &Basis, period: Option<f64>, ts: &[f64]| -> Vec<(usize, [f64; M])> {
            ts.iter()
                .map(|&t| {
                    let t = b.wrap(t, period);
                    let i = b.span(t);
                    (
                        i,
                        b.ders_then(i, t, 0, |n| {
                            std::array::from_fn(|j| n[0].get(j).map_or(0.0, |x| *x))
                        }),
                    )
                })
                .collect()
        };
        let nu = basis(&self.bu, self.u_period, us);
        let nv = basis(&self.bv, self.v_period, vs);
        let mut out = Vec::with_capacity(us.len() * vs.len());
        for (iu, fu) in &nu {
            for (iv, fv) in &nv {
                out.push(self.point_from(*iu, *iv, fu, fv));
            }
        }
        out
    }

    /// [`NurbsSurface::d1`] from the span polynomials where they evaluate the surface: a Horner
    /// pass in place of the basis recurrence, equal to rounding.
    pub(crate) fn d1_by_polys(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        if !self.by_polys() {
            return self.d1(u, v);
        }
        let (h, n) = self.poly_ders(u, v, 1);
        let [s, su, sv, ..] = rational(&h, n);
        (s, su, sv)
    }

    /// The surface along one parameter with the other held at `at`, as a curve of its own: along
    /// `v` with `u = at` when `along_v`, else along `u` with `v = at`. Each control point of the
    /// curve is a row of the net, in homogeneous coordinates, blended by the held parameter's basis
    /// there, so the curve is the surface exactly.
    pub fn iso_curve(&self, along_v: bool, at: f64) -> NurbsCurve {
        let (held, free, period) = if along_v {
            (&self.bu, &self.bv, self.u_period)
        } else {
            (&self.bv, &self.bu, self.v_period)
        };
        let at = held.wrap(at, period);
        let (p, i) = (held.degree, held.span(at));
        let funs = &held.ders::<M>(i, at, 0)[0];
        let index = |h: usize, f: usize| {
            if along_v {
                h * self.bv.n + f
            } else {
                f * self.bv.n + h
            }
        };
        let mut ctrl = Vec::with_capacity(free.n);
        let mut weights = Vec::with_capacity(free.n);
        for f in 0..free.n {
            let (mut sum, mut w) = (Vec3::ZERO, 0.0);
            for (k, n) in funs.iter().enumerate().take(p + 1) {
                let idx = index(i - p + k, f);
                let wk = self.weights.as_ref().map_or(1.0, |ws| ws[idx]) * n;
                sum += self.ctrl[idx] * wk;
                w += wk;
            }
            ctrl.push(sum / w);
            weights.push(w);
        }
        NurbsCurve::new(free.degree, free.knots.clone(), ctrl, Some(weights))
            .expect("a row of a valid net is a valid curve")
    }

    /// The surface along one parameter with the other held at `at`: [`IsoLine::d1`] gives what
    /// [`NurbsSurface::d1`] does there, the held parameter's basis worked out once.
    #[allow(clippy::needless_range_loop)]
    pub fn iso(&self, along_v: bool, at: f64) -> IsoLine<'_> {
        let (held, free) = if along_v {
            (&self.bu, &self.bv)
        } else {
            (&self.bv, &self.bu)
        };
        let at = held.wrap(
            at,
            if along_v {
                self.u_period
            } else {
                self.v_period
            },
        );
        let (p, i) = (held.degree, held.span(at));
        if self.by_polys() {
            let w = held.knots[i + 1] - held.knots[i];
            let x = (at - held.knots[i]) / w;
            let (x, end) = if x > 0.5 { (x - 1.0, 1) } else { (x, 0) };
            return IsoLine {
                surface: self,
                along_v,
                rows: Vec::new(),
                held: Some((i, x, end, w)),
                series: RefCell::new(Vec::new()),
            };
        }
        let rows = held.ders_then(i, at, 1, |nh| {
            let mut rows = Vec::with_capacity(free.n);
            for j in 0..free.n {
                let mut r = [Vec3::ZERO; 2];
                let mut w = [0.0; 2];
                for k in 0..=p {
                    let idx = if along_v {
                        (i - p + k) * self.bv.n + j
                    } else {
                        j * self.bv.n + i - p + k
                    };
                    let wt = self.weights.as_ref().map_or(1.0, |w| w[idx]);
                    let pw = self.ctrl[idx] * wt;
                    for m in 0..2 {
                        r[m] += pw * nh[m][k];
                        w[m] += wt * nh[m][k];
                    }
                }
                rows.push((r, w));
            }
            rows
        });
        IsoLine {
            surface: self,
            along_v,
            rows,
            held: None,
            series: RefCell::new(Vec::new()),
        }
    }

    /// The control points that move the point at `(u, v)`, as `(index, R)`: the point moves by
    /// `R` times a control point's displacement (the rational basis function).
    pub fn basis(&self, u: f64, v: f64) -> Vec<(usize, f64)> {
        let u = self.bu.wrap(u, self.u_period);
        let v = self.bv.wrap(v, self.v_period);
        let (pu, pv) = (self.bu.degree, self.bv.degree);
        let (iu, iv) = (self.bu.span(u), self.bv.span(v));
        let nu = self.bu.ders::<M>(iu, u, 0);
        let nv = self.bv.ders::<M>(iv, v, 0);
        let stride = self.bv.n;
        let mut out = Vec::with_capacity((pu + 1) * (pv + 1));
        let mut total = 0.0;
        for (k, bu) in nu[0].iter().enumerate().take(pu + 1) {
            for (l, bv) in nv[0].iter().enumerate().take(pv + 1) {
                let idx = (iu - pu + k) * stride + iv - pv + l;
                let r = self.weights.as_ref().map_or(1.0, |w| w[idx]) * bu * bv;
                total += r;
                out.push((idx, r));
            }
        }
        for o in &mut out {
            o.1 /= total;
        }
        out
    }

    /// The curve and steady motion the surface was swept by, where it keeps them: it does through
    /// a copy and a rigid move, which moves them with it, and not once its shape may have changed.
    pub(crate) fn swept(&self) -> Option<&Swept> {
        self.swept.as_deref()
    }

    /// The surface, known to be `generator` at its `v` carried through its `u` of `motion`.
    pub(crate) fn swept_by(mut self, swept: Swept) -> NurbsSurface {
        self.swept = Some(Arc::new(swept));
        self
    }

    /// The same surface with its control points moved to `ctrl`.
    pub fn with_control_points(&self, ctrl: Vec<Point3>) -> Result<NurbsSurface> {
        NurbsSurface::new(
            (self.bu.degree, self.bv.degree),
            (self.bu.knots.clone(), self.bv.knots.clone()),
            (self.bu.n, self.bv.n),
            ctrl,
            self.weights.clone(),
        )
    }

    pub fn transformed(&self, t: &Transform) -> NurbsSurface {
        // A rigid motion carries what is already worked out with it: the seeds are points of the
        // surface, and each span polynomial's homogeneous coefficients map as the control points
        // do, `(R·c + t·w, w)`.
        let seeds = OnceLock::new();
        if let Some(grid) = self.kept.seeds.get() {
            let _ = seeds.set(SeedGrid {
                seeds: grid.seeds.iter().map(|&(uv, p)| (uv, t.point(p))).collect(),
                ..*grid
            });
        }
        let polys = OnceLock::new();
        if let Some(sp) = self.kept.polys.get() {
            let cells = sp
                .cells
                .iter()
                .map(|cell| {
                    let out = OnceLock::new();
                    if let Some(c) = cell.get() {
                        let moved: Box<[[f64; 4]]> = c
                            .iter()
                            .map(|&[x, y, z, w]| {
                                let v = t.vector(Vec3::new(x, y, z)) + t.t * w;
                                [v.x, v.y, v.z, w]
                            })
                            .collect();
                        let _ = out.set(moved);
                    }
                    out
                })
                .collect();
            let _ = polys.set(SpanPolys {
                u_slot: sp.u_slot.clone(),
                v_slot: sp.v_slot.clone(),
                u_basis: sp.u_basis.clone(),
                v_basis: sp.v_basis.clone(),
                cells,
            });
        }
        NurbsSurface {
            bu: self.bu.clone(),
            bv: self.bv.clone(),
            ctrl: self.ctrl.iter().map(|&p| t.point(p)).collect(),
            weights: self.weights.clone(),
            u_period: self.u_period,
            v_period: self.v_period,
            swept: self.swept.as_ref().map(|s| Arc::new(s.transformed(t))),
            kept: Arc::new(Kept {
                seeds,
                polys,
                ..Kept::default()
            }),
        }
    }

    /// Works out now what projection and evaluation otherwise build on first use, so that copies
    /// moved rigidly from this one (see [`NurbsSurface::transformed`]) start with it.
    pub fn prepare(&self) {
        self.seed_grid();
        if self.by_polys() {
            let polys = self.polys();
            for (iu, &su) in polys.u_slot.iter().enumerate() {
                for (iv, &sv) in polys.v_slot.iter().enumerate() {
                    if su != usize::MAX && sv != usize::MAX {
                        for e in 0..4 {
                            self.cell(iu, iv, e / 2, e % 2);
                        }
                    }
                }
            }
        }
    }

    /// The box of the control points, which holds the surface (its weights are positive).
    pub fn hull_box(&self) -> crate::math::Aabb {
        let mut bb = crate::math::Aabb::EMPTY;
        for &p in &self.ctrl {
            bb.add(p);
        }
        bb
    }

    /// Boxes whose union holds the surface over its whole domain: the knot spans, each cut into
    /// pieces about a sixteenth of the surface's longer side, and each piece's box the hull of its
    /// own control points after knot insertion (weights are positive, so the hull holds it).
    pub fn patches(&self) -> &[Patch] {
        self.kept.patches.get_or_init(|| {
            let lengths = |along_u: bool| -> f64 {
                let (a0, a1) = if along_u {
                    self.bu.domain()
                } else {
                    self.bv.domain()
                };
                let (b0, b1) = if along_u {
                    self.bv.domain()
                } else {
                    self.bu.domain()
                };
                (0..=4)
                    .map(|j| {
                        let y = b0 + (b1 - b0) * j as f64 / 4.0;
                        let at = |x: f64| {
                            if along_u {
                                self.eval(x, y)
                            } else {
                                self.eval(y, x)
                            }
                        };
                        (0..32)
                            .map(|i| {
                                let x = |k: usize| a0 + (a1 - a0) * k as f64 / 32.0;
                                at(x(i)).dist(at(x(i + 1)))
                            })
                            .sum::<f64>()
                    })
                    .fold(0.0, f64::max)
            };
            let (lu, lv) = (lengths(true), lengths(false));
            let h = lu.max(lv) / 16.0;
            let splits = |b: &Basis, len: f64| -> Vec<f64> {
                let spans = b.spans();
                let (a, z) = b.domain();
                let mut out = Vec::new();
                for w in spans.windows(2) {
                    let share = len * (w[1] - w[0]) / (z - a);
                    let k = if h > 0.0 {
                        ((share / h).ceil() as usize).clamp(1, 16)
                    } else {
                        1
                    };
                    out.extend((1..k).map(|i| w[0] + (w[1] - w[0]) * i as f64 / k as f64));
                }
                out
            };
            let (nu, nv) = (self.bu.n, self.bv.n);
            let mut net: Vec<Vec<[f64; 4]>> = (0..nu)
                .map(|i| {
                    (0..nv)
                        .map(|j| {
                            let k = i * nv + j;
                            let w = self.weights.as_ref().map_or(1.0, |w| w[k]);
                            let p = self.ctrl[k];
                            [p.x * w, p.y * w, p.z * w, w]
                        })
                        .collect()
                })
                .collect();
            let mut kv = self.bv.knots.clone();
            let mut rows: Vec<Vec<[f64; 4]>> = Vec::new();
            for (r, row) in net.iter().enumerate() {
                let mut knots = self.bv.knots.clone();
                let mut pts = row.clone();
                for &t in &splits(&self.bv, lv) {
                    insert_knot(&mut knots, self.bv.degree, &mut pts, t);
                }
                if r == 0 {
                    kv = knots;
                }
                rows.push(pts);
            }
            net = rows;
            let mut ku = self.bu.knots.clone();
            let cols = net[0].len();
            let mut columns: Vec<Vec<[f64; 4]>> = Vec::new();
            for c in 0..cols {
                let mut knots = self.bu.knots.clone();
                let mut pts: Vec<[f64; 4]> = net.iter().map(|row| row[c]).collect();
                for &t in &splits(&self.bu, lu) {
                    insert_knot(&mut knots, self.bu.degree, &mut pts, t);
                }
                if c == 0 {
                    ku = knots;
                }
                columns.push(pts);
            }
            let (pu, pv) = (self.bu.degree, self.bv.degree);
            let mut out = Vec::new();
            for i in pu..ku.len() - pu - 1 {
                if ku[i + 1] <= ku[i] {
                    continue;
                }
                for j in pv..kv.len() - pv - 1 {
                    if kv[j + 1] <= kv[j] {
                        continue;
                    }
                    let mut bbox = Aabb::EMPTY;
                    for col in &columns[j - pv..=j] {
                        for q in &col[i - pu..=i] {
                            bbox.add(Vec3::new(q[0] / q[3], q[1] / q[3], q[2] / q[3]));
                        }
                    }
                    out.push(Patch {
                        bbox,
                        lo: Vec2::new(ku[i], kv[j]),
                        hi: Vec2::new(ku[i + 1], kv[j + 1]),
                    });
                }
            }
            out
        })
    }

    /// The surface at `(u, v)` in homogeneous coordinates: the weighted sum of its control points
    /// and the sum of the weights. A polynomial over each pair of spans, where the point is not.
    fn homogeneous(&self, u: f64, v: f64) -> [f64; 4] {
        let u = self.bu.wrap(u, self.u_period);
        let v = self.bv.wrap(v, self.v_period);
        let (pu, pv) = (self.bu.degree, self.bv.degree);
        let (iu, iv) = (self.bu.span(u), self.bv.span(v));
        let nu = self.bu.ders::<M>(iu, u, 0);
        let nv = self.bv.ders::<M>(iv, v, 0);
        let stride = self.bv.n;
        let mut out = [0.0; 4];
        for (k, bu) in nu[0].iter().enumerate().take(pu + 1) {
            for (l, bv) in nv[0].iter().enumerate().take(pv + 1) {
                let idx = (iu - pu + k) * stride + iv - pv + l;
                let w = self.weights.as_ref().map_or(1.0, |w| w[idx]) * bu * bv;
                let p = self.ctrl[idx];
                out = [
                    out[0] + w * p.x,
                    out[1] + w * p.y,
                    out[2] + w * p.z,
                    out[3] + w,
                ];
            }
        }
        out
    }

    /// A hierarchy over the boxes of [`NurbsSurface::patches`], whose indices its queries return.
    pub(crate) fn patch_tree(&self) -> &Bvh {
        self.kept
            .patch_tree
            .get_or_init(|| Bvh::of_boxes(self.patches().iter().map(|q| q.bbox).collect()))
    }

    /// Each of [`NurbsSurface::patches`] as a Bézier patch. A patch lies within one span of each
    /// parameter, where the surface is a polynomial in homogeneous coordinates, so its values at a
    /// grid of points fix its Bézier control points.
    pub(crate) fn beziers(&self) -> &[Bezier] {
        self.kept.beziers.get_or_init(|| {
            let (pu, pv) = (self.bu.degree, self.bv.degree);
            let (mu, mv) = (bernstein_from_values(pu), bernstein_from_values(pv));
            self.patches()
                .iter()
                .map(|patch| {
                    // Values a hair inside the patch's far sides, on its own spans' polynomial.
                    let node = |lo: f64, hi: f64, j: usize, n: usize| {
                        lo + (hi - lo) * (j as f64 / n.max(1) as f64).min(1.0 - 1e-12)
                    };
                    let values: Vec<Vec<[f64; 4]>> = (0..=pu)
                        .map(|i| {
                            let u = node(patch.lo.x, patch.hi.x, i, pu);
                            (0..=pv)
                                .map(|j| self.homogeneous(u, node(patch.lo.y, patch.hi.y, j, pv)))
                                .collect()
                        })
                        .collect();
                    // Along v in each row, then along u in each column.
                    let across: Vec<Vec<[f64; 4]>> = values
                        .iter()
                        .map(|row| {
                            (0..=pv)
                                .map(|k| {
                                    std::array::from_fn(|c| {
                                        (0..=pv).map(|j| mv[k][j] * row[j][c]).sum()
                                    })
                                })
                                .collect()
                        })
                        .collect();
                    let mut net = Vec::with_capacity((pu + 1) * (pv + 1));
                    for weights in &mu {
                        for col in 0..=pv {
                            net.push(std::array::from_fn(|c| {
                                weights
                                    .iter()
                                    .zip(&across)
                                    .map(|(w, row)| w * row[col][c])
                                    .sum()
                            }));
                        }
                    }
                    Bezier {
                        lo: patch.lo,
                        hi: patch.hi,
                        rows: pu + 1,
                        cols: pv + 1,
                        net,
                    }
                })
                .collect()
        })
    }

    fn seeds(&self) -> &[(Vec2, Point3)] {
        &self.seed_grid().seeds
    }

    /// The seeds, and the spacing they were laid out at along the surface.
    fn seed_grid(&self) -> &SeedGrid {
        self.kept.seeds.get_or_init(|| {
            let samples = |b: &Basis, total: usize| -> Vec<f64> {
                let spans = b.spans();
                let per = (24 / (spans.len() - 1)).clamp(2, 6);
                let per = per.max(total.div_ceil(spans.len() - 1)).min(256);
                let mut out = Vec::new();
                for w in spans.windows(2) {
                    for k in 0..per {
                        out.push(w[0] + (w[1] - w[0]) * k as f64 / per as f64);
                    }
                }
                out.push(spans[spans.len() - 1]);
                out
            };
            let (us, vs) = (samples(&self.bu, 0), samples(&self.bv, 0));
            let pts = self.grid(&us, &vs);
            // A long narrow patch (a twisted loft's side between close corners) needs seeds about
            // as close along its length as across it, or the nearest seed can lie in the basin of
            // a farther point.
            let (nu, nv) = (us.len(), vs.len());
            let at = |i: usize, j: usize| pts[i * nv + j];
            let lu = (0..nv)
                .map(|j| (1..nu).map(|i| at(i - 1, j).dist(at(i, j))).sum::<f64>())
                .fold(0.0, f64::max);
            let lv = (0..nu)
                .map(|i| (1..nv).map(|j| at(i, j - 1).dist(at(i, j))).sum::<f64>())
                .fold(0.0, f64::max);
            // At most about 2048 seeds, since a global projection measures many of them.
            let h = (lu.min(lv) / 4.0).max((lu * lv / 2048.0).sqrt());
            let (us, vs, pts) = if h > 0.0 {
                let n = |l: f64| ((l / h).ceil() as usize).min(256);
                let (us, vs) = (samples(&self.bu, n(lu)), samples(&self.bv, n(lv)));
                let pts = self.grid(&us, &vs);
                (us, vs, pts)
            } else {
                (us, vs, pts)
            };
            let params = us
                .iter()
                .flat_map(|&u| vs.iter().map(move |&v| Vec2::new(u, v)));
            SeedGrid {
                seeds: params.zip(pts).collect(),
                spacing: h,
                row: vs.len(),
            }
        })
    }

    fn seed_blocks(&self) -> &SeedBlocks {
        self.kept.seed_blocks.get_or_init(|| {
            let grid = self.seed_grid();
            let hull = |pts: &[(Vec2, Point3)]| {
                let mut bb = Aabb::EMPTY;
                pts.iter().for_each(|&(_, q)| bb.add(q));
                bb
            };
            let rows = grid.seeds.chunks(grid.row).map(hull).collect();
            let runs = grid
                .seeds
                .chunks(grid.row)
                .flat_map(|row| row.chunks(SEED_RUN).map(hull))
                .collect();
            // Seeds in order along the axis they spread widest on: one within `same` of another
            // is within its root of it along that axis, and so a neighbour or nearly one.
            let same = self.same_seed();
            let all = hull(&grid.seeds);
            let size = all.max - all.min;
            let axis = (0..3)
                .max_by(|&a, &b| size[a].total_cmp(&size[b]))
                .unwrap_or(0);
            let mut order: Vec<usize> = (0..grid.seeds.len()).collect();
            order.sort_by(|&a, &b| grid.seeds[a].1[axis].total_cmp(&grid.seeds[b].1[axis]));
            let reach = 2.0 * same.sqrt();
            let mut crowded = vec![false; grid.seeds.len()];
            for (k, &i) in order.iter().enumerate() {
                let q = grid.seeds[i].1;
                for &j in &order[k + 1..] {
                    let r = grid.seeds[j].1;
                    if r[axis] - q[axis] > reach {
                        break;
                    }
                    if (r - q).norm2() <= same {
                        crowded[i] = true;
                        crowded[j] = true;
                    }
                }
            }
            SeedBlocks {
                rows,
                runs,
                runs_per_row: grid.row.div_ceil(SEED_RUN),
                crowded,
            }
        })
    }

    /// The squared distance within which two seeds are one point: a pole's row of seeds is one
    /// point, and a mapped surface's pole row agrees only to rounding.
    fn same_seed(&self) -> f64 {
        (1e-6 * self.seed_grid().spacing).powi(2)
    }

    /// The four nearest seeds at distinct points, nearest first, each with its squared distance
    /// from `p`, found by measuring every seed in turn; fewer are left at infinity.
    fn nearest_seeds_in_turn(&self, p: Point3) -> [(f64, Vec2); 4] {
        let same = self.same_seed();
        let mut near = [(f64::INFINITY, Vec2::ZERO); 4];
        let mut at = [Point3::ZERO; 4];
        for &(uv, q) in self.seeds() {
            let d2 = (q - p).norm2();
            let taken = |i: usize| near[i].0.is_finite() && (at[i] - q).norm2() <= same;
            if d2 < near[3].0 && !(0..4).any(taken) {
                let mut i = 3;
                while i > 0 && near[i - 1].0 > d2 {
                    near[i] = near[i - 1];
                    at[i] = at[i - 1];
                    i -= 1;
                }
                near[i] = (d2, uv);
                at[i] = q;
            }
        }
        near
    }

    /// [`NurbsSurface::nearest_seeds_in_turn`] measuring only the seeds whose boxes come as near
    /// as the fourth nearest so far. Where no seed that near has another at its own point, the
    /// scan in turn keeps the four seeds least in distance, the earlier of two equally far, in
    /// that order, since it never passes one over for being at a kept seed's point: those four
    /// are found here in any order. `None` where a seed that near does have another at its
    /// point, and the scan in turn must choose between them.
    fn nearest_seeds_by_blocks(&self, p: Point3) -> Option<[(f64, Vec2); 4]> {
        let grid = self.seed_grid();
        let blocks = self.seed_blocks();
        // Squared distance and seed number, least first; an unfilled place loses to every seed.
        let mut near = [(f64::INFINITY, 0usize); 4];
        let mut crowded = f64::INFINITY;
        let mut row = |r: usize, near: &mut [(f64, usize); 4]| {
            let runs = &blocks.runs[r * blocks.runs_per_row..(r + 1) * blocks.runs_per_row];
            for (k, run) in runs.iter().enumerate() {
                if box_gap2(run, p) > near[3].0 {
                    continue;
                }
                let first = r * grid.row + k * SEED_RUN;
                let last = (first + SEED_RUN).min((r + 1) * grid.row);
                for n in first..last {
                    let d2 = (grid.seeds[n].1 - p).norm2();
                    if blocks.crowded[n] {
                        crowded = crowded.min(d2);
                    }
                    if d2 < near[3].0 || (d2 == near[3].0 && n < near[3].1) {
                        let mut i = 3;
                        while i > 0 && (near[i - 1].0, near[i - 1].1) > (d2, n) {
                            near[i] = near[i - 1];
                            i -= 1;
                        }
                        near[i] = (d2, n);
                    }
                }
            }
        };
        // The row whose box is nearest first, so the rest are measured against a near fourth.
        let gaps: Vec<f64> = blocks.rows.iter().map(|b| box_gap2(b, p)).collect();
        let first = (0..gaps.len()).min_by(|&a, &b| gaps[a].total_cmp(&gaps[b]))?;
        row(first, &mut near);
        for (r, &gap) in gaps.iter().enumerate() {
            if r != first && gap <= near[3].0 {
                row(r, &mut near);
            }
        }
        if crowded <= near[3].0 {
            return None;
        }
        let uv = |d2: f64, n: usize| {
            if d2.is_finite() {
                grid.seeds[n].0
            } else {
                Vec2::ZERO
            }
        };
        Some(near.map(|(d2, n)| (d2, uv(d2, n))))
    }

    /// The four nearest seeds at distinct points, nearest first, each with its squared distance
    /// from `p`; fewer are left at infinity.
    fn nearest_seeds(&self, p: Point3) -> [(f64, Vec2); 4] {
        if self.kept.seed_blocks.get().is_none()
            && self.kept.seed_scans.fetch_add(1, Ordering::Relaxed) < SCANS_BEFORE_BLOCKS
        {
            return self.nearest_seeds_in_turn(p);
        }
        self.nearest_seeds_by_blocks(p)
            .unwrap_or_else(|| self.nearest_seeds_in_turn(p))
    }

    /// The point, its first partial derivatives and its second, `(S, Su, Sv, Suu, Suv, Svv)`.
    pub(crate) fn d2(&self, u: f64, v: f64) -> [Vec3; 6] {
        if !self.by_polys() {
            return if self.bu.degree.max(self.bv.degree) < 8 {
                self.d2_n::<8>(u, v)
            } else {
                self.d2_n::<M>(u, v)
            };
        }
        let (h, n) = self.poly_ders(u, v, 2);
        rational(&h, n)
    }

    /// Whether the span polynomials evaluate this surface: past degree 6 their monomial
    /// coefficients cancel by more than the precision the kernel keeps, and the basis recurrence
    /// does instead.
    fn by_polys(&self) -> bool {
        self.bu.degree.max(self.bv.degree) <= MAX_POLY_DEGREE
    }

    fn polys(&self) -> &SpanPolys {
        self.kept.polys.get_or_init(|| {
            let slots = |b: &Basis| {
                let mut slot = vec![usize::MAX; b.n];
                let mut basis = Vec::new();
                for (i, s) in slot.iter_mut().enumerate().skip(b.degree) {
                    if b.knots[i + 1] > b.knots[i] {
                        *s = basis.len();
                        basis.push([b.span_poly(i, false), b.span_poly(i, true)]);
                    }
                }
                (slot, basis)
            };
            let (u_slot, u_basis) = slots(&self.bu);
            let (v_slot, v_basis) = slots(&self.bv);
            let cells = (0..4 * u_basis.len() * v_basis.len())
                .map(|_| OnceLock::new())
                .collect();
            SpanPolys {
                u_slot,
                v_slot,
                u_basis,
                v_basis,
                cells,
            }
        })
    }

    /// The homogeneous coefficients over spans `iu` and `iv`, about the ends `eu` and `ev` (each
    /// 0 for the start, 1 for the end), row by row along u.
    #[allow(clippy::needless_range_loop)]
    fn cell(&self, iu: usize, iv: usize, eu: usize, ev: usize) -> &[[f64; 4]] {
        let polys = self.polys();
        let (su, sv) = (polys.u_slot[iu], polys.v_slot[iv]);
        let at = (2 * su + eu) * 2 * polys.v_basis.len() + 2 * sv + ev;
        polys.cells[at].get_or_init(|| {
            let (pu, pv) = (self.bu.degree, self.bv.degree);
            let (cu, cv) = (&polys.u_basis[su][eu], &polys.v_basis[sv][ev]);
            let stride = self.bv.n;
            // Each control row summed against the v polynomials, then the rows against u's.
            let mut rows = vec![[0.0; 4]; (pu + 1) * (pv + 1)];
            for a in 0..=pu {
                for b in 0..=pv {
                    let idx = (iu - pu + a) * stride + iv - pv + b;
                    let w = self.weights.as_ref().map_or(1.0, |w| w[idx]);
                    let p = self.ctrl[idx];
                    let pw = [p.x * w, p.y * w, p.z * w, w];
                    for l in 0..=pv {
                        let c = cv[b][l];
                        for (r, x) in rows[a * (pv + 1) + l].iter_mut().zip(pw) {
                            *r += x * c;
                        }
                    }
                }
            }
            let mut out = vec![[0.0; 4]; (pu + 1) * (pv + 1)];
            for k in 0..=pu {
                for a in 0..=pu {
                    let c = cu[a][k];
                    for l in 0..=pv {
                        let row = rows[a * (pv + 1) + l];
                        for (o, x) in out[k * (pv + 1) + l].iter_mut().zip(row) {
                            *o += x * c;
                        }
                    }
                }
            }
            out.into_boxed_slice()
        })
    }

    /// The homogeneous point and its derivatives up to `order` (at most 2), in the order
    /// `(0,0), (1,0), (0,1), (2,0), (1,1), (0,2)`, and how many of them are filled.
    fn poly_ders(&self, u: f64, v: f64, order: usize) -> ([[f64; 4]; 6], usize) {
        let u = self.bu.wrap(u, self.u_period);
        let v = self.bv.wrap(v, self.v_period);
        let (pu, pv) = (self.bu.degree, self.bv.degree);
        let (iu, iv) = (self.bu.span(u), self.bv.span(v));
        let (ku, kv) = (&self.bu.knots, &self.bv.knots);
        let (hu, hv) = (ku[iu + 1] - ku[iu], kv[iv + 1] - kv[iv]);
        // The series about the nearer end of each span: exact at both ends, and |s| ≤ ½.
        let local = |x: f64| if x > 0.5 { (x - 1.0, 1) } else { (x, 0) };
        let ((s, eu), (t, ev)) = (local((u - ku[iu]) / hu), local((v - kv[iv]) / hv));
        let c = self.cell(iu, iv, eu, ev);
        let add = |a: [f64; 4], b: [f64; 4], k: f64| std::array::from_fn(|i| a[i] * k + b[i]);
        // Each row along u evaluated in t with its first two t derivatives (the second halved).
        let mut r = [[[0.0; 4]; 3]; MAX_POLY_DEGREE + 1];
        for k in 0..=pu {
            let row = &c[k * (pv + 1)..(k + 1) * (pv + 1)];
            let (mut p0, mut p1, mut p2) = (row[pv], [0.0; 4], [0.0; 4]);
            for l in (0..pv).rev() {
                if order >= 2 {
                    p2 = add(p2, p1, t);
                }
                if order >= 1 {
                    p1 = add(p1, p0, t);
                }
                p0 = add(p0, row[l], t);
            }
            r[k] = [p0, p1, std::array::from_fn(|i| 2.0 * p2[i])];
        }
        // The rows' values, t derivatives and second t derivatives, each in s.
        let horner = |j: usize, want: usize| {
            let (mut q0, mut q1, mut q2) = (r[pu][j], [0.0; 4], [0.0; 4]);
            for k in (0..pu).rev() {
                if want >= 2 {
                    q2 = add(q2, q1, s);
                }
                if want >= 1 {
                    q1 = add(q1, q0, s);
                }
                q0 = add(q0, r[k][j], s);
            }
            (q0, q1, std::array::from_fn(|i| 2.0 * q2[i]))
        };
        let scale = |a: [f64; 4], k: f64| std::array::from_fn(|i| a[i] * k);
        let mut out = [[0.0; 4]; 6];
        let (a, a_s, a_ss) = horner(0, order);
        out[0] = a;
        if order == 0 {
            return (out, 1);
        }
        let (a_t, a_st, _) = horner(1, order - 1);
        out[1] = scale(a_s, 1.0 / hu);
        out[2] = scale(a_t, 1.0 / hv);
        if order == 1 {
            return (out, 3);
        }
        let (a_tt, _, _) = horner(2, 0);
        out[3] = scale(a_ss, 1.0 / (hu * hu));
        out[4] = scale(a_st, 1.0 / (hu * hv));
        out[5] = scale(a_tt, 1.0 / (hv * hv));
        (out, 6)
    }

    #[allow(clippy::needless_range_loop)]
    fn d2_n<const N: usize>(&self, u: f64, v: f64) -> [Vec3; 6] {
        let u = self.bu.wrap(u, self.u_period);
        let v = self.bv.wrap(v, self.v_period);
        let (pu, pv) = (self.bu.degree, self.bv.degree);
        let (iu, iv) = (self.bu.span(u), self.bv.span(v));
        let nu = self.bu.ders::<N>(iu, u, 2);
        let nv = self.bv.ders::<N>(iv, v, 2);
        // Homogeneous sums for the orders (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), each row of
        // the net summed along v first for the three v orders.
        const ORD: [(usize, usize); 6] = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)];
        let mut a = [Vec3::ZERO; 6];
        let mut w = [0.0; 6];
        let stride = self.bv.n;
        for k in 0..=pu {
            let row = (iu - pu + k) * stride + iv - pv;
            let mut ra = [Vec3::ZERO; 3];
            let mut rw = [0.0; 3];
            for l in 0..=pv {
                let idx = row + l;
                let wt = self.weights.as_ref().map_or(1.0, |w| w[idx]);
                let pw = self.ctrl[idx] * wt;
                for j in 0..3 {
                    ra[j] += pw * nv[j][l];
                    rw[j] += wt * nv[j][l];
                }
            }
            for m in 0..6 {
                let (du, dv) = ORD[m];
                a[m] += ra[dv] * nu[du][k];
                w[m] += rw[dv] * nu[du][k];
            }
        }
        let s = a[0] / w[0];
        let su = (a[1] - s * w[1]) / w[0];
        let sv = (a[2] - s * w[2]) / w[0];
        let suu = (a[3] - su * (2.0 * w[1]) - s * w[3]) / w[0];
        let suv = (a[4] - su * w[2] - sv * w[1] - s * w[4]) / w[0];
        let svv = (a[5] - sv * (2.0 * w[2]) - s * w[5]) / w[0];
        [s, su, sv, suu, suv, svv]
    }

    /// Newton steps towards the surface point nearest `p` from `uv`, held inside the domain in a
    /// parameter that is not periodic, given the second derivatives at `uv` when the caller has
    /// them. Each step is the full Newton step of the squared distance where its Hessian is
    /// positive definite, which converges quadratically even off the surface, and the
    /// Gauss–Newton step elsewhere.
    ///
    /// Returns the parameters, their point's distance from `p`, and whether the search settled on
    /// a foot of `p`: a point from which no way along the surface inside the domain comes nearer.
    /// A search that stops short, where no fraction of its step comes nearer, has found a point
    /// and not a foot.
    fn search(&self, p: Point3, mut uv: Vec2, start: Option<[Vec3; 6]>) -> (Vec2, f64, bool) {
        let ((u0, u1), (v0, v1)) = self.domain();
        let clamp = |uv: Vec2| {
            Vec2::new(
                if self.u_period.is_some() {
                    uv.x
                } else {
                    uv.x.clamp(u0, u1)
                },
                if self.v_period.is_some() {
                    uv.y
                } else {
                    uv.y.clamp(v0, v1)
                },
            )
        };
        let inside = clamp(uv);
        let mut at = match start {
            Some(at) if inside == uv => at,
            _ => self.d2(inside.x, inside.y),
        };
        uv = inside;
        // At a pole one parameter moves nothing, and no step solved there says which of the lines
        // out of it leads to the foot: the line is chosen first, and the search starts down it.
        if let Some(free) = self.pole_at(uv, &at) {
            match self.off_pole(p, uv, free) {
                Some(off) => {
                    uv = off;
                    at = self.d2(uv.x, uv.y);
                }
                None => return (self.wrapped(uv), at[0].dist(p), true),
            }
        }
        let mut settled = false;
        for _ in 0..60 {
            let [q, su, sv, suu, suv, svv] = at;
            let r = q - p;
            let (g0, g1) = (su.dot(r), sv.dot(r));
            let (a, b, d) = (su.dot(su), su.dot(sv), sv.dot(sv));
            let (ha, hb, hd) = (a + suu.dot(r), b + suv.dot(r), d + svv.dot(r));
            // A parameter at an end of its range, with the nearer point beyond it, stays there: the
            // foot is on that boundary, and the search goes on along it alone. A step solved in
            // both and then cut back to the range is not the nearest point of the boundary.
            let held = |x: f64, lo: f64, hi: f64, g: f64, free: bool| {
                !free && ((x <= lo && g > 0.0) || (x >= hi && g < 0.0))
            };
            let hold_u = held(uv.x, u0, u1, g0, self.u_period.is_some());
            let hold_v = held(uv.y, v0, v1, g1, self.v_period.is_some());
            let one = |g: f64, h: f64, gauss: f64| {
                let by = if h > 0.0 { h } else { gauss };
                (by > 0.0).then(|| -g / by)
            };
            let step = match (hold_u, hold_v) {
                // A corner of the domain, the nearer points beyond both its sides.
                (true, true) => return (self.wrapped(uv), r.norm(), true),
                (true, false) => one(g1, hd, d).map(|dv| (0.0, dv)),
                (false, true) => one(g0, ha, a).map(|du| (du, 0.0)),
                (false, false) => (ha > 0.0 && ha * hd - hb * hb > 0.0)
                    .then(|| solve2(ha, hb, hb, hd, -g0, -g1))
                    .flatten()
                    .or_else(|| solve2(a, b, b, d, -g0, -g1))
                    .or_else(|| {
                        (self.pole_at(uv, &at).is_none())
                            .then(|| eigen_step(ha, hb, hd, g0, g1))
                            .flatten()
                    }),
            };
            let Some((du, dv)) = step else {
                break;
            };
            // A step this short moves the point by rounding alone: take it and stop, rather than
            // evaluate the surface again to confirm it. Where the two tangents lie along one line,
            // at a cusp, a step that moves the point by nothing still moves the parameters, and
            // the foot is found by going on: the distance there falls only with a power of the
            // step, so Newton nears it by a constant fraction each time.
            if (su * du + sv * dv).norm() <= 1e-6 * tol::LINEAR {
                let size = 1.0 + (u1 - u0) + (v1 - v0);
                uv = clamp(Vec2::new(uv.x + du, uv.y + dv));
                if du.abs() + dv.abs() <= 1e-12 * size || self.pole_at(uv, &at).is_some() {
                    return (self.wrapped(uv), self.eval(uv.x, uv.y).dist(p), true);
                }
                at = self.d2(uv.x, uv.y);
                settled = true;
                continue;
            }
            // A full step from far off can overshoot on a strongly curved or twisted patch; halve
            // it until it gets nearer, or stop where no fraction of it does.
            let mut f = 1.0;
            let next = loop {
                let next = clamp(Vec2::new(uv.x + f * du, uv.y + f * dv));
                let there = self.d2(next.x, next.y);
                if there[0].dist(p) <= r.norm() * (1.0 + 1e-12) + 1e-15 {
                    break Some((next, there));
                }
                f *= 0.5;
                if f < 1e-3 {
                    break None;
                }
            };
            let Some((next, there)) = next else {
                break;
            };
            let moved = (next.x - uv.x).abs() + (next.y - uv.y).abs();
            uv = next;
            at = there;
            if moved <= 1e-15 * (1.0 + uv.x.abs() + uv.y.abs()) {
                break;
            }
        }
        (self.wrapped(uv), at[0].dist(p), settled)
    }

    /// Where a search from the pole at `uv` starts: a point a way down the line out of the pole
    /// that heads for `p` most directly, at the nearest point of that line to `p`. `free` is the
    /// parameter that moves nothing at the pole (0 for u, 1 for v); the other is at an end of its
    /// range there. `None` where every line out of the pole leads away from `p`: the pole is the
    /// foot.
    fn off_pole(&self, p: Point3, uv: Vec2, free: usize) -> Option<Vec2> {
        let ((u0, u1), (v0, v1)) = self.domain();
        let ranges = [(u0, u1), (v0, v1)];
        let (lo, hi) = ranges[free];
        let (w0, w1) = ranges[1 - free];
        let pinned = [uv.x, uv.y][1 - free];
        // Into the domain from the end the pole is on; a pole inside a range has no such way.
        let inward = if pinned <= w0 {
            1.0
        } else if pinned >= w1 {
            -1.0
        } else {
            return Some(uv);
        };
        let at = |t: f64, w: f64| {
            if free == 0 {
                Vec2::new(t, w)
            } else {
                Vec2::new(w, t)
            }
        };
        // How fast the distance to `p` falls on the way in along the line at `t`.
        let fall = |t: f64| {
            let x = at(t, pinned);
            let [q, su, sv, ..] = self.d2(x.x, x.y);
            let line = if free == 0 { sv } else { su };
            let n = line.norm();
            if n > 0.0 {
                -inward * line.dot(q - p) / n
            } else {
                f64::NEG_INFINITY
            }
        };
        // The best of the lines a grid names, then the best between its neighbours: the fall is
        // one smooth wave round a pole, so its peak lies between the grid's best and one of them.
        const LINES: usize = 16;
        let grid = |k: usize| lo + (hi - lo) * k as f64 / LINES as f64;
        let best = (0..=LINES)
            .map(|k| (fall(grid(k)), k))
            .max_by(|x, y| x.0.total_cmp(&y.0))?;
        let (mut a, mut b) = (
            grid(best.1.saturating_sub(1)),
            grid((best.1 + 1).min(LINES)),
        );
        let golden = 0.5 * (5f64.sqrt() - 1.0);
        let (mut c, mut d) = (b - golden * (b - a), a + golden * (b - a));
        let (mut fc, mut fd) = (fall(c), fall(d));
        for _ in 0..40 {
            if fc > fd {
                (b, d, fd) = (d, c, fc);
                c = b - golden * (b - a);
                fc = fall(c);
            } else {
                (a, c, fc) = (c, d, fd);
                d = a + golden * (b - a);
                fd = fall(d);
            }
        }
        let t = 0.5 * (a + b);
        if fall(t) <= 0.0 {
            return None;
        }
        // The nearest point of that line to `p`: Newton in the one parameter, from the pole.
        let mut w = pinned;
        for _ in 0..40 {
            let x = at(t, w);
            let [q, su, sv, suu, _, svv] = self.d2(x.x, x.y);
            let (line, bend) = if free == 0 { (sv, svv) } else { (su, suu) };
            let r = q - p;
            let (g, gauss) = (line.dot(r), line.dot(line));
            let h = gauss + bend.dot(r);
            let by = if h > 0.0 { h } else { gauss };
            if by <= 0.0 {
                break;
            }
            let mut step = -g / by;
            let mut next = (w + step).clamp(w0, w1);
            // No further than gets nearer, as the search over both parameters does.
            while self.eval(at(t, next).x, at(t, next).y).dist(p) > r.norm() * (1.0 + 1e-12)
                && step.abs() > 1e-15 * (w1 - w0)
            {
                step *= 0.5;
                next = (w + step).clamp(w0, w1);
            }
            let moved = (next - w).abs();
            w = next;
            if moved <= 1e-15 * (1.0 + w.abs()) {
                break;
            }
        }
        // Still on the pole: the line fell at first and no step down it came nearer.
        if w == pinned {
            return None;
        }
        Some(at(t, w))
    }

    /// [`pole_of`] at `uv`, where the line of the parameter that moves nothing is one point: a
    /// vanishing derivative alone, along an end line whose control points coincide, is a line on
    /// which the surface stands still, not a pole, and no line out of it need be chosen.
    fn pole_at(&self, uv: Vec2, at: &[Vec3; 6]) -> Option<usize> {
        let free = pole_of(at)?;
        let ((u0, u1), (v0, v1)) = self.domain();
        let (lo, hi) = if free == 0 { (u0, u1) } else { (v0, v1) };
        let point = (0..=8).all(|k| {
            let x = lo + (hi - lo) * k as f64 / 8.0;
            let q = if free == 0 {
                self.eval(x, uv.y)
            } else {
                self.eval(uv.x, x)
            };
            q.dist(at[0]) <= 2.0 * tol::LINEAR
        });
        point.then_some(free)
    }

    fn wrapped(&self, uv: Vec2) -> Vec2 {
        Vec2::new(
            self.bu.wrap(uv.x, self.u_period),
            self.bv.wrap(uv.y, self.v_period),
        )
    }

    /// The parameters of the surface point nearest `p`, periodic parameters in the knot domain.
    pub fn project(&self, p: Point3) -> Vec2 {
        // The four nearest seeds at distinct points, nearest first: a pole's row of seeds is one
        // point, and four of them would be one search.
        let near = self.nearest_seeds(p);
        let mut best = (f64::INFINITY, Vec2::ZERO);
        for &(d2, seed) in &near {
            if d2.is_infinite() {
                break;
            }
            let (uv, d, settled) = self.search(p, seed, None);
            if d < best.0 {
                best = (d, uv);
            }
            // A foot this near the surface is the nearest, as `project_near` takes it: the seed
            // grid resolves the surface, so no other sheet comes that near it. A search that
            // stopped short of a foot says nothing of the kind.
            if settled && best.0 <= self.near_radius() {
                break;
            }
        }
        best.1
    }

    /// How near the surface a point must lie for [`NurbsSurface::project_near`] to take the point
    /// Newton finds from its hint as the nearest.
    pub fn near_radius(&self) -> f64 {
        tol::LINEAR.max(NEAR_SHARE * self.seed_grid().spacing)
    }

    /// The parameters of the surface point nearest `p`, found from `hint` when the surface passes
    /// within tolerance there, globally otherwise.
    pub fn project_near(&self, p: Point3, hint: Vec2) -> Vec2 {
        // A hint on a pole says nothing of the free parameter: Newton from it finds the point
        // nearest along that one line out of the pole, which can be a boundary's, not the patch's.
        let at = self.d2(hint.x, hint.y);
        if self.pole_at(hint, &at).is_some() {
            return self.project(p);
        }
        let (local, d, settled) = self.search(p, hint, Some(at));
        // Global projection trusts the seed grid to resolve the surface: no two of its sheets come
        // within a fraction of the seed spacing of one point. A foot that near the surface is
        // therefore the nearest.
        if settled && d <= self.near_radius() {
            return local;
        }
        let global = self.project(p);
        if self.eval(global.x, global.y).dist(p) < d {
            global
        } else {
            local
        }
    }
}

/// The Newton step `−H⁻¹g` for the symmetric `H = [[a, b], [b, d]]`, solved in the frame of its
/// eigenvectors, where both eigenvalues are positive: at a cusp, where the surface's two tangents
/// lie along one line, `H` is too near singular for a determinant to tell, yet the step along the
/// line the surface barely moves on is still the one that nears the foot.
fn eigen_step(a: f64, b: f64, d: f64, g0: f64, g1: f64) -> Option<(f64, f64)> {
    let mean = 0.5 * (a + d);
    let half = (0.5 * (a - d)).hypot(b);
    let (big, small) = (mean + half, mean - half);
    if small.is_nan() || small <= 0.0 {
        return None;
    }
    // The eigenvector of the larger eigenvalue, from whichever row of `H − big·I` is longer.
    let e = if (a - big).abs() >= (d - big).abs() {
        Vec2::new(b, big - a)
    } else {
        Vec2::new(big - d, b)
    };
    let e = if e.norm() > 0.0 {
        e * (1.0 / e.norm())
    } else {
        Vec2::new(1.0, 0.0)
    };
    let f = Vec2::new(-e.y, e.x);
    let g = Vec2::new(g0, g1);
    let step = e * (-g.dot(e) / big) + f * (-g.dot(f) / small);
    Some((step.x, step.y))
}

/// The parameter that moves nothing at a point with these derivatives, 0 for u and 1 for v: a
/// pole, where a whole row of the surface is one point.
fn pole_of(at: &[Vec3; 6]) -> Option<usize> {
    let (a, d) = (at[1].norm(), at[2].norm());
    if a <= 1e-9 * d {
        Some(0)
    } else if d <= 1e-9 * a {
        Some(1)
    } else {
        None
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use std::f64::consts::{FRAC_1_SQRT_2, TAU};

    /// The unit circle as a rational quadratic in four quarter arcs.
    fn circle() -> NurbsCurve {
        let w = FRAC_1_SQRT_2;
        let pts = [
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (-1.0, 1.0),
            (-1.0, 0.0),
            (-1.0, -1.0),
            (0.0, -1.0),
            (1.0, -1.0),
            (1.0, 0.0),
        ];
        let ctrl = pts.iter().map(|&(x, y)| Vec3::new(x, y, 0.0)).collect();
        let weights = (0..9).map(|i| if i % 2 == 0 { 1.0 } else { w }).collect();
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0];
        NurbsCurve::new(2, knots, ctrl, Some(weights)).unwrap()
    }

    #[test]
    fn rational_circle_is_exact_and_periodic() {
        let c = circle();
        assert_eq!(c.period(), Some(4.0));
        for i in 0..=40 {
            let t = -2.0 + 8.0 * i as f64 / 40.0;
            let (p, d1, d2) = c.d2(t);
            assert!((p.norm() - 1.0).abs() < 1e-15, "{t}");
            assert!(p.dot(d1).abs() < 1e-14);
            let h = 1e-5;
            let fd = (c.eval(t + h) - c.eval(t - h)) / (2.0 * h);
            let fd2 = (c.eval(t + h) - p * 2.0 + c.eval(t - h)) / (h * h);
            if (t - t.round()).abs() > 2.0 * h {
                assert!(d1.dist(fd) < 1e-8, "{t}");
                assert!(d2.dist(fd2) < 1e-4, "{t}");
            }
        }
    }

    #[test]
    fn surface_derivatives_and_projection() {
        // A quarter of a cylinder of radius 1, height 2, as a rational quadratic by linear patch.
        let w = FRAC_1_SQRT_2;
        let arc = [(1.0, 0.0, 1.0), (1.0, 1.0, w), (0.0, 1.0, 1.0)];
        let mut ctrl = Vec::new();
        let mut weights = Vec::new();
        for &(x, y, wt) in &arc {
            for z in [0.0, 2.0] {
                ctrl.push(Vec3::new(x, y, z));
                weights.push(wt);
            }
        }
        let s = NurbsSurface::new(
            (2, 1),
            (vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], vec![0.0, 0.0, 1.0, 1.0]),
            (3, 2),
            ctrl,
            Some(weights),
        )
        .unwrap();
        assert_eq!(s.u_period(), None);
        for &(u, v) in &[(0.3, 0.4), (0.9, 0.1), (0.5, 0.5)] {
            let (p, su, sv) = s.d1(u, v);
            assert!((p.x.hypot(p.y) - 1.0).abs() < 1e-15);
            let h = 1e-6;
            let fu = (s.eval(u + h, v) - s.eval(u - h, v)) / (2.0 * h);
            let fv = (s.eval(u, v + h) - s.eval(u, v - h)) / (2.0 * h);
            assert!(su.dist(fu) < 1e-8 && sv.dist(fv) < 1e-8);
            let uv = s.project(p + su.cross(sv).unit().unwrap() * 0.5);
            assert!(
                (uv.x - u).abs() < 1e-12 && (uv.y - v).abs() < 1e-12,
                "{uv:?}"
            );
            let near = s.project_near(p, Vec2::new(0.5, 0.5));
            assert!(s.eval(near.x, near.y).dist(p) < 1e-12);
        }
        let angle = TAU / 4.0;
        assert!(
            s.eval(1.0, 0.0)
                .dist(Vec3::new(angle.cos(), angle.sin(), 0.0))
                < 1e-15
        );
    }

    /// A rational bicubic with a twist and uneven weights, two spans along `u`.
    fn warped() -> NurbsSurface {
        let (nu, nv) = (5, 4);
        let mut ctrl = Vec::new();
        let mut weights = Vec::new();
        for i in 0..nu {
            for j in 0..nv {
                let (x, y) = (i as f64, j as f64);
                let a = 0.3 * x * y;
                ctrl.push(Vec3::new(
                    x * a.cos(),
                    y + x * a.sin(),
                    (x - 2.0) * (y - 1.5) * 0.4,
                ));
                weights.push(1.0 + 0.1 * ((i + 2 * j) % 3) as f64);
            }
        }
        NurbsSurface::new(
            (3, 3),
            (
                vec![0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0],
                vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            ),
            (nu, nv),
            ctrl,
            Some(weights),
        )
        .unwrap()
    }

    /// The least distance from `p` to a grid of the patch's points, its boundary among them: no
    /// point of the patch is much nearer than this, and its nearest point is no farther.
    fn nearest_sampled(s: &NurbsSurface, p: Point3, n: usize) -> f64 {
        let mut best = f64::INFINITY;
        for i in 0..=n {
            for j in 0..=n {
                let q = s.eval(i as f64 / n as f64, j as f64 / n as f64);
                best = best.min(q.dist(p));
            }
        }
        best
    }

    #[test]
    fn a_point_beyond_an_edge_projects_onto_the_nearest_point_of_that_edge() {
        let s = warped();
        // Points off each of the four edges, a way along and a way out, and off a corner.
        let mut tried = 0;
        for k in 1..12 {
            let t = k as f64 / 12.0;
            for (on, inside) in [
                (Vec2::new(t, 1.0), Vec2::new(t, 0.9)),
                (Vec2::new(t, 0.0), Vec2::new(t, 0.1)),
                (Vec2::new(1.0, t), Vec2::new(0.9, t)),
                (Vec2::new(0.0, t), Vec2::new(0.1, t)),
            ] {
                let (edge, back) = (s.eval(on.x, on.y), s.eval(inside.x, inside.y));
                let out = (edge - back).unit().unwrap();
                let [_, su, sv, ..] = s.d2(on.x, on.y);
                let lift = su.cross(sv).unit().unwrap() * 0.05;
                for reach in [0.02, 0.2, 0.6] {
                    let p = edge + out * reach + lift;
                    let uv = s.project(p);
                    let got = s.eval(uv.x, uv.y).dist(p);
                    let sampled = nearest_sampled(&s, p, 600);
                    assert!(
                        got <= sampled + 1e-9,
                        "off {on:?} by {reach}: {got} against {sampled} at {uv:?}"
                    );
                    tried += 1;
                }
            }
        }
        assert_eq!(tried, 132);
    }

    /// The point of a Bézier patch at `(s, t)` in its own unit square, by de Casteljau.
    fn bezier_point(b: &Bezier, s: f64, t: f64) -> Point3 {
        let fold = |mut level: Vec<[f64; 4]>, x: f64| {
            while level.len() > 1 {
                level = level
                    .windows(2)
                    .map(|w| std::array::from_fn(|c| (1.0 - x) * w[0][c] + x * w[1][c]))
                    .collect();
            }
            level[0]
        };
        let rows: Vec<[f64; 4]> = (0..b.rows)
            .map(|i| fold(b.net[i * b.cols..(i + 1) * b.cols].to_vec(), t))
            .collect();
        let q = fold(rows, s);
        Vec3::new(q[0] / q[3], q[1] / q[3], q[2] / q[3])
    }

    #[test]
    fn bezier_patches_are_the_surface_and_their_halves_are_too() {
        for s in uneven().iter().chain([&warped()]) {
            let mut pieces: Vec<Bezier> = s.beziers().to_vec();
            assert_eq!(pieces.len(), s.patches().len());
            // Each patch, and the quarters two cuts leave of the first few.
            let quarters: Vec<Bezier> = pieces
                .iter()
                .take(6)
                .flat_map(|b| b.halves(true))
                .flat_map(|b| b.halves(false))
                .collect();
            pieces.extend(quarters);
            for b in &pieces {
                let bb = b.bbox();
                for (x, y) in [(0.0, 0.0), (1.0, 1.0), (0.3, 0.8), (0.5, 0.5), (0.9, 0.1)] {
                    let (u, v) = (
                        b.lo.x + (b.hi.x - b.lo.x) * x,
                        b.lo.y + (b.hi.y - b.lo.y) * y,
                    );
                    // A hair inside the far sides, where the patch's own spans end.
                    let on = s.eval(u.min(b.hi.x - 1e-12), v.min(b.hi.y - 1e-12));
                    let got = bezier_point(b, x, y);
                    assert!(got.dist(on) < 1e-9, "{got:?} against {on:?}");
                    let inside = (0..3).all(|k| {
                        on.to_array()[k] >= bb.min.to_array()[k] - 1e-9
                            && on.to_array()[k] <= bb.max.to_array()[k] + 1e-9
                    });
                    assert!(inside, "{on:?} outside {bb:?}");
                }
            }
        }
    }

    #[test]
    fn an_iso_curve_is_the_surface_along_one_parameter() {
        let mut surfaces = vec![warped()];
        surfaces.extend(uneven());
        for s in surfaces {
            let ((u0, u1), (v0, v1)) = s.domain();
            let size = s.eval(u0, v0).norm() + 1.0;
            for along_v in [true, false] {
                for a in [0.0, 0.3, 0.5, 0.75, 1.0] {
                    let at = if along_v {
                        u0 + (u1 - u0) * a
                    } else {
                        v0 + (v1 - v0) * a
                    };
                    let c = s.iso_curve(along_v, at);
                    let (t0, t1) = c.domain();
                    assert_eq!((t0, t1), if along_v { (v0, v1) } else { (u0, u1) });
                    for b in [0.0, 0.2, 0.45, 0.5, 0.9, 1.0] {
                        let t = t0 + (t1 - t0) * b;
                        let want = if along_v {
                            s.eval(at, t)
                        } else {
                            s.eval(t, at)
                        };
                        let off = c.eval(t).dist(want);
                        assert!(off < 1e-12 * size, "{off} at {at}, {t}, along v {along_v}");
                    }
                }
            }
        }
    }

    #[test]
    fn iso_lines_agree_with_the_surface() {
        let s = warped();
        for &(at, t) in &[(0.2, 0.3), (0.55, 0.8), (0.4, 0.4), (1.0, 0.0)] {
            for along_v in [true, false] {
                let (u, v) = if along_v { (at, t) } else { (t, at) };
                let (p, su, sv) = s.iso(along_v, at).d1(t);
                let (q, qu, qv) = s.d1(u, v);
                assert!(p.dist(q) < 1e-14 && su.dist(qu) < 1e-13 && sv.dist(qv) < 1e-13);
            }
        }
        // The span polynomials evaluate up to degree 6, across a double knot and on both sides of
        // each span's middle, where the series switches ends; past it the basis recurrence does.
        for s in uneven() {
            let ((u0, u1), (v0, v1)) = s.domain();
            let one = s.iso(true, u0).d1(v0);
            let size = one.0.norm() + one.1.norm() + one.2.norm();
            for (a, b) in [
                (0.0, 0.0),
                (0.15, 0.9),
                (0.55, 0.45),
                (0.5, 0.5),
                (1.0, 1.0),
            ] {
                for u in [u0 + (u1 - u0) * a, 1.1] {
                    let v = v0 + (v1 - v0) * b;
                    let (q, qu, qv) = s.d1(u, v);
                    for (p, su, sv) in [s.iso(true, u).d1(v), s.iso(false, v).d1(u)] {
                        let off = p.dist(q).max(su.dist(qu)).max(sv.dist(qv));
                        assert!(off < 1e-12 * size, "{off} at ({u}, {v})");
                    }
                }
            }
        }
    }

    /// A degree-(3, 6) polynomial surface with a double interior knot along `u`, a degree-(6, 6)
    /// rational one, the highest the span polynomials take, and a degree-12 rational one: spans of
    /// unequal width, a knot where the surface is only C¹, a high degree.
    fn uneven() -> [NurbsSurface; 3] {
        let net = |nu: usize, nv: usize, rational: bool| {
            let mut ctrl = Vec::new();
            let mut weights = Vec::new();
            for i in 0..nu {
                for j in 0..nv {
                    let (x, y) = (i as f64, j as f64);
                    ctrl.push(Vec3::new(x + 0.2 * y, y - 0.1 * x * x, (x * 0.7).sin() * y));
                    weights.push(0.6 + 0.15 * ((3 * i + j) % 5) as f64);
                }
            }
            (ctrl, rational.then_some(weights))
        };
        let (c1, _) = net(7, 8, false);
        let a = NurbsSurface::new(
            (3, 6),
            (
                vec![0.0, 0.0, 0.0, 0.0, 0.3, 1.1, 1.1, 2.0, 2.0, 2.0, 2.0],
                vec![
                    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.5, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
                    2.0,
                ],
            ),
            (7, 8),
            c1,
            None,
        )
        .unwrap();
        let (c3, w3) = net(9, 8, true);
        let c = NurbsSurface::new(
            (6, 6),
            (
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 1.9, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
                ],
                vec![1.0; 7]
                    .into_iter()
                    .chain([2.5])
                    .chain(vec![4.0; 7])
                    .collect(),
            ),
            (9, 8),
            c3,
            w3,
        )
        .unwrap();
        let (c2, w2) = net(13, 4, true);
        let mut ku = vec![0.0; 13];
        ku.extend(vec![1.0; 13]);
        let b = NurbsSurface::new(
            (12, 2),
            (ku, vec![0.0, 0.0, 0.0, 0.25, 1.0, 1.0, 1.0]),
            (13, 4),
            c2,
            w2,
        )
        .unwrap();
        [a, c, b]
    }

    #[test]
    fn span_polynomials_agree_with_the_basis_recurrence() {
        // Each span's polynomial must give what the Cox–de Boor recurrence does: the point and
        // every derivative up to the second, at span ends, inside spans and on a double knot.
        let mut surfaces = uneven().to_vec();
        surfaces.push(warped());
        // Interior knots of full multiplicity: C⁰ joints, as a corner patch of three arcs has.
        let mut ctrl = Vec::new();
        let mut weights = Vec::new();
        for i in 0..7 {
            for j in 0..3 {
                let a = 0.9 * i as f64;
                ctrl.push(Vec3::new(
                    a.cos() * (2.0 - j as f64),
                    a.sin() * (2.0 - j as f64),
                    j as f64,
                ));
                weights.push(if i % 2 == 1 { FRAC_1_SQRT_2 } else { 1.0 });
            }
        }
        surfaces.push(
            NurbsSurface::new(
                (2, 2),
                (
                    vec![0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0],
                    vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                ),
                (7, 3),
                ctrl,
                Some(weights),
            )
            .unwrap(),
        );
        for s in &surfaces {
            let ((u0, u1), (v0, v1)) = s.domain();
            let mut knots: Vec<f64> = s.bu.knots.clone();
            knots.extend((0..=37).map(|k| u0 + (u1 - u0) * k as f64 / 37.0));
            let mut vk: Vec<f64> = s.bv.knots.clone();
            vk.extend((0..=29).map(|k| v0 + (v1 - v0) * k as f64 / 29.0));
            for &u in knots.iter().filter(|&&u| u >= u0 && u <= u1) {
                for &v in vk.iter().filter(|&&v| v >= v0 && v <= v1) {
                    let got = s.d2(u, v);
                    let want = s.d2_n::<M>(u, v);
                    // Rounding in the rational quotient grows with each derivative taken.
                    let tol = [1e-12, 1e-11, 1e-11, 1e-10, 1e-10, 1e-10];
                    for (k, (g, w)) in got.iter().zip(&want).enumerate() {
                        let size = 1.0 + w.norm();
                        assert!(
                            g.dist(*w) <= tol[k] * size,
                            "order {k} at ({u}, {v}): {g:?} against {w:?}"
                        );
                    }
                    let (p, su, sv) = s.d1(u, v);
                    let (q, qu, qv) = s.d1_n::<M>(u, v);
                    assert!(p.dist(q) <= 1e-12 * (1.0 + q.norm()));
                    assert!(su.dist(qu) <= 1e-11 * (1.0 + qu.norm()));
                    assert!(sv.dist(qv) <= 1e-11 * (1.0 + qv.norm()));
                    assert!(s.eval(u, v).dist(q) <= 1e-12 * (1.0 + q.norm()));
                }
            }
        }
    }

    #[test]
    fn second_derivatives_match_differences() {
        let s = warped();
        let h = 1e-5;
        for &(u, v) in &[(0.2, 0.3), (0.55, 0.8), (0.9, 0.1)] {
            let [_, su, sv, suu, suv, svv] = s.d2(u, v);
            let (_, su0, sv0) = s.d1(u, v);
            assert!(su.dist(su0) < 1e-12 && sv.dist(sv0) < 1e-12);
            let du = |u: f64, v: f64| s.d1(u, v).1;
            let dv = |u: f64, v: f64| s.d1(u, v).2;
            assert!(suu.dist((du(u + h, v) - du(u - h, v)) / (2.0 * h)) < 1e-6);
            assert!(suv.dist((du(u, v + h) - du(u, v - h)) / (2.0 * h)) < 1e-6);
            assert!(svv.dist((dv(u, v + h) - dv(u, v - h)) / (2.0 * h)) < 1e-6);
        }
    }

    #[test]
    fn patch_boxes_hold_the_surface_over_their_parameters() {
        let s = warped();
        let patches = s.patches();
        assert!(patches.len() > 2);
        for i in 0..=60 {
            for j in 0..=60 {
                let uv = Vec2::new(i as f64 / 60.0, j as f64 / 60.0);
                let p = s.eval(uv.x, uv.y);
                let holding = patches.iter().filter(|q| {
                    uv.x >= q.lo.x && uv.x <= q.hi.x && uv.y >= q.lo.y && uv.y <= q.hi.y
                });
                let mut any = false;
                for q in holding {
                    any = true;
                    let pad = 1e-12;
                    assert!(
                        (0..3).all(|k| p[k] >= q.bbox.min[k] - pad && p[k] <= q.bbox.max[k] + pad),
                        "{uv:?}"
                    );
                }
                assert!(any, "{uv:?}");
            }
        }
    }

    #[test]
    fn projection_from_off_the_surface_lands_on_the_nearest_point() {
        let s = warped();
        for &(u, v, off) in &[(0.3, 0.6, 0.05), (0.7, 0.2, -0.08), (0.45, 0.45, 0.02)] {
            let (p, su, sv) = s.d1(u, v);
            let q = p + su.cross(sv).unit().unwrap() * off;
            let uv = s.project(q);
            assert!(
                (uv.x - u).abs() < 1e-10 && (uv.y - v).abs() < 1e-10,
                "{uv:?}"
            );
            let near = s.project_near(q, Vec2::new(u + 0.02, v - 0.02));
            assert!(near.dist(uv) < 1e-10, "{near:?}");
        }
    }

    /// The octant of the sphere of radius `r` about the origin, as the rational biquadratic of
    /// Piegl and Tiller: `u` round the equator from +x to +y, `v` up to the pole at +z, `v = 1`.
    pub(crate) fn octant(r: f64) -> crate::geom::Surface {
        let w = FRAC_1_SQRT_2;
        let equator = [(r, 0.0, 1.0), (r, r, w), (0.0, r, 1.0)];
        let mut ctrl = Vec::new();
        let mut weights = Vec::new();
        for &(x, y, wu) in &equator {
            for (j, wv) in [1.0, w, 1.0].into_iter().enumerate() {
                ctrl.push(match j {
                    0 => Vec3::new(x, y, 0.0),
                    1 => Vec3::new(x, y, r),
                    _ => Vec3::new(0.0, 0.0, r),
                });
                weights.push(wu * wv);
            }
        }
        let k = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let n = NurbsSurface::new((2, 2), (k.clone(), k), (3, 3), ctrl, Some(weights)).unwrap();
        crate::geom::Surface::Nurbs(Box::new(n))
    }

    /// The unit cylinder round z, from z = 0 to 3: [`circle`] along `u`, closed, by a line.
    fn tube() -> NurbsSurface {
        let c = circle();
        let mut ctrl = Vec::new();
        let mut weights = Vec::new();
        for (p, w) in c.control_points().iter().zip(c.weights().unwrap()) {
            ctrl.extend([*p, *p + Vec3::Z * 3.0]);
            weights.extend([*w, *w]);
        }
        let knots = (c.knots().to_vec(), vec![0.0, 0.0, 1.0, 1.0]);
        NurbsSurface::new((2, 1), knots, (9, 2), ctrl, Some(weights)).unwrap()
    }

    fn bits(p: Point3) -> [u64; 3] {
        p.to_array().map(f64::to_bits)
    }

    #[test]
    fn a_point_alone_is_the_point_its_derivatives_come_with() {
        let crate::geom::Surface::Nurbs(pole) = octant(2.0) else {
            unreachable!()
        };
        let surfaces = uneven().into_iter().chain([warped(), tube(), *pole]);
        for s in surfaces {
            let ((u0, u1), (v0, v1)) = s.domain();
            // Past the ends too, where a closed parameter wraps and an open one runs on.
            let along = |a: f64, b: f64| -> Vec<f64> {
                (0..=23)
                    .map(|i| a - 0.2 * (b - a) + 1.4 * (b - a) * i as f64 / 23.0)
                    .chain([a, b])
                    .collect()
            };
            let (us, vs) = (along(u0, u1), along(v0, v1));
            let grid = s.grid(&us, &vs);
            for (i, &u) in us.iter().enumerate() {
                for (j, &v) in vs.iter().enumerate() {
                    let want = bits(s.d1(u, v).0);
                    assert_eq!(bits(s.eval(u, v)), want);
                    assert_eq!(bits(grid[i * vs.len() + j]), want);
                }
            }
        }
        let c = circle();
        for i in 0..=200 {
            let t = -2.0 + 8.0 * i as f64 / 200.0;
            let (p, d1, _) = c.d2(t);
            assert_eq!(bits(c.eval(t)), bits(p));
            assert_eq!(bits(c.d1(t).0), bits(p));
            assert_eq!(bits(c.d1(t).1), bits(d1));
        }
    }

    #[test]
    fn the_seeds_are_the_surface_at_an_even_spacing() {
        // The spacing is from the longest lines of a first, coarser grid along each parameter,
        // measured here from the surface point by point.
        for s in uneven().into_iter().chain([warped(), tube()]) {
            let grid = s.seed_grid();
            assert_eq!(grid.seeds.len() % grid.row, 0);
            for &(uv, q) in &grid.seeds {
                assert_eq!(bits(q), bits(s.d1(uv.x, uv.y).0));
            }
            let coarse = |b: &Basis| -> Vec<f64> {
                let spans = b.spans();
                let per = (24 / (spans.len() - 1)).clamp(2, 6);
                let mut out = Vec::new();
                for w in spans.windows(2) {
                    out.extend((0..per).map(|k| w[0] + (w[1] - w[0]) * k as f64 / per as f64));
                }
                out.push(spans[spans.len() - 1]);
                out
            };
            let (us, vs) = (coarse(&s.bu), coarse(&s.bv));
            let line = |pts: Vec<Point3>| pts.windows(2).map(|w| w[0].dist(w[1])).sum::<f64>();
            let lu = vs
                .iter()
                .map(|&v| line(us.iter().map(|&u| s.d1(u, v).0).collect()))
                .fold(0.0, f64::max);
            let lv = us
                .iter()
                .map(|&u| line(vs.iter().map(|&v| s.d1(u, v).0).collect()))
                .fold(0.0, f64::max);
            let h = (lu.min(lv) / 4.0).max((lu * lv / 2048.0).sqrt());
            assert_eq!(grid.spacing.to_bits(), h.to_bits());
        }
    }

    /// A sheet folded back over itself: out from near z = 0 to the fold at x = 2.25 and back near
    /// z = 1, the two sheets about 1 apart, 2 long in y. `u` runs round the fold and `v` along y.
    fn hairpin() -> NurbsSurface {
        let rows = [(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (0.0, 1.0)];
        let ctrl = rows
            .iter()
            .flat_map(|&(x, z)| [Vec3::new(x, 0.0, z), Vec3::new(x, 2.0, z)])
            .collect();
        let knots = (
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
        );
        NurbsSurface::new((3, 1), knots, (4, 2), ctrl, None).unwrap()
    }

    #[test]
    fn a_hint_is_kept_only_for_a_foot_beside_the_surface() {
        let s = hairpin();
        let near = s.near_radius();
        assert!(near < 0.05, "{near}");
        let dist = |uv: Vec2, p: Point3| s.eval(uv.x, uv.y).dist(p);
        let (low, high) = (Vec2::new(0.1, 0.5), Vec2::new(0.9, 0.5));
        let (below, above) = (s.eval(low.x, low.y), s.eval(high.x, high.y));
        assert!(below.z < 0.1 && above.z > 0.9 && (below.x - above.x).abs() < 1e-9);
        // The foot Newton finds from a hint alone, and whether it settled there.
        let local = |p: Point3, hint: Vec2| s.search(p, hint, None);

        // Beside the surface, from a hint beside the foot: the hint's foot is kept, and it is
        // the global search's.
        let p = above + Vec3::new(0.02, 0.1, -0.5 * near);
        let (kept, d, settled) = local(p, high);
        assert!(settled && d <= near);
        let got = s.project_near(p, high);
        assert_eq!(
            (got.x.to_bits(), got.y.to_bits()),
            (kept.x.to_bits(), kept.y.to_bits())
        );
        assert!((dist(got, p) - dist(s.project(p), p)).abs() < 1e-12);

        // The same point from a hint on the other sheet: Newton settles on that sheet's foot, a
        // whole gap away, which is no foot beside the surface, and the global search answers.
        let (wrong, d, settled) = local(p, low);
        assert!(settled && d > 0.9 && s.eval(wrong.x, wrong.y).z < 0.1);
        let got = s.project_near(p, low);
        assert!(dist(got, p) <= near && s.eval(got.x, got.y).z > 0.9);

        // Between the sheets the nearest point jumps from one to the other as the point crosses
        // the middle: whichever sheet the hint is on, the answer is the nearer sheet's.
        for (z, hint) in [(0.4, high), (0.6, low), (0.4, low), (0.6, high)] {
            let p = Vec3::new(below.x, 1.0, z);
            let got = s.project_near(p, hint);
            assert_eq!(s.eval(got.x, got.y).z > 0.5, z > 0.5);
            let sampled = nearest_sampled(&s, p, 400);
            assert!(dist(got, p) <= sampled + 1e-9 && dist(got, p) > sampled - 1e-3);
        }

        // Beyond the sheet's edge a step from the hint leaves the domain: the search is held to
        // the edge, and the foot there, beside the surface or not, is the nearest point.
        for out in [0.5 * near, 0.3] {
            let p = s.eval(low.x, 1.0) + Vec3::Y * out;
            let got = s.project_near(p, low);
            assert!(
                got.y == 1.0 && (dist(got, p) - out).abs() < 1e-12,
                "{got:?}"
            );
            assert!((dist(got, p) - dist(s.project(p), p)).abs() < 1e-12);
        }

        // A hint outside the domain is brought into it before the search starts.
        let got = s.project_near(below, Vec2::new(-3.0, 7.0));
        assert!(dist(got, below) < 1e-12);

        // A search that stops short of a foot is not kept, however near it stopped: at the fold
        // a point on the mid-plane inside the bend has the fold's line as a farthest point, not
        // a nearest, and the answer is the nearer of what the hint and the global search find.
        let fold = s.eval(0.5, 0.5);
        let p = Vec3::new(fold.x - 0.3, 1.0, 0.5);
        let got = s.project_near(p, Vec2::new(0.5, 0.5));
        assert!(
            dist(got, p) <= nearest_sampled(&s, p, 800) + 1e-9,
            "{}",
            dist(got, p)
        );
        assert!(dist(got, p) <= dist(s.project(p), p));
    }

    // A swept surface is its curve carried by its motion, within the fit of a screw's columns.
    #[test]
    fn a_swept_surface_keeps_its_curve_and_motion_until_its_shape_may_change() {
        use crate::{geom::Circle, geom::Curve, loft::rail::Rail, math::Frame, sweep};
        let fit = 2.0 * sweep::SCREW_FIT;
        let follows = |s: &NurbsSurface| {
            let swept = s.swept().expect("a swept surface");
            let ((u0, u1), (v0, v1)) = s.domain();
            (0..=20).all(|i| {
                (0..=8).all(|j| {
                    let (u, v) = (
                        u0 + (u1 - u0) * i as f64 / 20.0,
                        v0 + (v1 - v0) * j as f64 / 8.0,
                    );
                    let carried = swept.motion.pose(u).point(swept.generator.eval(v));
                    carried.dist(s.eval(u, v)) <= fit
                })
            })
        };
        let arc = Curve::Circle(Circle {
            frame: Frame::new(
                Point3::new(3.0, 1.0, 0.5),
                Vec3::new(0.0, 1.0, 0.2),
                Vec3::X,
            )
            .unwrap(),
            radius: 1.5,
        });
        let rail = Rail::of_curve(&arc, 0.2, 2.9).unwrap();
        for (rise, sign) in [(0.4, 1.0), (-0.7, -1.0), (0.0, 1.0)] {
            let screw = sweep::Helix {
                origin: Point3::new(0.5, -0.5, 1.0),
                axis: Vec3::new(0.6, 0.0, 0.8),
                rise,
                sign,
            };
            let s = sweep::screw_sheet(&rail, &screw, 3.0).unwrap();
            assert!(follows(&s));
            // A copy keeps them, and a rigid move moves them with the surface.
            assert!(follows(&s.clone()));
            let t = Transform::rotation(Point3::new(1.0, 2.0, 3.0), Vec3::new(0.0, 0.6, 0.8), 1.1)
                .then(&Transform::translation(Vec3::new(4.0, -2.0, 1.0)));
            let moved = s.transformed(&t);
            assert!(follows(&moved));
            let axis = moved.swept().unwrap().motion.axis;
            assert!(axis.dist(t.vector(screw.axis)) <= 1e-15);
            // A surface made from parts or control points may be any shape.
            let same = s.with_control_points(s.control_points().to_vec()).unwrap();
            assert!(same == s && same.swept().is_none());
            let (ku, kv) = s.knots();
            let rebuilt = NurbsSurface::new(
                s.degrees(),
                (ku.to_vec(), kv.to_vec()),
                s.counts(),
                s.control_points().to_vec(),
                s.weights().map(|w| w.to_vec()),
            )
            .unwrap();
            assert!(rebuilt == s && rebuilt.swept().is_none());
        }
    }

    #[test]
    fn a_copy_shares_what_its_surface_keeps_and_a_changed_surface_starts_afresh() {
        let s = warped();
        let copy = s.clone();
        // Seeds built through the copy are the original's too.
        copy.project(Vec3::new(1.0, 1.0, 3.0));
        assert!(s.kept.seeds.get().is_some());
        assert!(std::ptr::eq(s.seeds().as_ptr(), copy.seeds().as_ptr()));
        // Other control points: nothing is carried, and projection is onto the new surface.
        let mut ctrl = s.control_points().to_vec();
        ctrl[9] += Vec3::Z * 0.5;
        let bent = s.with_control_points(ctrl).unwrap();
        assert!(bent.kept.seeds.get().is_none());
        let q = bent.eval(0.5, 0.4);
        let on = |x: &NurbsSurface| {
            let uv = x.project(q);
            x.eval(uv.x, uv.y).dist(q)
        };
        assert!(
            on(&bent) < 1e-12 && on(&s) > 1e-3,
            "{} {}",
            on(&bent),
            on(&s)
        );
        assert!(
            bent.seeds()
                .iter()
                .all(|&(uv, p)| p == bent.eval(uv.x, uv.y))
        );
        // Moved: its seeds are the original's moved with it, on the moved surface, and what it
        // builds from them is its own.
        let turn = Transform::rotation(Vec3::new(0.3, 1.0, -2.0), Vec3::new(2.0, -1.0, 2.0), 1.1);
        let moved = s.transformed(&turn);
        assert!(!Arc::ptr_eq(&moved.kept, &s.kept));
        assert!(moved.kept.seed_blocks.get().is_none() && moved.kept.patches.get().is_none());
        for (&(uv, p), &(uv0, p0)) in moved.seeds().iter().zip(s.seeds()) {
            assert!(uv == uv0 && p == turn.point(p0));
            assert!(moved.eval(uv.x, uv.y).dist(p) < 1e-12);
        }
        let far = turn.point(q);
        let uv = moved.project(far);
        assert!((moved.eval(uv.x, uv.y).dist(far) - on(&s)).abs() < 1e-12);
    }

    #[test]
    fn the_nearest_seeds_by_their_boxes_are_the_ones_found_in_turn() {
        let crate::geom::Surface::Nurbs(pole) = octant(2.0) else {
            unreachable!()
        };
        // Open patches, a tube whose seam's seeds lie on each other, and a patch with a pole;
        // each as built and moved, its seeds carried with it.
        let turn = Transform::rotation(Vec3::new(0.3, 1.0, -2.0), Vec3::new(2.0, -1.0, 2.0), 1.1);
        let mut surfaces: Vec<NurbsSurface> = uneven().into_iter().collect();
        surfaces.extend([warped(), tube(), (*pole).clone()]);
        let moved: Vec<NurbsSurface> = surfaces
            .iter()
            .map(|s| {
                s.seed_grid();
                s.transformed(&turn)
            })
            .collect();
        surfaces.extend(moved);
        let (mut by_blocks, mut in_turn) = (0, 0);
        for s in &surfaces {
            let bb = s.hull_box();
            let (c, d) = (bb.center(), bb.diagonal());
            let n = 14usize;
            let mut points: Vec<Point3> = Vec::new();
            for i in 0..n * n * n {
                let f = |k: usize| 1.6 * ((i / n.pow(k as u32)) % n) as f64 / (n - 1) as f64 - 0.8;
                points.push(c + Vec3::new(f(0), f(1), f(2)) * d);
            }
            // On the surface, on its seeds, and just off them.
            let seeds = s.seeds();
            points.extend(seeds.iter().step_by(7).map(|x| x.1));
            points.extend(
                seeds
                    .iter()
                    .step_by(11)
                    .map(|x| x.1 + Vec3::new(1e-9, -2e-9, 1e-9)),
            );
            points.push(Vec3::new(f64::NAN, 0.0, 0.0));
            for &p in &points {
                let want = s.nearest_seeds_in_turn(p);
                match s.nearest_seeds_by_blocks(p) {
                    Some(got) => {
                        by_blocks += 1;
                        for (g, w) in got.iter().zip(&want) {
                            assert_eq!(g.0.to_bits(), w.0.to_bits(), "{p:?}");
                            assert_eq!(
                                (g.1.x.to_bits(), g.1.y.to_bits()),
                                (w.1.x.to_bits(), w.1.y.to_bits())
                            );
                        }
                    }
                    None => in_turn += 1,
                }
                // Whichever way the seeds are found, the projection is the same.
                let uv = s.project(p);
                let again = s.project(p);
                assert_eq!(
                    (uv.x.to_bits(), uv.y.to_bits()),
                    (again.x.to_bits(), again.y.to_bits())
                );
            }
        }
        // Most points are settled by the boxes; those near a seam or a pole are not.
        assert!(
            by_blocks > 20 * in_turn && in_turn > 0,
            "{by_blocks} {in_turn}"
        );
        // Beside the tube's seam and the patch's pole the scan in turn decides.
        let seam = Vec3::new(1.05, 0.0, 1.5);
        assert!(tube().nearest_seeds_by_blocks(seam).is_none());
        assert!(
            pole.nearest_seeds_by_blocks(Vec3::new(0.0, 0.0, 2.1))
                .is_none()
        );
        assert!(warped().nearest_seeds_by_blocks(seam).is_some());
    }

    #[test]
    fn projection_beside_a_pole_lands_on_the_nearest_point() {
        // Every seed of the pole's row is the pole: Newton from there never moves, so the point
        // must be found from the seeds beside it. A point at polar angle θ and azimuth φ, lifted
        // off the sphere by `off`, has its nearest point on the sphere at the same angles.
        let r = 2.0;
        let s = octant(r);
        for theta in [3e-8f64, 1e-4, 1e-3, 1e-2, 5e-2] {
            for phi in [0.2f64, 0.7, 1.3] {
                for off in [0.0, 1e-3, -1e-3] {
                    let dir = Vec3::new(
                        theta.sin() * phi.cos(),
                        theta.sin() * phi.sin(),
                        theta.cos(),
                    );
                    let uv = s.project(dir * (r + off));
                    let got = s.eval(uv.x, uv.y);
                    assert!(got.dist(dir * r) < 1e-11, "θ {theta} φ {phi}: {got:?}");
                    // From a hint on the pole itself, where every `u` is the same point.
                    let uv = s.project_near(dir * (r + off), Vec2::new(0.5, 1.0));
                    let got = s.eval(uv.x, uv.y);
                    assert!(
                        got.dist(dir * r) < 1e-11,
                        "near: θ {theta} φ {phi}: {got:?}"
                    );
                    // From a hint just off the pole on the far boundary line, `u = 1` (φ = π/2):
                    // `su` there is under the tolerance, yet the hint's `u` is not the same point.
                    let uv = s.project_near(dir * (r + off), Vec2::new(1.0, 1.0 - 1e-8));
                    let got = s.eval(uv.x, uv.y);
                    assert!(
                        got.dist(dir * r) < 1e-11,
                        "beside: θ {theta} φ {phi}: {got:?}"
                    );
                }
            }
        }
    }

    /// A cubic's control point for three knots in a row `a`, `b`, `c` is its polar form there:
    /// `(a + b + c) / 3` for `t`, `(ab + bc + ca) / 3` for `t²` and `abc` for `t³`.
    #[test]
    fn a_piece_of_a_cubic_has_the_polar_form_of_its_own_knots() {
        let polar = |k: &[f64]| -> Vec<[f64; 4]> {
            k[1..k.len() - 1]
                .windows(3)
                .map(|w| {
                    let (a, b, c) = (w[0], w[1], w[2]);
                    [
                        (a + b + c) / 3.0,
                        (a * b + b * c + c * a) / 3.0,
                        a * b * c,
                        1.0,
                    ]
                })
                .collect()
        };
        let knots = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0];
        let (got_knots, got) = piece(3, &knots, &polar(&knots), (0.5, 2.25));
        let want_knots = [0.5, 0.5, 0.5, 0.5, 1.0, 2.0, 2.25, 2.25, 2.25, 2.25];
        assert_eq!(got_knots, want_knots);
        let same = |got: &[[f64; 4]], want: &[[f64; 4]]| {
            assert_eq!(got.len(), want.len());
            for (g, w) in got.iter().zip(want) {
                assert!(
                    (0..4).all(|i| (g[i] - w[i]).abs() < 1e-14),
                    "{g:?} vs {w:?}"
                );
            }
        };
        same(&got, &polar(&want_knots));
        same(&got[..1], &[[0.5, 0.25, 0.125, 1.0]]);
        // From one knot to another, with nothing to put in but the ends' own multiplicity.
        let (got_knots, got) = piece(3, &knots, &polar(&knots), (1.0, 2.0));
        assert_eq!(got_knots, [1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0]);
        same(&got, &polar(&got_knots));
    }
}
