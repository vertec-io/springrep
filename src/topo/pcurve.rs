use crate::{
    geom::{BSpline, Curve, Surface},
    math::{Point3, Vec2, Vec3},
    numeric::solve2,
    tol,
};

/// Where a coedge runs in its face's parameter plane, as a function of its edge's parameter `t`.
#[derive(Clone, Debug, PartialEq)]
pub enum PCurve {
    /// `origin + t·dir`: every edge a sweep makes along or across its surface's parameter lines.
    Affine { origin: Vec2, dir: Vec2 },
    /// The edge's point projected onto the surface. A periodic parameter continues from `seed`,
    /// the coedge's parameters where it starts, so a seam's two coedges stay apart. Where the edge
    /// meets a pole, `u` is the one its direction there leaves the pole by.
    Projected { seed: Vec2 },
    /// `uv` read off a B-spline in `t` whose control points lie in the plane `z = 0`: a curve
    /// whose parameters are known along it, as a swept sheet's own, where projection lands on
    /// the wrong foot because the surface's map is nearly singular.
    Spline { uv: BSpline },
}

impl PCurve {
    /// The pcurve with its parameter plane mapped by `f`, an affine map.
    pub(crate) fn mapped(&self, f: &dyn Fn(Vec2) -> Vec2) -> PCurve {
        match self {
            PCurve::Affine { origin, dir } => PCurve::Affine {
                origin: f(*origin),
                dir: f(*origin + *dir) - f(*origin),
            },
            PCurve::Projected { seed } => PCurve::Projected { seed: f(*seed) },
            PCurve::Spline { uv } => PCurve::Spline {
                uv: map_spline(uv, f),
            },
        }
    }
}

/// `s` with its control points mapped by `f`, an affine map of the parameter plane, which a
/// B-spline follows exactly.
pub(crate) fn map_spline(s: &BSpline, f: &dyn Fn(Vec2) -> Vec2) -> BSpline {
    let ctrl = s
        .control_points()
        .iter()
        .map(|p| {
            let q = f(Vec2::new(p.x, p.y));
            Vec3::new(q.x, q.y, 0.0)
        })
        .collect();
    BSpline::new(s.degree(), s.knots().to_vec(), ctrl).expect("the same knots and degree")
}

/// `s` in the parameter `k·t`, for an edge whose parameter is scaled by `k > 0`.
pub(crate) fn stretch_spline(s: &BSpline, k: f64) -> BSpline {
    let knots = s.knots().iter().map(|x| x * k).collect();
    BSpline::new(s.degree(), knots, s.control_points().to_vec()).expect("knots scaled by k > 0")
}

/// The parameter-plane point of a [`PCurve::Spline`]'s spline at `t`.
pub(crate) fn spline_uv(s: &BSpline, t: f64) -> Vec2 {
    let p = s.eval(t);
    Vec2::new(p.x, p.y)
}

/// The rate of a [`PCurve::Spline`]'s spline at `t`.
pub(crate) fn spline_rate(s: &BSpline, t: f64) -> Vec2 {
    let d = s.eval_ders(t, 1)[1];
    Vec2::new(d.x, d.y)
}

/// How far a spline pcurve's point may lie from its curve's.
const SPLINE_FIT: f64 = tol::LINEAR / 4.0;

/// The cubic B-spline in `t` through `pts`, parameter-plane points of `s` along `curve` from
/// `range.0` to `range.1`, each at the parameter of `curve` its point is nearest, so that `s` at
/// the spline's point lies within [`SPLINE_FIT`] of `curve` all along. Each span whose middle
/// strays further gets a point there, its parameters moved from the spline's by the least change
/// that brings it onto the curve, so along a direction the map nearly drops it stays put (Nocedal
/// and Wright 2006, ch. 10). `None` where no such spline is found.
pub(crate) fn fit_spline(
    s: &Surface,
    curve: &Curve,
    range: (f64, f64),
    pts: &[Vec2],
) -> Option<BSpline> {
    let (lo, hi) = range;
    if pts.len() < 2 || hi <= lo {
        return None;
    }
    // A link of few points (a straight one, along a line of the plane) gets points between them,
    // each at the parameter its share of the chords in space gives as a first guess.
    const LEAST: usize = 32;
    let per = LEAST.div_ceil(pts.len() - 1).max(1);
    let mut dense = Vec::with_capacity((pts.len() - 1) * per + 1);
    for w in pts.windows(2) {
        for j in 0..per {
            dense.push(w[0].lerp(w[1], j as f64 / per as f64));
        }
    }
    dense.push(pts[pts.len() - 1]);
    let at: Vec<Point3> = dense.iter().map(|q| s.eval(q.x, q.y)).collect();
    let mut along = vec![0.0];
    for w in at.windows(2) {
        along.push(along[along.len() - 1] + w[0].dist(w[1]));
    }
    let total = along[along.len() - 1].max(1e-300);
    // Points all on one parameter line keep to it: there the pcurve is that line exactly.
    let held = if pts.iter().all(|q| q.x == pts[0].x) {
        Some(true)
    } else if pts.iter().all(|q| q.y == pts[0].y) {
        Some(false)
    } else {
        None
    };
    let mut ts = vec![lo];
    let mut uvs = vec![pts[0]];
    let spacing = (hi - lo) * 1e-9;
    for (k, q) in dense.iter().enumerate().take(dense.len() - 1).skip(1) {
        let prev = *ts.last().expect("starts at lo");
        let guess = lo + (hi - lo) * along[k] / total;
        let t = nearest_param(curve, at[k], guess.max(prev), (prev, hi));
        if t <= prev + spacing || t >= hi - spacing {
            continue;
        }
        let p = curve.eval(t);
        let q = if s.eval(q.x, q.y).dist(p) <= SPLINE_FIT {
            Some(*q)
        } else {
            settle(s, p, *q, held)
        };
        if let Some(q) = q {
            ts.push(t);
            uvs.push(q);
        }
    }
    ts.push(hi);
    uvs.push(pts[pts.len() - 1]);
    for _ in 0..12 {
        let at: Vec<Vec3> = uvs.iter().map(|q| Vec3::new(q.x, q.y, 0.0)).collect();
        let spline = BSpline::interpolate_at(&at, &ts).ok()?;
        let mut next_t = vec![ts[0]];
        let mut next_uv = vec![uvs[0]];
        let mut fits = true;
        for i in 0..ts.len() - 1 {
            for f in [0.25, 0.5, 0.75] {
                let t = ts[i] + (ts[i + 1] - ts[i]) * f;
                let uv = spline_uv(&spline, t);
                let p = curve.eval(t);
                if s.eval(uv.x, uv.y).dist(p) <= SPLINE_FIT {
                    continue;
                }
                fits = false;
                let mid = 0.5 * (ts[i] + ts[i + 1]);
                if mid > ts[i] + spacing && mid < ts[i + 1] - spacing {
                    let guess = spline_uv(&spline, mid);
                    next_t.push(mid);
                    next_uv.push(settle(s, curve.eval(mid), guess, held)?);
                }
                break;
            }
            next_t.push(ts[i + 1]);
            next_uv.push(uvs[i + 1]);
        }
        if fits {
            return Some(spline);
        }
        if next_t.len() == ts.len() {
            return None;
        }
        (ts, uvs) = (next_t, next_uv);
    }
    None
}

/// The parameter of `curve` in `(lo, hi)` nearest `p`, by Newton steps from `guess`, or the
/// nearest over the range where those leave it.
fn nearest_param(curve: &Curve, p: Point3, guess: f64, (lo, hi): (f64, f64)) -> f64 {
    let mut t = guess.clamp(lo, hi);
    for _ in 0..32 {
        let (c, d, dd) = curve.d2(t);
        let r = c - p;
        let slope = d.dot(d) + r.dot(dd);
        if slope <= 0.0 {
            break;
        }
        let step = r.dot(d) / slope;
        t = (t - step).clamp(lo, hi);
        if step.abs() <= 1e-15 * (1.0 + t.abs()) {
            return t;
        }
    }
    if curve.eval(t).dist(p) <= SPLINE_FIT {
        t
    } else {
        curve.project(p, lo, hi)
    }
}

/// Parameters near `uv` where `s` reaches `p` to within a tenth of [`SPLINE_FIT`], by damped
/// Gauss–Newton steps, each the least change that meets the linearised map; with `held`, `u`
/// (`Some(true)`) or `v` stays as it is.
fn settle(s: &Surface, p: Point3, mut uv: Vec2, held: Option<bool>) -> Option<Vec2> {
    for _ in 0..64 {
        let (at, su, sv) = s.d1(uv.x, uv.y);
        let r = p - at;
        if r.norm() <= 0.1 * SPLINE_FIT {
            return Some(uv);
        }
        let (a, b, d) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let damp = 1e-10 * (a + d);
        let (x, y) = match held {
            Some(true) => (0.0, sv.dot(r) / (d + damp)),
            Some(false) => (su.dot(r) / (a + damp), 0.0),
            None => solve2(a + damp, b, b, d + damp, su.dot(r), sv.dot(r))?,
        };
        if !(x.is_finite() && y.is_finite()) {
            return None;
        }
        uv = Vec2::new(uv.x + x, uv.y + y);
    }
    let at = s.eval(uv.x, uv.y);
    (at.dist(p) <= SPLINE_FIT).then_some(uv)
}

/// `curve(t)` projected onto `s` near `hint`. At a pole, where every `u` is the same point, `u` is
/// the angle the curve leaves the pole at towards `toward`, a parameter beside `t`: the limit of
/// its neighbours' `u`, which every caller finds alike whatever its sampling.
pub(crate) fn project_on(s: &Surface, curve: &Curve, t: f64, toward: f64, hint: Vec2) -> Vec2 {
    let (p, dp, ddp) = curve.d2(t);
    let uv = s.project_near(p, hint);
    let away = dp * (toward - t).signum();
    // The way the curve bends: its second derivative less its part along the curve.
    let bend = dp.unit().map(|d| ddp - d * ddp.dot(d));
    let (_, su, sv) = s.d1(uv.x, uv.y);
    match pole_param(s, p, uv, (su, sv), away, bend, hint, false) {
        Some(u) => Vec2::new(u, uv.y),
        None => match pole_param(s, p, uv, (su, sv), away, None, hint, true) {
            Some(v) => Vec2::new(uv.x, v),
            None => uv,
        },
    }
}

/// `p` projected onto `s` near `hint`, a point of an edge between its ends. Where it lands on a
/// pole along which `v` is free (an imported edge that misses its face by its tolerance can run
/// past one), `v` stays at the hint's, as the analytic surfaces keep a free `u`.
pub(crate) fn project_past(s: &Surface, p: Point3, hint: Vec2) -> Vec2 {
    let uv = s.project_near(p, hint);
    let (_, _, sv) = s.d1(uv.x, uv.y);
    if sv.norm() <= tol::LINEAR && s.line_is_point(uv, true) {
        Vec2::new(uv.x, hint.y)
    } else {
        uv
    }
}

/// [`project_past`] of `curve(t1)`, continued from `hint` at `curve(t0)`. A periodic parameter
/// continues to the image nearest the last one, which is the curve's own only while each step turns
/// it by well under half a period; beside a pole a short step can turn it by nearly half, either
/// way, so such a step is halved until each part turns it by under a quarter.
pub(crate) fn follow(s: &Surface, curve: &Curve, (t0, t1): (f64, f64), hint: Vec2) -> Vec2 {
    let uv = project_past(s, curve.eval(t1), hint);
    let jumps = |per: Option<f64>, d: f64| per.is_some_and(|p| d.abs() > 0.25 * p);
    let turned = jumps(s.u_period(), uv.x - hint.x) || jumps(s.v_period(), uv.y - hint.y);
    let mid = 0.5 * (t0 + t1);
    if !turned || mid == t0 || mid == t1 {
        return uv;
    }
    let h = follow(s, curve, (t0, mid), hint);
    follow(s, curve, (mid, t1), h)
}

/// At a pole of `s` where `u` (or, with `in_v`, `v`) is free, the value of it at which the
/// surface leaves the pole along `away`; `None` away from such a pole. A curve leaving along the
/// surface's axis, as every curve through the pole of a torus whose tube touches its axis does,
/// leaves at the angle it bends toward, `bend`. `su` and `sv` are the derivatives at `uv`.
#[allow(clippy::too_many_arguments)]
fn pole_param(
    s: &Surface,
    p: Point3,
    uv: Vec2,
    (su, sv): (Vec3, Vec3),
    away: Vec3,
    bend: Option<Vec3>,
    hint: Vec2,
    in_v: bool,
) -> Option<f64> {
    let d = away.unit_or_zero();
    if (if in_v { sv } else { su }).norm() > tol::LINEAR
        || d.norm() == 0.0
        || !s.line_is_point(uv, in_v)
    {
        return None;
    }
    let (x, period, near) = if in_v {
        let v = leaving_v(s, uv.x, d).unwrap_or_else(|| s.project(p + d * 1e-3).y);
        (v, s.v_period(), hint.y)
    } else {
        // Round a pole of an analytic surface any step along a direction in the tangent plane has
        // the same angle about it; a spline's pole need not be round, and a point beside it
        // projects onto a line leaving it at another angle.
        let u = match s {
            Surface::Nurbs(n) => leaving_u(n, uv.y, d),
            _ => None,
        };
        let axis = match s {
            Surface::Sphere(x) => Some(x.frame.z),
            Surface::Cone(x) => Some(x.frame.z),
            Surface::Torus(x) => Some(x.frame.z),
            _ => None,
        };
        let across = |v: Vec3, z: Vec3| v - z * v.dot(z);
        let step = match (axis, bend) {
            (Some(z), Some(k)) if across(d, z).norm() <= 1e-6 => across(k, z).unit().unwrap_or(d),
            _ => d,
        };
        let u = u.unwrap_or_else(|| s.project(p + step * 1e-3).x);
        (u, s.u_period(), hint.x)
    };
    Some(match period {
        Some(per) => x + ((near - x) / per).round() * per,
        None => x,
    })
}

/// At a pole of `s` along the line `u = u0`, where `v` is free, the `v` whose `u`-line leaves the
/// pole along `away`, by the surface's own derivative rather than by projecting a point near the
/// pole, which on a thin patch lands at any `v`. `None` where `v` has no known range.
fn leaving_v(s: &Surface, u0: f64, away: Vec3) -> Option<f64> {
    let (lo, hi) = match (s, s.v_period()) {
        (Surface::Nurbs(n), _) => n.domain().1,
        (_, Some(per)) => (0.0, per),
        _ => return None,
    };
    let inward = match s {
        Surface::Nurbs(n) if u0 > 0.5 * (n.domain().0.0 + n.domain().0.1) => -1.0,
        _ => 1.0,
    };
    let leaving = |v: f64| (s.d1(u0, v).1 * inward).unit();
    best_of(&leaving, away, (lo, hi))
}

/// At a pole of `n` along the line `v = v0`, where `u` is free, the `u` whose `v`-line leaves the
/// pole along `away`, by the surface's own derivative.
fn leaving_u(n: &crate::geom::NurbsSurface, v0: f64, away: Vec3) -> Option<f64> {
    let ((lo, hi), (v_lo, v_hi)) = n.domain();
    let inward = if v0 > 0.5 * (v_lo + v_hi) { -1.0 } else { 1.0 };
    let leaving = |u: f64| (n.d1(u, v0).2 * inward).unit();
    best_of(&leaving, away, (lo, hi))
}

/// Where `leaving` points along `away`, on `[lo, hi]`: the nearest of a run of samples, refined
/// round it.
fn best_of(leaving: &dyn Fn(f64) -> Option<Vec3>, away: Vec3, (lo, hi): (f64, f64)) -> Option<f64> {
    const SAMPLES: usize = 128;
    let cos = |x: f64| leaving(x).map_or(-2.0, |l| l.dot(away));
    let at = |k: usize| lo + (hi - lo) * k as f64 / SAMPLES as f64;
    let best = (0..=SAMPLES).max_by(|&a, &b| cos(at(a)).total_cmp(&cos(at(b))))?;
    // The cosine is flat at its peak, so rounding of 1e-16 moves the peak by 1e-8. The sine
    // between the directions, signed about the pole's normal, crosses zero there at full slope.
    for (j, k) in [
        (best.saturating_sub(1), best),
        (best, (best + 1).min(SAMPLES)),
    ] {
        let (mut a, mut b) = (at(j), at(k));
        let (Some(la), Some(lb)) = (leaving(a), leaving(b)) else {
            continue;
        };
        let Some(normal) = la.cross(lb).unit() else {
            continue;
        };
        let sine = |x: f64| leaving(x).map_or(0.0, |l| l.cross(away).dot(normal));
        let (sa, sb) = (sine(a), sine(b));
        for (x, s) in [(a, sa), (b, sb)] {
            if s.abs() <= 1e-15 && cos(x) > 0.0 {
                return Some(x);
            }
        }
        if (sa > 0.0) == (sb > 0.0) {
            // Where the direction leaves along the domain's end, rounding can put the crossing a
            // hair beyond it.
            for (x, s) in [(a, sa), (b, sb)] {
                if (x == lo || x == hi) && s.abs() <= 1e-9 && cos(x) > 0.0 {
                    return Some(x);
                }
            }
            continue;
        }
        for _ in 0..64 {
            let m = 0.5 * (a + b);
            if (sine(m) > 0.0) == (sa > 0.0) {
                a = m;
            } else {
                b = m;
            }
        }
        return Some(0.5 * (a + b));
    }
    let (mut a, mut b) = (at(best.saturating_sub(1)), at((best + 1).min(SAMPLES)));
    // Golden-section search on the bracket round the best sample.
    let g = 0.5 * (5f64.sqrt() - 1.0);
    for _ in 0..60 {
        let (c, d) = (b - g * (b - a), a + g * (b - a));
        if cos(c) >= cos(d) {
            b = d;
        } else {
            a = c;
        }
    }
    Some(0.5 * (a + b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        edit::{Affine, affine},
        geom::NurbsSurface,
        math::{Frame, Transform},
        primitive::sphere,
    };

    // The sheet S(s, t) = (s, t·(s² + δ), s²) on [-1, 1]² moves δ along t at s = 0: there its map
    // nearly drops rank. The curve S(x, x) crosses that line; the spline through its plane points
    // (x, x) at the curve's own parameter x is that line of the plane, and S over it is the curve.
    #[test]
    fn a_spline_pcurve_across_a_nearly_singular_line_keeps_its_parameters() {
        let delta = 1e-6;
        let g = [1.0 + delta, delta - 1.0, 1.0 + delta];
        let ctrl: Vec<Point3> = (0..6)
            .map(|i| {
                let (k, t) = (i / 2, if i % 2 == 0 { -1.0 } else { 1.0 });
                let s = [-1.0, 0.0, 1.0][k];
                Vec3::new(s, t * g[k], [1.0, -1.0, 1.0][k])
            })
            .collect();
        let knots = (
            vec![-1.0, -1.0, -1.0, 1.0, 1.0, 1.0],
            vec![-1.0, -1.0, 1.0, 1.0],
        );
        let sheet = Surface::Nurbs(Box::new(
            NurbsSurface::new((2, 1), knots, (3, 2), ctrl, None).unwrap(),
        ));
        let on = |x: f64| Vec3::new(x, x * (x * x + delta), x * x);
        let xs = [-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0];
        let pts: Vec<Vec3> = xs.iter().map(|&x| on(x)).collect();
        let curve = Curve::BSpline(BSpline::interpolate_at(&pts, &xs).unwrap());
        let plane: Vec<Vec2> = (0..=8)
            .map(|i| {
                let x = -1.0 + 0.25 * i as f64;
                Vec2::new(x, x)
            })
            .collect();
        let spline = fit_spline(&sheet, &curve, (-1.0, 1.0), &plane).unwrap();
        for i in 0..=200 {
            let t = -1.0 + 0.01 * i as f64;
            let uv = spline_uv(&spline, t);
            assert!(uv.dist(Vec2::new(t, t)) <= 1e-9, "{t}: {uv:?}");
            let off = sheet.eval(uv.x, uv.y).dist(curve.eval(t));
            assert!(off <= SPLINE_FIT, "{t}: {off:e}");
        }
    }

    // A direction turning once round on [-π, π] and one a hair past its end, as rounding leaves a
    // pole's line along a spline's seam: the end is the answer, not the flat peak's 1e-8 guess.
    #[test]
    fn a_direction_a_hair_past_the_domains_end_is_that_end() {
        use std::f64::consts::PI;
        let leaving = |x: f64| Some(Vec3::new(x.cos(), x.sin(), 0.0));
        for past in [1e-13, -1e-13, 0.0] {
            for (end, sign) in [(PI, 1.0), (-PI, -1.0)] {
                let a = end + sign * past;
                let found = best_of(&leaving, Vec3::new(a.cos(), a.sin(), 0.0), (-PI, PI)).unwrap();
                let off = (found - end).abs().min((found + end).abs());
                assert!(off <= 1e-10, "{past} past {end}: {found}");
            }
        }
    }

    // A ball stretched and sheared has a spline skin with a pole at each end of `v`. Each `u`'s
    // line leaves a pole along its own direction, and that direction gives the same `u` back
    // wherever the ball is placed: a coedge ending at the pole takes it, and the pole's own coedge
    // keeps the one found where the part was built. At the seam it is the domain's end exactly.
    #[test]
    fn a_direction_leaving_a_pole_names_its_parameter_to_a_hair_wherever_the_surface_is() {
        let map = Affine {
            m: [
                [0.1028, -0.9280, 0.5227],
                [-0.6650, 0.2041, 0.6977],
                [1.2920, 0.5654, 0.3662],
            ],
            t: Vec3::new(20.147, -35.624, 7.546),
        };
        let f = Frame::new(Vec3::new(-2.9, -1.3, -1.2), Vec3::Z, Vec3::X).unwrap();
        let ball = affine(&sphere(&f, 1.8, "sph").unwrap(), &map).unwrap();
        let s = &ball.face(ball.face_ids().next().unwrap()).surface;
        let Surface::Nurbs(n) = s else {
            panic!("a stretched ball's skin is a spline")
        };
        let ((u0, u1), (v0, v1)) = n.domain();
        let moves = [
            Vec3::ZERO,
            Vec3::new(100.0, -50.0, 300.0),
            Vec3::new(-731.0, 0.4, 0.0),
        ];
        let mut placed = Vec::new();
        for (k, moved) in moves.into_iter().enumerate() {
            for turn in 0..8 {
                let axis = Vec3::new(1.0, k as f64 - 1.0, 0.3 * turn as f64)
                    .unit()
                    .unwrap();
                let t = Transform::rotation(Vec3::ZERO, axis, 0.77 * turn as f64)
                    .then(&Transform::translation(moved));
                placed.push((k, t));
            }
        }
        for (k, t) in placed {
            let Surface::Nurbs(n) = s.transformed(&t) else {
                unreachable!()
            };
            for (pole, inward) in [(v0, 1.0), (v1, -1.0)] {
                assert!(
                    n.d1(0.3, pole).1.norm() <= tol::LINEAR,
                    "v {pole} is a pole"
                );
                for j in 0..=12 {
                    let off = if j == 0 || j == 12 { 0.0 } else { 0.01 };
                    let u = u0 + (u1 - u0) * j as f64 / 12.0 + off;
                    let away = (n.d1(u, pole).2 * inward).unit().unwrap();
                    let found = leaving_u(&n, pole, away).unwrap();
                    // The seam's two ends are one line.
                    let per = u1 - u0;
                    let off = (found - u) - ((found - u) / per).round() * per;
                    assert!(off.abs() <= 1e-10, "move {k}, pole {pole}: {found} for {u}");
                }
            }
        }
    }
}
