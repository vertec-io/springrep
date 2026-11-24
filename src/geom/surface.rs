//! Surfaces: the analytic five in closed form, and the swept and ruled surfaces built on curves.
//!
//! Each surface's natural normal is `∂S/∂u × ∂S/∂v`; a face whose material lies on that side is
//! not reversed. On every surface swept about an axis `u` is the angle about it.

use std::f64::consts::{PI, TAU};

use super::{Curve, NurbsSurface};
use crate::{
    math::{Frame, Point3, Transform, Vec2, Vec3},
    numeric::solve2,
    tol,
};

#[derive(Clone, Debug, PartialEq)]
pub struct Plane {
    pub frame: Frame,
}

/// `origin + radius·radial(u) + v·z`.
#[derive(Clone, Debug, PartialEq)]
pub struct Cylinder {
    pub frame: Frame,
    pub radius: f64,
}

/// `origin + (radius + v·sin a)·radial(u) + v·cos a·z`, `a` the signed half angle, so `v` is arc
/// length along a generator and `radius` the radius at `v = 0`.
#[derive(Clone, Debug, PartialEq)]
pub struct Cone {
    pub frame: Frame,
    pub radius: f64,
    pub half_angle: f64,
}

/// `origin + radius·(cos v·radial(u) + sin v·z)`, `v` the latitude.
#[derive(Clone, Debug, PartialEq)]
pub struct Sphere {
    pub frame: Frame,
    pub radius: f64,
}

/// `origin + (major + minor·cos v)·radial(u) + minor·sin v·z`.
#[derive(Clone, Debug, PartialEq)]
pub struct Torus {
    pub frame: Frame,
    pub major: f64,
    pub minor: f64,
}

/// `curve(u) + v·dir`, `dir` a unit vector.
#[derive(Clone, Debug, PartialEq)]
pub struct Extrusion {
    pub curve: Curve,
    pub dir: Vec3,
}

/// `curve(v)` turned by `u` about the line through `origin` along the unit vector `axis`.
#[derive(Clone, Debug, PartialEq)]
pub struct Revolution {
    pub curve: Curve,
    pub origin: Point3,
    pub axis: Vec3,
}

/// `(1 − v)·a(u) + v·b(u)`, each curve's range mapped onto `u ∈ [0, 1]`.
#[derive(Clone, Debug, PartialEq)]
pub struct Ruled {
    pub a: Curve,
    pub a_range: (f64, f64),
    pub b: Curve,
    pub b_range: (f64, f64),
}

#[derive(Clone, Debug, PartialEq)]
pub enum Surface {
    Plane(Plane),
    Cylinder(Cylinder),
    Cone(Cone),
    Sphere(Sphere),
    Torus(Torus),
    Extrusion(Box<Extrusion>),
    Revolution(Box<Revolution>),
    Ruled(Box<Ruled>),
    /// A rational B-spline surface of any degree, as STEP files carry them.
    Nurbs(Box<NurbsSurface>),
}

fn d_radial(f: &Frame, u: f64) -> Vec3 {
    f.y * u.cos() - f.x * u.sin()
}

fn angle_of(local: Vec3) -> f64 {
    if local.x.abs() < 1e-300 && local.y.abs() < 1e-300 {
        0.0
    } else {
        local.y.atan2(local.x).rem_euclid(TAU)
    }
}

impl Surface {
    pub fn kind(&self) -> &'static str {
        match self {
            Surface::Plane(_) => "plane",
            Surface::Cylinder(_) => "cylinder",
            Surface::Cone(_) => "cone",
            Surface::Sphere(_) => "sphere",
            Surface::Torus(_) => "torus",
            Surface::Extrusion(_) => "extrusion",
            Surface::Revolution(_) => "revolution",
            Surface::Ruled(_) => "ruled",
            Surface::Nurbs(_) => "nurbs",
        }
    }

    pub fn eval(&self, u: f64, v: f64) -> Point3 {
        self.d1(u, v).0
    }

    /// The point and its two first partial derivatives.
    pub fn d1(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        match self {
            Surface::Plane(s) => (s.frame.at(u, v), s.frame.x, s.frame.y),
            Surface::Cylinder(s) => {
                let f = &s.frame;
                let e = f.radial(u);
                (
                    f.origin + e * s.radius + f.z * v,
                    d_radial(f, u) * s.radius,
                    f.z,
                )
            }
            Surface::Cone(s) => {
                let f = &s.frame;
                let (sa, ca) = s.half_angle.sin_cos();
                let rho = s.radius + v * sa;
                let e = f.radial(u);
                (
                    f.origin + e * rho + f.z * (v * ca),
                    d_radial(f, u) * rho,
                    e * sa + f.z * ca,
                )
            }
            Surface::Sphere(s) => {
                let f = &s.frame;
                let (sv, cv) = v.sin_cos();
                let e = f.radial(u);
                let r = s.radius;
                (
                    f.origin + (e * cv + f.z * sv) * r,
                    d_radial(f, u) * (r * cv),
                    (f.z * cv - e * sv) * r,
                )
            }
            Surface::Torus(s) => {
                let f = &s.frame;
                let (sv, cv) = v.sin_cos();
                let e = f.radial(u);
                let rho = s.major + s.minor * cv;
                (
                    f.origin + e * rho + f.z * (s.minor * sv),
                    d_radial(f, u) * rho,
                    (f.z * cv - e * sv) * s.minor,
                )
            }
            Surface::Extrusion(s) => {
                let (p, d) = s.curve.d1(u);
                (p + s.dir * v, d, s.dir)
            }
            Surface::Revolution(s) => {
                let (q, dq) = s.curve.d1(v);
                let rot = Transform::rotation(s.origin, s.axis, u);
                let p = rot.point(q);
                (p, s.axis.cross(p - s.origin), rot.vector(dq))
            }
            Surface::Ruled(s) => {
                let (ta, sa) = (
                    s.a_range.0 + u * (s.a_range.1 - s.a_range.0),
                    s.a_range.1 - s.a_range.0,
                );
                let (tb, sb) = (
                    s.b_range.0 + u * (s.b_range.1 - s.b_range.0),
                    s.b_range.1 - s.b_range.0,
                );
                let (pa, da) = s.a.d1(ta);
                let (pb, db) = s.b.d1(tb);
                (
                    pa.lerp(pb, v),
                    da * (sa * (1.0 - v)) + db * (sb * v),
                    pb - pa,
                )
            }
            Surface::Nurbs(s) => s.d1(u, v),
        }
    }

    /// The unit natural normal. At a pole, where the derivatives' cross product vanishes, it is
    /// the limit from inside the parameter domain.
    pub fn normal(&self, u: f64, v: f64) -> Vec3 {
        match self {
            Surface::Plane(s) => s.frame.z,
            Surface::Cylinder(s) => s.frame.radial(u),
            Surface::Cone(s) => {
                let (sa, ca) = s.half_angle.sin_cos();
                let rho = s.radius + v * sa;
                let n = s.frame.radial(u) * ca - s.frame.z * sa;
                if rho < 0.0 { -n } else { n }
            }
            Surface::Sphere(s) => {
                let (sv, cv) = v.sin_cos();
                s.frame.radial(u) * cv + s.frame.z * sv
            }
            Surface::Torus(s) => {
                let (sv, cv) = v.sin_cos();
                let n = s.frame.radial(u) * cv + s.frame.z * sv;
                if s.major + s.minor * cv < 0.0 { -n } else { n }
            }
            _ => {
                let (_, su, sv) = self.d1(u, v);
                if let Some(n) = su.cross(sv).unit().filter(|_| su.cross(sv).norm() > 1e-12) {
                    return n;
                }
                // At a pole where `su` vanishes it grows as `suv·(v − v_pole)`, so the normal is
                // the limit of `suv × sv` taken from inside the domain (and alike for `sv`).
                if let Surface::Nurbs(n) = self {
                    let ((u0, u1), (v0, v1)) = n.domain();
                    let [_, su, sv, _, suv, _] = n.d2(u, v);
                    let (limit, top) = if su.norm() <= sv.norm() {
                        (suv.cross(sv), v - v0 > v1 - v)
                    } else {
                        (su.cross(suv), u - u0 > u1 - u)
                    };
                    if let Some(l) = limit.unit() {
                        return if top { -l } else { l };
                    }
                }
                for (du, dv) in [(0.0, 1e-6), (0.0, -1e-6), (1e-6, 0.0), (-1e-6, 0.0)] {
                    let (_, su, sv) = self.d1(u + du, v + dv);
                    if let Some(n) = su.cross(sv).unit() {
                        return n;
                    }
                }
                Vec3::ZERO
            }
        }
    }

    /// Whether the line of `u` (of `v`, `along_v`) through `uv` is a single point, as at a pole. A
    /// vanishing derivative alone is not one: a spline whose end control points coincide stands still
    /// along its end line, which is still a curve, and the loop runs along it as along any other.
    /// The line is sampled across its parameter's period or domain, a pole's points all within the
    /// tolerance of each other.
    pub(crate) fn line_is_point(&self, uv: Vec2, along_v: bool) -> bool {
        let s = self;
        let (period, domain) = match (along_v, s) {
            (false, Surface::Nurbs(n)) => (s.u_period(), Some(n.domain().0)),
            (true, Surface::Nurbs(n)) => (s.v_period(), Some(n.domain().1)),
            (false, _) => (s.u_period(), None),
            (true, _) => (s.v_period(), None),
        };
        let at = |x: f64| {
            if along_v {
                s.eval(uv.x, x)
            } else {
                s.eval(x, uv.y)
            }
        };
        let here = s.eval(uv.x, uv.y);
        let own = if along_v { uv.y } else { uv.x };
        let samples: Vec<f64> = match (period, domain) {
            (Some(p), _) => (1..8).map(|k| own + p * k as f64 / 8.0).collect(),
            (None, Some((lo, hi))) => (0..=8).map(|k| lo + (hi - lo) * k as f64 / 8.0).collect(),
            // No known extent: a step either way, far beyond where a stationary point stays put.
            (None, None) => [-1e-2, 1e-2].map(|d| own + d * (1.0 + own.abs())).to_vec(),
        };
        samples
            .into_iter()
            .all(|x| at(x).dist(here) <= 2.0 * tol::LINEAR)
    }

    pub fn u_period(&self) -> Option<f64> {
        match self {
            Surface::Cylinder(_)
            | Surface::Cone(_)
            | Surface::Sphere(_)
            | Surface::Torus(_)
            | Surface::Revolution(_) => Some(TAU),
            Surface::Extrusion(s) => s.curve.period(),
            Surface::Nurbs(s) => s.u_period(),
            Surface::Plane(_) | Surface::Ruled(_) => None,
        }
    }

    pub fn v_period(&self) -> Option<f64> {
        match self {
            Surface::Torus(_) => Some(TAU),
            Surface::Revolution(s) => s.curve.period(),
            Surface::Nurbs(s) => s.v_period(),
            _ => None,
        }
    }

    /// Parameters strictly inside (u0, u1) where the surface may lose smoothness along `u`.
    pub fn u_breakpoints(&self, u0: f64, u1: f64) -> Vec<f64> {
        match self {
            Surface::Extrusion(s) => s.curve.breakpoints(u0, u1),
            Surface::Ruled(s) => {
                let map = |r: (f64, f64), c: &Curve| -> Vec<f64> {
                    let (lo, hi) = (r.0 + u0 * (r.1 - r.0), r.0 + u1 * (r.1 - r.0));
                    let (lo, hi) = (lo.min(hi), lo.max(hi));
                    c.breakpoints(lo, hi)
                        .into_iter()
                        .map(|t| (t - r.0) / (r.1 - r.0))
                        .collect()
                };
                let mut out = map(s.a_range, &s.a);
                out.extend(map(s.b_range, &s.b));
                out.sort_by(f64::total_cmp);
                out.dedup();
                out
            }
            Surface::Nurbs(s) => s.u_breakpoints(u0, u1),
            _ => Vec::new(),
        }
    }

    /// Parameters strictly inside (v0, v1) where the surface may lose smoothness along `v`.
    pub fn v_breakpoints(&self, v0: f64, v1: f64) -> Vec<f64> {
        match self {
            Surface::Revolution(s) => s.curve.breakpoints(v0, v1),
            Surface::Nurbs(s) => s.v_breakpoints(v0, v1),
            _ => Vec::new(),
        }
    }

    pub fn transformed(&self, t: &Transform) -> Surface {
        match self {
            Surface::Plane(s) => Surface::Plane(Plane {
                frame: s.frame.transformed(t),
            }),
            Surface::Cylinder(s) => Surface::Cylinder(Cylinder {
                frame: s.frame.transformed(t),
                ..*s
            }),
            Surface::Cone(s) => Surface::Cone(Cone {
                frame: s.frame.transformed(t),
                ..*s
            }),
            Surface::Sphere(s) => Surface::Sphere(Sphere {
                frame: s.frame.transformed(t),
                ..*s
            }),
            Surface::Torus(s) => Surface::Torus(Torus {
                frame: s.frame.transformed(t),
                ..*s
            }),
            Surface::Extrusion(s) => Surface::Extrusion(Box::new(Extrusion {
                curve: s.curve.transformed(t),
                dir: t.vector(s.dir),
            })),
            Surface::Revolution(s) => Surface::Revolution(Box::new(Revolution {
                curve: s.curve.transformed(t),
                origin: t.point(s.origin),
                axis: t.vector(s.axis),
            })),
            Surface::Ruled(s) => Surface::Ruled(Box::new(Ruled {
                a: s.a.transformed(t),
                b: s.b.transformed(t),
                ..(**s).clone()
            })),
            Surface::Nurbs(s) => Surface::Nurbs(Box::new(s.transformed(t))),
        }
    }

    /// `p`'s parameters on the generator of a cone, or the tube circle of a torus, in the half
    /// plane through `p` and the axis: the foot on that line or circle, though the other one in the
    /// same meridian plane may be nearer. Its `u` is `p`'s own angle, so it turns with `p` along a
    /// path, as a pcurve carried across a spindle torus's axis needs. Any other surface's
    /// [`Surface::project`].
    pub(crate) fn project_in_half_plane(&self, p: Point3) -> Vec2 {
        match self {
            Surface::Cone(s) => {
                let l = s.frame.to_local(p);
                let (sa, ca) = s.half_angle.sin_cos();
                Vec2::new(angle_of(l), (l.x.hypot(l.y) - s.radius) * sa + l.z * ca)
            }
            Surface::Torus(s) => {
                let l = s.frame.to_local(p);
                let rho = l.x.hypot(l.y);
                Vec2::new(angle_of(l), l.z.atan2(rho - s.major).rem_euclid(TAU))
            }
            _ => self.project(p),
        }
    }

    /// The parameters of the surface point nearest `p`. A periodic parameter comes back in its
    /// first period.
    pub fn project(&self, p: Point3) -> Vec2 {
        match self {
            Surface::Plane(s) => {
                let l = s.frame.to_local(p);
                Vec2::new(l.x, l.y)
            }
            Surface::Cylinder(s) => {
                let l = s.frame.to_local(p);
                Vec2::new(angle_of(l), l.z)
            }
            // The meridian plane through `p` holds two of a cone's generators, the one at `p`'s
            // own angle and the one half a turn round, which crosses the apex onto `p`'s side; so
            // it holds two of a torus's tube circles, and a spindle torus's far one reaches past
            // the axis to `p`'s side. The nearest point is on the nearer of the two.
            Surface::Cone(s) => {
                let l = s.frame.to_local(p);
                let rho = l.x.hypot(l.y);
                let (sa, ca) = s.half_angle.sin_cos();
                let u = angle_of(l);
                let near = Vec2::new(u, (rho - s.radius) * sa + l.z * ca);
                let far = Vec2::new((u + PI).rem_euclid(TAU), (-rho - s.radius) * sa + l.z * ca);
                self.nearer(p, near, far)
            }
            Surface::Sphere(s) => {
                let l = s.frame.to_local(p);
                Vec2::new(angle_of(l), l.z.atan2(l.x.hypot(l.y)))
            }
            Surface::Torus(s) => {
                let l = s.frame.to_local(p);
                let rho = l.x.hypot(l.y);
                let u = angle_of(l);
                let near = Vec2::new(u, l.z.atan2(rho - s.major).rem_euclid(TAU));
                if s.minor <= s.major {
                    return near;
                }
                let far = Vec2::new(
                    (u + PI).rem_euclid(TAU),
                    l.z.atan2(-rho - s.major).rem_euclid(TAU),
                );
                self.nearer(p, near, far)
            }
            Surface::Extrusion(s) if s.curve.spline().is_some() => {
                let (a, b) = s.curve.domain();
                let spline = s.curve.spline().expect("checked above");
                // The nearest point lies on the generator through the curve point nearest `p`
                // seen along the extrusion.
                let u = spline.nearest(p, Some(s.dir), a, b);
                let v = (p - s.curve.eval(u)).dot(s.dir);
                self.refine_projection(p, std::iter::once(Vec2::new(u, v)))
            }
            Surface::Extrusion(s) => {
                let (a, b) = finite_domain(&s.curve);
                let seeds = curve_seeds(&s.curve, a, b)
                    .into_iter()
                    .map(|u| Vec2::new(u, (p - s.curve.eval(u)).dot(s.dir)));
                self.refine_projection(p, seeds)
            }
            Surface::Revolution(s) => {
                let m = Meridian::of(s);
                let seeds = meridian_seeds(s, &m, p);
                let seeds = seeds.into_iter().map(|v| m.seed(s, p, v));
                self.refine_projection(p, seeds)
            }
            Surface::Ruled(_) => {
                let seeds = (0..=16)
                    .flat_map(|i| (0..=4).map(move |j| Vec2::new(i as f64 / 16.0, j as f64 / 4.0)));
                self.refine_projection(p, seeds)
            }
            Surface::Nurbs(s) => s.project(p),
        }
    }

    /// [`Surface::project`] onto a surface of revolution whose meridian `m` holds, for many points
    /// in turn: the same seeds and the same point. A seed turned to the angle of `p` lies as far
    /// from `p` as its place in the meridian plane lies from that of `p`, so only the seeds within
    /// rounding of the nearest there are placed on the surface to pick the nearest from.
    pub(crate) fn project_on(&self, m: &Meridian, p: Point3) -> Vec2 {
        let Surface::Revolution(s) = self else {
            return self.project(p);
        };
        let target = m.to_meridian(p);
        let seeds = meridian_seeds(s, m, p);
        let gaps: Vec<f64> = (0..seeds.len())
            .map(|k| match m.at.get(k) {
                Some(q) => q.dist(target),
                None => m.to_meridian(s.curve.eval(seeds[k])).dist(target),
            })
            .collect();
        let least = gaps.iter().copied().fold(f64::INFINITY, f64::min);
        let near = least + 1e-9 * (1.0 + least + target.norm());
        let nearest = (0..seeds.len())
            .filter(|&k| gaps[k] <= near)
            .map(|k| m.seed(s, p, seeds[k]));
        self.refine_projection(p, nearest)
    }

    /// As [`Surface::project`], with each periodic parameter moved by whole periods to lie
    /// nearest `hint`. At a pole, where every `u` is the same point, `u` stays at the hint's, so a
    /// curve through the pole keeps the `u` it arrived with.
    pub fn project_near(&self, p: Point3, hint: Vec2) -> Vec2 {
        self.project_near_at(p, hint).0
    }

    /// [`Surface::project_near`] with the surface point there and its unit natural normal, from
    /// the one evaluation that places the parameters.
    pub fn project_near_at(&self, p: Point3, hint: Vec2) -> (Vec2, Point3, Vec3) {
        let uv = self.project_near_raw(p, hint);
        let (q, su, sv) = self.d1(uv.x, uv.y);
        let out = self.toward(uv, hint, su);
        let n = match self {
            Surface::Plane(_)
            | Surface::Cylinder(_)
            | Surface::Cone(_)
            | Surface::Sphere(_)
            | Surface::Torus(_) => self.normal(out.x, out.y),
            _ => su
                .cross(sv)
                .unit()
                .filter(|_| su.cross(sv).norm() > 1e-12)
                .unwrap_or_else(|| self.normal(out.x, out.y)),
        };
        (out, q, n)
    }

    /// [`Surface::project_near`] with the partial derivatives there, a spline's from its span
    /// polynomials.
    pub(crate) fn project_near_d1(&self, p: Point3, hint: Vec2) -> (Vec2, Vec3, Vec3) {
        let uv = self.project_near_raw(p, hint);
        let (_, su, sv) = match self {
            Surface::Nurbs(n) => n.d1_by_polys(uv.x, uv.y),
            _ => self.d1(uv.x, uv.y),
        };
        (self.toward(uv, hint, su), su, sv)
    }

    /// The parameters [`Surface::project_near`] places `p` at, before they are moved toward `hint`.
    fn project_near_raw(&self, p: Point3, hint: Vec2) -> Vec2 {
        match self {
            // A point on the surface near the hint is found by Newton from the hint alone; the
            // global seed search is needed only when that lands off the surface.
            Surface::Extrusion(e) if e.curve.spline().is_some() => {
                match extrusion_near(e, p, hint.x) {
                    Some(uv) if self.eval(uv.x, uv.y).dist(p) <= tol::LINEAR => uv,
                    _ => self.project(p),
                }
            }
            Surface::Extrusion(_) | Surface::Revolution(_) | Surface::Ruled(_) => {
                let local = self.refine_projection(p, std::iter::once(hint));
                if self.eval(local.x, local.y).dist(p) <= tol::LINEAR {
                    local
                } else {
                    self.project(p)
                }
            }
            Surface::Nurbs(s) => s.project_near(p, hint),
            // A point stepped off a face, as a pcurve followed past its edge is, can lie nearer
            // the other nappe or the other tube circle, a part of the surface the face does not
            // reach; continued from the hint, it keeps its own half plane's foot.
            Surface::Cone(_) | Surface::Torus(_) => self.project_in_half_plane(p),
            _ => self.project(p),
        }
    }

    /// `uv` with each periodic parameter moved by whole periods to lie nearest `hint`, and `u` the
    /// hint's where the point is a pole: where the hint's `u` names the same point.
    fn toward(&self, uv: Vec2, hint: Vec2, su: Vec3) -> Vec2 {
        let shift = |x: f64, h: f64, period: Option<f64>| match period {
            Some(t) => x + ((h - x) / t).round() * t,
            None => x,
        };
        // A small `su` alone is not a pole: just beside one, the hint's `u` can be a whole
        // tolerance away, which a fit to half the tolerance cannot absorb.
        let pole = || self.eval(hint.x, uv.y).dist(self.eval(uv.x, uv.y)) <= 1e-3 * tol::LINEAR;
        let u = if su.norm() <= tol::LINEAR && pole() {
            hint.x
        } else {
            shift(uv.x, hint.x, self.u_period())
        };
        Vec2::new(u, shift(uv.y, hint.y, self.v_period()))
    }

    /// Of the parameters `near` and `far`, those whose point is nearer `p`, as the surface itself
    /// gives it; `near` where they are as near, as at a point on both sheets, where a spindle torus
    /// crosses its axis or at a cone's apex.
    fn nearer(&self, p: Point3, near: Vec2, far: Vec2) -> Vec2 {
        if self.eval(far.x, far.y).dist(p) < self.eval(near.x, near.y).dist(p) {
            far
        } else {
            near
        }
    }

    fn refine_projection(&self, p: Point3, seeds: impl Iterator<Item = Vec2>) -> Vec2 {
        let mut best = (f64::INFINITY, Vec2::ZERO);
        for s in seeds {
            let d = self.eval(s.x, s.y).dist(p);
            if d < best.0 {
                best = (d, s);
            }
        }
        let mut uv = best.1;
        for _ in 0..60 {
            let (q, su, sv) = self.d1(uv.x, uv.y);
            let r = q - p;
            // At a pole one partial vanishes, and the step along the other alone leaves it: a
            // seed on a revolution's axis would otherwise hold a point beside the pole there.
            let along = |s: Vec3| (s.dot(s) > 0.0).then(|| -s.dot(r) / s.dot(s));
            let Some((du, dv)) = solve2(
                su.dot(su),
                su.dot(sv),
                su.dot(sv),
                sv.dot(sv),
                -su.dot(r),
                -sv.dot(r),
            )
            .or_else(|| match (su.norm() <= sv.norm(), along(su), along(sv)) {
                (true, _, Some(dv)) => Some((0.0, dv)),
                (false, Some(du), _) => Some((du, 0.0)),
                _ => None,
            }) else {
                break;
            };
            let next = Vec2::new(uv.x + du, uv.y + dv);
            let done = du.abs() + dv.abs() <= 1e-15 * (1.0 + uv.x.abs() + uv.y.abs());
            if self.eval(next.x, next.y).dist(p) > r.norm() * (1.0 + 1e-12) + 1e-15 {
                break;
            }
            uv = next;
            if done {
                break;
            }
        }
        let wrap = |x: f64, period: Option<f64>, start: f64| match period {
            Some(t) => start + (x - start).rem_euclid(t),
            None => x,
        };
        let u0 = match self {
            Surface::Extrusion(s) => s.curve.domain().0,
            _ => 0.0,
        };
        let v0 = match self {
            Surface::Revolution(s) => s.curve.domain().0,
            _ => 0.0,
        };
        Vec2::new(
            wrap(uv.x, self.u_period(), u0),
            wrap(uv.y, self.v_period(), v0),
        )
    }
}

/// The nearest point of an extrusion by Newton's method on the curve parameter alone, from `u`:
/// across the extrusion's direction the surface is its curve, and `v` follows in closed form.
fn extrusion_near(e: &Extrusion, p: Point3, mut u: f64) -> Option<Vec2> {
    let flat = |w: Vec3| w - e.dir * w.dot(e.dir);
    let mut step = f64::INFINITY;
    for _ in 0..30 {
        let (c, d1, d2) = e.curve.d2(u);
        let (r, d1) = (flat(c - p), flat(d1));
        let f = r.dot(d1);
        let df = d1.norm2() + r.dot(flat(d2));
        if df <= 0.0 {
            return None;
        }
        step = f / df;
        u -= step;
        if step.abs() <= 1e-15 * (1.0 + u.abs()) {
            break;
        }
    }
    // Rounding may keep the last steps from vanishing; the caller checks the point is on `e`.
    (step.abs() <= 1e-9 * (1.0 + u.abs())).then(|| Vec2::new(u, (p - e.curve.eval(u)).dot(e.dir)))
}

fn finite_domain(c: &Curve) -> (f64, f64) {
    let (a, b) = c.domain();
    (a.max(-1e4), b.min(1e4))
}

/// A surface of revolution's curve read once for the projections onto it: its samples, and where
/// each lies in the meridian plane, by its distance from the axis and its place along it.
pub(crate) struct Meridian {
    frame: Frame,
    /// The end of the curve's domain the samples stop at.
    end: f64,
    samples: Vec<f64>,
    at: Vec<Vec2>,
}

impl Meridian {
    pub(crate) fn of(s: &Revolution) -> Meridian {
        let (a, b) = finite_domain(&s.curve);
        let frame = Frame::from_axis(s.origin, s.axis).unwrap_or(Frame::WORLD);
        let samples = curve_seeds(&s.curve, a, b);
        let mut m = Meridian {
            frame,
            end: b,
            samples,
            at: Vec::new(),
        };
        m.at = (m.samples.iter())
            .map(|&v| m.to_meridian(s.curve.eval(v)))
            .collect();
        m
    }

    fn to_meridian(&self, q: Point3) -> Vec2 {
        let l = self.frame.to_local(q);
        Vec2::new(l.x.hypot(l.y), l.z)
    }

    /// The parameters of the curve's point at `v`, turned to the angle of `p` about the axis.
    fn seed(&self, s: &Revolution, p: Point3, v: f64) -> Vec2 {
        let frame = &self.frame;
        let ang = angle_of(frame.to_local(p));
        let (q, dq) = s.curve.d1(v);
        let q = frame.to_local(q);
        // A point on the axis has no angle of its own: the meridian leaves it along its tangent,
        // which is the way the curve's points beside it lie.
        let at = if q.x.hypot(q.y) <= tol::LINEAR {
            let dq = frame.to_local(frame.origin + dq);
            angle_of(if v < self.end { dq } else { -dq })
        } else {
            angle_of(q)
        };
        Vec2::new((ang - at).rem_euclid(TAU), v)
    }
}

/// The curve's sample parameters, each local least of the distance in the meridian plane (from the
/// axis and along it) from `p` refined by Newton's method on that distance alone. The point of a
/// surface of revolution nearest `p` lies in the half-plane through `p` and the axis, so this is
/// the whole problem; the two-parameter steps stall where the curve crosses its own meridian
/// square to it, where the surface folds back on itself and its partials are parallel.
fn meridian_seeds(s: &Revolution, m: &Meridian, p: Point3) -> Vec<f64> {
    let frame = &m.frame;
    let target = m.to_meridian(p);
    let mut seeds = m.samples.clone();
    let gaps: Vec<f64> = m.at.iter().map(|q| q.dist(target)).collect();
    let n = seeds.len();
    for k in 0..n {
        let below = |j: usize| gaps[k] <= gaps[j];
        if (k == 0 || below(k - 1)) && (k + 1 == n || below(k + 1)) {
            // Where the meridian turns back the sample can sit on a stationary point of the
            // distance with the least either side of it, so each side is started on its own.
            let (lo, v, hi) = (
                seeds[k.saturating_sub(1)],
                seeds[k],
                seeds[(k + 1).min(n - 1)],
            );
            for (start, range) in [
                (v, (lo, hi)),
                (0.5 * (lo + v), (lo, v)),
                (0.5 * (v + hi), (v, hi)),
            ] {
                seeds.push(meridian_newton(s, frame, target, start, range));
            }
        }
    }
    seeds
}

/// Newton's method on half the squared meridian distance from `target`, from `v` within `range`,
/// each step taken only where it brings the curve nearer.
fn meridian_newton(
    s: &Revolution,
    frame: &Frame,
    target: Vec2,
    mut v: f64,
    range: (f64, f64),
) -> f64 {
    let state = |v: f64| {
        let (q, d1, d2) = s.curve.d2(v);
        let (l, l1, l2) = (
            frame.to_local(q),
            frame.dir_to_local(d1),
            frame.dir_to_local(d2),
        );
        let rho = l.x.hypot(l.y);
        (rho, l, l1, l2)
    };
    let gap2 = |v: f64| {
        let (rho, l, ..) = state(v);
        (rho - target.x).powi(2) + (l.z - target.y).powi(2)
    };
    let mut g = gap2(v);
    for _ in 0..30 {
        let (rho, l, l1, l2) = state(v);
        if rho <= 1e-12 {
            break;
        }
        let r1 = (l.x * l1.x + l.y * l1.y) / rho;
        let r2 = (l1.x * l1.x + l1.y * l1.y + l.x * l2.x + l.y * l2.y) / rho - r1 * r1 / rho;
        let (dr, dh) = (rho - target.x, l.z - target.y);
        let grad = dr * r1 + dh * l1.z;
        let hess = r1 * r1 + dr * r2 + l1.z * l1.z + dh * l2.z;
        if grad == 0.0 || hess <= 0.0 {
            break;
        }
        let mut step = -grad / hess;
        let mut taken = false;
        for _ in 0..8 {
            let next = (v + step).clamp(range.0, range.1);
            let gn = gap2(next);
            if gn < g {
                taken = true;
                let done = (next - v).abs() <= 1e-15 * (1.0 + v.abs());
                (v, g) = (next, gn);
                if done {
                    return v;
                }
                break;
            }
            step *= 0.5;
        }
        if !taken {
            break;
        }
    }
    v
}

fn curve_seeds(c: &Curve, a: f64, b: f64) -> Vec<f64> {
    let mut cuts = vec![a];
    cuts.extend(c.breakpoints(a, b));
    cuts.push(b);
    let mut out = Vec::new();
    for w in cuts.windows(2) {
        for i in 0..8 {
            out.push(w[0] + (w[1] - w[0]) * i as f64 / 8.0);
        }
    }
    out.push(b);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{BSpline, Circle};

    fn frame() -> Frame {
        Frame::new(Vec3::new(1.0, -2.0, 0.5), Vec3::new(0.2, 0.3, 1.0), Vec3::X).unwrap()
    }

    fn all() -> Vec<Surface> {
        let f = frame();
        let spline = BSpline::interpolate(&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(2.0, 1.0, 0.0),
            Vec3::new(4.0, 0.5, 0.0),
            Vec3::new(6.0, 2.0, 0.0),
        ])
        .unwrap();
        let profile = Curve::BSpline(
            BSpline::interpolate(&[
                Vec3::new(5.0, 0.0, 0.0),
                Vec3::new(6.0, 0.0, 2.0),
                Vec3::new(5.5, 0.0, 4.0),
                Vec3::new(7.0, 0.0, 6.0),
            ])
            .unwrap(),
        );
        let circle = Curve::Circle(Circle {
            frame: Frame::WORLD,
            radius: 3.0,
        });
        vec![
            Surface::Plane(Plane { frame: f }),
            Surface::Cylinder(Cylinder {
                frame: f,
                radius: 4.0,
            }),
            Surface::Cone(Cone {
                frame: f,
                radius: 4.0,
                half_angle: 0.4,
            }),
            Surface::Sphere(Sphere {
                frame: f,
                radius: 4.0,
            }),
            Surface::Torus(Torus {
                frame: f,
                major: 6.0,
                minor: 2.0,
            }),
            Surface::Extrusion(Box::new(Extrusion {
                curve: Curve::BSpline(spline),
                dir: Vec3::new(0.0, 0.6, 0.8),
            })),
            Surface::Revolution(Box::new(Revolution {
                curve: profile,
                origin: Vec3::ZERO,
                axis: Vec3::Z,
            })),
            Surface::Ruled(Box::new(Ruled {
                a: circle.clone(),
                a_range: (0.0, TAU),
                b: circle.transformed(&Transform::translation(Vec3::new(0.0, 0.0, 5.0))),
                b_range: (0.0, TAU),
            })),
        ]
    }

    #[test]
    fn derivatives_and_normals_match_differences() {
        for s in all() {
            let (u, v) = (0.7, 0.4);
            let h = 1e-6;
            let (_, su, sv) = s.d1(u, v);
            let fu = (s.eval(u + h, v) - s.eval(u - h, v)) / (2.0 * h);
            let fv = (s.eval(u, v + h) - s.eval(u, v - h)) / (2.0 * h);
            assert!(su.dist(fu) < 1e-7, "{} su", s.kind());
            assert!(sv.dist(fv) < 1e-7, "{} sv", s.kind());
            let n = s.normal(u, v);
            assert!(
                n.dist(su.cross(sv).unit().unwrap()) < 1e-12,
                "{} normal",
                s.kind()
            );
        }
    }

    #[test]
    fn a_spline_normal_at_a_pole_is_its_limit() {
        // At the octant's pole every `u` is the point (0, 0, r), where the sphere's normal is +z.
        let s = crate::geom::octant(2.0);
        for u in [0.0, 0.3, 0.5, 1.0] {
            let n = s.normal(u, 1.0);
            assert!(n.dist(Vec3::Z) < 1e-14, "u {u}: {n:?}");
        }
    }

    #[test]
    fn projection_beside_a_revolutions_pole_lands_on_the_point() {
        // The meridian starts on the axis, so the seed nearest a point a hair round from it is
        // the pole, where every `u` is one point and the partial along `u` vanishes.
        let mut pts: Vec<Vec3> = (0..4)
            .map(|i| Vec3::new(i as f64, 0.0, 0.125 * (i * i) as f64))
            .collect();
        for _ in 0..2 {
            let meridian = BSpline::interpolate(&pts).unwrap();
            pts.reverse();
            let (a, b) = meridian.domain();
            let s = Surface::Revolution(Box::new(Revolution {
                curve: Curve::BSpline(meridian),
                origin: Vec3::ZERO,
                axis: Vec3::Z,
            }));
            let near_pole = if s.eval(0.0, a).norm() < 1e-12 {
                a + 1e-3 * (b - a)
            } else {
                b - 1e-3 * (b - a)
            };
            for u in [0.0, 1.0, 2.5, 4.0, 5.9] {
                let p = s.eval(u, near_pole);
                let uv = s.project(p);
                assert!(s.eval(uv.x, uv.y).dist(p) < 1e-9, "at {u}: {uv:?}");
            }
        }
    }

    // Three meridians: a spline that wanders, a circle out of every meridian plane, whose image
    // there runs out and turns back on itself, and a line that starts on the axis.
    #[test]
    fn projection_from_a_meridian_read_once_is_projection() {
        let wander = Curve::BSpline(
            BSpline::interpolate(&[
                Vec3::new(5.0, 0.0, 0.0),
                Vec3::new(6.0, 0.5, 2.0),
                Vec3::new(5.5, -0.5, 4.0),
                Vec3::new(7.0, 0.0, 6.0),
                Vec3::new(6.0, 1.0, 7.0),
            ])
            .unwrap(),
        );
        let askew = Curve::Circle(Circle {
            frame: Frame::new(Vec3::new(5.0, 0.0, 0.0), Vec3::X, Vec3::Y).unwrap(),
            radius: 0.5,
        });
        let spoke = Curve::BSpline(
            BSpline::interpolate(&[
                Vec3::ZERO,
                Vec3::new(1.0, 0.0, 0.125),
                Vec3::new(2.0, 0.0, 0.5),
                Vec3::new(3.0, 0.0, 1.125),
            ])
            .unwrap(),
        );
        for (curve, span) in [(wander, 8.0), (askew, 6.0), (spoke, 4.0)] {
            let r = Revolution {
                curve,
                origin: Vec3::ZERO,
                axis: Vec3::Z,
            };
            let m = Meridian::of(&r);
            let (a, b) = finite_domain(&r.curve);
            let s = Surface::Revolution(Box::new(r));
            for k in 0..3000 {
                let x = k as f64 * 0.618;
                // Points on the surface, a hair off it, and anywhere about it.
                let p = match k % 3 {
                    0 => s.eval(x, a + (b - a) * (x * 0.37).sin().abs()),
                    1 => {
                        let (u, v) = (x, a + (b - a) * (x * 0.37).sin().abs());
                        s.eval(u, v) + s.normal(u, v) * (1e-4 * (x * 1.7).sin())
                    }
                    _ => Vec3::new(
                        span * (x * 1.3).sin(),
                        span * (x * 0.7).cos(),
                        span * (x * 2.1).sin(),
                    ),
                };
                assert_eq!(s.project(p), s.project_on(&m, p), "{k}");
            }
        }
    }

    #[test]
    fn projection_inverts_evaluation() {
        for s in all() {
            for &(u, v) in &[(0.7, 0.4), (2.5, 0.9), (5.9, 0.1)] {
                let (u, v) = match s {
                    Surface::Ruled(_) => (u / 6.0, v),
                    Surface::Extrusion(_) | Surface::Revolution(_) => (u, v * 4.0),
                    _ => (u, v),
                };
                let p = s.eval(u, v);
                let uv = s.project_near(p, Vec2::new(u, v));
                assert!(
                    s.eval(uv.x, uv.y).dist(p) < 1e-9,
                    "{} at {u},{v}: {uv:?}",
                    s.kind()
                );
                assert!(
                    (uv.x - u).abs() < 1e-7 && (uv.y - v).abs() < 1e-7,
                    "{}",
                    s.kind()
                );
            }
        }
    }

    /// A tilted cone, and a tilted spindle torus whose tube is wider than its reach from the axis,
    /// with points in a box round each, below the cone's apex and inside the torus's spindle too.
    fn far_sheets() -> Vec<(Surface, Vec<Point3>)> {
        let f = frame();
        let cone = Surface::Cone(Cone {
            frame: f,
            radius: 2.0,
            half_angle: 0.4,
        });
        let spindle = Surface::Torus(Torus {
            frame: f,
            major: 2.0,
            minor: 3.0,
        });
        let grid = |half: f64, below: f64| -> Vec<Point3> {
            let mut pts = Vec::new();
            for i in 0..9 {
                for j in 0..9 {
                    for k in 0..9 {
                        let at = |n: usize, lo: f64, hi: f64| lo + (hi - lo) * n as f64 / 8.0;
                        let l = Vec3::new(
                            at(i, -half, half) + 0.013,
                            at(j, -half, half) - 0.007,
                            at(k, below, half),
                        );
                        pts.push(f.to_world(l));
                    }
                }
            }
            pts
        };
        vec![(cone, grid(12.0, -16.0)), (spindle, grid(7.0, -7.0))]
    }

    /// The distance from a point at `rho` from the axis and `z` along it to a cone or a torus, in
    /// the meridian plane through the point: two lines through the cone's apex, two circles of the
    /// torus's tube, one each side of the axis, the nearer of the two the nearest point's.
    fn meridian_gap(s: &Surface, rho: f64, z: f64) -> f64 {
        match s {
            Surface::Cone(c) => {
                let (sa, ca) = c.half_angle.sin_cos();
                let near = ((rho - c.radius) * ca - z * sa).abs();
                let far = ((rho + c.radius) * ca + z * sa).abs();
                near.min(far)
            }
            Surface::Torus(t) => {
                let near = ((rho - t.major).hypot(z) - t.minor).abs();
                let far = ((rho + t.major).hypot(z) - t.minor).abs();
                near.min(far)
            }
            _ => unreachable!("cones and tori"),
        }
    }

    fn frame_of(s: &Surface) -> Frame {
        match s {
            Surface::Cone(c) => c.frame,
            Surface::Torus(t) => t.frame,
            _ => unreachable!("cones and tori"),
        }
    }

    #[test]
    fn a_point_projects_onto_the_far_sheet_of_a_cone_or_a_spindle_torus_where_that_is_nearer() {
        for (s, pts) in far_sheets() {
            let f = frame_of(&s);
            for p in pts {
                let l = f.to_local(p);
                let want = meridian_gap(&s, l.x.hypot(l.y), l.z);
                let uv = s.project(p);
                let got = s.eval(uv.x, uv.y).dist(p);
                assert!(
                    (got - want).abs() <= 1e-9,
                    "{} at {l:?}: {got} vs {want}",
                    s.kind()
                );
            }
        }
        // The cone and point of the report: the nearest point is on the other nappe, past the apex.
        let cone = Surface::Cone(Cone {
            frame: Frame::new(Vec3::ZERO, Vec3::Z, Vec3::X).unwrap(),
            radius: 2.0,
            half_angle: 0.4,
        });
        let p = Vec3::new(3.0, 0.0, -10.0);
        let uv = cone.project(p);
        assert!((cone.eval(uv.x, uv.y).dist(p) - meridian_gap(&cone, 3.0, -10.0)).abs() < 1e-12);
        assert!((uv.x - std::f64::consts::PI).abs() < 1e-12);
    }
}
