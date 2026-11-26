//! Area, volume, mass properties and tight bounding box, from the exact geometry rather than a
//! mesh.
//!
//! A face integral `∬ F du dv` over the trimmed parameter domain becomes a boundary integral by
//! Green's theorem: `∮ G dv` with `G(u, v) = ∫ F(s, v) ds` from a fixed `u`, or `−∮ H du` with the
//! inner integral along `v`. The inner integral runs along the parameter in which the integrand is
//! a low-order polynomial or trigonometric polynomial, so a fixed Gauss rule is exact to rounding;
//! the boundary integral is adaptive.

mod sample;
mod section;
mod slice;

pub use sample::{
    edge_length, extent, face_extent, face_samples, face_uv_bounds, solid_area, solid_bbox,
};
pub use section::{SectionEdge, section, section_area, section_curves};
pub use slice::{SliceCurve, slice, slices};

use std::{
    f64::consts::TAU,
    sync::atomic::{AtomicUsize, Ordering},
};

use crate::{
    check::coedge_params,
    geom::{IsoLine, NurbsSurface, Patch, Surface},
    math::{Aabb, Point3, Vec2, Vec3},
    numeric::{KRONROD_7, KRONROD_15, gauss, solve2},
    threads,
    topo::{
        Brep, CoedgeId, FaceId, FaceVolume, IntegralKey, LoopId, PCurve, SolidId, VertexId,
        pcurve::{spline_rate, spline_uv},
    },
};

/// Whether a face's inner integral runs along v, and the lowest value its loops reach along that
/// parameter, where every inner integral starts; `None` for a face without loops.
fn inner_start(b: &Brep, f: FaceId) -> Option<(bool, f64)> {
    let face = b.face(f);
    let inner_v = match &face.surface {
        Surface::Extrusion(_) | Surface::Ruled(_) => true,
        // The inner rule runs once per knot span, so it runs along the parameter with fewer.
        Surface::Nurbs(n) => n.span_counts().1 < n.span_counts().0,
        _ => false,
    };
    // Any start serves, since the inner integral from one start differs from that from another by
    // a function of the outer parameter alone, which integrates to nothing round a closed loop.
    // The lowest the loops reach keeps the inner integrals short.
    let (mut lo, mut has) = (f64::INFINITY, false);
    let tables = pcurve_tables(b, f);
    for (&lp, tables) in face.loops.iter().zip(tables) {
        for (&co, table) in b.lp(lp).coedges.iter().zip(tables) {
            let uvs = if table.is_empty() {
                b.coedge_uvs(co, &coedge_params(b, co, 4))
            } else {
                table.clone()
            };
            for uv in uvs {
                lo = lo.min(if inner_v { uv.y } else { uv.x });
                has = true;
            }
        }
    }
    has.then_some((inner_v, lo))
}

/// `∫ F` along the inner parameter from `lo` to `to` at the other parameter `at`, by a Gauss rule
/// of `order` on each smooth piece no longer than a turn, `F` given the point and `su × sv`.
#[allow(clippy::too_many_arguments)]
fn inner_rule<const K: usize>(
    s: &Surface,
    iso: Option<&IsoLine<'_>>,
    inner_v: bool,
    order: usize,
    (lo, at, to): (f64, f64, f64),
    density: &dyn Fn(Point3, Vec3) -> [f64; K],
) -> [f64; K] {
    let mut cuts = vec![lo];
    let (a, z) = (lo.min(to), lo.max(to));
    let bp = if inner_v {
        s.v_breakpoints(a, z)
    } else {
        s.u_breakpoints(a, z)
    };
    if to >= lo {
        cuts.extend(bp);
    } else {
        cuts.extend(bp.into_iter().rev());
    }
    cuts.push(to);
    // Only an angle needs pieces no longer than a turn. A plane's parameters, an extrusion's or a
    // ruled surface's v and a spline's are lengths or arbitrary, and the integrand is one
    // polynomial piece between cuts however long that is: a part metres across is not thousands
    // of turns.
    let turns = match s {
        Surface::Plane(_) | Surface::Nurbs(_) => false,
        Surface::Extrusion(_) | Surface::Ruled(_) => !inner_v,
        _ => true,
    };
    let mut sum = [0.0; K];
    for w in cuts.windows(2) {
        let n_chunks = if turns {
            ((w[1] - w[0]).abs() / TAU).ceil().max(1.0) as usize
        } else {
            1
        };
        for k in 0..n_chunks {
            let x0 = w[0] + (w[1] - w[0]) * k as f64 / n_chunks as f64;
            let x1 = w[0] + (w[1] - w[0]) * (k + 1) as f64 / n_chunks as f64;
            let (m, h) = ((x0 + x1) * 0.5, (x1 - x0) * 0.5);
            for &(x, wt) in gauss(order) {
                let q = m + h * x;
                let (p, su, sv) = match iso {
                    Some(line) => line.d1(q),
                    None if inner_v => s.d1(at, q),
                    None => s.d1(q, at),
                };
                let d = density(p, su.cross(sv));
                for i in 0..K {
                    sum[i] += wt * h * d[i];
                }
            }
        }
    }
    sum
}

/// `∬ F du dv` over the face's trimmed parameter domain, through Green's theorem, each value to
/// its own tolerance in `tol`. `inner` is the integral along the inner parameter, from `lo` to
/// `to` at `at`.
fn face_sum<const K: usize>(
    b: &Brep,
    f: FaceId,
    inner_v: bool,
    inner: &dyn Fn(f64, f64) -> [f64; K],
    tol: [f64; K],
    rel: f64,
) -> [f64; K] {
    let mut total = [0.0; K];
    let tables = pcurve_tables(b, f);
    for (&lp, tables) in b.face(f).loops.iter().zip(tables) {
        for (&co, table) in b.lp(lp).coedges.iter().zip(tables) {
            let v = coedge_integral(b, co, table, inner, inner_v, tol, rel);
            for i in 0..K {
                total[i] += v[i];
            }
        }
    }
    total
}

/// The face's area and its contribution `⅓∬ (p − c)·n dA` to the volume of its solid.
fn face_integrals(b: &Brep, f: FaceId, c: Point3, scale: f64, rel: f64) -> (f64, f64) {
    let key = IntegralKey::new(c, scale, rel, b.face(f).reversed);
    if let Some((k, v)) = *b.memo(f).integrals.lock().unwrap()
        && k == key
    {
        return v;
    }
    let v = face_integrals_anew(b, f, c, scale, rel);
    *b.memo(f).integrals.lock().unwrap() = Some((key, v));
    v
}

fn face_integrals_anew(b: &Brep, f: FaceId, c: Point3, scale: f64, rel: f64) -> (f64, f64) {
    let face = b.face(f);
    let s = &face.surface;
    let sigma = if face.reversed { -1.0 } else { 1.0 };
    let Some((inner_v, lo)) = inner_start(b, f) else {
        return (0.0, 0.0);
    };
    // Along the inner parameter a plane's and an extrusion's integrands are at most linear.
    // A rough integral settles only a sign, which a shorter rule along a smooth piece cannot flip.
    let order = if matches!(s, Surface::Plane(_) | Surface::Extrusion(_)) {
        8
    } else if rel > EXACT {
        12
    } else {
        24
    };
    let density = |p: Point3, n: Vec3| [n.norm(), (p - c).dot(n) / 3.0];
    let inner = |at: f64, to: f64| -> [f64; 2] {
        let h = to - lo;
        match (s, inner_v) {
            // Along a plane's u and an extrusion's generators the normal is constant and the
            // point linear, so the inner integral has a closed form.
            (Surface::Plane(_), false) | (Surface::Extrusion(_), true) => {
                let (p, su, sv) = if inner_v { s.d1(at, lo) } else { s.d1(lo, at) };
                let n = su.cross(sv);
                let along = if inner_v { sv } else { su };
                let vol = ((p - c).dot(n) * h + along.dot(n) * h * h * 0.5) / 3.0;
                return [n.norm() * h, vol];
            }
            // Round a cylinder the normal is r·(a cos u + b sin u), a trigonometric polynomial
            // of degree one, as is the point.
            (Surface::Cylinder(cy), false) => {
                let f = &cy.frame;
                let r = cy.radius;
                let (a, b) = (f.y.cross(f.z), f.z.cross(f.x));
                let w = f.origin + f.z * at - c;
                let m = 0.5 * (lo + to);
                let (sm, cm) = m.sin_cos();
                let (s2m, c2m) = (2.0 * m).sin_cos();
                let sh2 = (0.5 * h).sin();
                let sh = h.sin();
                let int_cos = 2.0 * cm * sh2;
                let int_sin = 2.0 * sm * sh2;
                let int_cos2 = 0.5 * h + 0.5 * c2m * sh;
                let int_sin2 = 0.5 * h - 0.5 * c2m * sh;
                let int_sc = 0.5 * s2m * sh;
                let vol = r
                    * (w.dot(a) * int_cos
                        + w.dot(b) * int_sin
                        + r * (f.x.dot(a) * int_cos2
                            + (f.x.dot(b) + f.y.dot(a)) * int_sc
                            + f.y.dot(b) * int_sin2))
                    / 3.0;
                return [r * h, vol];
            }
            (Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_), false) => {
                if let Some(turned) = round_axis(s, at, (lo, to), c) {
                    return turned;
                }
            }
            _ => {}
        }
        let line = match s {
            Surface::Nurbs(n) => Some(n.iso(inner_v, at)),
            _ => None,
        };
        inner_rule(s, line.as_ref(), inner_v, order, (lo, at, to), &density)
    };
    let rel = rel.max(noise_floor(c, scale));
    let tol = [rel * scale * scale, rel * scale * scale * scale];
    let [a, v] = face_sum(b, f, inner_v, &inner, tol, rel);
    (sigma * a, v)
}

/// `∫ |n| du` and `⅓∫ (p − c)·n du` round the axis of a cone, sphere or torus from `lo` to `to`
/// at the parameter `v`, `n = su × sv`, in closed form. Each is a circle `ρ(v)` from the axis at
/// the height `ζ(v)` turned through `u`: `p = o + ρ·e(u) + ζ·z` with `e = x cos u + y sin u`, so
/// `n = ρ·(ζ′·(a cos u + b sin u) − ρ′·x × y)` with `a = y × z` and `b = z × x`. The point and the
/// normal are trigonometric polynomials of degree one in `u`, and their product of degree two.
fn round_axis(s: &Surface, v: f64, (lo, to): (f64, f64), c: Point3) -> Option<[f64; 2]> {
    let (f, [rho, zeta, d_rho, d_zeta]) = match s {
        Surface::Cone(k) => {
            let (sa, ca) = k.half_angle.sin_cos();
            (&k.frame, [k.radius + v * sa, v * ca, sa, ca])
        }
        Surface::Sphere(k) => {
            let (sv, cv) = v.sin_cos();
            let r = k.radius;
            (&k.frame, [r * cv, r * sv, -r * sv, r * cv])
        }
        Surface::Torus(k) => {
            let (sv, cv) = v.sin_cos();
            let r = k.minor;
            (&k.frame, [k.major + r * cv, r * sv, -r * sv, r * cv])
        }
        _ => return None,
    };
    let (a, b, w) = (f.y.cross(f.z), f.z.cross(f.x), f.x.cross(f.y));
    let from = f.origin + f.z * zeta - c;
    let h = to - lo;
    let m = 0.5 * (lo + to);
    let (sm, cm) = m.sin_cos();
    let (s2m, c2m) = (2.0 * m).sin_cos();
    let (sh2, sh) = ((0.5 * h).sin(), h.sin());
    let (int_cos, int_sin) = (2.0 * cm * sh2, 2.0 * sm * sh2);
    let int_cos2 = 0.5 * h + 0.5 * c2m * sh;
    let int_sin2 = 0.5 * h - 0.5 * c2m * sh;
    let int_sc = 0.5 * s2m * sh;
    let along = d_zeta
        * (from.dot(a) * int_cos
            + from.dot(b) * int_sin
            + rho
                * (f.x.dot(a) * int_cos2
                    + (f.x.dot(b) + f.y.dot(a)) * int_sc
                    + f.y.dot(b) * int_sin2));
    let across = d_rho * (from.dot(w) * h + rho * (f.x.dot(w) * int_cos + f.y.dot(w) * int_sin));
    Some([
        rho.abs() * d_rho.hypot(d_zeta) * h,
        rho * (along - across) / 3.0,
    ])
}

/// The face's contributions to the volume, first and second moments of its solid about `c`:
/// `⅓∬ x·n`, `¼∬ xᵢ (x·n)` and `⅕∬ xᵢxⱼ (x·n) dA` with `x = p − c`, by the divergence theorem,
/// in the order 1, x, y, z, xx, yy, zz, xy, yz, zx.
fn face_moments(b: &Brep, f: FaceId, c: Point3, scale: f64) -> [f64; 10] {
    let m = face_moments_and_area(b, f, c, scale);
    std::array::from_fn(|i| m[i])
}

/// [`face_moments`] with the face's area after them. One pass gives both, so a solid measured for
/// its volume and area and then its centre integrates each of these faces once.
fn face_moments_and_area(b: &Brep, f: FaceId, c: Point3, scale: f64) -> [f64; 11] {
    let key = IntegralKey::new(c, scale, EXACT, b.face(f).reversed);
    if let Some((k, v)) = *b.memo(f).moments.lock().unwrap()
        && k == key
    {
        return v;
    }
    let v = face_moments_anew(b, f, c, scale);
    *b.memo(f).moments.lock().unwrap() = Some((key, v));
    v
}

/// The face's area and volume share about `c` to the exact tolerance: by [`face_integrals`] where
/// its inner integral has a closed form, which the moments have not, and with the moments
/// otherwise, whose pass costs the same.
fn face_exact(b: &Brep, f: FaceId, c: Point3, scale: f64) -> (f64, f64) {
    let face = b.face(f);
    if matches!(
        face.surface,
        Surface::Plane(_)
            | Surface::Extrusion(_)
            | Surface::Cylinder(_)
            | Surface::Cone(_)
            | Surface::Sphere(_)
            | Surface::Torus(_)
    ) {
        return face_integrals(b, f, c, scale, EXACT);
    }
    let m = face_moments_and_area(b, f, c, scale);
    (if face.reversed { -m[10] } else { m[10] }, m[0])
}

fn face_moments_anew(b: &Brep, f: FaceId, c: Point3, scale: f64) -> [f64; 11] {
    let s = &b.face(f).surface;
    let Some((inner_v, lo)) = inner_start(b, f) else {
        return [0.0; 11];
    };
    // Along a plane's u the integrands are polynomials of degree at most three.
    let order = if matches!(s, Surface::Plane(_) | Surface::Extrusion(_)) {
        8
    } else {
        24
    };
    let density = |p: Point3, n: Vec3| {
        let x = p - c;
        let flux = x.dot(n);
        [
            flux / 3.0,
            x.x * flux / 4.0,
            x.y * flux / 4.0,
            x.z * flux / 4.0,
            x.x * x.x * flux / 5.0,
            x.y * x.y * flux / 5.0,
            x.z * x.z * flux / 5.0,
            x.x * x.y * flux / 5.0,
            x.y * x.z * flux / 5.0,
            x.z * x.x * flux / 5.0,
            n.norm(),
        ]
    };
    let inner = |at: f64, to: f64| {
        let line = match s {
            Surface::Nurbs(n) => Some(n.iso(inner_v, at)),
            _ => None,
        };
        inner_rule(s, line.as_ref(), inner_v, order, (lo, at, to), &density)
    };
    let (s3, s4, s5) = (scale.powi(3), scale.powi(4), scale.powi(5));
    let rel = EXACT.max(noise_floor(c, scale));
    let s2 = scale * scale;
    let tol = [s3, s4, s4, s4, s5, s5, s5, s5, s5, s5, s2].map(|x| rel * x);
    face_sum(b, f, inner_v, &inner, tol, rel)
}

/// How many even steps along a coedge its table of parameters takes.
const TABLE: usize = 64;

/// Each projected coedge's surface parameters at [`TABLE`] + 1 even steps along its edge in its
/// own direction, loop by loop; empty for a straight pcurve. They seed every projection the face's
/// integrals make, and its memo keeps them for each later integral.
fn pcurve_tables(b: &Brep, f: FaceId) -> &[Vec<Vec<Vec2>>] {
    b.memo(f).pcurves.get_or_init(|| {
        b.face(f)
            .loops
            .iter()
            .map(|&lp| {
                b.lp(lp)
                    .coedges
                    .iter()
                    .map(|&co| pcurve_table(b, co))
                    .collect()
            })
            .collect()
    })
}

/// One coedge's entry of [`pcurve_tables`].
fn pcurve_table(b: &Brep, co: CoedgeId) -> Vec<Vec2> {
    let coedge = b.coedge(co);
    if matches!(coedge.pcurve, PCurve::Affine { .. }) || b.edge(coedge.edge).curve.is_none() {
        return Vec::new();
    }
    b.coedge_uvs(co, &coedge_params(b, co, TABLE))
}

/// How fast the parameters move along a curve with derivative `dp` on a surface with partial
/// derivatives `su` and `sv`.
fn uv_rate(su: Vec3, sv: Vec3, dp: Vec3) -> Vec2 {
    solve2(
        su.dot(su),
        su.dot(sv),
        su.dot(sv),
        sv.dot(sv),
        su.dot(dp),
        sv.dot(dp),
    )
    .map(|(a, b)| Vec2::new(a, b))
    .unwrap_or(Vec2::ZERO)
}

/// `∮ G dv` (or `−∮ H du`) along one coedge, in its direction, each value to its own tolerance.
/// `table` is its entry of [`pcurve_tables`].
fn coedge_integral<const K: usize>(
    b: &Brep,
    co: CoedgeId,
    table: &[Vec2],
    inner: &dyn Fn(f64, f64) -> [f64; K],
    inner_v: bool,
    tol: [f64; K],
    rel: f64,
) -> [f64; K] {
    let coedge = b.coedge(co);
    let edge = b.edge(coedge.edge);
    let surface = &b.face(b.coedge_face(co)).surface;
    let (t0, t1) = edge.range;
    let (a, z) = if coedge.reversed { (t1, t0) } else { (t0, t1) };
    // Between the two entries either side of `t`, Newton starts nearer its foot than from either.
    // A straight blend stays between them, so a periodic parameter keeps their period; a cubic
    // through their rates can overshoot by one beside a pole.
    let hint = |t: f64| -> Vec2 {
        let x = ((t - a) / (z - a) * TABLE as f64).clamp(0.0, TABLE as f64);
        let k = (x.floor() as usize).min(TABLE - 1);
        let s = x - k as f64;
        table[k] * (1.0 - s) + table[k + 1] * s
    };
    let uv_and_rate = |t: f64| -> (Vec2, Vec2) {
        match &coedge.pcurve {
            PCurve::Affine { origin, dir } => (*origin + *dir * t, *dir),
            PCurve::Spline { uv } => (spline_uv(uv, t), spline_rate(uv, t)),
            PCurve::Projected { .. } => {
                let curve = edge.curve.as_ref().expect("projected pcurves have curves");
                let (p, dp) = curve.d1(t);
                let (uv, su, sv) = surface.project_near_d1(p, hint(t));
                (uv, uv_rate(su, sv, dp))
            }
        }
    };
    let f = |t: f64| -> [f64; K] {
        let (uv, rate) = uv_and_rate(t);
        if inner_v {
            inner(uv.x, uv.y).map(|g| -g * rate.x)
        } else {
            inner(uv.y, uv.x).map(|g| g * rate.y)
        }
    };
    let mut cuts = vec![t0];
    // A fitted intersection curve has hundreds of breakpoints. A short rule settles the smooth
    // piece between two of them, where the halving below takes far more points to find each one.
    if let Some(c) = &edge.curve {
        cuts.extend(c.breakpoints(t0, t1));
    }
    if let PCurve::Spline { uv } = &coedge.pcurve {
        cuts.extend(uv.breakpoints(t0, t1));
        cuts.sort_by(f64::total_cmp);
    }
    // Where a straight pcurve crosses the surface's knot lines the inner integral has a kink too:
    // along a pole, which has no curve of its own, the adaptive rule would have to find them.
    if let PCurve::Affine { origin, dir } = &coedge.pcurve {
        let (a, b) = (*origin + *dir * t0, *origin + *dir * t1);
        let knots = |lo: f64, hi: f64, u: bool| {
            if u {
                surface.u_breakpoints(lo.min(hi), lo.max(hi))
            } else {
                surface.v_breakpoints(lo.min(hi), lo.max(hi))
            }
        };
        for (u, d, o) in [(true, dir.x, origin.x), (false, dir.y, origin.y)] {
            if d.abs() <= f64::EPSILON {
                continue;
            }
            let (lo, hi) = if u { (a.x, b.x) } else { (a.y, b.y) };
            for k in knots(lo, hi, u) {
                let t = (k - o) / d;
                if t > t0 && t < t1 {
                    cuts.push(t);
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
    }
    cuts.push(t1);
    let sum = match iso_line(table) {
        Some((v_fixed, at)) => {
            // Along a parameter line of the surface the integral runs along that parameter
            // between the edge's two ends, and no point of the edge is projected. Beside a cusp of
            // the line the edge's points have two feet within its tolerance, on either side of
            // the cusp, and a projection's pick between them is noise no rule can settle.
            let (start, end) = if coedge.reversed {
                (table[TABLE], table[0])
            } else {
                (table[0], table[TABLE])
            };
            let (s0, s1) = if v_fixed {
                (start.x, end.x)
            } else {
                (start.y, end.y)
            };
            let g = |s: f64| -> [f64; K] {
                match (inner_v, v_fixed) {
                    (true, true) => inner(s, at).map(|g| -g),
                    (false, false) => inner(s, at),
                    _ => [0.0; K],
                }
            };
            let mut cuts = vec![s0];
            if inner_v == v_fixed {
                let (lo, hi) = (s0.min(s1), s0.max(s1));
                let mut knots = if v_fixed {
                    surface.u_breakpoints(lo, hi)
                } else {
                    surface.v_breakpoints(lo, hi)
                };
                knots.retain(|&k| k > lo && k < hi);
                if s1 < s0 {
                    knots.reverse();
                }
                cuts.extend(knots);
            }
            cuts.push(s1);
            integrate_pieces(&g, &cuts, tol, rel)
        }
        None => integrate_pieces(&f, &cuts, tol, rel),
    };
    if coedge.reversed {
        sum.map(|x| -x)
    } else {
        sum
    }
}

/// Whether a projected coedge's table runs along one parameter line of its surface: `v` the same
/// at every entry (`true`) or `u`, and that value.
fn iso_line(table: &[Vec2]) -> Option<(bool, f64)> {
    let first = *table.first()?;
    let same = |x: f64, y: f64| (x - y).abs() <= 1e-12 * (1.0 + y.abs());
    if table.iter().all(|p| same(p.y, first.y)) {
        Some((true, first.y))
    } else if table.iter().all(|p| same(p.x, first.x)) {
        Some((false, first.x))
    } else {
        None
    }
}

/// `∫ f` over the pieces between successive `cuts`, each value to its own tolerance in `tol`.
fn integrate_pieces<const K: usize>(
    f: &dyn Fn(f64) -> [f64; K],
    cuts: &[f64],
    tol: [f64; K],
    rel: f64,
) -> [f64; K] {
    let order = if rel > EXACT { 8 } else { 16 };
    let mut sum = [0.0; K];
    for w in cuts.windows(2) {
        // Most pieces are short and smooth, and a pair of nested low rules settles them with a
        // third of the points the halving below takes to accept one. To a rough tolerance the
        // shorter pair settles most of them with half as many again.
        let short = (rel > EXACT)
            .then(|| kronrod_k(f, w[0], w[1], tol, &KRONROD_7))
            .flatten();
        let r = match short.or_else(|| kronrod_k(f, w[0], w[1], tol, &KRONROD_15)) {
            Some(r) => r,
            None => {
                let (whole, _) = gauss_k(f, w[0], w[1], order);
                adapt_k(f, (w[0], w[1]), whole, tol, order, 0, f64::INFINITY)
            }
        };
        for i in 0..K {
            sum[i] += r[i];
        }
    }
    sum
}

/// How many rounding errors of a rule's own sum two rules may differ by and still agree: halving
/// cannot settle a difference that small, since each half's sum is rounded as coarsely.
const ROUNDING: f64 = 64.0 * f64::EPSILON;

/// The share of its last miss a halving must cut an adaptive rule's miss to, below which the miss
/// is the rule's error and above which it is noise in the integrand.
const STALL: f64 = 0.75;

/// The largest share of the size of an integrand's values that its noise is taken to reach. A face
/// integral's boundary integrand is itself an inner integral, by a fixed rule on each knot span
/// the inner parameter crosses, so it jumps by that rule's error, near 1e-8 of it on a fitted
/// sheet's `|su × sv|`, each time the boundary crosses a knot.
const NOISE: f64 = 1e-6;

/// `∫ f` over `[a, b]` by a Kronrod rule where the Gauss rule inside it agrees to `tol`: the Gauss
/// rule is then that accurate itself, and the Kronrod rule more so.
fn kronrod_k<const K: usize>(
    f: &dyn Fn(f64) -> [f64; K],
    a: f64,
    b: f64,
    tol: [f64; K],
    rule: &[(f64, f64, f64)],
) -> Option<[f64; K]> {
    let (m, h) = ((a + b) * 0.5, (b - a) * 0.5);
    let (mut kr, mut ga, mut size) = ([0.0; K], [0.0; K], [0.0; K]);
    for &(x, wk, wg) in rule {
        let (l, r) = if x == 0.0 {
            (f(m), [0.0; K])
        } else {
            (f(m - h * x), f(m + h * x))
        };
        for i in 0..K {
            kr[i] += wk * (l[i] + r[i]);
            ga[i] += wg * (l[i] + r[i]);
            size[i] += wk * (l[i].abs() + r[i].abs());
        }
    }
    let ok = (0..K).all(|i| {
        let floor = (1e-14 * kr[i].abs()).max(ROUNDING * size[i]);
        ((kr[i] - ga[i]) * h).abs() <= tol[i].max(floor * h.abs())
    });
    ok.then(|| kr.map(|x| x * h))
}

/// The rule's sum, and the sum of its terms' sizes, which its rounding error is a share of.
fn gauss_k<const K: usize>(
    f: &dyn Fn(f64) -> [f64; K],
    a: f64,
    b: f64,
    order: usize,
) -> ([f64; K], [f64; K]) {
    let (m, h) = ((a + b) * 0.5, (b - a) * 0.5);
    let (mut s, mut size) = ([0.0; K], [0.0; K]);
    for &(x, w) in gauss(order) {
        let v = f(m + h * x);
        for i in 0..K {
            s[i] += w * v[i];
            size[i] += w * v[i].abs();
        }
    }
    (s.map(|x| x * h), size.map(|x| x * h.abs()))
}

/// `∫ f` over `[a, b]`, halved until the rule over the halves agrees with `whole`, the rule over
/// the whole, to `tol`. `last` is how many times its tolerance the halves of the piece this one was
/// cut from missed by. The tolerance halves with the piece, so where the rules differ by noise in
/// `f` (small jumps where an inner rule changes) and not by the rule's error, which shrinks by
/// orders with each halving, they miss by about as many times at every depth: once a halving no
/// longer cuts the miss to below `STALL` of the last, and the miss is within `NOISE` of the size
/// of `f` there, the halves are as near as `f` allows and are taken (the round-off test of
/// Piessens et al. 1983). Across a kink the miss halves each time.
fn adapt_k<const K: usize>(
    f: &dyn Fn(f64) -> [f64; K],
    (a, b): (f64, f64),
    whole: [f64; K],
    tol: [f64; K],
    order: usize,
    depth: u32,
    last: f64,
) -> [f64; K] {
    let m = 0.5 * (a + b);
    let ((l, ls), (r, rs)) = (gauss_k(f, a, m, order), gauss_k(f, m, b, order));
    let mut both = [0.0; K];
    for i in 0..K {
        both[i] = l[i] + r[i];
    }
    let miss = (0..K)
        .map(|i| {
            let floor = (1e-14 * both[i].abs()).max(ROUNDING * (ls[i] + rs[i]));
            (both[i] - whole[i]).abs() / tol[i].max(floor)
        })
        .fold(0.0, f64::max);
    // A miss larger than `NOISE` of the size of `f` is the rule's error, however slowly it shrinks.
    let small = (0..K).all(|i| (both[i] - whole[i]).abs() <= NOISE * (ls[i] + rs[i]));
    if miss <= 1.0 || depth >= 20 || (small && miss >= STALL * last) {
        return both;
    }
    let half = tol.map(|t| t * 0.5);
    let x = adapt_k(f, (a, m), l, half, order, depth + 1, miss);
    let y = adapt_k(f, (m, b), r, half, order, depth + 1, miss);
    let mut out = [0.0; K];
    for i in 0..K {
        out[i] = x[i] + y[i];
    }
    out
}

/// A centre and size for conditioning the integrals. The vertices' box serves: it costs no
/// surface search, which matters when every solid of a large arena is measured in turn. Fewer than
/// four vertices span no volume (a ring has one), and a tolerance sized to them is unreachable.
fn reference(b: &Brep) -> (Point3, f64) {
    if b.vertex_ids().nth(3).is_none() {
        return moment_reference(b);
    }
    let mut bb = Aabb::EMPTY;
    for v in b.vertex_ids() {
        bb.add(b.vertex(v).point);
    }
    if bb.is_empty() {
        (Vec3::ZERO, 1.0)
    } else {
        (bb.center(), bb.diagonal().max(1.0))
    }
}

/// The signed area loop `lp` encloses in its face's parameter plane, `∮ u dv` integrated on its
/// coedges' exact parameter curves. The checker's sampled outline cannot settle the sign of a sliver
/// narrower than its chords' sag, and this can.
pub(crate) fn loop_area(b: &Brep, lp: LoopId) -> f64 {
    // About the loop's own corner and to its own size: about the plane's origin a sliver's area
    // is the small difference of large terms, below a tolerance sized to the whole model. A
    // billionth of its box settles a sliver's sign; tighter, the quadrature splits to its limit.
    let uvs: Vec<Vec2> = b
        .lp(lp)
        .coedges
        .iter()
        .flat_map(|&co| b.coedge_uvs(co, &coedge_params(b, co, 8)))
        .collect();
    let Some(&corner) = uvs.first() else {
        return 0.0;
    };
    let mut bb = Aabb::EMPTY;
    for uv in &uvs {
        bb.add(Point3::new(uv.x, uv.y, 0.0));
    }
    let size = bb.diagonal();
    let tol = [ROUGH * size * size];
    b.lp(lp)
        .coedges
        .iter()
        .map(|&co| {
            let table = pcurve_table(b, co);
            coedge_integral(b, co, &table, &|_, to| [to - corner.x], false, tol, EXACT)[0]
        })
        .sum()
}

pub fn face_area(b: &Brep, f: FaceId) -> f64 {
    let (c, scale) = reference(b);
    face_integrals(b, f, c, scale, EXACT).0
}

/// The area of every face in the arena.
pub fn area(b: &Brep) -> f64 {
    let (c, scale) = moment_reference(b);
    let faces: Vec<FaceId> = b.face_ids().collect();
    each_face(b, &faces, |f| face_exact(b, f, c, scale).0)
        .into_iter()
        .sum()
}

/// `work` on each of `faces`, in their order, spread over the machine's threads where there is
/// enough of it to pay for them, heaviest first so that no long face starts last.
///
/// Each face's integrals are its own and are summed in the caller's order, so the totals are the
/// same to the last bit on any number of threads. A face's cost is about its boundary integral's
/// smooth pieces, one per span of each edge's curve: an assembly's few fitted sheets can hold
/// tens of thousands of them, and its hundreds of planes and cylinders a handful each.
fn each_face<T: Send>(b: &Brep, faces: &[FaceId], work: impl Fn(FaceId) -> T + Sync) -> Vec<T> {
    let order: Vec<(usize, usize)> = faces
        .iter()
        .enumerate()
        .map(|(k, &f)| (pieces(b, f), k))
        .collect();
    // A small solid is measured on the calling thread alone: its few pieces take less time
    // than starting threads would.
    let most = (order.iter().map(|w| w.0).sum::<usize>() / PIECES_PER_THREAD).min(faces.len());
    each_face_on(most, faces, order, work)
}

/// How many smooth pieces face `f`'s boundary integrals take at the least: one per span of each
/// of its edges' curves.
fn pieces(b: &Brep, f: FaceId) -> usize {
    b.face(f)
        .loops
        .iter()
        .flat_map(|&lp| &b.lp(lp).coedges)
        .map(|&co| {
            let edge = b.edge(b.coedge(co).edge);
            let (t0, t1) = edge.range;
            1 + edge
                .curve
                .as_ref()
                .map_or(0, |c| c.breakpoints(t0, t1).len())
        })
        .sum()
}

/// [`each_face`] on at most `most` threads, `order` each face's weight and place in `faces`.
fn each_face_on<T: Send>(
    most: usize,
    faces: &[FaceId],
    mut order: Vec<(usize, usize)>,
    work: impl Fn(FaceId) -> T + Sync,
) -> Vec<T> {
    if most <= 1 {
        return faces.iter().map(|&f| work(f)).collect();
    }
    order.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    let next = AtomicUsize::new(0);
    let runs = threads::spread(most, || {
        let mut mine = Vec::new();
        while let Some(&(_, k)) = order.get(next.fetch_add(1, Ordering::Relaxed)) {
            mine.push((k, work(faces[k])));
        }
        mine
    });
    let mut done: Vec<(usize, T)> = runs.into_iter().flatten().collect();
    done.sort_by_key(|d| d.0);
    done.into_iter().map(|d| d.1).collect()
}

/// How many smooth boundary pieces, summed over an arena's faces, are worth a thread of their own
/// in [`each_face`].
const PIECES_PER_THREAD: usize = 512;

/// The boundary integrals' tolerance relative to the model's size: area to `EXACT·size²`, volume
/// to `EXACT·size³`.
const EXACT: f64 = 1e-14;
/// The tolerance that settles only a volume's sign.
const ROUGH: f64 = 1e-9;

/// The least relative tolerance a model `scale` across about `c` can be integrated to: every point
/// it evaluates is rounded to the last digits of its coordinates, `ε·|c|`, so far from the origin
/// [`EXACT`] is below that noise and refinement never stops.
fn noise_floor(c: Point3, scale: f64) -> f64 {
    64.0 * f64::EPSILON * c.norm() / scale
}

fn solid_volume_to(b: &Brep, s: SolidId, rel: f64) -> f64 {
    let (c, scale) = reference(b);
    b.solid_faces(s)
        .into_iter()
        .map(|f| face_integrals(b, f, c, scale, rel).1)
        .sum()
}

/// The enclosed volume: positive for an outward-oriented solid.
pub fn solid_volume(b: &Brep, s: SolidId) -> f64 {
    solid_volume_to(b, s, EXACT)
}

/// The enclosed volume where its sign is all that matters: a rough integral when that is far
/// from zero, the exact one otherwise. It costs a fraction of [`solid_volume`] on a solid with
/// many spline faces, and each face's rough share is worked out once for its geometry, so a
/// face a boolean passes through costs nothing the next time.
pub fn solid_volume_sign(b: &Brep, s: SolidId) -> f64 {
    let (c, scale) = reference(b);
    // Over a closed shell the vector areas sum to zero, so each face's share may be taken about
    // its own anchor and moved to `c` by its vector area.
    let rough: f64 = b
        .solid_faces(s)
        .into_iter()
        .map(|f| {
            let fv = face_volume(b, f);
            fv.volume + (fv.anchor - c).dot(fv.vector_area) / 3.0
        })
        .sum();
    if rough.abs() > 1e4 * ROUGH * scale * scale * scale {
        rough
    } else {
        solid_volume(b, s)
    }
}

/// The volume the closed surface of `faces` encloses, outward normals positive, where its sign is
/// all that matters: [`solid_volume_sign`] for any set of faces.
pub(crate) fn enclosed_volume_sign(b: &Brep, faces: &[FaceId]) -> f64 {
    let (c, scale) = reference(b);
    let rough: f64 = faces
        .iter()
        .map(|&f| {
            let fv = face_volume(b, f);
            fv.volume + (fv.anchor - c).dot(fv.vector_area) / 3.0
        })
        .sum();
    if rough.abs() > 1e4 * ROUGH * scale * scale * scale {
        rough
    } else {
        let rel = EXACT;
        faces
            .iter()
            .map(|&f| face_integrals(b, f, c, scale, rel).1)
            .sum()
    }
}

/// Works out face `f`'s rough volume share now, in its memo: a copy moved rigidly keeps it.
pub(crate) fn prime_face_volume(b: &Brep, f: FaceId) {
    face_volume(b, f);
}

fn face_volume(b: &Brep, f: FaceId) -> FaceVolume {
    *b.memo(f).volume.get_or_init(|| {
        let (anchor, scale) = face_reference(b, f);
        FaceVolume {
            anchor,
            volume: face_integrals(b, f, anchor, scale, ROUGH).1,
            vector_area: vector_area(b, f, anchor, scale),
        }
    })
}

/// [`reference`] for one face, from the vertices its edges end at and three points along each
/// edge: a disc's rim is one edge from a vertex back to it, and a tolerance sized to that one
/// point is unreachable however large the disc.
fn face_reference(b: &Brep, f: FaceId) -> (Point3, f64) {
    let mut bb = Aabb::EMPTY;
    for &lp in &b.face(f).loops {
        for &co in &b.lp(lp).coedges {
            let e = b.edge(b.coedge(co).edge);
            bb.add(b.vertex(e.start).point);
            bb.add(b.vertex(e.end).point);
            if let Some(c) = &e.curve {
                let (t0, t1) = e.range;
                for k in 1..4 {
                    bb.add(c.eval(t0 + (t1 - t0) * k as f64 / 4.0));
                }
            }
        }
    }
    if bb.is_empty() {
        (Vec3::ZERO, 1.0)
    } else {
        (bb.center(), bb.diagonal().max(1.0))
    }
}

/// `∬ n dA` over the face, as `½∮ (p − a) × dp` round its loops (Stokes), to the rough tolerance.
fn vector_area(b: &Brep, f: FaceId, a: Point3, scale: f64) -> Vec3 {
    let tol = ROUGH * scale * scale;
    let mut sum = Vec3::ZERO;
    for &lp in &b.face(f).loops {
        for &co in &b.lp(lp).coedges {
            let coedge = b.coedge(co);
            let edge = b.edge(coedge.edge);
            let Some(c) = &edge.curve else { continue };
            let g = |t: f64| {
                let (p, dp) = c.d1(t);
                (p - a).cross(dp) * 0.5
            };
            let (t0, t1) = edge.range;
            let mut cuts = vec![t0];
            cuts.extend(c.breakpoints(t0, t1));
            cuts.push(t1);
            let mut along = Vec3::ZERO;
            for w in cuts.windows(2) {
                along += integrate_vec(&g, w[0], w[1], tol);
            }
            sum = if coedge.reversed {
                sum - along
            } else {
                sum + along
            };
        }
    }
    sum
}

fn gauss_vec(g: &dyn Fn(f64) -> Vec3, a: f64, b: f64) -> Vec3 {
    let (m, h) = ((a + b) * 0.5, (b - a) * 0.5);
    let mut s = Vec3::ZERO;
    for &(x, w) in gauss(8) {
        s += g(m + h * x) * w;
    }
    s * h
}

/// `∫ g` over `[a, b]`, halved until two levels agree to `tol`.
fn integrate_vec(g: &dyn Fn(f64) -> Vec3, a: f64, b: f64, tol: f64) -> Vec3 {
    fn rec(g: &dyn Fn(f64) -> Vec3, a: f64, b: f64, whole: Vec3, tol: f64, depth: u32) -> Vec3 {
        let m = 0.5 * (a + b);
        let (l, r) = (gauss_vec(g, a, m), gauss_vec(g, m, b));
        let both = l + r;
        if (both - whole).norm() <= tol || depth >= 20 {
            return both;
        }
        rec(g, a, m, l, 0.5 * tol, depth + 1) + rec(g, m, b, r, 0.5 * tol, depth + 1)
    }
    rec(g, a, b, gauss_vec(g, a, b), tol, 0)
}

/// The volume of every solid in the arena.
pub fn volume(b: &Brep) -> f64 {
    let (c, scale) = moment_reference(b);
    let faces: Vec<FaceId> = b.solid_ids().flat_map(|s| b.solid_faces(s)).collect();
    each_face(b, &faces, |f| face_exact(b, f, c, scale).1)
        .into_iter()
        .sum()
}

/// [`volume`] to `rel` of the model's size cubed rather than to the last digits: each face's
/// integral stops refining once two levels agree that well.
pub fn volume_to(b: &Brep, rel: f64) -> f64 {
    b.solid_ids().map(|s| solid_volume_to(b, s, rel)).sum()
}

/// [`area`] and [`volume`] together, integrating each face once.
pub fn area_and_volume(b: &Brep) -> (f64, f64) {
    let (c, scale) = moment_reference(b);
    let faces: Vec<FaceId> = b.face_ids().collect();
    let per = each_face(b, &faces, |f| face_exact(b, f, c, scale));
    let area = per.iter().map(|x| x.0).sum();
    let volume = b
        .solid_ids()
        .map(|s| {
            b.solid_faces(s)
                .into_iter()
                .map(|f| per[f.index()].1)
                .sum::<f64>()
        })
        .sum();
    (area, volume)
}

/// Volume, centre and inertia of solids at unit density, in mm: `inertia` is the tensor
/// `∭ (|r|²·I − r rᵀ) dV` about the centre, so its diagonal holds the moments about axes through
/// the centre and its off-diagonal the negated products.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct MassProperties {
    pub volume: f64,
    pub centre: Point3,
    pub inertia: [[f64; 3]; 3],
}

/// The mass properties of every solid in the arena together.
pub fn mass_properties(b: &Brep) -> MassProperties {
    let (c, scale) = moment_reference(b);
    let faces: Vec<FaceId> = b.solid_ids().flat_map(|s| b.solid_faces(s)).collect();
    let mut m = [0.0; 10];
    for fm in each_face(b, &faces, |f| face_moments(b, f, c, scale)) {
        for i in 0..10 {
            m[i] += fm[i];
        }
    }
    from_moments(m, c)
}

/// The mass properties of one solid.
pub fn solid_mass_properties(b: &Brep, s: SolidId) -> MassProperties {
    let (c, scale) = moment_reference(b);
    let mut m = [0.0; 10];
    for f in b.solid_faces(s) {
        let fm = face_moments(b, f, c, scale);
        for i in 0..10 {
            m[i] += fm[i];
        }
    }
    from_moments(m, c)
}

/// [`reference`] from the geometry's loose box rather than its vertices: a moment that symmetry
/// makes zero settles only to a tolerance in the model's size, and a torus or a ball has one
/// vertex or none to give it.
fn moment_reference(b: &Brep) -> (Point3, f64) {
    let bb = loose_bbox(b);
    if bb.is_empty() {
        (Vec3::ZERO, 1.0)
    } else {
        (bb.center(), bb.diagonal().max(1.0))
    }
}

/// Mass properties from the volume and the first and second moments about `c`, in
/// [`face_moments`]' order: the second moments move to the centre by `S − V d dᵀ`.
fn from_moments(m: [f64; 10], c: Point3) -> MassProperties {
    let v = m[0];
    let d = if v == 0.0 {
        [0.0; 3]
    } else {
        [m[1] / v, m[2] / v, m[3] / v]
    };
    let about_c = [[m[4], m[7], m[9]], [m[7], m[5], m[8]], [m[9], m[8], m[6]]];
    let mut s = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            s[i][j] = about_c[i][j] - v * d[i] * d[j];
        }
    }
    let trace = s[0][0] + s[1][1] + s[2][2];
    let mut inertia = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            inertia[i][j] = if i == j { trace } else { 0.0 } - s[i][j];
        }
    }
    MassProperties {
        volume: v,
        centre: c + Vec3::new(d[0], d[1], d[2]),
        inertia,
    }
}

/// The Betti numbers of the arena's point set: its separate pieces, the tunnels through them and
/// the voids they enclose, where pieces of a solid that touch are one piece.
///
/// The point set is its interior and its boundary complex `B`, so its Euler characteristic is
/// `χ(B) − χ(M)/2`, `M` the closed surface `B` becomes with every contact pulled apart: an edge
/// that `k` pairs of faces use becomes `k` edges, and a vertex `m` fans of faces meet at becomes
/// `m`. Each solid is one piece. Faces' outer sides joined across each edge's empty wedges bound
/// the regions of empty space; a group whose faces, with their outward normals, enclose negative
/// volume faces into the region it bounds, a void, and one enclosing positive volume is the
/// outside of the material it wraps. Tunnels are what is left: `b1 = b0 + b2 − χ`.
/// How many regions of material the fans of faces `fans` at `v`, a point where pieces touch and no
/// edge they share, bound round it: each fan's faces cut a small sphere round `v` along a circle,
/// and the circles part the sphere into regions, material and empty by turns. Two fans bound one
/// region of material when each lies on the other's material side and no fan parts them, read
/// from a point of each fan's faces `scale` from `v`.
fn material_regions(b: &Brep, v: VertexId, fans: &[Vec<FaceId>], scale: f64) -> usize {
    let p = b.vertex(v).point;
    let near: Vec<Option<Point3>> = fans.iter().map(|fs| point_near(b, fs, p, scale)).collect();
    let n = fans.len();
    // side[i][j]: whether fan j lies on fan i's material side.
    let side: Vec<Vec<Option<bool>>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| near[j].and_then(|q| material_side(b, &fans[i], v, q)))
                .collect()
        })
        .collect();
    let mut parent: Vec<usize> = (0..n).collect();
    fn root(p: &mut [usize], mut i: usize) -> usize {
        while p[i] != i {
            i = p[i];
        }
        i
    }
    for i in 0..n {
        for j in i + 1..n {
            let mutual = side[i][j] == Some(true) && side[j][i] == Some(true);
            let parted = (0..n).any(|k| k != i && k != j && side[k][i] != side[k][j]);
            if mutual && !parted {
                let (x, y) = (root(&mut parent, i), root(&mut parent, j));
                parent[x.max(y)] = x.min(y);
            }
        }
    }
    (0..n).filter(|&i| root(&mut parent, i) == i).count()
}

/// A point inside one of the faces `faces`, about `step` from their vertex at `p`.
fn point_near(b: &Brep, faces: &[FaceId], p: Point3, step: f64) -> Option<Point3> {
    for &f in faces {
        let dom = crate::boolean::face_domain(b, f);
        let uv0 = dom.surface.project(p);
        let (_, su, sv) = dom.surface.d1(uv0.x, uv0.y);
        let speed = su.norm().max(sv.norm());
        if speed <= 0.0 {
            continue;
        }
        for k in 0..16 {
            let a = TAU * (k as f64 + 0.5) / 16.0;
            let uv = uv0 + Vec2::new(a.cos(), a.sin()) * (step / speed);
            let q = dom.surface.eval(uv.x, uv.y);
            if let crate::boolean::domain::Where::Inside(_) = dom.locate(q) {
                return Some(q);
            }
        }
    }
    None
}

/// Whether `q`, near vertex `v`, lies on the material side of the fan of `faces` there: behind
/// the nearest of the faces, or, nearest an edge of theirs at `v`, behind the faces either side.
fn material_side(b: &Brep, faces: &[FaceId], v: VertexId, q: Point3) -> Option<bool> {
    use crate::boolean::domain::Where;
    let mut best: Option<(f64, Point3, Vec3)> = None;
    fn keep(best: &mut Option<(f64, Point3, Vec3)>, d: f64, r: Point3, n: Vec3) {
        if best.is_none_or(|x| d < x.0) {
            *best = Some((d, r, n));
        }
    }
    for &f in faces {
        let dom = crate::boolean::face_domain(b, f);
        let uv = dom.surface.project(q);
        let r = dom.surface.eval(uv.x, uv.y);
        if dom.locate(r) != Where::Outside {
            keep(&mut best, r.dist(q), r, dom.outward(uv));
        }
    }
    if best.is_none() {
        let outward = |g: FaceId, r: Point3| {
            let face = b.face(g);
            let uv = face.surface.project(r);
            let n = face.surface.normal(uv.x, uv.y);
            if face.reversed { -n } else { n }
        };
        for &f in faces {
            for &lp in &b.face(f).loops {
                for &c in &b.lp(lp).coedges {
                    let id = b.coedge(c).edge;
                    let e = b.edge(id);
                    let Some(curve) = &e.curve else { continue };
                    if e.start != v && e.end != v {
                        continue;
                    }
                    let r = curve.eval(curve.project(q, e.range.0, e.range.1));
                    let n = b
                        .edge_uses(id)
                        .iter()
                        .map(|&d| b.coedge_face(d))
                        .filter(|g| faces.contains(g))
                        .fold(Vec3::ZERO, |a, g| a + outward(g, r));
                    keep(&mut best, r.dist(q), r, n);
                }
            }
        }
    }
    let (_, r, n) = best?;
    let d = (q - r).dot(n);
    (d.abs() > 1e-3 * crate::tol::LINEAR).then_some(d < 0.0)
}

pub fn betti(b: &Brep) -> [u64; 3] {
    use std::collections::{HashMap, HashSet};
    let mut verts: HashSet<VertexId> = HashSet::new();
    let (mut edges, mut pulled) = (0i64, 0i64);
    for e in b.edge_ids() {
        let edge = b.edge(e);
        let us = b.edge_uses(e);
        if us.is_empty() {
            continue;
        }
        verts.insert(edge.start);
        verts.insert(edge.end);
        if edge.is_degenerate() {
            continue;
        }
        edges += 1;
        pulled += (us.len() as i64 / 2 - 1).max(0);
    }
    let faces: i64 = b.face_ids().map(|f| 2 - b.face(f).loops.len() as i64).sum();
    let chi_b = verts.len() as i64 - edges + faces;
    // Resolving each vertex where fans meet into the manifold boundary of the interior replaces
    // it by its link, the material on a small sphere round it. Where the vertex ends an edge
    // pieces share, that is one piece of link per fan, and the edge's own term splits them.
    let on_radial: HashSet<VertexId> = b
        .radial_edges()
        .flat_map(|e| [b.edge(e).start, b.edge(e).end])
        .collect();
    let scale = 1e-3 * bbox(b).diagonal();
    let fans: i64 = crate::check::fans_at_vertices(b)
        .into_iter()
        .filter(|(_, fs)| fs.len() > 1)
        .map(|(v, fs)| {
            let f = fs.len() as i64;
            if on_radial.contains(&v) {
                return f - 1;
            }
            let faces: Vec<Vec<FaceId>> = fs
                .iter()
                .map(|fan| fan.iter().map(|&c| b.coedge_face(c)).collect())
                .collect();
            2 * material_regions(b, v, &faces, scale) as i64 - f - 1
        })
        .sum();
    let chi_m = chi_b + fans - pulled;
    let chi = chi_b - chi_m / 2;
    // Faces' outer sides joined across the empty wedges at each edge.
    let ids: Vec<FaceId> = b.face_ids().collect();
    let at: HashMap<FaceId, usize> = ids.iter().enumerate().map(|(i, &f)| (f, i)).collect();
    let mut parent: Vec<usize> = (0..ids.len()).collect();
    fn root(p: &mut [usize], mut i: usize) -> usize {
        while p[i] != i {
            p[i] = p[p[i]];
            i = p[i];
        }
        i
    }
    for c in b.coedge_ids() {
        if let Some(o) = b.beside(c) {
            let (x, y) = (
                root(&mut parent, at[&b.coedge_face(c)]),
                root(&mut parent, at[&b.coedge_face(o)]),
            );
            if x != y {
                parent[x.max(y)] = x.min(y);
            }
        }
    }
    let mut groups: HashMap<usize, Vec<FaceId>> = HashMap::new();
    for (i, &f) in ids.iter().enumerate() {
        groups.entry(root(&mut parent, i)).or_default().push(f);
    }
    let voids = groups
        .values()
        .filter(|fs| enclosed_volume_sign(b, fs) < 0.0)
        .count() as u64;
    let pieces = b.solid_ids().count() as u64;
    let tunnels = pieces as i64 + voids as i64 - chi;
    [pieces, u64::try_from(tunnels).unwrap_or(0), voids]
}

/// The face's loops as parameter-plane polygons, for point classification.
pub(crate) fn uv_polygons(b: &Brep, f: FaceId) -> Vec<Vec<Vec2>> {
    uv_polygons_at(b, f, 64)
}

fn uv_polygons_at(b: &Brep, f: FaceId, n: usize) -> Vec<Vec<Vec2>> {
    b.face(f)
        .loops
        .iter()
        .map(|&lp| {
            let mut pts = Vec::new();
            for &co in &b.lp(lp).coedges {
                let uv = b.coedge_uvs(co, &coedge_params(b, co, n));
                pts.extend_from_slice(&uv[..uv.len() - 1]);
            }
            pts
        })
        .collect()
}

/// How many times the loops turn round `p`: the edges that cross `p`'s level to its right,
/// counted up for one rising and down for one falling. An edge owns its lower end and not its
/// upper, so a vertex on the level is counted once.
fn winding(polys: &[Vec<Vec2>], p: Vec2) -> i32 {
    let mut turns = 0;
    for poly in polys {
        for (k, &from) in poly.iter().enumerate() {
            let to = poly[(k + 1) % poly.len()];
            let rising = from.y <= p.y && p.y < to.y;
            let falling = to.y <= p.y && p.y < from.y;
            // `p` lies left of a rising edge, or right of a falling one, when the edge passes
            // on its right.
            let left = (to - from).cross(p - from);
            if rising && left > 0.0 {
                turns += 1;
            } else if falling && left < 0.0 {
                turns -= 1;
            }
        }
    }
    turns
}

/// A face's outline, for telling whether a point of its surface lies on the face.
pub(crate) struct Outline {
    /// Its loops as polygons in its surface's parameters.
    polys: Vec<Vec<Vec2>>,
    /// How far those polygons stray from its edges.
    sag: f64,
    /// Points along each edge, and how far the chords between them stray from the edge.
    rim: Vec<Vec<Point3>>,
    rim_sag: f64,
}

impl Outline {
    pub(crate) fn of(b: &Brep, f: FaceId) -> Outline {
        let mut rim = Vec::new();
        let mut rim_sag: f64 = 0.0;
        for &lp in &b.face(f).loops {
            for &co in &b.lp(lp).coedges {
                let Some(c) = &b.edge(b.coedge(co).edge).curve else {
                    continue;
                };
                let ts = coedge_params(b, co, 64);
                let pts: Vec<Point3> = ts.iter().map(|&t| c.eval(t)).collect();
                for k in 0..ts.len() - 1 {
                    let m = c.eval(0.5 * (ts[k] + ts[k + 1]));
                    rim_sag = rim_sag.max(chord_distance(m, pts[k], pts[k + 1]));
                }
                rim.push(pts);
            }
        }
        Outline {
            polys: uv_polygons(b, f),
            sag: outline_sag(b, f),
            rim,
            rim_sag,
        }
    }

    /// Where `p`, the point at `uv` of face `f`'s surface, lies on the face: its parameters
    /// there. The polygon decides, except within its chord sag of the boundary, where a polygon
    /// cuts inside a curved edge: there the nearest edge decides by its own side, and a point
    /// within `within` of an edge is on each face it bounds. The edges are asked only where `p`
    /// is that near the chords along them.
    pub(crate) fn holds(
        &self,
        b: &Brep,
        f: FaceId,
        uv: Vec2,
        p: Point3,
        within: f64,
    ) -> Option<Vec2> {
        let polygon = in_face(&b.face(f).surface, &self.polys, uv);
        let near = within.max(2.0 * self.sag) + self.rim_sag;
        let by_rim = self
            .rim
            .iter()
            .flat_map(|pts| pts.windows(2))
            .any(|w| chord_distance(p, w[0], w[1]) <= near);
        if !by_rim {
            return polygon;
        }
        match boundary_side(b, f, p) {
            Some((d, _)) if d <= within => Some(polygon.unwrap_or(uv)),
            Some((d, Some(side))) if d <= 2.0 * self.sag => side.then_some(polygon.unwrap_or(uv)),
            _ => polygon,
        }
    }
}

/// The distance from `p` to the chord from `a` to `b`.
fn chord_distance(p: Point3, a: Point3, b: Point3) -> f64 {
    let ab = b - a;
    let l2 = ab.dot(ab);
    let t = if l2 > 0.0 {
        ((p - a).dot(ab) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    p.dist(a + ab * t)
}

/// How far, in mm, the outline polygons of [`uv_polygons`] stray from face `f`'s edges: a point
/// nearer its boundary than this may lie between a chord and its arc, where a polygon's answer is
/// not the face's.
pub(crate) fn outline_sag(b: &Brep, f: FaceId) -> f64 {
    let s = &b.face(f).surface;
    let mut sag: f64 = 0.0;
    for &lp in &b.face(f).loops {
        for &co in &b.lp(lp).coedges {
            let Some(c) = &b.edge(b.coedge(co).edge).curve else {
                continue;
            };
            let ts = coedge_params(b, co, 64);
            let uv = b.coedge_uvs(co, &ts);
            for k in 0..ts.len() - 1 {
                let m = (uv[k] + uv[k + 1]) * 0.5;
                sag = sag.max(s.eval(m.x, m.y).dist(c.eval(0.5 * (ts[k] + ts[k + 1]))));
            }
        }
    }
    sag
}

/// The distance from `p`, a point on face `f`'s surface, to the face's nearest edge, and which
/// side of that edge it lies on by the edge itself: a loop has the face's material on its left
/// seen from outside, so `p` is on the face where the outward normal crossed with the coedge's
/// direction points towards it. `None` for the side where the nearest point is an edge's end,
/// whose two edges need not agree.
pub(crate) fn boundary_side(b: &Brep, f: FaceId, p: Point3) -> Option<(f64, Option<bool>)> {
    let face = b.face(f);
    let uv = face.surface.project(p);
    let n = face.surface.normal(uv.x, uv.y);
    let n = if face.reversed { -n } else { n };
    let mut best: Option<(f64, Option<bool>)> = None;
    for &lp in &face.loops {
        for &co in &b.lp(lp).coedges {
            let coedge = b.coedge(co);
            let edge = b.edge(coedge.edge);
            let Some(c) = &edge.curve else {
                continue;
            };
            let (t0, t1) = edge.range;
            let t = c.project(p, t0, t1);
            let (q, d1) = c.d1(t);
            let dist = q.dist(p);
            if best.is_some_and(|(d, _)| d <= dist) {
                continue;
            }
            let end = (t - t0).abs() <= 1e-9 * (1.0 + t0.abs())
                || (t1 - t).abs() <= 1e-9 * (1.0 + t1.abs());
            let along = if coedge.reversed { -d1 } else { d1 };
            // A seam, or a slit the face runs along and back, has the face on both sides.
            let twice = b
                .edge_uses(coedge.edge)
                .iter()
                .filter(|&&u| b.coedge_face(u) == f)
                .count()
                > 1;
            let side = (!end).then(|| twice || n.cross(along).dot(p - q) > 0.0);
            best = Some((dist, side));
        }
    }
    best
}

/// Whether `uv`, or a whole period away from it, lies inside the face.
/// The edge nearest `p`, and how far: of edges equally near, one between two faces before a seam,
/// whose two sides are one face, as where a seam's end meets a rim.
pub fn edge_at(b: &Brep, p: Point3) -> Option<(crate::topo::EdgeId, f64)> {
    let seam = |e: crate::topo::EdgeId| {
        let fs: Vec<FaceId> = b
            .coedge_ids()
            .filter(|&c| b.coedge(c).edge == e)
            .map(|c| b.coedge_face(c))
            .collect();
        fs.len() == 2 && fs[0] == fs[1]
    };
    b.edge_ids()
        .filter_map(|e| {
            let x = b.edge(e);
            let c = x.curve.as_ref()?;
            let (t0, t1) = (x.range.0.min(x.range.1), x.range.0.max(x.range.1));
            Some((e, c.eval(c.project(p, t0, t1)).dist(p)))
        })
        .min_by(|x, y| {
            let key = |z: &(crate::topo::EdgeId, f64)| (z.1 > 10.0 * crate::tol::LINEAR, seam(z.0));
            key(x).cmp(&key(y)).then(x.1.total_cmp(&y.1))
        })
}

/// The outward normal of the boundary at `p`, a point on it: every face holding `p`, within
/// `within`, gives its normal, and where they agree (a face's inside, or a smooth seam such as a
/// fillet's) that is it; at a crease, where they do not, there is none.
pub fn normal_at(b: &Brep, p: Point3, within: f64) -> Option<Vec3> {
    Outlines::new(b).normal_at(p, within)
}

/// A solid's faces, each with its outline in its surface's parameters: drawn once, to ask about
/// many points of its boundary.
pub struct Outlines<'a> {
    b: &'a Brep,
    faces: Vec<(FaceId, Outline)>,
}

impl<'a> Outlines<'a> {
    pub fn new(b: &'a Brep) -> Outlines<'a> {
        Outlines {
            b,
            faces: b.face_ids().map(|f| (f, Outline::of(b, f))).collect(),
        }
    }

    /// Where `p`, a point of the surface of the face `held` names, lies on that face.
    fn on_face(&self, (f, outline): &(FaceId, Outline), p: Point3, within: f64) -> Option<Vec2> {
        let uv = self.b.face(*f).surface.project(p);
        outline.holds(self.b, *f, uv, p, within)
    }

    /// As [`normal_at`].
    pub fn normal_at(&self, p: Point3, within: f64) -> Option<Vec3> {
        let b = self.b;
        let found: Vec<Vec3> = self
            .faces
            .iter()
            .filter_map(|held| {
                let face = b.face(held.0);
                let s = &face.surface;
                let near = s.project(p);
                if s.eval(near.x, near.y).dist(p) > within {
                    return None;
                }
                let uv = self.on_face(held, p, within)?;
                let n = s.normal(uv.x, uv.y);
                Some(if face.reversed { -n } else { n })
            })
            .collect();
        let first = *found.first()?;
        found
            .iter()
            .all(|n| (*n - first).norm() <= 1e-6)
            .then_some(first)
    }

    /// As [`curvature`].
    pub fn curvature(&self, p: Point3) -> Option<[f64; 2]> {
        let b = self.b;
        // Each surface's point nearest `p`, nearest first: the first its face holds is the
        // nearest point of the faces.
        let mut near: Vec<(f64, usize, Point3)> = self
            .faces
            .iter()
            .enumerate()
            .map(|(k, held)| {
                let s = &b.face(held.0).surface;
                let uv = s.project(p);
                let q = s.eval(uv.x, uv.y);
                (q.dist(p), k, q)
            })
            .collect();
        near.sort_by(|x, y| x.0.total_cmp(&y.0));
        let face = near.into_iter().find_map(|(d, k, q)| {
            let held = &self.faces[k];
            self.on_face(held, q, 0.0).map(|uv| (d, held.0, uv))
        });
        let edge = b
            .edge_ids()
            .filter_map(|e| {
                let edge = b.edge(e);
                let c = edge.curve.as_ref()?;
                let (t0, t1) = (
                    edge.range.0.min(edge.range.1),
                    edge.range.0.max(edge.range.1),
                );
                let q = c.eval(c.project(p, t0, t1));
                Some((q.dist(p), e, q))
            })
            .min_by(|x, y| x.0.total_cmp(&y.0));
        let (f, uv) = match (face, edge) {
            (Some((df, f, uv)), Some((de, ..))) if df <= de + crate::tol::LINEAR => (f, uv),
            (_, Some((_, e, q))) => {
                // On an edge: a crease unless its faces share their tangent plane there, as a
                // fillet and the face it meets do, where either face's curvature answers.
                let faces: Vec<FaceId> = b
                    .coedge_ids()
                    .filter(|&c| b.coedge(c).edge == e)
                    .map(|c| b.coedge_face(c))
                    .collect();
                let normal = |f: FaceId| {
                    let s = &b.face(f).surface;
                    let uv = s.project(q);
                    let n = s.normal(uv.x, uv.y);
                    (if b.face(f).reversed { -n } else { n }, uv)
                };
                let (&f0, &f1) = (faces.first()?, faces.get(1)?);
                let ((n0, uv0), (n1, _)) = (normal(f0), normal(f1));
                if n0.cross(n1).norm() > 1e-9 || n0.dot(n1) < 0.0 {
                    return None;
                }
                (f0, uv0)
            }
            (Some((_, f, uv)), None) => (f, uv),
            (None, None) => return None,
        };
        let face = b.face(f);
        principal(&face.surface, uv, face.reversed)
    }
}

/// The face nearest `p` among those whose own region, not just their surface, holds the point of
/// their surface nearest `p`; with that distance.
pub fn face_at(b: &Brep, p: Point3) -> Option<(FaceId, f64)> {
    use crate::boolean::{domain::Where, face_domain};
    b.face_ids()
        .filter_map(|f| {
            let s = &b.face(f).surface;
            let uv = s.project(p);
            let foot = s.eval(uv.x, uv.y);
            // By the face's exact domain: a polygon of its outline cuts inside a curved edge, and
            // takes a point of a hole's opening for the face round it.
            let on = !matches!(face_domain(b, f).locate(foot), Where::Outside);
            on.then(|| (f, foot.dist(p)))
        })
        .min_by(|x, y| x.1.total_cmp(&y.1))
}

/// The principal curvatures `[k1, k2]`, `k1 ≥ k2`, 1/mm, of the boundary at its point nearest
/// `p`: positive where the boundary bends away from its outward normal, as on a ball. `None` where
/// that point lies on a crease, an edge whose faces meet at an angle there.
pub fn curvature(b: &Brep, p: Point3) -> Option<[f64; 2]> {
    Outlines::new(b).curvature(p)
}

/// The principal curvatures of `s` at `uv` against its outward normal, `k1 ≥ k2`: the shape
/// operator from the first and second fundamental forms, the second derivatives by central
/// differences of the exact first ones. Where the parameters degenerate (a ball's pole), a point
/// a hair away answers for it.
fn principal(s: &Surface, uv: Vec2, reversed: bool) -> Option<[f64; 2]> {
    let at = |u: f64, v: f64| {
        let (_, su, sv) = s.d1(u, v);
        let n = su.cross(sv).unit()?;
        let n = {
            let natural = s.normal(u, v);
            let n = if n.dot(natural) < 0.0 { -n } else { n };
            if reversed { -n } else { n }
        };
        let (hu, hv) = (1e-5 * (1.0 + u.abs()), 1e-5 * (1.0 + v.abs()));
        let (_, su1, sv1) = s.d1(u + hu, v);
        let (_, su0, sv0) = s.d1(u - hu, v);
        let (_, su3, sv3) = s.d1(u, v + hv);
        let (_, su2, sv2) = s.d1(u, v - hv);
        let suu = (su1 - su0) * (0.5 / hu);
        let svv = (sv3 - sv2) * (0.5 / hv);
        let suv = ((sv1 - sv0) * (0.5 / hu) + (su3 - su2) * (0.5 / hv)) * 0.5;
        let (e, f, g) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let det = e * g - f * f;
        if det <= 1e-24 * (e * g).max(1e-300) {
            return None;
        }
        let (l, m, nn) = (suu.dot(n), suv.dot(n), svv.dot(n));
        // Bending toward the outward normal is negative.
        let h = -(e * nn + g * l - 2.0 * f * m) / (2.0 * det);
        let k = (l * nn - m * m) / det;
        let r = (h * h - k).max(0.0).sqrt();
        Some([h + r, h - r])
    };
    at(uv.x, uv.y).or_else(|| {
        let nudge = 1e-7 * (1.0 + uv.y.abs());
        at(uv.x, uv.y + nudge).or_else(|| at(uv.x, uv.y - nudge))
    })
}

pub(crate) fn in_face(s: &Surface, polys: &[Vec<Vec2>], uv: Vec2) -> Option<Vec2> {
    // Every whole-period shift that lands within the outline's own extent: a face can wind round
    // its surface several times (a bore's wall left between a thread's turns), its outline
    // running on across the parameter plane.
    let (mut lo, mut hi) = (
        Vec2::new(f64::INFINITY, f64::INFINITY),
        Vec2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for q in polys.iter().flatten() {
        lo = Vec2::new(lo.x.min(q.x), lo.y.min(q.y));
        hi = Vec2::new(hi.x.max(q.x), hi.y.max(q.y));
    }
    let shifts = |p: Option<f64>, x: f64, lo: f64, hi: f64| match p {
        Some(t) if lo.is_finite() => {
            let (k0, k1) = (((lo - x) / t).floor() as i64, ((hi - x) / t).ceil() as i64);
            let mut out: Vec<f64> = (k0..=k1).map(|k| k as f64 * t).collect();
            out.sort_by(|a, b| a.abs().total_cmp(&b.abs()));
            out
        }
        _ => vec![0.0],
    };
    for du in shifts(s.u_period(), uv.x, lo.x, hi.x) {
        for dv in shifts(s.v_period(), uv.y, lo.y, hi.y) {
            let q = Vec2::new(uv.x + du, uv.y + dv);
            if winding(polys, q) != 0 {
                return Some(q);
            }
        }
    }
    None
}

/// Points of the face where its coordinate along `e` is stationary and that a boundary need not
/// reach: the critical sets of planes, cylinders, cones and extrusions are lines that always meet
/// the boundary or a seam, so only the remaining surfaces have candidates.
fn stationary(s: &Surface, polys: &[Vec<Vec2>], e: Vec3, starts: Option<Vec<Vec2>>) -> Vec<Vec2> {
    let mut cands = Vec::new();
    match s {
        Surface::Sphere(sp) => {
            for d in [e, -e] {
                cands.push(s.project(sp.frame.origin + d * sp.radius));
            }
        }
        Surface::Torus(t) => {
            let l = t.frame.dir_to_local(e);
            let perp = l.x.hypot(l.y);
            if perp > 1e-12 {
                let u0 = l.y.atan2(l.x);
                for (u, r) in [(u0, perp), (u0 + std::f64::consts::PI, -perp)] {
                    let v = l.z.atan2(r);
                    cands.push(Vec2::new(u, v));
                    cands.push(Vec2::new(u, v + std::f64::consts::PI));
                }
            }
        }
        Surface::Revolution(_) | Surface::Ruled(_) | Surface::Nurbs(_) => {
            let (lo, hi) = uv_bounds(polys);
            let n = 12;
            let starts = starts.unwrap_or_else(|| {
                (0..=n)
                    .flat_map(|i| {
                        (0..=n).map(move |j| {
                            Vec2::new(
                                lo.x + (hi.x - lo.x) * i as f64 / n as f64,
                                lo.y + (hi.y - lo.y) * j as f64 / n as f64,
                            )
                        })
                    })
                    .collect()
            });
            for start in starts {
                let mut uv = start;
                let grad = |uv: Vec2| {
                    let (_, su, sv) = s.d1(uv.x, uv.y);
                    Vec2::new(su.dot(e), sv.dot(e))
                };
                let mut ok = false;
                for _ in 0..40 {
                    let g = grad(uv);
                    let h = 1e-6;
                    let gu = (grad(Vec2::new(uv.x + h, uv.y)) - grad(Vec2::new(uv.x - h, uv.y)))
                        / (2.0 * h);
                    let gv = (grad(Vec2::new(uv.x, uv.y + h)) - grad(Vec2::new(uv.x, uv.y - h)))
                        / (2.0 * h);
                    let Some((du, dv)) = solve2(gu.x, gv.x, gu.y, gv.y, -g.x, -g.y) else {
                        break;
                    };
                    uv = Vec2::new(uv.x + du, uv.y + dv);
                    if du.abs() + dv.abs() < 1e-13 {
                        ok = true;
                        break;
                    }
                }
                if ok {
                    cands.push(uv);
                }
            }
        }
        _ => {}
    }
    cands
        .into_iter()
        .filter_map(|uv| in_face(s, polys, uv))
        .collect()
}

/// Where to look for a B-spline face's stationary points along axis `i`: a few points in each
/// patch that the face's parameter box reaches and whose box pokes out of `bb` on that axis.
fn nurbs_starts(n: &NurbsSurface, polys: &[Vec<Vec2>], bb: &Aabb, i: usize) -> Vec<Vec2> {
    let (lo, hi) = uv_bounds(polys);
    let mut out = Vec::new();
    for q in reached_patches(n, polys) {
        if q.bbox.min[i] >= bb.min[i] && q.bbox.max[i] <= bb.max[i] {
            continue;
        }
        let (u0, u1, v0, v1) = if n.u_period().is_none() && n.v_period().is_none() {
            (
                q.lo.x.max(lo.x),
                q.hi.x.min(hi.x),
                q.lo.y.max(lo.y),
                q.hi.y.min(hi.y),
            )
        } else {
            (q.lo.x, q.hi.x, q.lo.y, q.hi.y)
        };
        for a in 0..3 {
            for c in 0..3 {
                out.push(Vec2::new(
                    u0 + (u1 - u0) * (a as f64 + 0.5) / 3.0,
                    v0 + (v1 - v0) * (c as f64 + 0.5) / 3.0,
                ));
            }
        }
    }
    out
}

fn uv_bounds(polys: &[Vec<Vec2>]) -> (Vec2, Vec2) {
    let (mut lo, mut hi) = (
        Vec2::new(f64::INFINITY, f64::INFINITY),
        Vec2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for p in polys.iter().flatten() {
        lo = Vec2::new(lo.x.min(p.x), lo.y.min(p.y));
        hi = Vec2::new(hi.x.max(p.x), hi.y.max(p.y));
    }
    (lo, hi)
}

/// The patches of `n` that the parameter box of a face's loops `polys` reaches. A periodic
/// parameter's face may sit a period away from the knot domain; every patch along it counts.
fn reached_patches<'a>(
    n: &'a NurbsSurface,
    polys: &[Vec<Vec2>],
) -> impl Iterator<Item = &'a Patch> + 'a {
    let (lo, hi) = uv_bounds(polys);
    patches_in(n, lo, hi)
}

/// The patches of `n` that the parameter box `lo`–`hi` reaches.
fn patches_in<'a>(n: &'a NurbsSurface, lo: Vec2, hi: Vec2) -> impl Iterator<Item = &'a Patch> + 'a {
    let (pu, pv) = (n.u_period().is_some(), n.v_period().is_some());
    n.patches().iter().filter(move |q| {
        (pu || (q.lo.x <= hi.x && lo.x <= q.hi.x)) && (pv || (q.lo.y <= hi.y && lo.y <= q.hi.y))
    })
}

/// A box that holds every face in the arena, cheaper than [`bbox`] and at most a little larger:
/// a B-spline face counts with the boxes of the patches its parameter box reaches.
pub fn loose_bbox(b: &Brep) -> Aabb {
    tight_or_loose_bbox(b, false)
}

/// A box holding face `f`: the smallest for every surface but a B-spline, where it may be a
/// little larger.
pub fn face_bbox(b: &Brep, f: FaceId) -> Aabb {
    face_box(b, f)
}

/// One face's share of [`loose_bbox`], from its memo.
pub(crate) fn face_box(b: &Brep, f: FaceId) -> Aabb {
    *b.memo(f).loose_box.get_or_init(|| face_loose_box(b, f))
}

/// One face's share of [`loose_bbox`]: its edges' extremes, and for a B-spline face the boxes of
/// the patches it reaches, for any other curved surface its stationary points.
fn face_loose_box(b: &Brep, f: FaceId) -> Aabb {
    let mut bb = Aabb::EMPTY;
    let axes = [Vec3::X, Vec3::Y, Vec3::Z];
    for &lp in &b.face(f).loops {
        for &co in &b.lp(lp).coedges {
            let edge = b.edge(b.coedge(co).edge);
            bb.add(b.vertex(edge.start).point);
            bb.add(b.vertex(edge.end).point);
            if let Some(c) = &edge.curve {
                for ax in axes {
                    for t in c.extremes(ax, edge.range.0, edge.range.1) {
                        bb.add(c.eval(t));
                    }
                }
            }
        }
    }
    let s = &b.face(f).surface;
    match s {
        Surface::Nurbs(n) => {
            // A sparse outline, its parameter box widened by a quarter each way, still holds
            // the face's parameters: a pcurve does not bulge that far between samples.
            let (lo, hi) = uv_bounds(&uv_polygons_at(b, f, 8));
            let pad = (hi - lo) * 0.25;
            for q in patches_in(n, lo - pad, hi + pad) {
                bb.union(&q.bbox);
            }
        }
        Surface::Sphere(_) | Surface::Torus(_) | Surface::Revolution(_) | Surface::Ruled(_) => {
            let polys = uv_polygons(b, f);
            for ax in axes {
                for uv in stationary(s, &polys, ax, None) {
                    bb.add(s.eval(uv.x, uv.y));
                }
            }
        }
        _ => {}
    }
    bb
}

/// A box that holds [`face_box`]'s, from control points and closed forms alone: a spline edge or
/// face counts with its control points' box, a sphere or torus with its whole extent. It costs a
/// pass over the face's edges, where [`face_box`] finds extremes, so it rules out faces far from
/// something before their [`face_box`] is worked out.
pub(crate) fn face_hull_box(b: &Brep, f: FaceId) -> Aabb {
    if let Some(bb) = b.memo(f).loose_box.get() {
        return *bb;
    }
    let mut bb = Aabb::EMPTY;
    let axes = [Vec3::X, Vec3::Y, Vec3::Z];
    for &lp in &b.face(f).loops {
        for &co in &b.lp(lp).coedges {
            let edge = b.edge(b.coedge(co).edge);
            bb.add(b.vertex(edge.start).point);
            bb.add(b.vertex(edge.end).point);
            match &edge.curve {
                None | Some(crate::geom::Curve::Line(_)) => {}
                Some(crate::geom::Curve::BSpline(s)) => {
                    s.control_points().iter().for_each(|&p| bb.add(p))
                }
                Some(crate::geom::Curve::Nurbs(n)) => {
                    n.control_points().iter().for_each(|&p| bb.add(p))
                }
                Some(c) => {
                    for ax in axes {
                        for t in c.extremes(ax, edge.range.0, edge.range.1) {
                            bb.add(c.eval(t));
                        }
                    }
                }
            }
        }
    }
    match &b.face(f).surface {
        Surface::Nurbs(n) => bb.union(&n.hull_box()),
        Surface::Sphere(s) => {
            let r = Vec3::new(s.radius, s.radius, s.radius);
            bb.add(s.frame.origin - r);
            bb.add(s.frame.origin + r);
        }
        Surface::Torus(t) => {
            let r = t.major + t.minor;
            let r = Vec3::new(r, r, r);
            bb.add(t.frame.origin - r);
            bb.add(t.frame.origin + r);
        }
        Surface::Revolution(_) | Surface::Ruled(_) => bb.union(&face_box(b, f)),
        _ => {}
    }
    // Knot insertion and the rational division round; the box gives that rounding room.
    let pad = 1e-9 * (1.0 + bb.diagonal());
    Aabb {
        min: bb.min - Vec3::new(pad, pad, pad),
        max: bb.max + Vec3::new(pad, pad, pad),
    }
}

/// The tight axis-aligned box of every face in the arena.
pub fn bbox(b: &Brep) -> Aabb {
    tight_or_loose_bbox(b, true)
}

fn tight_or_loose_bbox(b: &Brep, tight: bool) -> Aabb {
    if !tight {
        // Every edge of a solid bounds a face, so the faces' own loose boxes cover the arena.
        let mut bb = Aabb::EMPTY;
        for f in b.face_ids() {
            bb.union(&face_box(b, f));
        }
        for e in b.edge_ids() {
            let edge = b.edge(e);
            bb.add(b.vertex(edge.start).point);
            bb.add(b.vertex(edge.end).point);
        }
        return bb;
    }
    let mut bb = Aabb::EMPTY;
    let axes = [Vec3::X, Vec3::Y, Vec3::Z];
    for e in b.edge_ids() {
        let edge = b.edge(e);
        bb.add(b.vertex(edge.start).point);
        bb.add(b.vertex(edge.end).point);
        if let Some(c) = &edge.curve {
            for ax in axes {
                for t in c.extremes(ax, edge.range.0, edge.range.1) {
                    bb.add(c.eval(t));
                }
            }
        }
    }
    for f in b.face_ids() {
        let s = &b.face(f).surface;
        if matches!(
            s,
            Surface::Sphere(_)
                | Surface::Torus(_)
                | Surface::Revolution(_)
                | Surface::Ruled(_)
                | Surface::Nurbs(_)
        ) {
            let polys = uv_polygons(b, f);
            for (i, ax) in axes.into_iter().enumerate() {
                // A B-spline lies inside its patches' boxes: only a patch the face reaches into
                // that pokes out of the box on this axis can widen it, and only there is a
                // stationary point searched for.
                let starts = match s {
                    Surface::Nurbs(n) => {
                        let starts = nurbs_starts(n, &polys, &bb, i);
                        if starts.is_empty() {
                            continue;
                        }
                        Some(starts)
                    }
                    _ => None,
                };
                for uv in stationary(s, &polys, ax, starts) {
                    bb.add(s.eval(uv.x, uv.y));
                }
            }
        }
    }
    bb
}

/// The face's area centroid, by Green's theorem over its parameter-plane loops, each edge taken as
/// 32 points: close, not exact, and the same for the same face.
pub fn face_centroid(b: &Brep, f: FaceId) -> Point3 {
    let face = b.face(f);
    let s = &face.surface;
    let mut polys: Vec<Vec<Vec2>> = Vec::new();
    for &lp in &face.loops {
        let mut pts = Vec::new();
        for &c in &b.lp(lp).coedges {
            let uv = b.coedge_uvs(c, &coedge_params(b, c, 32));
            pts.extend_from_slice(&uv[..uv.len() - 1]);
        }
        polys.push(pts);
    }
    let Some(u_ref) = polys.iter().flatten().next().map(|p| p.x) else {
        return Vec3::ZERO;
    };
    // ∬ g du dv = ∮ G dv, with G(u, v) = ∫ g from u_ref to u.
    let (mut area, mut moment) = (0.0, Vec3::ZERO);
    for poly in &polys {
        for i in 0..poly.len() {
            let (p, q) = (poly[i], poly[(i + 1) % poly.len()]);
            let dv = q.y - p.y;
            if dv == 0.0 {
                continue;
            }
            for &(x, w) in gauss(8) {
                let t = 0.5 * (x + 1.0);
                let (u, v) = (p.x + (q.x - p.x) * t, p.y + dv * t);
                let h = 0.5 * (u - u_ref);
                for &(y, wy) in gauss(8) {
                    let su = u_ref + h * (y + 1.0);
                    let (pt, du, dvv) = s.d1(su, v);
                    let jac = du.cross(dvv).norm();
                    let k = 0.5 * w * dv * h * wy * jac;
                    area += k;
                    moment += pt * k;
                }
            }
        }
    }
    if area == 0.0 {
        return polys[0].first().map_or(Vec3::ZERO, |uv| s.eval(uv.x, uv.y));
    }
    moment * (1.0 / area)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        math::Frame,
        profile::{NamedSegment, Profile, ProfileLoop, Segment},
        sweep::extrude,
    };

    // x² + 1e-9·sin(1e7·x) on [0, 1]: no rule settles the ripple until its pieces are shorter
    // than its wavelength, 23 halvings down, so to a tolerance of 1e-14 its rules differ by about
    // as much at every depth. The integral is 1/3 + 1e-16·(1 − cos 1e7), and it is found to within
    // the ripple's own size after a few halvings, not the millions to depth 20.
    #[test]
    fn an_integrand_noisier_than_its_tolerance_is_integrated_to_its_noise_in_few_steps() {
        let calls = std::cell::Cell::new(0usize);
        let f = |x: f64| {
            calls.set(calls.get() + 1);
            [x * x, x * x + 1e-9 * (1e7 * x).sin()]
        };
        let got = integrate_pieces(&f, &[0.0, 1.0], [1e-14, 1e-14], EXACT);
        let want = 1.0 / 3.0 + 1e-16 * (1.0 - 1e7f64.cos());
        assert!((got[0] - 1.0 / 3.0).abs() <= 1e-15, "{}", got[0]);
        assert!((got[1] - want).abs() <= 1e-9, "{}", got[1]);
        assert!(calls.get() <= 10_000, "{} calls", calls.get());
    }

    #[test]
    fn the_integral_round_an_axis_in_closed_form_is_the_rule_s_and_a_zone_s_by_hand() {
        use crate::geom::{Cone, Sphere, Torus};
        // A frame off the origin and tilted, the volume taken about a point off its axis.
        let frame = Frame::new(
            Vec3::new(3.0, -2.0, 5.0),
            Vec3::new(1.0, 2.0, 2.0),
            Vec3::new(2.0, 1.0, -2.0),
        )
        .unwrap();
        let c = Vec3::new(-7.0, 4.0, 1.5);
        let surfaces = [
            Surface::Cone(Cone {
                frame,
                radius: 4.0,
                half_angle: 0.3,
            }),
            Surface::Sphere(Sphere { frame, radius: 9.0 }),
            Surface::Torus(Torus {
                frame,
                major: 11.0,
                minor: 3.0,
            }),
        ];
        let density = |p: Point3, n: Vec3| [n.norm(), (p - c).dot(n) / 3.0];
        for s in &surfaces {
            for (v, lo, to) in [(0.4, -1.0, 2.5), (-1.1, 0.3, -4.0), (2.0, 0.0, TAU)] {
                let closed = round_axis(s, v, (lo, to), c).unwrap();
                let rule = inner_rule(s, None, false, 24, (lo, v, to), &density);
                for k in 0..2 {
                    let scale = rule[k].abs().max(1.0);
                    assert!(
                        (closed[k] - rule[k]).abs() <= 1e-13 * scale,
                        "{} at {v}: {} vs the rule's {}",
                        s.kind(),
                        closed[k],
                        rule[k]
                    );
                }
            }
        }
        // By hand: the zone of a sphere of radius 9 at the latitude v is a circle 2π·9·cos v
        // round, and `n` there is 81·cos v long, pointing out, so about the sphere's centre
        // the volume's share is a third of 9·(2π·81·cos v).
        let v = 0.4f64;
        let zone = round_axis(&surfaces[1], v, (0.0, TAU), frame.origin).unwrap();
        assert!((zone[0] - TAU * 81.0 * v.cos()).abs() <= 1e-12);
        assert!((zone[1] - TAU * 243.0 * v.cos()).abs() <= 1e-12);
    }

    #[test]
    fn a_face_whose_edges_run_through_a_cusp_of_its_parameter_lines_has_its_closed_form_area() {
        // S(u, v) = (v², v³, u) over [0, 1] × [−1, 1]: each line u = const is a cuspidal cubic,
        // and the two edges along such lines pass the cusp at v = 0, where an edge point near it
        // has a foot on either branch within the tolerance. The area is
        // ∫∫ |v| √(4 + 9v²) dv du = 2 (13√13 − 8) / 27.
        use crate::geom::{Curve, Line, NurbsCurve, NurbsSurface};
        use crate::topo::{CoedgeSpec, PCurve};
        let row = |u: f64| {
            [
                (1.0, -1.0),
                (-1.0 / 3.0, 1.0),
                (-1.0 / 3.0, -1.0),
                (1.0, 1.0),
            ]
            .map(|(x, y)| Point3::new(x, y, u))
        };
        let cubic = vec![-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0];
        let net: Vec<Point3> = row(0.0).into_iter().chain(row(1.0)).collect();
        let surface = NurbsSurface::new(
            (1, 3),
            (vec![0.0, 0.0, 1.0, 1.0], cubic.clone()),
            (2, 4),
            net,
            None,
        )
        .unwrap();
        let mut b = Brep::new();
        let corner = |u: f64, v: f64| Point3::new(v * v, v * v * v, u);
        let [a, c, d, e] = [(0.0, -1.0), (1.0, -1.0), (1.0, 1.0), (0.0, 1.0)]
            .map(|(u, v)| b.add_vertex(corner(u, v)));
        let line = |from: Point3, dir: Vec3| Curve::Line(Line { origin: from, dir });
        let cusp = |u: f64| {
            let c = NurbsCurve::new(3, cubic.clone(), row(u).to_vec(), None).unwrap();
            Curve::Nurbs(Box::new(c))
        };
        let e0 = b.add_edge(line(corner(0.0, -1.0), Vec3::Z), (0.0, 1.0), a, c);
        let e1 = b.add_edge(cusp(1.0), (-1.0, 1.0), c, d);
        let e2 = b.add_edge(line(corner(1.0, 1.0), -Vec3::Z), (0.0, 1.0), d, e);
        let e3 = b.add_edge(cusp(0.0), (-1.0, 1.0), a, e);
        let f = b.add_face(Surface::Nurbs(Box::new(surface)), false, []);
        let spec = |edge, reversed, pcurve| CoedgeSpec {
            edge,
            reversed,
            pcurve,
        };
        b.add_loop(
            f,
            vec![
                spec(
                    e0,
                    false,
                    PCurve::Affine {
                        origin: Vec2::new(0.0, -1.0),
                        dir: Vec2::new(1.0, 0.0),
                    },
                ),
                spec(
                    e1,
                    false,
                    PCurve::Projected {
                        seed: Vec2::new(1.0, -1.0),
                    },
                ),
                spec(
                    e2,
                    false,
                    PCurve::Affine {
                        origin: Vec2::new(1.0, 1.0),
                        dir: Vec2::new(-1.0, 0.0),
                    },
                ),
                spec(
                    e3,
                    true,
                    PCurve::Projected {
                        seed: Vec2::new(0.0, 1.0),
                    },
                ),
            ],
        );
        let (area, _) = face_integrals(&b, f, Point3::ZERO, 2.0, EXACT);
        let want = 2.0 * (13.0 * 13f64.sqrt() - 8.0) / 27.0;
        assert!((area - want).abs() < 1e-12, "{area} against {want}");
    }

    /// `k` prisms 5 tall side by side, each on the regular 48-gon of circumradius 3 drawn as a
    /// closed spline of degree one, so each of its edges has 48 spans.
    fn spline_prisms(k: usize) -> Brep {
        let mut arena = Brep::new();
        for i in 0..k {
            let corners: Vec<Vec3> = (0..48)
                .map(|j| {
                    let t = TAU * j as f64 / 48.0;
                    Vec3::new(3.0 * t.cos(), 3.0 * t.sin(), 0.0)
                })
                .collect();
            let curve = crate::geom::BSpline::periodic_uniform(1, &corners).unwrap();
            let side = Segment::Spline {
                curve,
                reversed: false,
            };
            let at = Frame::new(Vec3::new(7.0 * i as f64, 0.0, 0.0), Vec3::Z, Vec3::X).unwrap();
            let outline = Profile::outline(vec![NamedSegment::new(side, "side")]);
            arena.append(&extrude(&outline, &at, 5.0, "prism").unwrap());
        }
        arena
    }

    // Each face's integrals are its own and summed in the caller's order, so the faces measured
    // on many threads, heaviest first, give the same values to the bit as one after another. Two
    // arenas built alike keep the second from reading the first's memos.
    #[test]
    fn faces_measured_on_many_threads_match_those_measured_on_one_to_the_bit() {
        let (one, many) = (spline_prisms(12), spline_prisms(12));
        let faces: Vec<FaceId> = one.face_ids().collect();
        assert_eq!(faces, many.face_ids().collect::<Vec<_>>());
        let (c, scale) = moment_reference(&one);
        let order: Vec<(usize, usize)> = faces
            .iter()
            .enumerate()
            .map(|(k, &f)| (pieces(&one, f), k))
            .collect();
        // Enough pieces that the arena's own measures are spread too.
        let total: usize = order.iter().map(|w| w.0).sum();
        assert!(total >= 2 * PIECES_PER_THREAD, "{total} pieces");
        let bits = |b: &Brep, most: usize| -> Vec<[u64; 2]> {
            each_face_on(most, &faces, order.clone(), |f| {
                let (a, v) = face_exact(b, f, c, scale);
                [a.to_bits(), v.to_bits()]
            })
        };
        assert_eq!(bits(&one, 1), bits(&many, 8));
    }

    #[test]
    fn a_face_moved_keeps_its_box_moved() {
        // An assembly measures a part against another moved before each boolean: a slid part's
        // faces keep their boxes, slid, rather than finding each spline edge's extremes again.
        let rod = extrude(
            &Profile::circle(Vec2::ZERO, 20.0),
            &Frame::new(Vec3::ZERO, Vec3::Z, Vec3::X).unwrap(),
            20.0,
            "r",
        )
        .unwrap();
        let boxes: Vec<Aabb> = rod.face_ids().map(|f| face_box(&rod, f)).collect();
        let shift = Vec3::new(3.5, -120.0, 7.25);
        let mut slid = rod.clone();
        slid.transform(&crate::math::Transform::translation(shift));
        for (f, before) in slid.face_ids().zip(boxes) {
            let kept = slid.memo(f).loose_box.get().copied();
            assert_eq!(
                kept.map(|b| (b.min, b.max)),
                Some((before.min + shift, before.max + shift))
            );
        }
        // Turned, a box is found again: the turned box's own would be looser.
        let mut turned = rod.clone();
        turned.transform(&crate::math::Transform::rotation(Vec3::ZERO, Vec3::X, 0.3));
        assert!(
            turned
                .face_ids()
                .all(|f| turned.memo(f).loose_box.get().is_none())
        );
    }

    #[test]
    fn a_rods_end_just_inside_its_rim_is_flat() {
        // A rod of radius 10: its end's outline, drawn as chords, cuts up to 0.012 inside the
        // rim, and a point of the end between a chord and the rim is still a point of the end.
        let rod = extrude(
            &Profile::circle(Vec2::ZERO, 20.0),
            &Frame::new(Vec3::ZERO, Vec3::Z, Vec3::X).unwrap(),
            20.0,
            "r",
        )
        .unwrap();
        for k in 0..256 {
            let a = std::f64::consts::TAU * k as f64 / 256.0;
            let p = Point3::new(9.999 * a.cos(), 9.999 * a.sin(), 20.0);
            let n = normal_at(&rod, p, 1e-6).unwrap_or_else(|| panic!("no normal at {p:?}"));
            assert!((n - Vec3::Z).norm() < 1e-12, "{n:?} at {p:?}");
            let bend = curvature(&rod, p).unwrap_or_else(|| panic!("no curvature at {p:?}"));
            assert!(
                bend[0].abs() < 1e-9 && bend[1].abs() < 1e-9,
                "{bend:?} at {p:?}"
            );
        }
    }

    #[test]
    fn a_drilled_rods_side_bends_away_its_hole_toward_and_its_rim_is_a_crease() {
        // A rod of radius 10, 20 long, drilled through along its axis at radius 4.
        let circle = |radius: f64, ccw: bool, name: &str| ProfileLoop {
            segments: vec![NamedSegment::new(
                Segment::Circle {
                    center: Vec2::ZERO,
                    radius,
                    ccw,
                },
                name,
            )],
        };
        let p = Profile {
            loops: vec![circle(10.0, true, "out"), circle(4.0, false, "in")],
        };
        let rod = extrude(
            &p,
            &Frame::new(Vec3::ZERO, Vec3::Z, Vec3::X).unwrap(),
            20.0,
            "r",
        )
        .unwrap();
        let near =
            |a: [f64; 2], b: [f64; 2]| (a[0] - b[0]).abs() < 1e-9 && (a[1] - b[1]).abs() < 1e-9;
        let side = curvature(&rod, Point3::new(11.0, 1.0, 7.0)).unwrap();
        assert!(near(side, [0.1, 0.0]), "{side:?}");
        let hole = curvature(&rod, Point3::new(0.5, 3.0, 13.0)).unwrap();
        assert!(near(hole, [0.0, -0.25]), "{hole:?}");
        let end = curvature(&rod, Point3::new(6.0, 2.0, 20.5)).unwrap();
        assert!(near(end, [0.0, 0.0]), "{end:?}");
        assert!(curvature(&rod, Point3::new(10.5, 0.0, 20.5)).is_none());
    }
}
