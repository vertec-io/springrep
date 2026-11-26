//! Mass properties against closed forms: volume, centre and the inertia tensor about the centre.

use std::f64::consts::{PI, TAU};

use springrep::{
    boolean::intersection,
    codec,
    geom::{Curve, Surface},
    loft::{Section, loft},
    math::{Frame, Transform, Vec2, Vec3},
    measure::{self, MassProperties},
    profile::{NamedSegment, Profile, ProfileLoop, Segment},
    sweep::{extrude, revolve},
    topo::{Brep, PCurve, SolidId},
};

/// The diagonal tensor with these moments about the axes.
fn diagonal(i: [f64; 3]) -> [[f64; 3]; 3] {
    [[i[0], 0.0, 0.0], [0.0, i[1], 0.0], [0.0, 0.0, i[2]]]
}

fn assert_mass(b: &Brep, volume: f64, centre: [f64; 3], inertia: [[f64; 3]; 3]) {
    let MassProperties {
        volume: v,
        centre: c,
        inertia: i,
        ..
    } = measure::mass_properties(b);
    assert!(
        (v - volume).abs() <= 1e-12 * volume,
        "volume {v} vs {volume}"
    );
    let size = volume.cbrt();
    let d = (c - Vec3::from_array(centre)).max_abs();
    assert!(d <= 1e-12 * size, "centre {c:?} vs {centre:?}");
    let scale = inertia.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs()));
    for r in 0..3 {
        for k in 0..3 {
            assert!(
                (i[r][k] - inertia[r][k]).abs() <= 1e-11 * scale,
                "inertia {i:?} vs {inertia:?}"
            );
        }
    }
}

fn xz() -> Frame {
    Frame::new(Vec3::ZERO, Vec3::new(0.0, -1.0, 0.0), Vec3::X).unwrap()
}

fn seg(segment: Segment, name: &str) -> NamedSegment {
    NamedSegment::new(segment, name)
}

#[test]
fn a_box_off_the_origin_has_the_box_formulas() {
    // x 1..11, y 2..8, z 3..6: a 10 × 6 × 3 box of 180 about (6, 5, 4.5), each axis's moment
    // V(b² + c²)/12 in the other two sides, no products.
    let frame = Frame::new(Vec3::new(0.0, 0.0, 3.0), Vec3::Z, Vec3::X).unwrap();
    let b = extrude(
        &Profile::rect(Vec2::new(6.0, 5.0), 10.0, 6.0),
        &frame,
        3.0,
        "b",
    )
    .unwrap();
    let v = 180.0;
    let i = [v * 45.0 / 12.0, v * 109.0 / 12.0, v * 136.0 / 12.0];
    assert_mass(&b, v, [6.0, 5.0, 4.5], diagonal(i));
}

#[test]
fn a_cylinder_off_the_axis_has_the_cylinder_formulas() {
    // Radius 4 about (2, −1), height 5: V = 80π, across the axis V(3r² + h²)/12, along it Vr²/2.
    let b = extrude(
        &Profile::circle(Vec2::new(2.0, -1.0), 8.0),
        &Frame::WORLD,
        5.0,
        "c",
    )
    .unwrap();
    let v = 80.0 * PI;
    let across = v * (48.0 + 25.0) / 12.0;
    assert_mass(&b, v, [2.0, -1.0, 2.5], diagonal([across, across, v * 8.0]));
}

#[test]
fn a_ball_has_two_fifths_v_r_squared_about_every_axis() {
    let r = 3.0;
    let p = Profile {
        loops: vec![ProfileLoop {
            segments: vec![
                seg(
                    Segment::Arc {
                        center: Vec2::ZERO,
                        radius: r,
                        start: -PI / 2.0,
                        sweep: PI,
                    },
                    "surface",
                ),
                seg(
                    Segment::Line {
                        a: Vec2::new(0.0, r),
                        b: Vec2::new(0.0, -r),
                    },
                    "axis",
                ),
            ],
        }],
    };
    let b = revolve(&p, &xz(), Vec3::ZERO, Vec3::Z, TAU, "ball").unwrap();
    let v = 4.0 / 3.0 * PI * r * r * r;
    let i = 0.4 * v * r * r;
    assert_mass(&b, v, [0.0; 3], diagonal([i, i, i]));
}

#[test]
fn a_torus_has_the_torus_formulas() {
    // Tube radius a = 2 at R = 10 round z: V = 2π²Ra², along the axis V(R² + 3a²/4), across it
    // V(R²/2 + 5a²/8).
    let (big, a) = (10.0, 2.0);
    let p = Profile::circle(Vec2::new(big, 0.0), 2.0 * a);
    let b = revolve(&p, &xz(), Vec3::ZERO, Vec3::Z, TAU, "t").unwrap();
    let v = 2.0 * PI * PI * big * a * a;
    let across = v * (big * big / 2.0 + 5.0 * a * a / 8.0);
    let along = v * (big * big + 0.75 * a * a);
    assert_mass(&b, v, [0.0; 3], diagonal([across, across, along]));
}

#[test]
fn an_l_shaped_plate_has_products_of_inertia() {
    // The L is two 1-thick boxes: A, x 0..4 by y 0..2, and B, x 0..2 by y 2..6, each 8 about
    // (2, 1, ½) and (1, 4, ½). Their sum is 16 about (1.5, 2.5, ½); each box adds its own moments
    // V(b² + c²)/12 and, by the parallel-axis rule, V(d² − d dᵀ) for its offset d.
    let pts = [
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 2.0),
        (2.0, 2.0),
        (2.0, 6.0),
        (0.0, 6.0),
    ]
    .map(|(x, y)| Vec2::new(x, y));
    let names = ["s", "e", "n1", "e2", "n2", "w"];
    let b = extrude(&Profile::polygon(&pts, &names), &Frame::WORLD, 1.0, "l").unwrap();
    let centre = [1.5, 2.5, 0.5];
    let mut want = [[0.0; 3]; 3];
    for (size, at) in [
        ([4.0, 2.0, 1.0], [2.0, 1.0, 0.5]),
        ([2.0, 4.0, 1.0], [1.0, 4.0, 0.5]),
    ] {
        let v: f64 = size.iter().product();
        let s2 = size.map(|x| x * x);
        let own = [
            v * (s2[1] + s2[2]) / 12.0,
            v * (s2[0] + s2[2]) / 12.0,
            v * (s2[0] + s2[1]) / 12.0,
        ];
        let d = [at[0] - centre[0], at[1] - centre[1], at[2] - centre[2]];
        let dd: f64 = d.iter().map(|x| x * x).sum();
        for r in 0..3 {
            for k in 0..3 {
                let shift = v * (if r == k { dd } else { 0.0 } - d[r] * d[k]);
                want[r][k] += shift + if r == k { own[r] } else { 0.0 };
            }
        }
    }
    assert!(want[0][1] != 0.0);
    assert_mass(&b, 16.0, centre, want);
}

#[test]
fn a_tilted_cylinder_turns_its_tensor_with_it() {
    // A cylinder along a unit axis z has I = I⊥·1 + (I∥ − I⊥)·z zᵀ: every entry is nonzero when
    // z leans away from all three axes.
    let z = Vec3::new(1.0, 2.0, 2.0) * (1.0 / 3.0);
    let x = Vec3::new(2.0, -2.0, 1.0) * (1.0 / 3.0);
    let o = Vec3::new(1.0, -3.0, 2.0);
    let frame = Frame::new(o, z, x).unwrap();
    let (r, h) = (1.5, 6.0);
    let b = extrude(&Profile::circle(Vec2::ZERO, 2.0 * r), &frame, h, "c").unwrap();
    let v = PI * r * r * h;
    let (across, along) = (v * (3.0 * r * r + h * h) / 12.0, v * r * r / 2.0);
    let zs = [z.x, z.y, z.z];
    let mut want = [[0.0; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            want[i][k] = if i == k { across } else { 0.0 } + (along - across) * zs[i] * zs[k];
        }
    }
    let c = o + z * (h / 2.0);
    assert_mass(&b, v, [c.x, c.y, c.z], want);
}

/// Each way a solid's volume is measured, each to its own tolerance and on a fresh copy so no face
/// memo carries from one to the next, against `want`, to `rel` of it.
fn assert_volume_every_way(b: &Brep, want: f64, rel: f64) {
    let fresh = || codec::from_bytes(&codec::to_bytes(std::slice::from_ref(b))).unwrap();
    let per_solid = |f: fn(&Brep, SolidId) -> f64| -> f64 {
        let c = &fresh()[0];
        c.solid_ids().map(|s| f(c, s)).sum()
    };
    for (what, v) in [
        ("rough sign", per_solid(measure::solid_volume_sign)),
        ("to 1e-10", measure::volume_to(&fresh()[0], 1e-10)),
        ("exact", per_solid(measure::solid_volume)),
        ("moments", measure::volume(&fresh()[0])),
    ] {
        assert!((v - want).abs() <= rel * want, "{what}: {v} vs {want}");
    }
}

#[test]
fn a_lofted_cylinder_cut_on_a_slant_measures_its_closed_form_every_way() {
    // A smooth loft through equal circles is a cylinder on a rational spline surface. Below the
    // plane z = h + x·tan θ it holds πr²h, the slant cancelling across the axis, and its area is
    // the disc, the ellipse πr²/cos θ and the side 2πrh. The plane meets the spline side along a
    // fitted curve, so the boundary integrals there run through projections.
    let (r, h, theta) = (5.0, 7.0, 0.3);
    let sections: Vec<Section> = [0.0, 5.0, 10.0, 15.0]
        .iter()
        .map(|&z| {
            Section::new(
                Profile::circle(Vec2::ZERO, 2.0 * r),
                Frame::new(Vec3::new(0.0, 0.0, z), Vec3::Z, Vec3::X).unwrap(),
            )
        })
        .collect();
    let can = loft(&sections, false, "can").unwrap();
    assert!(
        can.face_ids()
            .any(|f| matches!(can.face(f).surface, Surface::Nurbs(_)))
    );
    let below = Frame::new(Vec3::new(0.0, 0.0, -100.0), Vec3::Z, Vec3::X).unwrap();
    let mut cut = extrude(
        &Profile::rect(Vec2::ZERO, 100.0, 100.0),
        &below,
        100.0,
        "cut",
    )
    .unwrap();
    cut.transform(
        &Transform::rotation(Vec3::ZERO, Vec3::Y, -theta)
            .then(&Transform::translation(Vec3::new(0.0, 0.0, h))),
    );
    let b = intersection(&can, &cut).unwrap();
    let fitted = b.face_ids().any(|f| {
        matches!(b.face(f).surface, Surface::Nurbs(_))
            && b.face(f).loops.iter().any(|&lp| {
                b.lp(lp).coedges.iter().any(|&co| {
                    matches!(b.coedge(co).pcurve, PCurve::Projected { .. })
                        && matches!(
                            b.edge(b.coedge(co).edge).curve,
                            Some(Curve::Intersection(_))
                        )
                })
            })
    });
    assert!(fitted);
    assert_volume_every_way(&b, PI * r * r * h, 1e-9);
    let area = PI * r * r * (1.0 + 1.0 / theta.cos()) + TAU * r * h;
    let got = measure::area(&b);
    assert!((got - area).abs() <= 1e-9 * area, "area {got} vs {area}");
}

/// An assembly measured face by face over the machine's threads: sixteen prisms 7 apart along x,
/// each 5 tall on the regular 48-gon of circumradius 3 drawn as a closed spline of degree one,
/// whose edges' 48 spans each are pieces of the boundary integrals, several thousand in all.
///
/// The 48-gon is 48 triangles from its centre, each `½R² sin θ` with `θ = 2π/48`, so its area is
/// `A = 24R² sin θ` and its outline `96R sin(θ/2)`. A triangle from the origin to `a` and `b`
/// has `∫|r|² = (Aₜ/6)(|a|² + a·b + |b|²)`, so the 48-gon's polar moment is `A R²(2 + cos θ)/6`,
/// half of it about each axis in its plane by symmetry, and a prism's moments about its centre
/// are those times its height, plus `V h²/12` about the axes across it.
#[test]
fn an_assembly_of_spline_prisms_has_the_sums_of_its_prisms_closed_forms() {
    let (n, r, h, pitch, k) = (48usize, 3.0, 5.0, 7.0, 16usize);
    let theta = TAU / n as f64;
    let mut arena = Brep::new();
    for i in 0..k {
        let corners: Vec<Vec3> = (0..n)
            .map(|j| {
                let t = theta * j as f64;
                Vec3::new(r * t.cos(), r * t.sin(), 0.0)
            })
            .collect();
        let curve = springrep::geom::BSpline::periodic_uniform(1, &corners).unwrap();
        let side = Segment::Spline {
            curve,
            reversed: false,
        };
        let outline = Profile::outline(vec![seg(side, "side")]);
        let at = Frame::new(Vec3::new(pitch * i as f64, 0.0, 0.0), Vec3::Z, Vec3::X).unwrap();
        arena.append(&extrude(&outline, &at, h, "prism").unwrap());
    }
    let section = 0.5 * n as f64 * r * r * theta.sin();
    let outline = 2.0 * n as f64 * r * (0.5 * theta).sin();
    let (volume, skin) = (
        k as f64 * section * h,
        k as f64 * (2.0 * section + outline * h),
    );

    let (a, v) = measure::area_and_volume(&arena);
    assert!(
        (v - volume).abs() <= 1e-12 * volume,
        "volume {v} vs {volume}"
    );
    assert!((a - skin).abs() <= 1e-12 * skin, "area {a} vs {skin}");
    let (a, v) = (measure::area(&arena), measure::volume(&arena));
    assert!(
        (v - volume).abs() <= 1e-12 * volume,
        "volume {v} vs {volume}"
    );
    assert!((a - skin).abs() <= 1e-12 * skin, "area {a} vs {skin}");

    let one = section * h;
    let mid = pitch * (k - 1) as f64 / 2.0;
    let spread: f64 = (0..k).map(|i| (pitch * i as f64 - mid).powi(2)).sum();
    let across = k as f64 * one * (r * r * (2.0 + theta.cos()) / 12.0 + h * h / 12.0);
    let about_axis = k as f64 * one * r * r * (2.0 + theta.cos()) / 6.0;
    let inertia = diagonal([across, across + one * spread, about_axis + one * spread]);
    assert_mass(&arena, volume, [mid, 0.0, h / 2.0], inertia);
}
