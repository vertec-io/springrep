use std::f64::consts::{PI, TAU};

use super::{ReadOptions, WriteOptions, read, read_part_with, read_with, write, write_with};
use crate::{
    boolean, check,
    geom::{Circle, Cone, Curve, Line, Ruled, Surface},
    math::{Frame, Vec2, Vec3},
    measure, mesh,
    profile::{NamedSegment, PointsInterp, Profile, ProfileLoop, Segment},
    sweep::{extrude, revolve},
    topo::Brep,
};

fn seg(segment: Segment, name: &str) -> NamedSegment {
    NamedSegment::new(segment, name)
}

fn block(id: &str, lo: [f64; 3], hi: [f64; 3]) -> Brep {
    let c = Vec2::new(0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]));
    let frame = Frame::new(Vec3::new(0.0, 0.0, lo[2]), Vec3::Z, Vec3::X).unwrap();
    let p = Profile::rect(c, hi[0] - lo[0], hi[1] - lo[1]);
    extrude(&p, &frame, hi[2] - lo[2], id).unwrap()
}

fn xz() -> Frame {
    Frame::new(Vec3::ZERO, Vec3::new(0.0, -1.0, 0.0), Vec3::X).unwrap()
}

fn ball(r: f64) -> Brep {
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
    revolve(&p, &xz(), Vec3::ZERO, Vec3::Z, TAU, "ball").unwrap()
}

/// Writes `b`, reads it back, and checks the copy is valid with the same measures and names.
/// `b` written as STEP, with pcurves and without, and read back as the same solid each way.
fn round_trip(b: &Brep) -> Brep {
    assert_within_step_domains(b);
    let bare = write_with(b, "part", &WriteOptions::default().pcurves(false)).unwrap();
    for kind in [
        "PCURVE",
        "SURFACE_CURVE",
        "SEAM_CURVE",
        "DEFINITIONAL_REPRESENTATION",
    ] {
        assert_eq!(instances(&bare, kind), 0, "{kind} without pcurves");
    }
    read_back(b, &bare);
    let text = write(b, "part").unwrap();
    assert_eq!(
        write_with(b, "part", &WriteOptions::default()).unwrap(),
        text
    );
    read_back(b, &text)
}

fn read_back(b: &Brep, text: &str) -> Brep {
    assert!(check::check(b).is_empty());
    let back = read(text.as_bytes(), "imp").unwrap();
    assert!(back.rejected.is_empty(), "{:?}", back.rejected);
    let c = back.brep;
    assert!(check::check(&c).is_empty(), "{:?}", check::check(&c));
    assert_eq!(c.solid_ids().count(), b.solid_ids().count());
    assert_eq!(c.face_count(), b.face_count());
    for (got, want, what) in [
        (measure::volume(&c), measure::volume(b), "volume"),
        (measure::area(&c), measure::area(b), "area"),
    ] {
        assert!((got - want).abs() <= 1e-9 * want, "{what} {got} vs {want}");
    }
    let (bb, cb) = (measure::bbox(b), measure::bbox(&c));
    assert!((bb.min - cb.min).max_abs().max((bb.max - cb.max).max_abs()) < 1e-9);
    assert_eq!(c.name_counts(), b.name_counts());
    c
}

#[test]
fn planar_solids_round_trip() {
    round_trip(&block("blk", [0.0, 0.0, 0.0], [10.0, 4.0, 3.0]));
    let hole = Profile::circle(Vec2::new(2.0, 0.0), 2.0).loops.remove(0);
    let p = Profile::rect(Vec2::ZERO, 10.0, 6.0).with_holes([hole]);
    round_trip(&extrude(&p, &Frame::WORLD, 2.0, "plate").unwrap());
}

#[test]
fn spheres_and_tori_round_trip_through_their_poles_and_seams() {
    round_trip(&ball(5.0));
    let ring = Profile::circle(Vec2::new(8.0, 0.0), 4.0);
    round_trip(&revolve(&ring, &xz(), Vec3::ZERO, Vec3::Z, TAU, "ring").unwrap());
    let half = revolve(&ring, &xz(), Vec3::ZERO, Vec3::Z, PI, "arc").unwrap();
    round_trip(&half);
}

/// Cones of all four kinds of half angle, each written with a STEP semi-angle in `[0, π/2)`.
#[test]
fn cones_round_trip_both_ways() {
    for pts in [
        [(0.0, 0.0), (6.0, 0.0), (6.0, 2.0), (3.0, 5.0), (0.0, 5.0)],
        [(0.0, 0.0), (3.0, 0.0), (6.0, 3.0), (6.0, 5.0), (0.0, 5.0)],
    ] {
        for (flip, axis) in [
            (1.0, Vec3::Z),
            (-1.0, Vec3::Z),
            (1.0, -Vec3::Z),
            (-1.0, -Vec3::Z),
        ] {
            let pts: Vec<Vec2> = pts.iter().map(|&(x, y)| Vec2::new(x, flip * y)).collect();
            let p = Profile::polygon(&pts, &["a", "b", "c", "d", "e"]);
            let b = revolve(&p, &xz(), Vec3::ZERO, axis, TAU, "t").unwrap();
            let text = write(&b, "t").unwrap();
            for l in text.lines().filter(|l| l.contains("=CONICAL_SURFACE(")) {
                let a: f64 = l
                    .rsplit(',')
                    .next()
                    .unwrap()
                    .trim_end_matches(");")
                    .parse()
                    .unwrap();
                assert!((0.0..PI / 2.0).contains(&a), "{l}");
            }
            round_trip(&b);
        }
    }
    // Each quadrant of half angle reads back as the same points, with the normal flip the face's
    // sense undoes.
    for a in [0.4, -0.4, PI - 0.4, 0.4 - PI] {
        let cone = Surface::Cone(Cone {
            frame: Frame::new(Vec3::new(1.0, 2.0, 3.0), Vec3::new(0.0, 0.6, 0.8), Vec3::X).unwrap(),
            radius: 4.0,
            half_angle: a,
        });
        let mut w = super::write::Writer::new();
        let id = w.surface(&cone).unwrap();
        let text = format!(
            "ISO-10303-21;\nDATA;\n{}ENDSEC;\nEND-ISO-10303-21;\n",
            w.out
        );
        let file = super::part21::parse(text.as_bytes()).unwrap();
        let back = super::geometry::Geometry::new(&file)
            .surface(id as u64)
            .unwrap();
        let flips = a.abs() > PI / 2.0;
        for (u, v) in [(0.3, 0.5), (2.0, -1.0), (4.0, 2.5)] {
            let p = cone.eval(u, v);
            let q = back.project(p);
            assert!(back.eval(q.x, q.y).dist(p) < 1e-9, "{a}");
            let dot = back.normal(q.x, q.y).dot(cone.normal(u, v));
            assert!(
                (dot - if flips { -1.0 } else { 1.0 }).abs() < 1e-9,
                "{a}: {dot}"
            );
        }
    }
}

#[test]
fn spline_extrusions_and_revolutions_round_trip() {
    let pts: Vec<Vec2> = (0..7)
        .map(|i| {
            let a = TAU * i as f64 / 7.0;
            let r = 10.0 + 2.0 * (3.0 * a).cos();
            Vec2::new(r * a.cos(), r * a.sin())
        })
        .collect();
    let p = Profile::points(&pts, PointsInterp::PeriodicSpline, "side").unwrap();
    round_trip(&extrude(&p, &Frame::WORLD, 4.0, "cam").unwrap());
    let open = [
        (0.0, 0.0),
        (4.0, 0.0),
        (5.0, 2.0),
        (4.5, 4.0),
        (5.5, 6.0),
        (0.0, 6.0),
    ];
    let open: Vec<Vec2> = open.iter().map(|&(x, y)| Vec2::new(x, y)).collect();
    let prof = Profile::points(&open, PointsInterp::Smooth { corner: 0.5 }, "s").unwrap();
    round_trip(&revolve(&prof, &xz(), Vec3::ZERO, Vec3::Z, TAU, "vase").unwrap());
}

#[test]
fn boolean_results_with_voids_round_trip() {
    let outer = block("o", [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
    let inner = block("i", [3.0, 3.0, 3.0], [6.0, 6.0, 6.0]);
    let hollow = boolean::difference(&outer, &inner).unwrap();
    assert_eq!(
        hollow
            .solid(hollow.solid_ids().next().unwrap())
            .shells
            .len(),
        2
    );
    round_trip(&hollow);
}

/// Some systems write a void's faces already facing into it and turn its shell as well, so the
/// shell's own orientation and the solid's say opposite things: the void is read by the volume it
/// encloses, which is negative.
#[test]
fn a_void_turned_twice_is_read_as_a_void() {
    let outer = block("o", [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
    let inner = block("i", [3.0, 3.0, 3.0], [6.0, 6.0, 6.0]);
    let hollow = boolean::difference(&outer, &inner).unwrap();
    let text = write(&hollow, "part").unwrap();
    // The void's faces: those of the shell its ORIENTED_CLOSED_SHELL names.
    let entity = |id: &str| {
        let at = text.find(&format!("\n#{id}=")).unwrap() + 1;
        text[at..at + text[at..].find(";\n").unwrap()].to_string()
    };
    let refs_in = |s: &str| -> Vec<String> {
        s.split('#')
            .skip(1)
            .map(|r| r.chars().take_while(char::is_ascii_digit).collect())
            .collect()
    };
    let at = text.find("=ORIENTED_CLOSED_SHELL(").unwrap();
    let line = &text[text[..at].rfind('\n').unwrap() + 1..at + text[at..].find(";\n").unwrap()];
    let shell = refs_in(line.split_once('=').unwrap().1)[0].clone();
    // Each face's sense and each of its bounds' orientation turned: the same face, facing the
    // other way, its loops still keeping it on their left.
    let mut turned = text.clone();
    let mut turn = |id: &str| {
        let line = entity(id);
        let flipped = match line.strip_suffix(".T.)") {
            Some(head) => format!("{head}.F.)"),
            None => format!("{}.T.)", line.strip_suffix(".F.)").unwrap()),
        };
        turned = turned.replace(&format!("\n{line};"), &format!("\n{flipped};"));
    };
    for f in refs_in(entity(&shell).split_once('=').unwrap().1) {
        let face = entity(&f);
        let bounds = face
            .split_once("(")
            .unwrap()
            .1
            .split_once(")")
            .unwrap()
            .0
            .to_string();
        for bd in refs_in(&bounds) {
            turn(&bd);
        }
        turn(&f);
    }
    assert_ne!(turned, text);
    let back = read(turned.as_bytes(), "imp").unwrap();
    assert!(back.rejected.is_empty(), "{:?}", back.rejected);
    assert!(check::check(&back.brep).is_empty());
    assert!((measure::volume(&back.brep) - 973.0).abs() < 1e-9);
}

/// How many instances of each of `kinds` the file has.
fn instances(text: &str, kind: &str) -> usize {
    text.matches(&format!("={kind}(")).count()
}

/// Writes a solid whose pieces touch, checks the file holds them as manifold solids, as systems
/// that write the manifold boundary of a solid's interior do: each piece its own solid, the
/// contact's vertices and edges written once for each side, no `VERTEX_LOOP`. Returns the file and
/// what reading it gives.
fn written_as_pieces(b: &Brep, solids: usize) -> (String, Brep) {
    assert!(check::check(b).is_empty());
    let text = write(b, "part").unwrap();
    let written = instances(&text, "MANIFOLD_SOLID_BREP") + instances(&text, "BREP_WITH_VOIDS");
    assert_eq!(written, solids, "{text}");
    assert_eq!(instances(&text, "VERTEX_LOOP"), 0);
    let back = read(text.as_bytes(), "imp").unwrap();
    assert!(back.rejected.is_empty(), "{:?}", back.rejected);
    let c = back.brep;
    assert!(
        check::check_each(&c).is_empty(),
        "{:?}",
        check::check_each(&c)
    );
    let (vb, vc) = (measure::volume(b), measure::volume(&c));
    assert!((vb - vc).abs() <= 1e-9 * vb, "volume {vc} vs {vb}");
    (text, c)
}

/// The vertices with a point, and the edges with a curve, of `b`.
fn vertices_and_edges(b: &Brep) -> (usize, usize) {
    let edges = b.edge_ids().filter(|&e| !b.edge(e).is_degenerate()).count();
    (b.vertex_ids().count(), edges)
}

/// Two cubes sharing an edge are two solids in the file, each with all 8 corners and 12 edges.
#[test]
fn pieces_touching_along_an_edge_are_written_as_two_solids() {
    let a = block("a", [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
    let b = block("b", [10.0, 10.0, 0.0], [20.0, 20.0, 10.0]);
    let u = boolean::union(&a, &b).unwrap();
    let (text, back) = written_as_pieces(&u, 2);
    assert_eq!(instances(&text, "VERTEX_POINT"), 16);
    assert_eq!(instances(&text, "EDGE_CURVE"), 24);
    assert_eq!(back.solid_ids().count(), 2);
}

/// A rod lying on a box's top face along a line: in the result the line is a slit in the top face,
/// which runs along it and back, and splits the rod's side in two. In the box's own boundary it is
/// no edge at all, so the file has a plain box; the rod keeps its two halves.
#[test]
fn a_slit_where_a_piece_touches_a_face_is_not_written() {
    let a = block("a", [0.0, 0.0, 0.0], [20.0, 10.0, 10.0]);
    let frame = Frame::new(Vec3::new(5.0, 5.0, 12.0), Vec3::X, Vec3::Y).unwrap();
    let rod = extrude(&Profile::circle(Vec2::ZERO, 4.0), &frame, 10.0, "c").unwrap();
    let u = boolean::union(&a, &rod).unwrap();
    assert_eq!(u.radial_edges().count(), 1);
    let (text, back) = written_as_pieces(&u, 2);
    // The box's 12 edges, and the rod's: the line between its halves, its seam and the two arcs
    // of each end.
    assert_eq!(instances(&text, "EDGE_CURVE"), 18);
    assert_eq!(back.solid_ids().count(), 2);
}

/// Two balls touching at a point are two balls in the file; the point is in neither.
#[test]
fn balls_touching_at_a_point_are_written_as_two_solids() {
    let mut b = ball(5.0);
    b.transform(&crate::math::Transform::translation(Vec3::new(
        0.0, 10.0, 0.0,
    )));
    let u = boolean::union(&ball(5.0), &b).unwrap();
    let (_, back) = written_as_pieces(&u, 2);
    assert_eq!(back.solid_ids().count(), 2);
}

/// A block less a bore whose top line lies on its top face is one solid touching itself along
/// that line: the file gives each side of the contact its own edge, and reading it joins them into
/// the edge four faces use again.
#[test]
fn a_solid_touching_itself_along_a_line_is_one_solid_read_back_joined() {
    let a = block("a", [0.0, 0.0, 0.0], [12.0, 10.0, 8.0]);
    let frame = Frame::new(Vec3::new(-1.0, 5.0, 5.0), Vec3::X, Vec3::Y).unwrap();
    let bore = extrude(&Profile::circle(Vec2::ZERO, 6.0), &frame, 14.0, "c").unwrap();
    let d = boolean::difference(&a, &bore).unwrap();
    assert_eq!(d.radial_edges().count(), 1);
    let (text, back) = written_as_pieces(&d, 1);
    let (vertices, edges) = vertices_and_edges(&d);
    // The contact line and its two ends, once more each.
    assert_eq!(instances(&text, "EDGE_CURVE"), edges + 1);
    assert_eq!(instances(&text, "VERTEX_POINT"), vertices + 2);
    assert!(check::check(&back).is_empty());
    assert_eq!(back.radial_edges().count(), 1);
    assert_eq!(vertices_and_edges(&back), (vertices, edges));
    assert_eq!(measure::betti(&back), [1, 1, 0]);
}

/// A cube less the ball touching its six faces: the void and the skin are written apart, the
/// six points where they touch in neither, so the file's solid has its void but not the contacts.
#[test]
fn a_void_touching_the_skin_at_points_is_written_apart_from_it() {
    let mut s = ball(5.0);
    s.transform(&crate::math::Transform::translation(Vec3::new(
        5.0, 5.0, 5.0,
    )));
    let d = boolean::difference(&block("c", [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]), &s).unwrap();
    assert!(d.loop_ids().any(|lp| d.is_vertex_loop(lp)));
    let (text, back) = written_as_pieces(&d, 1);
    assert_eq!(instances(&text, "BREP_WITH_VOIDS"), 1);
    assert_eq!(back.solid_ids().count(), 1);
    assert_eq!(measure::betti(&back), [1, 0, 1]);
}

#[test]
fn intersection_curves_are_written_as_their_splines() {
    let a = Surface::Cylinder(crate::geom::Cylinder {
        frame: Frame::WORLD,
        radius: 5.0,
    });
    let b = Surface::Cylinder(crate::geom::Cylinder {
        frame: Frame::new(Vec3::ZERO, Vec3::X, Vec3::Y).unwrap(),
        radius: 3.0,
    });
    let pts: Vec<Vec3> = (0..12)
        .map(|i| {
            let t = TAU * i as f64 / 12.0;
            let (y, z) = (3.0 * t.cos(), 3.0 * t.sin());
            Vec3::new((25.0 - y * y).sqrt(), y, z)
        })
        .collect();
    let fit = crate::geom::IntersectionCurve::fit(&pts, a.clone(), b.clone(), true).unwrap();
    let curve = Curve::Intersection(Box::new(fit));
    let mut w = super::write::Writer::new();
    let id = w.curve(&curve).unwrap();
    let text = format!(
        "ISO-10303-21;\nDATA;\n{}ENDSEC;\nEND-ISO-10303-21;\n",
        w.out
    );
    let file = super::part21::parse(text.as_bytes()).unwrap();
    let back = super::geometry::Geometry::new(&file)
        .curve(id as u64)
        .unwrap();
    assert!(back.period().is_some());
    let (t0, t1) = curve.domain();
    for i in 0..=50 {
        let t = t0 + (t1 - t0) * i as f64 / 50.0;
        assert!(back.eval(t).dist(curve.eval(t)) < 1e-12);
        assert!(crate::geom::distance_to(&a, back.eval(t)) < crate::tol::LINEAR);
        assert!(crate::geom::distance_to(&b, back.eval(t)) < crate::tol::LINEAR);
    }
}

#[test]
fn ruled_surfaces_are_written_as_the_same_surface() {
    let circle = |z: f64, r: f64| {
        Curve::Circle(Circle {
            frame: Frame::WORLD.translated(Vec3::new(0.0, 0.0, z)),
            radius: r,
        })
    };
    let line = |a: Vec3, b: Vec3| Curve::Line(Line::through(a, b).unwrap().0);
    let cases = [
        Ruled {
            a: circle(0.0, 3.0),
            a_range: (0.0, TAU),
            b: circle(5.0, 2.0),
            b_range: (1.0, 1.0 + TAU),
        },
        Ruled {
            a: line(Vec3::ZERO, Vec3::X),
            a_range: (0.0, 1.0),
            b: Curve::BSpline(
                crate::geom::BSpline::interpolate(&[
                    Vec3::new(0.0, 2.0, 0.0),
                    Vec3::new(0.3, 2.5, 0.0),
                    Vec3::new(0.7, 2.2, 1.0),
                    Vec3::new(1.0, 2.0, 0.0),
                ])
                .unwrap(),
            ),
            b_range: (0.0, 1.0),
        },
    ];
    for r in cases {
        let ruled = Surface::Ruled(Box::new(r));
        let mut w = super::write::Writer::new();
        let id = w.surface(&ruled).unwrap();
        let text = format!(
            "ISO-10303-21;\nDATA;\n{}ENDSEC;\nEND-ISO-10303-21;\n",
            w.out
        );
        let file = super::part21::parse(text.as_bytes()).unwrap();
        let mut geo = super::geometry::Geometry::new(&file);
        let Surface::Nurbs(n) = geo.surface(id as u64).unwrap() else {
            panic!("expected a B-spline surface")
        };
        // The B-spline's parameter along the rails need not be the ruled surface's, but each
        // of its rulings must be one of the ruled surface's, run at the same rate.
        let ((u0, u1), _) = n.domain();
        for i in 0..=10 {
            let u = u0 + (u1 - u0) * i as f64 / 10.0;
            let at = ruled.project(n.eval(u, 0.0));
            assert!(at.y.abs() < 1e-9);
            for j in 0..=4 {
                let v = j as f64 / 4.0;
                let d = n.eval(u, v).dist(ruled.eval(at.x, v));
                assert!(d < 1e-9, "{d} at {u},{v}");
            }
        }
    }
}

/// A part of two blocks placed by an assembly, in inches, as other systems write them. With
/// `cycle`, the part also uses itself.
fn assembly_text(cycle: bool) -> String {
    let part = write(&block("b", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), "blk").unwrap();
    let data = part
        .split("DATA;\n")
        .nth(1)
        .unwrap()
        .split("ENDSEC;")
        .next()
        .unwrap()
        .replace(
            "(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))",
            "(CONVERSION_BASED_UNIT('INCH',#9001) LENGTH_UNIT() NAMED_UNIT(*))",
        );
    let pd = data
        .lines()
        .find(|l| l.contains("=PRODUCT_DEFINITION('design'"))
        .and_then(|l| l[1..].split('=').next())
        .unwrap()
        .to_string();
    let rep = data
        .lines()
        .find(|l| l.contains("ADVANCED_BREP_SHAPE_REPRESENTATION"))
        .and_then(|l| l[1..].split('=').next())
        .unwrap()
        .to_string();
    let extra = format!(
        "#9000=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));
#9001=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#9000);
#9010=PRODUCT('asm','asm','',(#3));
#9011=PRODUCT_DEFINITION_FORMATION('','',#9010);
#9012=PRODUCT_DEFINITION('design','',#9011,#6);
#9013=PRODUCT_DEFINITION_SHAPE('','',#9012);
#9014=SHAPE_REPRESENTATION('asm',(#9020),#9099);
#9015=SHAPE_DEFINITION_REPRESENTATION(#9013,#9014);
#9020=AXIS2_PLACEMENT_3D('',#9021,#9022,#9023);
#9021=CARTESIAN_POINT('',(0.,0.,0.));
#9022=DIRECTION('',(0.,0.,1.));
#9023=DIRECTION('',(1.,0.,0.));
#9030=AXIS2_PLACEMENT_3D('',#9031,#9022,#9023);
#9031=CARTESIAN_POINT('',(10.,0.,0.));
#9040=AXIS2_PLACEMENT_3D('',#9041,#9042,#9023);
#9041=CARTESIAN_POINT('',(0.,0.,5.));
#9042=DIRECTION('',(0.,-1.,0.));
#9099=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#9000)) \
       REPRESENTATION_CONTEXT('','3D'));
{}",
        [("9030", "1", "9012"), ("9040", "2", "9012"), ("9030", "3", &pd)]
            .iter()
            .take(if cycle { 3 } else { 2 })
            .map(|(place, k, parent)| format!(
                "#91{k}0=NEXT_ASSEMBLY_USAGE_OCCURRENCE('{k}','','',#{parent},#{pd},$);
#91{k}1=ITEM_DEFINED_TRANSFORMATION('','',#9020,#{place});
#91{k}2=(REPRESENTATION_RELATIONSHIP('','',#{rep},#9014) \
       REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#91{k}1) SHAPE_REPRESENTATION_RELATIONSHIP());
#91{k}3=PRODUCT_DEFINITION_SHAPE('','',#91{k}0);
#91{k}4=CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#91{k}2,#91{k}3);
"
            ))
            .collect::<String>()
    );
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_NAME('a','',(''),(''),'','other','');\nENDSEC;\nDATA;\n\
         {data}{extra}ENDSEC;\nEND-ISO-10303-21;\n"
    )
}

#[test]
fn assemblies_place_each_instance_in_file_units() {
    let text = assembly_text(false);
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    let b = imp.brep;
    assert_eq!(b.solid_ids().count(), 2);
    let inch: f64 = 25.4;
    assert!((measure::volume(&b) - 2.0 * inch.powi(3)).abs() < 1e-9 * inch.powi(3));
    let bb = measure::bbox(&b);
    // The first block moves 10 mm along x; the second turns y→z, z→−y and rises 5 mm.
    assert!(
        (bb.min - Vec3::new(0.0, -inch, 0.0)).max_abs() < 1e-9,
        "{bb:?}"
    );
    assert!(
        (bb.max - Vec3::new(10.0 + inch, inch, 5.0 + inch)).max_abs() < 1e-9,
        "{bb:?}"
    );
    assert_eq!(b.name_counts().get("imp.face_0"), Some(&1));
    assert_eq!(b.name_counts().len(), 12);
}

/// The text of `b` written as STEP with the point of its first vertex moved by `d`.
fn with_vertex_moved(b: &Brep, d: Vec3) -> String {
    let text = write(b, "blk").unwrap();
    let vp = text.lines().find(|l| l.contains("=VERTEX_POINT(")).unwrap();
    with_point_moved(&text, vp, d)
}

/// `text`, a file this kernel wrote, with the point of its vertex `vp`, a line of it, moved by `d`,
/// and written by another system, so that a read heals it rather than keeping it as it was written.
fn with_point_moved(text: &str, vp: &str, d: Vec3) -> String {
    let pid = vp.split('#').nth(2).unwrap().trim_end_matches(");");
    let head = format!("#{pid}=CARTESIAN_POINT('',(");
    let line = text.lines().find(|l| l.starts_with(&head)).unwrap();
    let xyz: Vec<f64> = line[head.len()..line.len() - 3]
        .split(',')
        .map(|x| x.parse().unwrap())
        .collect();
    // A new point, so geometry that shared the vertex's point stays where it was.
    let moved = format!(
        "#99999=CARTESIAN_POINT('',({:?},{:?},{:?}));\n",
        xyz[0] + d.x,
        xyz[1] + d.y,
        xyz[2] + d.z
    );
    let vp_moved = vp.replace(&format!("#{pid})"), "#99999)");
    text.replace(vp, &vp_moved)
        .replace("ENDSEC;\nEND-ISO", &format!("{moved}ENDSEC;\nEND-ISO"))
        .replace("'springrep'", "'other'")
}

#[test]
fn a_vertex_off_its_faces_is_moved_back_onto_them() {
    let blk = block("b", [0.0, 0.0, 0.0], [10.0, 4.0, 3.0]);
    let text = with_vertex_moved(&blk, Vec3::new(2e-5, -1e-5, 3e-5));
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(check::check(&imp.brep).is_empty());
    assert!((measure::volume(&imp.brep) - 120.0).abs() < 1e-9);
    let bb = measure::bbox(&imp.brep);
    assert!(
        bb.min
            .max_abs()
            .max((bb.max - Vec3::new(10.0, 4.0, 3.0)).max_abs())
            < 1e-9
    );
}

#[test]
fn faces_that_do_not_meet_are_rejected_not_widened() {
    // A vertex a millimetre off its corner is far past anything healing may move.
    let blk = block("b", [0.0, 0.0, 0.0], [10.0, 4.0, 3.0]);
    let text = with_vertex_moved(&blk, Vec3::new(1.0, 0.0, 0.0));
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert_eq!(imp.brep.solid_ids().count(), 0);
    assert_eq!(imp.rejected.len(), 1);
}

#[test]
fn reference_cycles_and_short_points_are_errors() {
    assert!(read(assembly_text(true).as_bytes(), "imp").is_err());
    let blk = block("b", [0.0, 0.0, 0.0], [10.0, 4.0, 3.0]);
    let text = write(&blk, "blk").unwrap();
    let line = |kind: &str| {
        let l = text.lines().find(|l| l.contains(kind)).unwrap();
        (l.to_string(), l.split('=').next().unwrap().to_string())
    };
    let (vp, vid) = line("=VERTEX_POINT(");
    let (pl, pid) = line("=PLANE(");
    let (cp, _) = line("=CARTESIAN_POINT(");
    let mutants = [
        text.replace(&vp, &format!("{vid}=VERTEX_POINT('',{vid});")),
        text.replace(
            &pl,
            &format!("{pid}=RECTANGULAR_TRIMMED_SURFACE('',{pid},0.,1.,0.,1.,.T.,.T.);"),
        ),
        text.replace(
            &cp,
            &format!(
                "{}=CARTESIAN_POINT('',(1.0));",
                cp.split('=').next().unwrap()
            ),
        ),
    ];
    for m in mutants {
        let imp = read(m.as_bytes(), "imp").unwrap();
        assert_eq!(imp.brep.solid_ids().count(), 0);
        assert_eq!(imp.rejected.len(), 1, "{:?}", imp.rejected);
    }
}

/// A block with one vertical edge rounded: the round is a cylinder tangent to two planes.
fn rounded_block() -> Brep {
    let (w, h, r) = (10.0, 6.0, 1.0);
    let line = |a: (f64, f64), b: (f64, f64), name: &str| {
        seg(
            Segment::Line {
                a: Vec2::new(a.0, a.1),
                b: Vec2::new(b.0, b.1),
            },
            name,
        )
    };
    let round = seg(
        Segment::Arc {
            center: Vec2::new(w - r, h - r),
            radius: r,
            start: 0.0,
            sweep: PI / 2.0,
        },
        "round",
    );
    let p = Profile {
        loops: vec![ProfileLoop {
            segments: vec![
                line((0.0, 0.0), (w, 0.0), "s"),
                line((w, 0.0), (w, h - r), "e"),
                round,
                line((w - r, h), (0.0, h), "n"),
                line((0.0, h), (0.0, 0.0), "w"),
            ],
        }],
    };
    extrude(&p, &Frame::WORLD, 3.0, "rb").unwrap()
}

/// `b` written as another system's STEP with every cylinder's radius grown by `dr`, so the
/// rounds no longer touch the sides they are tangent to.
fn with_rounds_grown(b: &Brep, dr: f64) -> String {
    let text = write(b, "rb").unwrap().replace("'springrep'", "'other'");
    text.lines()
        .map(|l| match l.find("=CYLINDRICAL_SURFACE(") {
            Some(_) => {
                let (head, r) = l.rsplit_once(',').unwrap();
                let r: f64 = r.trim_end_matches(");").parse().unwrap();
                format!("{head},{:?});", r + dr)
            }
            None => l.to_string(),
        })
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn faces_held_apart_are_closed_by_moving_their_surfaces() {
    let b = rounded_block();
    let dr = 3e-5;
    let imp = read(with_rounds_grown(&b, dr).as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(check::check(&imp.brep).is_empty());
    assert_eq!(imp.healed.len(), 1);
    let h = &imp.healed[0];
    assert!(h.moved >= 0.5 * dr && h.moved <= 0.05, "{h:?}");
    assert!(h.surfaces > 0.0 && h.surfaces <= h.moved, "{h:?}");
    // Moving surfaces by at most `surfaces` changes the volume by at most that times the area.
    let (v0, v1) = (measure::volume(&b), measure::volume(&imp.brep));
    assert!(
        (v1 - v0).abs() <= measure::area(&b) * h.surfaces,
        "{v1} vs {v0}"
    );
}

#[test]
fn geometry_healing_may_not_close_is_kept_with_the_tolerances_it_measures() {
    // The rounds grown 0.2 mm: healing would move a surface 0.15 mm, more than it may, so the
    // file's surfaces stay and the round's four edges and four vertices hold the gap, whatever
    // the file states. An end arc of the round takes a curve between its faces instead, out on
    // the grown round mid-way along it: that moves it 0.2 mm, and the report says so.
    let imp = read(with_rounds_grown(&rounded_block(), 0.2).as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert_eq!(imp.healed.len(), 1, "{:?}", imp.healed);
    let h = &imp.healed[0];
    assert!(
        h.what.starts_with("edge #") && (h.moved - 0.2).abs() < 0.005 && h.surfaces == 0.0,
        "{h:?}"
    );
    assert!(check::check(&imp.brep).is_empty());
    assert_eq!((imp.tolerant[0].vertices, imp.tolerant[0].edges), (4, 4));
    assert!(
        (imp.tolerant[0].largest - 0.3).abs() < 1e-6,
        "{:?}",
        imp.tolerant
    );
}

#[test]
fn a_gap_too_wide_to_tell_an_edge_s_ends_apart_is_rejected() {
    // Grown 0.5 mm, the round's edges need tolerances of 0.75 mm, and vertices 1.4 mm apart are
    // then one point within them: the solid has no definite shape, whatever the file states.
    let text = with_uncertainty(&with_rounds_grown(&rounded_block(), 0.5), 5.0);
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert_eq!(imp.brep.solid_ids().count(), 0);
    let why = &imp.rejected[0];
    assert!(
        why.contains("within the tolerance, so the same point")
            || why.contains("too far for an edge"),
        "{why}"
    );
}

/// `text` with its file's stated uncertainty set to `mm`.
fn with_uncertainty(text: &str, mm: f64) -> String {
    let head = "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(";
    let at = text.find(head).unwrap() + head.len();
    let end = at + text[at..].find(')').unwrap();
    format!("{}{mm:?}{}", &text[..at], &text[end..])
}

/// The rounded block's rounds grown by 0.2 mm, past what healing may move (0.05 mm): the round
/// misses the two sides it was tangent to, and its two end arcs, by 0.2 mm.
fn rounds_held_apart(uncertainty: f64) -> super::Import {
    let text = with_uncertainty(&with_rounds_grown(&rounded_block(), 0.2), uncertainty);
    read(text.as_bytes(), "imp").unwrap()
}

#[test]
fn faces_held_apart_within_the_files_uncertainty_keep_local_tolerances() {
    let imp = rounds_held_apart(0.4);
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(imp.healed.is_empty());
    assert!(check::check(&imp.brep).is_empty());
    // The round's two lines and two arcs, and the four vertices they end at, each 0.2 mm off it,
    // get 1.5 times that; every other edge and vertex lies on its faces.
    let (vs, es) = imp.brep.tolerant();
    assert_eq!((vs.len(), es.len()), (4, 4));
    for t in vs.iter().map(|x| x.1).chain(es.iter().map(|x| x.1)) {
        assert!((t - 0.3).abs() < 1e-6, "{t}");
    }
    assert_eq!(imp.tolerant.len(), 1);
    assert_eq!((imp.tolerant[0].vertices, imp.tolerant[0].edges), (4, 4));
    // Faces bounded by edges that miss them by at most 0.3 enclose the block to within its
    // area times that.
    let b = rounded_block();
    let (v0, v1) = (measure::volume(&b), measure::volume(&imp.brep));
    assert!((v1 - v0).abs() <= 0.3 * measure::area(&b), "{v1} vs {v0}");
    // Clearing the tolerances leaves edges that miss their faces.
    let mut exact = imp.brep.clone();
    for (e, _) in es {
        exact.set_edge_tolerance(e, crate::tol::LINEAR);
    }
    for (v, _) in vs {
        exact.set_vertex_tolerance(v, crate::tol::LINEAR);
    }
    assert!(!check::check(&exact).is_empty());
}

#[test]
fn gaps_wider_than_the_file_states_are_kept_and_counted() {
    // The file states 0.15 mm and its round misses by 0.2: the eight tolerances of 0.3 mm are
    // all past what it states, and the report says so.
    let imp = rounds_held_apart(0.15);
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    let t = &imp.tolerant[0];
    assert_eq!((t.vertices, t.edges, t.past_stated), (4, 4, 8));
    assert_eq!(t.stated, Some(0.15));
    // Stated at 0.4 mm, none is past it.
    let within = rounds_held_apart(0.4);
    assert_eq!(within.tolerant[0].past_stated, 0);
}

/// `b` written as another system's STEP with its first straight edge's curve bent into a spline
/// through its two vertices whose middle lies `d` off the edge, along the unit vector (1, 1, 1)
/// leaves across it.
fn with_edge_bent(b: &Brep, d: f64) -> String {
    let text = write(b, "b").unwrap().replace("'springrep'", "'other'");
    let ents: std::collections::HashMap<&str, &str> =
        text.lines().filter_map(|l| l.split_once('=')).collect();
    let refs = |s: &str| -> Vec<String> {
        s.split('#')
            .skip(1)
            .map(|r| {
                format!(
                    "#{}",
                    r.split(|c: char| !c.is_ascii_digit()).next().unwrap()
                )
            })
            .collect()
    };
    let point = |id: &str| -> Vec3 {
        let body = ents[id];
        let inner = &body[body.rfind('(').unwrap() + 1..body.find(')').unwrap()];
        let c: Vec<f64> = inner
            .split(',')
            .map(|x| x.trim().parse().unwrap())
            .collect();
        Vec3::new(c[0], c[1], c[2])
    };
    let (curve, pa, pb) = ents
        .values()
        .filter(|body| body.starts_with("EDGE_CURVE("))
        .find_map(|body| {
            let r = refs(body);
            let geom = ents[r[2].as_str()];
            let curve = match geom.starts_with("SURFACE_CURVE(") {
                true => refs(geom)[0].clone(),
                false => r[2].clone(),
            };
            let vertex = |v: &str| refs(ents[v])[0].clone();
            ents[curve.as_str()]
                .starts_with("LINE(")
                .then(|| (curve, vertex(&r[0]), vertex(&r[1])))
        })
        .unwrap();
    let (a, b) = (point(&pa), point(&pb));
    let along = (b - a).unit_or_zero();
    let diag = Vec3::new(1.0, 1.0, 1.0);
    let off = (diag - along * diag.dot(along)).unit_or_zero();
    // A quadratic's middle lies halfway to its middle control point.
    let m = (a + b) * 0.5 + off * (2.0 * d);
    let old = format!("{curve}={}", ents[curve.as_str()]);
    text.replace(
        &old,
        &format!(
            "{curve}=B_SPLINE_CURVE_WITH_KNOTS('',2,({pa},#99990,{pb}),.UNSPECIFIED.,.F.,.F.,\
             (3,3),(0.,1.),.UNSPECIFIED.);"
        ),
    )
    .replace(
        "ENDSEC;\nEND-ISO",
        &format!(
            "#99990=CARTESIAN_POINT('',({:?},{:?},{:?}));\nENDSEC;\nEND-ISO",
            m.x, m.y, m.z
        ),
    )
}

#[test]
fn an_edge_off_its_faces_past_the_files_uncertainty_takes_a_curve_between_them() {
    // The bent edge's middle lies 0.1 mm off the block's edge along a unit vector across both
    // faces' normals, so off each plane by 0.1/√2, more than the file's 0.05 mm, and healing may
    // move it only 0.05 mm. Between two planes that cross, the curve midway is their line: the
    // block comes back exact, with no local tolerance. Its middle moved the 0.1 mm the bend put
    // it off, square to the line, and the report says so whether or not the file states how
    // closely its geometry meets.
    let b = block("b", [0.0; 3], [10.0, 6.0, 3.0]);
    let bent = with_edge_bent(&b, 0.1);
    for text in [with_uncertainty(&bent, 0.05), bent] {
        let imp = read(text.as_bytes(), "imp").unwrap();
        assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
        assert!(check::check(&imp.brep).is_empty());
        assert_eq!(imp.brep.tolerant(), (vec![], vec![]));
        assert!((measure::volume(&imp.brep) - 180.0).abs() < 1e-9);
        assert_eq!(imp.healed.len(), 1, "{:?}", imp.healed);
        let h = &imp.healed[0];
        assert!(h.what.starts_with("edge #"), "{h:?}");
        assert!((h.moved - 0.1).abs() < 1e-6, "{h:?}");
    }
}

#[test]
fn local_tolerances_are_written_as_the_files_uncertainty_and_read_back() {
    let imp = rounds_held_apart(0.4);
    let text = write(&imp.brep, "rb").unwrap();
    assert!(text.contains("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(0.3"));
    let back = read(text.as_bytes(), "imp").unwrap();
    assert!(back.rejected.is_empty(), "{:?}", back.rejected);
    assert_eq!(back.brep.tolerant(), imp.brep.tolerant());
    let (v0, v1) = (measure::volume(&imp.brep), measure::volume(&back.brep));
    assert!((v1 - v0).abs() <= 1e-9 * v0, "{v1} vs {v0}");
}

#[test]
fn a_solid_kept_with_local_tolerances_reads_back_as_it_was_written() {
    // A vertex 0.03 mm off its corner, in a file that states 0.01: healing would move it further
    // than that, so it is kept where it is. Written, the file states its 0.042 mm tolerance, within
    // which healing could move it; read back, it is still the solid that was written.
    let blk = block("b", [0.0, 0.0, 0.0], [10.0, 4.0, 3.0]);
    let moved = with_vertex_moved(&blk, Vec3::new(0.01, -0.02, 0.02));
    let imp = read(with_uncertainty(&moved, 0.01).as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert_eq!(imp.tolerant.len(), 1, "{:?}", imp.healed);
    let back = read(write(&imp.brep, "b").unwrap().as_bytes(), "imp").unwrap();
    assert!(back.rejected.is_empty(), "{:?}", back.rejected);
    assert!(back.healed.is_empty(), "{:?}", back.healed);
    assert_eq!(back.brep.tolerant(), imp.brep.tolerant());
    assert_eq!(measure::bbox(&back.brep), measure::bbox(&imp.brep));
    for (got, want) in [
        (measure::volume(&back.brep), measure::volume(&imp.brep)),
        (measure::area(&back.brep), measure::area(&imp.brep)),
    ] {
        assert!((got - want).abs() <= 1e-9 * want, "{got} vs {want}");
    }
}

#[test]
fn booleans_pass_tolerant_solids_they_do_not_reach() {
    let imp = rounds_held_apart(0.4).brep;
    let far = block("far", [20.0, 0.0, 0.0], [21.0, 1.0, 1.0]);
    let both = boolean::union(&imp, &far).unwrap();
    assert!(check::check(&both).is_empty());
    assert_eq!(both.tolerant().1.len(), 4);
}

#[test]
fn a_cut_away_from_tolerant_entities_stays_exact() {
    // A pocket in the block's top, 4 mm from the round: its edges and vertices are kernel-made
    // and touch nothing tolerant, so only the round's 4 edges and 4 vertices keep 0.3 mm.
    let imp = rounds_held_apart(0.4).brep;
    let pocket = block("pocket", [2.0, 2.0, 2.0], [4.0, 3.0, 4.0]);
    let cut = boolean::difference(&imp, &pocket).unwrap();
    assert!(check::check(&cut).is_empty());
    let (vs, es) = cut.tolerant();
    assert_eq!((vs.len(), es.len()), (4, 4));
    // The pocket takes 2 x 1 x 1 mm³ of exact material.
    let lost = measure::volume(&imp) - measure::volume(&cut);
    assert!((lost - 2.0).abs() < 1e-9, "{lost}");
}

#[test]
fn a_cut_through_tolerant_edges_gives_its_vertices_their_tolerance() {
    // A slab z 1..2 across the whole block cuts the round's two tolerant lines at z = 1 and 2: the
    // four new vertices on them inherit the lines' 0.3 mm (D. J. Jackson, "Boundary
    // representation modelling with local tolerances", Solid Modeling '95), the lines split into
    // four tolerant pieces, and the two arcs keep theirs. What remains is two solids.
    let imp = rounds_held_apart(0.4).brep;
    let slab = block("slab", [-1.0, -1.0, 1.0], [12.0, 8.0, 2.0]);
    let cut = boolean::difference(&imp, &slab).unwrap();
    assert!(check::check(&cut).is_empty());
    assert_eq!(cut.solid_ids().count(), 2);
    let (vs, es) = cut.tolerant();
    assert_eq!((vs.len(), es.len()), (8, 6));
    for t in vs.iter().map(|x| x.1).chain(es.iter().map(|x| x.1)) {
        assert!((t - 0.3).abs() < 1e-6, "{t}");
    }
    // The slab's part and the rest make up the whole, to within the area its faces bound by
    // tolerant edges can differ by: the slab's two faces cross the 1.2 mm round's 0.2 mm gap.
    let common = boolean::intersection(&imp, &slab).unwrap();
    assert!(check::check(&common).is_empty());
    let whole = measure::volume(&imp);
    let parts = measure::volume(&cut) + measure::volume(&common);
    assert!(
        (parts - whole).abs() <= 2.0 * 0.3 * 1.2,
        "{parts} vs {whole}"
    );
}

#[test]
fn a_tolerant_solid_and_its_copy_within_tolerance_are_one_solid() {
    let imp = rounds_held_apart(0.4).brep;
    let mut copy = imp.clone();
    copy.transform(&crate::math::Transform::translation(Vec3::new(
        1e-9, 0.0, 0.0,
    )));
    let u = boolean::union(&imp, &copy).unwrap();
    assert_eq!(u.face_count(), imp.face_count());
    assert_eq!(u.tolerant(), imp.tolerant());
    let i = boolean::intersection(&imp, &copy).unwrap();
    assert_eq!(measure::volume(&i), measure::volume(&imp));
    let d = boolean::difference(&imp, &copy).unwrap();
    assert_eq!(d.solid_ids().count(), 0);
}

/// A file of `data` lines, from another system.
fn step_file(data: &str) -> String {
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_NAME('a','',(''),(''),'','other','');\nENDSEC;\nDATA;\n\
         {data}\nENDSEC;\nEND-ISO-10303-21;\n"
    )
}

#[test]
fn faceted_solids_read_their_polygons() {
    let corners = [
        (0.0, 0.0, 0.0),
        (2.0, 0.0, 0.0),
        (2.0, 3.0, 0.0),
        (0.0, 3.0, 0.0),
        (0.0, 0.0, 4.0),
        (2.0, 0.0, 4.0),
        (2.0, 3.0, 4.0),
        (0.0, 3.0, 4.0),
    ];
    // Each side counter-clockwise seen from outside, with its outward normal.
    let sides: [([usize; 4], (f64, f64, f64)); 6] = [
        ([0, 3, 2, 1], (0.0, 0.0, -1.0)),
        ([4, 5, 6, 7], (0.0, 0.0, 1.0)),
        ([0, 1, 5, 4], (0.0, -1.0, 0.0)),
        ([2, 3, 7, 6], (0.0, 1.0, 0.0)),
        ([1, 2, 6, 5], (1.0, 0.0, 0.0)),
        ([0, 4, 7, 3], (-1.0, 0.0, 0.0)),
    ];
    let mut d = String::new();
    for (i, c) in corners.iter().enumerate() {
        d += &format!(
            "#{}=CARTESIAN_POINT('',({:?},{:?},{:?}));\n",
            1 + i,
            c.0,
            c.1,
            c.2
        );
    }
    let mut faces = Vec::new();
    for (k, (idx, n)) in sides.iter().enumerate() {
        let b = 100 + 10 * k;
        let o = corners[idx[0]];
        let pts: Vec<String> = idx.iter().map(|i| format!("#{}", 1 + i)).collect();
        let x = if n.0 == 0.0 {
            "(1.,0.,0.)"
        } else {
            "(0.,1.,0.)"
        };
        d += &format!("#{b}=POLY_LOOP('',({}));\n", pts.join(","));
        d += &format!("#{}=FACE_OUTER_BOUND('',#{b},.T.);\n", b + 1);
        d += &format!(
            "#{}=CARTESIAN_POINT('',({:?},{:?},{:?}));\n",
            b + 2,
            o.0,
            o.1,
            o.2
        );
        d += &format!(
            "#{}=DIRECTION('',({:?},{:?},{:?}));\n",
            b + 3,
            n.0,
            n.1,
            n.2
        );
        d += &format!("#{}=DIRECTION('',{x});\n", b + 4);
        d += &format!(
            "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});\n",
            b + 5,
            b + 2,
            b + 3,
            b + 4
        );
        d += &format!("#{}=PLANE('',#{});\n", b + 6, b + 5);
        d += &format!("#{}=FACE_SURFACE('',(#{}),#{},.T.);\n", b + 7, b + 1, b + 6);
        faces.push(format!("#{}", b + 7));
    }
    d += &format!(
        "#200=CLOSED_SHELL('',({}));\n#201=FACETED_BREP('',#200);",
        faces.join(",")
    );
    let imp = read(step_file(&d).as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(check::check(&imp.brep).is_empty());
    assert_eq!(imp.brep.face_count(), 6);
    assert!((measure::volume(&imp.brep) - 24.0).abs() < 1e-12);
}

/// A cone of base radius 6 and height 5 whose sloped face closes on its apex by a vertex loop,
/// as many systems write it, on `surface`.
fn cone_with_vertex_loop(surface: &str) -> String {
    step_file(&format!(
        "#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=CIRCLE('',#4,6.);
#6=CARTESIAN_POINT('',(6.,0.,0.));
#7=VERTEX_POINT('',#6);
#8=EDGE_CURVE('',#7,#7,#5,.T.);
#9=CARTESIAN_POINT('',(0.,0.,5.));
#10=VERTEX_POINT('',#9);
#11=DIRECTION('',(0.,0.,-1.));
#12=AXIS2_PLACEMENT_3D('',#1,#11,#3);
#13={surface};
#14=PLANE('',#12);
#20=ORIENTED_EDGE('',*,*,#8,.F.);
#21=EDGE_LOOP('',(#20));
#22=FACE_OUTER_BOUND('',#21,.T.);
#23=ADVANCED_FACE('',(#22),#14,.T.);
#30=ORIENTED_EDGE('',*,*,#8,.T.);
#31=EDGE_LOOP('',(#30));
#32=FACE_OUTER_BOUND('',#31,.T.);
#33=VERTEX_LOOP('',#10);
#34=FACE_BOUND('',#33,.T.);
#35=ADVANCED_FACE('',(#32,#34),#13,.T.);
#40=CLOSED_SHELL('',(#23,#35));
#41=MANIFOLD_SOLID_BREP('',#40);"
    ))
}

/// A block whose top face has a further bound: a point (a `VERTEX_LOOP`), or a slit, a line inside
/// it the bound runs along and back. Each marks where something touched the face, and is kept.
#[test]
fn a_point_or_a_slit_bounding_a_face_is_kept() {
    let text = write(&block("b", [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]), "part").unwrap();
    let top = text
        .lines()
        .find(|l| l.contains("ADVANCED_FACE('b.end'"))
        .unwrap()
        .to_string();
    let n = text.matches("\n#").count() + 100;
    let with = |extra: String, bound: usize| {
        let (head, tail) = top.split_once("),").unwrap();
        let face = format!("{head},#{bound}),{tail}");
        text.replace(&top, &face)
            .replace("ENDSEC;\nEND-ISO", &format!("{extra}ENDSEC;\nEND-ISO"))
    };
    let point = format!(
        "#{n}=CARTESIAN_POINT('',(5.,5.,10.));\n#{}=VERTEX_POINT('',#{n});\n\
         #{}=VERTEX_LOOP('',#{});\n#{}=FACE_BOUND('',#{},.T.);\n",
        n + 1,
        n + 2,
        n + 1,
        n + 3,
        n + 2
    );
    let slit = format!(
        "#{n}=CARTESIAN_POINT('',(3.,5.,10.));\n#{a}=VERTEX_POINT('',#{n});\n\
         #{b}=CARTESIAN_POINT('',(7.,5.,10.));\n#{c}=VERTEX_POINT('',#{b});\n\
         #{d}=DIRECTION('',(1.,0.,0.));\n#{e}=VECTOR('',#{d},1.);\n#{f}=LINE('',#{n},#{e});\n\
         #{g}=EDGE_CURVE('',#{a},#{c},#{f},.T.);\n#{h}=ORIENTED_EDGE('',*,*,#{g},.T.);\n\
         #{i}=ORIENTED_EDGE('',*,*,#{g},.F.);\n#{j}=EDGE_LOOP('',(#{h},#{i}));\n\
         #{k}=FACE_BOUND('',#{j},.T.);\n",
        a = n + 1,
        b = n + 2,
        c = n + 3,
        d = n + 4,
        e = n + 5,
        f = n + 6,
        g = n + 7,
        h = n + 8,
        i = n + 9,
        j = n + 10,
        k = n + 11
    );
    for (file, vertices, edges) in [(with(point, n + 3), 9, 13), (with(slit, n + 11), 10, 13)] {
        let back = read(file.as_bytes(), "imp").unwrap();
        assert!(back.rejected.is_empty(), "{:?}", back.rejected);
        let b = back.brep;
        assert!(check::check(&b).is_empty(), "{:?}", check::check(&b));
        assert_eq!(b.face_count(), 6);
        assert_eq!(
            (b.vertex_ids().count(), b.edge_ids().count()),
            (vertices, edges)
        );
        assert!((measure::volume(&b) - 1000.0).abs() < 1e-9);
        assert!((measure::area(&b) - 600.0).abs() < 1e-9);
    }
}

/// A ball of radius 5 written as one face bounded only by a `VERTEX_LOOP` at its south pole: a
/// closed surface needs no edge, and the loop names a point of it.
#[test]
fn a_ball_bounded_by_a_vertex_loop_alone_is_read_whole() {
    let text = step_file(
        "#1=CARTESIAN_POINT('',(1.,2.,3.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=SPHERICAL_SURFACE('',#4,5.);
#6=CARTESIAN_POINT('',(1.,2.,-2.));
#7=VERTEX_POINT('',#6);
#8=VERTEX_LOOP('',#7);
#9=FACE_BOUND('',#8,.T.);
#10=ADVANCED_FACE('',(#9),#5,.T.);
#11=CLOSED_SHELL('',(#10));
#12=MANIFOLD_SOLID_BREP('',#11);",
    );
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(check::check(&imp.brep).is_empty());
    let want = 4.0 / 3.0 * PI * 125.0;
    assert!((measure::volume(&imp.brep) - want).abs() < 1e-9 * want);
    assert!((measure::area(&imp.brep) - 100.0 * PI).abs() < 1e-9 * 100.0 * PI);
}

#[test]
fn vertex_loops_close_cones_on_their_apex() {
    let cone = cone_with_vertex_loop("CONICAL_SURFACE('',#12,6.,0.876058050598193)");
    let imp = read(cone.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(check::check(&imp.brep).is_empty());
    let want = 60.0 * PI;
    assert!((measure::volume(&imp.brep) - want).abs() < 1e-9 * want);
    // On a cylinder the vertex is no apex, so it would be a point the face holds, but it is off
    // the surface: refused by name.
    let bad = cone_with_vertex_loop("CYLINDRICAL_SURFACE('',#12,6.)");
    let imp = read(bad.as_bytes(), "imp").unwrap();
    assert_eq!(imp.rejected.len(), 1);
    assert!(
        imp.rejected[0].contains("VERTEX_LOOP"),
        "{:?}",
        imp.rejected
    );
}

/// A solid of revolution about z whose side is one face on `surface` bounded by two circles and no
/// seam, as some systems write a band: the circle of radius `rb` at height `zb` and of `rt` at
/// `zt`, each closed by a disc, with the top circle's vertex turned by `turn` from the bottom's.
fn seamless_band(surface: &str, (zb, rb): (f64, f64), (zt, rt): (f64, f64), turn: f64) -> String {
    let (c, s) = (turn.cos(), turn.sin());
    step_file(&format!(
        "#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=DIRECTION('',(0.,0.,-1.));
#10=CARTESIAN_POINT('',(0.,0.,{zb:?}));
#11=AXIS2_PLACEMENT_3D('',#10,#2,#3);
#12=CIRCLE('',#11,{rb:?});
#13=CARTESIAN_POINT('',({rb:?},0.,{zb:?}));
#14=VERTEX_POINT('',#13);
#15=EDGE_CURVE('',#14,#14,#12,.T.);
#16=AXIS2_PLACEMENT_3D('',#10,#5,#3);
#17=PLANE('',#16);
#20=CARTESIAN_POINT('',(0.,0.,{zt:?}));
#21=AXIS2_PLACEMENT_3D('',#20,#2,#3);
#22=CIRCLE('',#21,{rt:?});
#23=CARTESIAN_POINT('',({:?},{:?},{zt:?}));
#24=VERTEX_POINT('',#23);
#25=EDGE_CURVE('',#24,#24,#22,.T.);
#26=PLANE('',#21);
#30={surface};
#40=ORIENTED_EDGE('',*,*,#15,.F.);
#41=EDGE_LOOP('',(#40));
#42=FACE_OUTER_BOUND('',#41,.T.);
#43=ADVANCED_FACE('',(#42),#17,.T.);
#50=ORIENTED_EDGE('',*,*,#25,.T.);
#51=EDGE_LOOP('',(#50));
#52=FACE_OUTER_BOUND('',#51,.T.);
#53=ADVANCED_FACE('',(#52),#26,.T.);
#60=ORIENTED_EDGE('',*,*,#15,.T.);
#61=EDGE_LOOP('',(#60));
#62=FACE_BOUND('',#61,.T.);
#63=ORIENTED_EDGE('',*,*,#25,.F.);
#64=EDGE_LOOP('',(#63));
#65=FACE_BOUND('',#64,.T.);
#66=ADVANCED_FACE('',(#62,#65),#30,.T.);
#70=CLOSED_SHELL('',(#43,#53,#66));
#71=MANIFOLD_SOLID_BREP('',#70);",
        rt * c,
        rt * s,
    ))
}

/// `text` with a point `p` (a `VERTEX_LOOP`, as a file holds where something touches the face)
/// added to the bounds of the face whose bounds `bounds` lists, `#66=ADVANCED_FACE('',(…),#30,.T.)`.
fn holding(text: &str, bounds: &str, p: Vec3) -> String {
    let face = format!("#66=ADVANCED_FACE('',({bounds}),#30,.T.);");
    assert!(text.contains(&face));
    text.replace(
        &face,
        &format!(
            "#80=CARTESIAN_POINT('',({:?},{:?},{:?}));\n#81=VERTEX_POINT('',#80);\n\
             #82=VERTEX_LOOP('',#81);\n#83=FACE_BOUND('',#82,.T.);\n\
             #66=ADVANCED_FACE('',({bounds},#83),#30,.T.);",
            p.x, p.y, p.z
        ),
    )
}

/// Reads `text`, a solid whose face holds the point `p`, whole and of volume `want`, and meshes
/// it watertight with `corner`, `p` or the vertex `p` lies within the tolerance of, a corner of
/// the mesh. Every mesh point lies on the surfaces and every
/// triangle within the deflection of them, so the mesh's volume is within the area times the
/// deflection of the solid's.
fn meshes_through(text: &str, corner: Vec3, want: f64, name: &str) {
    let opt = crate::mesh::MeshOptions::within(0.01);
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{name}: {:?}", imp.rejected);
    let b = imp.brep;
    assert!(
        check::check(&b).is_empty(),
        "{name}: {:?}",
        check::check(&b)
    );
    assert!(b.loop_ids().any(|lp| b.is_vertex_loop(lp)), "{name}");
    let (area, volume) = measure::area_and_volume(&b);
    assert!((volume - want).abs() < 1e-9 * want, "{name}: {volume}");
    let m = crate::mesh::tessellate(&b, &opt).unwrap_or_else(|e| panic!("{name}: {e}"));
    let s = b.solid_ids().next().unwrap();
    let rep = crate::mesh::report(&b, &m, s);
    assert!(rep.watertight(), "{name}: {rep:?}");
    assert_eq!(rep.degenerate_triangles, 0, "{name}");
    assert!(
        (rep.volume.abs() - want).abs() <= area * opt.linear,
        "{name}: {}",
        rep.volume
    );
    let found = m
        .triangles
        .iter()
        .flatten()
        .any(|&i| m.positions[i as usize].dist(corner) <= 1e-12);
    assert!(found, "{name}: no triangle has {corner:?} as a corner");
}

/// A cylinder of radius 5 from z 0 to 10, and the outer half of a torus R 10 r 2, each with its
/// side bounded by its two circles and no seam, that holds a point at angle `at`. At the circles'
/// vertices' angle the point lies on the seam the side is given, and at z 0 on the cylinder's
/// bottom circle; elsewhere inside. Either way it is read whole, and meshes watertight through
/// the point.
#[test]
fn a_point_bound_on_a_seam_meshes_through_the_point() {
    let cylinder = 250.0 * PI;
    // As in `seamless_bands_of_revolution_are_given_a_seam`.
    let torus = 400.0 * PI + TAU * (10.0 + 8.0 / (3.0 * PI)) * TAU;
    for turn in [0.0, 1.0] {
        let band = seamless_band(
            "CYLINDRICAL_SURFACE('',#4,5.)",
            (0.0, 5.0),
            (10.0, 5.0),
            turn,
        );
        for (at, z) in [(0.0, 5.0), (1.0, 5.0), (PI, 5.0), (0.0, 0.0), (2.0, 0.0)] {
            let p = Vec3::new(5.0 * f64::cos(at), 5.0 * f64::sin(at), z);
            let name = format!("cylinder {turn} {at} {z}");
            meshes_through(&holding(&band, "#62,#65", p), p, cylinder, &name);
        }
        let band = seamless_band(
            "TOROIDAL_SURFACE('',#4,10.,2.)",
            (-2.0, 10.0),
            (2.0, 10.0),
            turn,
        );
        for at in [0.0, 1.0, PI] {
            let p = Vec3::new(12.0 * f64::cos(at), 12.0 * f64::sin(at), 0.0);
            let name = format!("torus {turn} {at}");
            meshes_through(&holding(&band, "#62,#65", p), p, torus, &name);
        }
    }
}

/// The cylinder of `a_point_bound_on_a_seam_meshes_through_the_point` with its seam written, a
/// line from the bottom circle's vertex to the top's, and the point on it, within the linear
/// tolerance of 1e-7 of the bottom vertex, or opposite.
#[test]
fn a_point_bound_on_a_written_seam_meshes_through_the_point() {
    let text = step_file(
        "#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=DIRECTION('',(0.,0.,-1.));
#11=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#12=CIRCLE('',#11,5.);
#13=CARTESIAN_POINT('',(5.,0.,0.));
#14=VERTEX_POINT('',#13);
#15=EDGE_CURVE('',#14,#14,#12,.T.);
#16=AXIS2_PLACEMENT_3D('',#1,#5,#3);
#17=PLANE('',#16);
#20=CARTESIAN_POINT('',(0.,0.,10.));
#21=AXIS2_PLACEMENT_3D('',#20,#2,#3);
#22=CIRCLE('',#21,5.);
#23=CARTESIAN_POINT('',(5.,0.,10.));
#24=VERTEX_POINT('',#23);
#25=EDGE_CURVE('',#24,#24,#22,.T.);
#26=PLANE('',#21);
#27=VECTOR('',#2,1.);
#28=LINE('',#13,#27);
#29=EDGE_CURVE('',#14,#24,#28,.T.);
#30=CYLINDRICAL_SURFACE('',#4,5.);
#40=ORIENTED_EDGE('',*,*,#15,.F.);
#41=EDGE_LOOP('',(#40));
#42=FACE_OUTER_BOUND('',#41,.T.);
#43=ADVANCED_FACE('',(#42),#17,.T.);
#50=ORIENTED_EDGE('',*,*,#25,.T.);
#51=EDGE_LOOP('',(#50));
#52=FACE_OUTER_BOUND('',#51,.T.);
#53=ADVANCED_FACE('',(#52),#26,.T.);
#60=ORIENTED_EDGE('',*,*,#15,.T.);
#61=ORIENTED_EDGE('',*,*,#29,.T.);
#62=ORIENTED_EDGE('',*,*,#25,.F.);
#63=ORIENTED_EDGE('',*,*,#29,.F.);
#64=EDGE_LOOP('',(#60,#61,#62,#63));
#65=FACE_OUTER_BOUND('',#64,.T.);
#66=ADVANCED_FACE('',(#65),#30,.T.);
#70=CLOSED_SHELL('',(#43,#53,#66));
#71=MANIFOLD_SOLID_BREP('',#70);",
    );
    for (at, z) in [(0.0, 5.0), (0.0, 1e-3), (0.0, 5e-8), (PI, 5.0)] {
        let p = Vec3::new(5.0 * f64::cos(at), 5.0 * f64::sin(at), z);
        let corner = if z < 1e-7 {
            Vec3::new(5.0, 0.0, 0.0)
        } else {
            p
        };
        let name = format!("written seam {at} {z}");
        meshes_through(&holding(&text, "#65", p), corner, 250.0 * PI, &name);
    }
}

#[test]
fn seamless_bands_of_revolution_are_given_a_seam() {
    // Cylinder r 4 from z 0 to 5: 80π. Sphere r 5 between z ±3: π∫(25 − z²)dz = 132π. The outer
    // half of a torus R 10 r 2 and the cylinder r 10 inside it: 400π, plus by Pappus the half
    // disc's 2π about its centroid at 10 + 8/3π.
    let torus = 400.0 * PI + TAU * (10.0 + 8.0 / (3.0 * PI)) * TAU;
    let cases = [
        (
            "CYLINDRICAL_SURFACE('',#4,4.)",
            (0.0, 4.0),
            (5.0, 4.0),
            80.0 * PI,
        ),
        (
            "SPHERICAL_SURFACE('',#4,5.)",
            (-3.0, 4.0),
            (3.0, 4.0),
            132.0 * PI,
        ),
        (
            "TOROIDAL_SURFACE('',#4,10.,2.)",
            (-2.0, 10.0),
            (2.0, 10.0),
            torus,
        ),
    ];
    for (surface, bottom, top, want) in cases {
        // Vertices on one line along the surface, and turned apart, which splits an edge.
        for turn in [0.0, 1.0] {
            let imp = read(seamless_band(surface, bottom, top, turn).as_bytes(), "imp").unwrap();
            assert!(
                imp.rejected.is_empty(),
                "{surface} {turn}: {:?}",
                imp.rejected
            );
            assert!(check::check(&imp.brep).is_empty(), "{surface} {turn}");
            let got = measure::volume(&imp.brep);
            assert!(
                (got - want).abs() < 1e-9 * want,
                "{surface} {turn}: {got} vs {want}"
            );
        }
    }
}

/// The part `above` or below the plane at height `z` of a solid on `surface`, as a disc and one
/// face bounded by the circle of radius √(25 − z²) there and no seam. `#4` frames along z at the
/// origin, `#6` along x, and `#11` along z at the circle.
fn cap(surface: &str, z: f64, above: bool) -> String {
    let r = (25.0 - z * z).sqrt();
    let (disc_axis, disc_use, cap_use) = if above {
        ("#16", ".F.", ".T.")
    } else {
        ("#11", ".T.", ".F.")
    };
    step_file(&format!(
        "#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=DIRECTION('',(0.,0.,-1.));
#6=AXIS2_PLACEMENT_3D('',#1,#3,#7);
#7=DIRECTION('',(0.,1.,0.));
#10=CARTESIAN_POINT('',(0.,0.,{z:?}));
#11=AXIS2_PLACEMENT_3D('',#10,#2,#3);
#12=CIRCLE('',#11,{r:?});
#13=CARTESIAN_POINT('',({r:?},0.,{z:?}));
#14=VERTEX_POINT('',#13);
#15=EDGE_CURVE('',#14,#14,#12,.T.);
#16=AXIS2_PLACEMENT_3D('',#10,#5,#3);
#17=PLANE('',{disc_axis});
#30={surface};
#40=ORIENTED_EDGE('',*,*,#15,{disc_use});
#41=EDGE_LOOP('',(#40));
#42=FACE_OUTER_BOUND('',#41,.T.);
#43=ADVANCED_FACE('',(#42),#17,.T.);
#60=ORIENTED_EDGE('',*,*,#15,{cap_use});
#61=EDGE_LOOP('',(#60));
#62=FACE_OUTER_BOUND('',#61,.T.);
#66=ADVANCED_FACE('',(#62),#30,.T.);
#70=CLOSED_SHELL('',(#43,#66));
#71=MANIFOLD_SOLID_BREP('',#70);"
    ))
}

#[test]
fn caps_without_a_seam_close_on_the_pole_or_apex_they_cover() {
    // A cap of height h on a sphere of radius 5: πh²(15 − h)/3. Cut at z 3, the small cap above
    // has h 2, and the rest of the ball, 500π/3, lies below. Cut at z 0, a hemisphere: 250π/3,
    // its circle through the poles of a sphere framed along x. A cone of radius 4 at z 3 and half
    // angle atan(4/3) has its apex at the origin, a drill's point: π·16·3/3 = 16π.
    let (along_z, along_x) = ("SPHERICAL_SURFACE('',#4,5.)", "SPHERICAL_SURFACE('',#6,5.)");
    let cone = format!("CONICAL_SURFACE('',#11,4.,{:?})", (4.0f64 / 3.0).atan());
    let cases = [
        (along_z, 3.0, true, 52.0 * PI / 3.0),
        (along_z, 3.0, false, 448.0 * PI / 3.0),
        (along_x, 0.0, true, 250.0 * PI / 3.0),
        (along_x, 0.0, false, 250.0 * PI / 3.0),
        (cone.as_str(), 3.0, false, 16.0 * PI),
    ];
    for (surface, z, above, want) in cases {
        let imp = read(cap(surface, z, above).as_bytes(), "imp").unwrap();
        assert!(
            imp.rejected.is_empty(),
            "{surface} {z} {above}: {:?}",
            imp.rejected
        );
        assert!(check::check(&imp.brep).is_empty(), "{surface} {z} {above}");
        let got = measure::volume(&imp.brep);
        assert!(
            (got - want).abs() < 1e-9 * want,
            "{surface} {z} {above}: {got} vs {want}"
        );
    }
}

/// The text of `b` with its first cylinder written as an offset by `d` of a cylinder `d`
/// thinner, or of `basis` where given.
fn with_offset_cylinder(b: &Brep, d: f64, basis: Option<&str>) -> String {
    let text = write(b, "c").unwrap().replace("'springrep'", "'other'");
    let l = text
        .lines()
        .find(|l| l.contains("=CYLINDRICAL_SURFACE("))
        .unwrap();
    let (id, rest) = l.split_once('=').unwrap();
    let (head, r) = rest.rsplit_once(',').unwrap();
    let r: f64 = r.trim_end_matches(");").parse().unwrap();
    let inner = match basis {
        Some(s) => s.to_string(),
        None => format!("{head},{:?})", r - d),
    };
    text.replace(l, &format!("{id}=OFFSET_SURFACE('',#99990,{d:?},.F.);"))
        .replace(
            "ENDSEC;\nEND-ISO",
            &format!("#99990={inner};\nENDSEC;\nEND-ISO"),
        )
}

#[test]
fn offset_surfaces_read_as_the_surface_they_are() {
    let rod = extrude(&Profile::circle(Vec2::ZERO, 4.0), &Frame::WORLD, 5.0, "c").unwrap();
    let imp = read(with_offset_cylinder(&rod, 0.5, None).as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(imp.healed.is_empty());
    let want = measure::volume(&rod);
    assert!((measure::volume(&imp.brep) - want).abs() < 1e-9 * want);
    // An offset with no closed form is refused by name, never approximated.
    let ext = "SURFACE_OF_LINEAR_EXTRUSION('',#99991,#99992);\n\
               #99991=LINE('',#99993,#99992);\n#99992=VECTOR('',#99994,1.);\n\
               #99993=CARTESIAN_POINT('',(0.,0.,0.));\n#99994=DIRECTION('',(0.,0.,1.))";
    let imp = read(with_offset_cylinder(&rod, 0.5, Some(ext)).as_bytes(), "imp").unwrap();
    assert_eq!(imp.rejected.len(), 1);
    assert!(
        imp.rejected[0].contains("OFFSET_SURFACE"),
        "{:?}",
        imp.rejected
    );
}

fn mapped_part(part: &Brep, mapping_origin: u64, mapping: &str) -> String {
    let part = write(part, "part").unwrap();
    let data = part
        .split("DATA;\n")
        .nth(1)
        .unwrap()
        .split("ENDSEC;")
        .next()
        .unwrap();
    let rep = data
        .lines()
        .find(|l| l.contains("ADVANCED_BREP_SHAPE_REPRESENTATION"))
        .and_then(|l| l[1..].split('=').next())
        .unwrap();
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_NAME('a','',(''),(''),'','other','');
ENDSEC;
DATA;
         {data}#9000=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));
         #9010=PRODUCT('asm','asm','',(#3));
         #9011=PRODUCT_DEFINITION_FORMATION('','',#9010);
         #9012=PRODUCT_DEFINITION('design','',#9011,#6);
         #9013=PRODUCT_DEFINITION_SHAPE('','',#9012);
         #9014=SHAPE_REPRESENTATION('asm',(#9016),#9099);
         #9015=SHAPE_DEFINITION_REPRESENTATION(#9013,#9014);
         #9096=REPRESENTATION_MAP(#{mapping_origin},#{rep});
         #9020=AXIS2_PLACEMENT_3D('',#9021,#9022,#9023);
         #9021=CARTESIAN_POINT('',(0.,0.,0.));
         #9022=DIRECTION('',(0.,0.,1.));
         #9023=DIRECTION('',(1.,0.,0.));
         #9024=DIRECTION('',(0.,1.,0.));
         #9025=DIRECTION('',(-1.,0.,0.));
         #9026=DIRECTION('',(0.,-1.,0.));
         #9031=CARTESIAN_POINT('',(10.,0.,0.));
         #9099=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#9000))
                REPRESENTATION_CONTEXT('','3D'));
         {mapping}ENDSEC;
END-ISO-10303-21;
"
    )
}

#[test]
fn scaled_transformations_place_scaled_copies() {
    let part = block("b", [0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
    let text = mapped_part(
        &part,
        9020,
        "#9016=MAPPED_ITEM('',#9096,#9018);
         #9018=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',#9023,#9024,#9031,2.,#9022);
",
    );
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    let b = imp.brep;
    // The part itself, and a copy twice its size whose corner sits at x = 10.
    assert_eq!(b.solid_ids().count(), 2);
    assert!((measure::volume(&b) - 9.0 * 6.0).abs() < 1e-9);
    let bb = measure::bbox(&b);
    assert!(
        (bb.max - Vec3::new(12.0, 4.0, 6.0)).max_abs() < 1e-9,
        "{bb:?}"
    );
}

#[test]
fn a_mapped_solid_with_coplanar_axes_is_reported_as_rejected() {
    let part = block("b", [0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
    let text = mapped_part(
        &part,
        9020,
        "#9016=MAPPED_ITEM('',#9096,#9018);
         #9018=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',#9023,#9023,#9031,1.,#9022);
",
    );
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert_eq!(imp.bodies.len(), 1);
    assert_eq!(imp.rejected.len(), 1, "{:?}", imp.rejected);
    let why = &imp.rejected[0];
    assert!(why.starts_with("#9016:"), "{why}");
    assert!(
        why.contains("STEP #9018: a transformation's axes are coplanar"),
        "{why}"
    );
}

fn block_with_off_centre_bore() -> Brep {
    let plate = Profile::rect(Vec2::new(0.5, 1.0), 1.0, 2.0);
    let bore = Profile::circle(Vec2::new(0.25, 0.5), 0.25);
    extrude(&plate.with_holes(bore.loops), &Frame::WORLD, 3.0, "plate").unwrap()
}

fn assert_bored_block(
    body: &Brep,
    solid: crate::topo::SolidId,
    min: Vec3,
    max: Vec3,
    bore_centre: Vec2,
    scale: f64,
) {
    let bounds = measure::solid_bbox(body, solid);
    assert!(bounds.min.dist(min) <= crate::tol::LINEAR, "{bounds:?}");
    assert!(bounds.max.dist(max) <= crate::tol::LINEAR, "{bounds:?}");

    let bore_radius = 0.125;
    let cylinders: Vec<_> = body
        .solid_faces(solid)
        .into_iter()
        .filter_map(|f| match &body.face(f).surface {
            Surface::Cylinder(cylinder) => Some(cylinder),
            _ => None,
        })
        .collect();
    assert_eq!(cylinders.len(), 1);
    let bore = cylinders[0];
    assert!((bore.frame.origin.x - bore_centre.x).abs() <= crate::tol::LINEAR);
    assert!((bore.frame.origin.y - bore_centre.y).abs() <= crate::tol::LINEAR);
    assert!((bore.radius - bore_radius * scale).abs() <= crate::tol::LINEAR);
    let volume = measure::solid_volume(body, solid);
    let expected = 3.0 * (2.0 - PI * bore_radius * bore_radius) * scale.powi(3);
    assert!(
        (volume - expected).abs() <= crate::tol::LINEAR * 22.0 * scale.powi(2),
        "{volume}"
    );
}

fn assert_reflected_bore(axis3: &str) {
    let part = block_with_off_centre_bore();
    let text = mapped_part(
        &part,
        9020,
        &format!(
            "#9016=MAPPED_ITEM('',#9096,#9018);
             #9018=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',#9025,#9024,#9031,2.,{axis3});
"
        ),
    );
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert_eq!(imp.bodies.len(), 2);
    let b = imp.brep;
    // The target maps (x, y, z) to (10 − 2x, 2y, 2z).
    assert_bored_block(
        &b,
        imp.bodies[1].solid,
        Vec3::new(8.0, 0.0, 0.0),
        Vec3::new(10.0, 4.0, 6.0),
        Vec2::new(9.5, 1.0),
        2.0,
    );
    assert!(check::check(&b).is_empty());
}

#[test]
fn a_mapped_solid_reflects_its_off_centre_bore() {
    assert_reflected_bore("#9022");
}

#[test]
fn a_reflected_mapping_defaults_its_third_axis_to_positive_z() {
    assert_reflected_bore("$");
}

#[test]
fn nested_reflections_compose_around_a_rotated_mapping_origin() {
    let part = block_with_off_centre_bore();
    let text = mapped_part(
        &part,
        9047,
        "#9016=MAPPED_ITEM('',#9040,#9041);
         #9040=REPRESENTATION_MAP(#9020,#9042);
         #9041=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',#9023,#9026,#9043,0.5,#9022);
         #9042=SHAPE_REPRESENTATION('nested',(#9044),#9099);
         #9043=CARTESIAN_POINT('',(0.,10.,2.));
         #9044=MAPPED_ITEM('',#9096,#9046);
         #9046=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',#9025,#9024,#9031,2.,#9022);
         #9047=AXIS2_PLACEMENT_3D('',#9048,#9022,#9024);
         #9048=CARTESIAN_POINT('',(1.,2.,3.));
",
    );
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert_eq!(imp.bodies.len(), 2);
    let b = imp.brep;
    assert_bored_block(
        &b,
        imp.bodies[0].solid,
        Vec3::ZERO,
        Vec3::new(1.0, 2.0, 3.0),
        Vec2::new(0.25, 0.5),
        1.0,
    );
    // The inner map is (14 − 2y, 2 − 2x, 2z − 6). The outer map scales by 1/2,
    // reflects Y, and translates by (0, 10, 2), giving (7 − y, 9 + x, z − 1).
    assert_bored_block(
        &b,
        imp.bodies[1].solid,
        Vec3::new(5.0, 9.0, -1.0),
        Vec3::new(7.0, 10.0, 2.0),
        Vec2::new(6.5, 9.25),
        1.0,
    );
    assert!(check::check(&b).is_empty());
}

#[test]
fn a_scaled_reflection_and_its_inverse_compose_to_identity() {
    let file = super::part21::parse(
        b"ISO-10303-21; DATA;
        #1=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',#2,#3,#5,2.,#4);
        #2=DIRECTION('',(0.,1.,0.));
        #3=DIRECTION('',(1.,0.,0.));
        #4=DIRECTION('',(0.,0.,1.));
        #5=CARTESIAN_POINT('',(1.,2.,3.));
        ENDSEC; END-ISO-10303-21;",
    )
    .unwrap();
    let placing = super::geometry::Geometry::new(&file)
        .transform_of(1)
        .unwrap();
    let identity = super::geometry::Placing::IDENTITY;
    assert_eq!(placing.then(&placing.inverse()), identity);
    assert_eq!(placing.inverse().then(&placing), identity);
}

#[test]
fn a_transformation_operator_rejects_a_second_axis_in_the_first_and_third_axes_plane() {
    let file = super::part21::parse(
        b"ISO-10303-21; DATA;
        #1=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',#2,#3,#5,1.,#4);
        #2=DIRECTION('',(1.,0.,0.));
        #3=DIRECTION('',(1.,0.,1.));
        #4=DIRECTION('',(0.,0.,1.));
        #5=CARTESIAN_POINT('',(0.,0.,0.));
        ENDSEC; END-ISO-10303-21;",
    )
    .unwrap();
    assert!(matches!(
        super::geometry::Geometry::new(&file).transform_of(1),
        Err(crate::Error::InvalidInput(message)) if message.contains("STEP #1:")
    ));
}

#[test]
fn a_transformation_operator_defaults_its_first_axis_to_y_when_its_third_axis_is_x() {
    let file = super::part21::parse(
        b"ISO-10303-21; DATA;
        #1=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',$,#2,#4,1.,#3);
        #2=DIRECTION('',(0.,0.,1.));
        #3=DIRECTION('',(1.,0.,0.));
        #4=CARTESIAN_POINT('',(0.,0.,0.));
        ENDSEC; END-ISO-10303-21;",
    )
    .unwrap();
    assert!(
        super::geometry::Geometry::new(&file)
            .transform_of(1)
            .is_ok()
    );
}

#[test]
fn touching_placed_solids_check_valid_each_but_not_as_one() {
    let mut b = block("a", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    b.append(&block("b", [1.0, 0.0, 0.0], [2.0, 1.0, 1.0]));
    let text = write(&b, "board").unwrap();
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert_eq!(imp.brep.solid_ids().count(), 2);
    assert!(!check::check(&imp.brep).is_empty());
    assert!(check::check_each(&imp.brep).is_empty());
}

/// Where a pcurve puts the edge at its parameter `t`, on the face's surface as STEP parameterizes
/// it: a cone by height along its STEP axis, every other surface as the kernel does.
fn on_step_surface(s: &Surface, c2: &super::write::Curve2, t: f64) -> Vec3 {
    let uv = match c2 {
        super::write::Curve2::Line { origin, dir } => *origin + *dir * t,
        super::write::Curve2::Conic {
            center,
            x,
            major,
            minor,
        } => *center + *x * (major * t.cos()) + x.perp() * (minor * t.sin()),
        super::write::Curve2::Spline {
            degree,
            knots,
            ctrl,
            weights,
        } => {
            let pts: Vec<Vec3> = ctrl.iter().map(|p| Vec3::new(p.x, p.y, 0.0)).collect();
            let c = crate::geom::NurbsCurve::new(*degree, knots.clone(), pts, weights.clone());
            let p = c.unwrap().eval(t);
            Vec2::new(p.x, p.y)
        }
    };
    match s {
        Surface::Cone(c) => {
            let (fr, a, _) = super::write::step_cone(c).unwrap();
            let (u, h) = (uv.x, uv.y);
            fr.origin + (fr.x * u.cos() + fr.y * u.sin()) * (c.radius + h * a.tan()) + fr.z * h
        }
        _ => s.eval(uv.x, uv.y),
    }
}

/// Every coedge's pcurve, read on its face's surface as STEP parameterizes it, follows its edge
/// within the tolerance, and each edge is written with one per face, a seam with both.
fn assert_pcurves_follow(b: &Brep) {
    for c in b.coedge_ids() {
        let edge = b.edge(b.coedge(c).edge);
        let Some(curve) = &edge.curve else { continue };
        let s = &b.face(b.coedge_face(c)).surface;
        let c2 = super::write::pcurve(b, c).unwrap_or_else(|| panic!("{} has none", s.kind()));
        let (t0, t1) = edge.range;
        for i in 0..=200 {
            let t = t0 + (t1 - t0) * i as f64 / 200.0;
            let d = on_step_surface(s, &c2, t).dist(curve.eval(t));
            assert!(d <= crate::tol::LINEAR, "{} {d:e}", s.kind());
        }
    }
    let text = write(b, "part").unwrap();
    let edges = b.edge_ids().filter(|&e| !b.edge(e).is_degenerate()).count();
    let with = text.matches("=SURFACE_CURVE(").count() + text.matches("=SEAM_CURVE(").count();
    assert_eq!(with, edges);
    assert_eq!(text.matches("=PCURVE(").count(), 2 * edges);
}

#[test]
fn pcurves_follow_their_edges_on_every_surface() {
    // A frustum narrowing up its axis (written about the reversed axis) and one widening, each cut
    // by a tilted block, so their cones carry traced conic sections as well as their own circles.
    for pts in [
        [(0.0, 0.0), (6.0, 0.0), (3.0, 5.0), (0.0, 5.0)],
        [(0.0, 0.0), (3.0, 0.0), (6.0, 5.0), (0.0, 5.0)],
    ] {
        let pts: Vec<Vec2> = pts.iter().map(|&(x, y)| Vec2::new(x, y)).collect();
        let p = Profile::polygon(&pts, &["a", "b", "c", "d"]);
        let cone = revolve(&p, &xz(), Vec3::ZERO, Vec3::Z, TAU, "t").unwrap();
        let mut cut = block("k", [0.0, -10.0, -10.0], [20.0, 10.0, 10.0]);
        cut.transform(&crate::math::Transform::rotation(
            Vec3::new(2.0, 0.0, 2.5),
            Vec3::new(0.3, 1.0, 0.2).unit().unwrap(),
            0.4,
        ));
        let d = boolean::difference(&cone, &cut).unwrap();
        assert_pcurves_follow(&d);
        round_trip(&d);
    }
    // A twisted loft's sides crossed by a rod: traced curves on rational B-spline sides, and the
    // loft's own edges on its plane caps.
    let sq = Profile::rect(Vec2::ZERO, 6.0, 6.0);
    let sections: Vec<crate::loft::Section> = (0..2)
        .map(|i| {
            let th = 0.6 * i as f64;
            crate::loft::Section {
                profile: sq.clone(),
                frame: Frame::new(
                    Vec3::new(0.0, 0.0, 8.0 * i as f64),
                    Vec3::Z,
                    Vec3::new(th.cos(), th.sin(), 0.0),
                )
                .unwrap(),
            }
        })
        .collect();
    let loft = crate::loft::loft(&sections, false, "l").unwrap();
    let rod = extrude(
        &Profile::circle(Vec2::new(4.0, 1.0), 3.0),
        &Frame::new(Vec3::new(-10.0, 0.0, 0.0), Vec3::X, Vec3::Y).unwrap(),
        20.0,
        "r",
    )
    .unwrap();
    assert_pcurves_follow(&boolean::difference(&loft, &rod).unwrap());
}

/// A written 2D line, circle or ellipse as a function of its parameter; `None` for any other curve.
fn flat_curve(file: &super::part21::File, id: u64) -> Option<Box<dyn Fn(f64) -> Vec2>> {
    let args = |id: u64, kind: &str| file.get(id).unwrap().part(kind).unwrap().to_vec();
    let pair = |id: u64, kind: &str| {
        let xs = args(id, kind)[1].as_list().unwrap().to_vec();
        assert_eq!(xs.len(), 2);
        Vec2::new(xs[0].as_f64().unwrap(), xs[1].as_f64().unwrap())
    };
    let placed = |id: u64| {
        let a = args(id, "AXIS2_PLACEMENT_2D");
        let x = pair(a[2].as_ref().unwrap(), "DIRECTION");
        assert!((x.norm() - 1.0).abs() < 1e-12);
        (pair(a[1].as_ref().unwrap(), "CARTESIAN_POINT"), x)
    };
    let e = file.get(id).unwrap();
    match e.kind() {
        "LINE" => {
            let a = args(id, "LINE");
            let o = pair(a[1].as_ref().unwrap(), "CARTESIAN_POINT");
            let v = args(a[2].as_ref().unwrap(), "VECTOR");
            let d = pair(v[1].as_ref().unwrap(), "DIRECTION") * v[2].as_f64().unwrap();
            Some(Box::new(move |t| o + d * t))
        }
        "CIRCLE" | "ELLIPSE" => {
            let a = args(id, e.kind());
            let (c, x) = placed(a[1].as_ref().unwrap());
            let r1 = a[2].as_f64().unwrap();
            let r2 = a.get(3).map_or(r1, |p| p.as_f64().unwrap());
            Some(Box::new(move |t| {
                c + x * (r1 * t.cos()) + x.perp() * (r2 * t.sin())
            }))
        }
        _ => None,
    }
}

/// How many pcurves of each kind `text` holds, a spline under a conic that turns against its
/// plane's parameters counted as `against`. Each that is a line or a conic on an elementary
/// surface, under a line or a conic in space, is read back from the text and lifted onto its
/// surface as STEP parameterizes it, where it must lie on the curve in space at equal parameters.
fn written_pcurves(text: &str) -> std::collections::BTreeMap<String, usize> {
    let file = super::part21::parse(text.as_bytes()).unwrap();
    let mut g = super::geometry::Geometry::new(&file);
    let mut kinds = std::collections::BTreeMap::new();
    for (&id, e) in file.ids.iter().zip(&file.entities) {
        if !matches!(e.kind(), "SURFACE_CURVE" | "SEAM_CURVE") {
            continue;
        }
        let a = e.part(e.kind()).unwrap();
        let space = g.curve(id).unwrap();
        for pc in a[2].as_list().unwrap() {
            let pc = file.get(pc.as_ref().unwrap()).unwrap();
            let pc = pc.part("PCURVE").unwrap();
            let surface = g.surface(pc[1].as_ref().unwrap()).unwrap();
            let rep = file.get(pc[2].as_ref().unwrap()).unwrap();
            let rep = rep.part("DEFINITIONAL_REPRESENTATION").unwrap();
            let c2 = rep[1].as_list().unwrap()[0].as_ref().unwrap();
            let against = match (&space, &surface) {
                (Curve::Circle(c), Surface::Plane(p)) => c.frame.z.dot(p.frame.z) < 0.0,
                (Curve::Ellipse(c), Surface::Plane(p)) => c.frame.z.dot(p.frame.z) < 0.0,
                _ => false,
            };
            let kind = match file.get(c2).unwrap().kind() {
                "" => "rational",
                "B_SPLINE_CURVE_WITH_KNOTS" if against => "against",
                k => k,
            };
            *kinds.entry(kind.to_string()).or_insert(0) += 1;
            let lift: Box<dyn Fn(Vec2) -> Vec3> = match &surface {
                Surface::Cone(c) => {
                    let c = c.clone();
                    Box::new(move |uv| {
                        let fr = c.frame;
                        let r = c.radius + uv.y * c.half_angle.tan();
                        fr.origin + (fr.x * uv.x.cos() + fr.y * uv.x.sin()) * r + fr.z * uv.y
                    })
                }
                Surface::Plane(_)
                | Surface::Cylinder(_)
                | Surface::Sphere(_)
                | Surface::Torus(_) => {
                    let s = surface.clone();
                    Box::new(move |uv| s.eval(uv.x, uv.y))
                }
                _ => continue,
            };
            if !matches!(space, Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_)) {
                continue;
            }
            let Some(flat) = flat_curve(&file, c2) else {
                continue;
            };
            for i in -4..=12 {
                let t = 0.7 * i as f64;
                let d = lift(flat(t)).dist(space.eval(t));
                assert!(
                    d <= crate::tol::LINEAR,
                    "{kind} on {} {d:e}",
                    surface.kind()
                );
            }
        }
    }
    kinds
}

/// `b`'s pcurves by kind as written, after checking that each follows its edge and that the
/// solid reads back.
fn pcurve_kinds(b: &Brep) -> Vec<(String, usize)> {
    assert_pcurves_follow(b);
    round_trip(b);
    written_pcurves(&write(b, "part").unwrap())
        .into_iter()
        .collect()
}

fn kinds(want: &[(&str, usize)]) -> Vec<(String, usize)> {
    want.iter().map(|&(k, n)| (k.to_string(), n)).collect()
}

/// Every pcurve of `b` is written as a line or a circle, but a circle's that turns against its
/// plane's parameters, which has no circle there with the edge's own parameter.
fn lines_and_circles_only(b: &Brep) {
    let got = pcurve_kinds(b);
    let exact = |k: &str| matches!(k, "CIRCLE" | "LINE" | "against");
    assert!(got.iter().all(|(k, _)| exact(k)), "{got:?}");
    assert!(got.iter().any(|(k, _)| k == "LINE"), "{got:?}");
}

fn frustum(pts: [(f64, f64); 4]) -> Brep {
    let pts: Vec<Vec2> = pts.iter().map(|&(x, y)| Vec2::new(x, y)).collect();
    let p = Profile::polygon(&pts, &["a", "b", "c", "d"]);
    revolve(&p, &xz(), Vec3::ZERO, Vec3::Z, TAU, "t").unwrap()
}

#[test]
fn circles_and_rulings_of_swept_solids_are_written_as_exact_pcurves() {
    // Two rims, each a circle on its cap and a line on the side, and the seam, a line twice.
    let rod = extrude(
        &Profile::circle(Vec2::new(1.0, 2.0), 10.0),
        &Frame::WORLD,
        7.0,
        "rod",
    );
    let rims_and_seam = kinds(&[("CIRCLE", 2), ("LINE", 4)]);
    assert_eq!(pcurve_kinds(&rod.unwrap()), rims_and_seam);
    for pts in [
        [(0.0, 0.0), (6.0, 0.0), (3.0, 5.0), (0.0, 5.0)],
        [(0.0, 0.0), (3.0, 0.0), (6.0, 5.0), (0.0, 5.0)],
    ] {
        assert_eq!(pcurve_kinds(&frustum(pts)), rims_and_seam);
    }
}

#[test]
fn circles_and_lines_a_boolean_leaves_on_surfaces_of_revolution_are_exact_pcurves() {
    // Two balls on one axis meet in a circle of latitude of each.
    let mut small = ball(4.0);
    small.transform(&crate::math::Transform::translation(Vec3::new(
        0.0, 0.0, 6.0,
    )));
    lines_and_circles_only(&boolean::union(&ball(5.0), &small).unwrap());
    // A ring halved across its axis, leaving circles of latitude that are circles on the cut too,
    // and one halved through its axis, leaving circles of its tube.
    let ring = Profile::circle(Vec2::new(8.0, 0.0), 4.0);
    let ring = revolve(&ring, &xz(), Vec3::ZERO, Vec3::Z, TAU, "ring").unwrap();
    let (above, _) = boolean::split(&ring, Vec3::new(0.0, 0.0, 0.5), Vec3::Z).unwrap();
    lines_and_circles_only(&above);
    let (half, _) = boolean::split(&ring, Vec3::ZERO, Vec3::Y).unwrap();
    lines_and_circles_only(&half);
    // A block bored through, the bore's rims circles on the block's faces.
    let bore = extrude(
        &Profile::circle(Vec2::new(1.0, 0.5), 3.0),
        &Frame::new(Vec3::new(0.0, 0.0, -1.0), Vec3::Z, Vec3::X).unwrap(),
        6.0,
        "bore",
    )
    .unwrap();
    let plate = block("b", [-5.0, -5.0, 0.0], [5.0, 5.0, 4.0]);
    lines_and_circles_only(&boolean::difference(&plate, &bore).unwrap());
    // A cone and a rod, each cut square across its axis and along it.
    let rod = extrude(&Profile::circle(Vec2::ZERO, 4.0), &Frame::WORLD, 5.0, "rod").unwrap();
    for b in [
        &frustum([(0.0, 0.0), (6.0, 0.0), (3.0, 5.0), (0.0, 5.0)]),
        &rod,
    ] {
        let (_, below) = boolean::split(b, Vec3::new(0.0, 0.0, 2.0), Vec3::Z).unwrap();
        lines_and_circles_only(&below);
        let (_, side) = boolean::split(b, Vec3::ZERO, Vec3::new(0.6, 0.8, 0.0)).unwrap();
        lines_and_circles_only(&side);
    }
}

/// The rate [`super::write::steady_rate`] reads for `curve` on `s` from its point at `t = 0.3`, to
/// which the surface's own projection gives the parameters.
fn rate_on(s: &Surface, curve: &Curve) -> Option<Vec2> {
    let uv = s.project(curve.eval(0.3));
    super::write::steady_rate(s, curve, (0.3, uv), 0.9)
}

// Each curve refused lies close to a parameter line of its surface but strays from it by more
// than the 0.45 of the tolerance a pcurve keeps to, by its tilt across its size or by its centre.
#[test]
fn a_curve_close_to_a_parameter_line_but_not_one_is_not_written_as_a_line() {
    let (lin, ang) = (crate::tol::LINEAR, crate::tol::ANGULAR);
    let circle = |origin: Vec3, z: Vec3, radius: f64| {
        let frame = Frame::from_axis(origin, z.unit().unwrap()).unwrap();
        Curve::Circle(Circle { frame, radius })
    };
    let line = |origin: Vec3, dir: Vec3| {
        Curve::Line(Line {
            origin,
            dir: dir.unit().unwrap(),
        })
    };
    let up = Vec3::new(0.0, 0.0, 2.0);
    // In radians: 5e-7 mm across a radius of 5 mm, and 6e-8 mm along the 0.6 mm of a line read.
    let lean = 1e-7;
    let tilt = Vec3::new(lean, 0.0, 1.0);
    let aside = Vec3::new(4.0 * lin, 0.0, 2.0);

    let rod = Surface::Cylinder(crate::geom::Cylinder {
        frame: Frame::WORLD,
        radius: 5.0,
    });
    let along = Some(Vec2::new(1.0, 0.0));
    assert_eq!(rate_on(&rod, &circle(up, Vec3::Z, 5.0)), along);
    assert_eq!(rate_on(&rod, &circle(up, -Vec3::Z, 5.0)), along.map(|d| -d));
    assert_eq!(rate_on(&rod, &circle(up, tilt, 5.0)), None);
    // A tilt of the angular tolerance strays 5e-10 mm across this radius and 1e-6 mm across 10 m.
    let slight = Vec3::new(ang, 0.0, 1.0);
    assert_eq!(rate_on(&rod, &circle(up, slight, 5.0)), along);
    let mast = Surface::Cylinder(crate::geom::Cylinder {
        frame: Frame::WORLD,
        radius: 1e4,
    });
    assert_eq!(rate_on(&mast, &circle(up, Vec3::Z, 1e4)), along);
    assert_eq!(rate_on(&mast, &circle(up, slight, 1e4)), None);
    assert_eq!(rate_on(&rod, &circle(aside, Vec3::Z, 5.0)), None);
    let foot = Vec3::new(5.0, 0.0, 0.0);
    assert_eq!(
        rate_on(&rod, &line(foot, Vec3::Z)),
        Some(Vec2::new(0.0, 1.0))
    );
    assert_eq!(
        rate_on(&rod, &line(foot, -Vec3::Z)),
        Some(Vec2::new(0.0, -1.0))
    );
    assert_eq!(rate_on(&rod, &line(foot, Vec3::new(0.0, lean, 1.0))), None);

    let cone = Surface::Cone(Cone {
        frame: Frame::WORLD,
        radius: 5.0,
        half_angle: 0.5,
    });
    let slant = Vec3::new(0.5f64.sin(), 0.0, 0.5f64.cos());
    let rim = 5.0 + 2.0 * 0.5f64.tan();
    assert_eq!(rate_on(&cone, &circle(up, Vec3::Z, rim)), along);
    assert_eq!(rate_on(&cone, &circle(up, tilt, rim)), None);
    assert_eq!(
        rate_on(&cone, &line(foot, slant)),
        Some(Vec2::new(0.0, 1.0))
    );
    let off_slant = slant + Vec3::new(0.0, lean, 0.0);
    assert_eq!(rate_on(&cone, &line(foot, off_slant)), None);

    let globe = Surface::Sphere(crate::geom::Sphere {
        frame: Frame::WORLD,
        radius: 5.0,
    });
    assert_eq!(rate_on(&globe, &circle(up, Vec3::Z, 21f64.sqrt())), along);
    assert_eq!(rate_on(&globe, &circle(up, tilt, 21f64.sqrt())), None);
    assert_eq!(rate_on(&globe, &circle(aside, Vec3::Z, 21f64.sqrt())), None);
    let meridian = Frame::new(Vec3::ZERO, -Vec3::Y, Vec3::X).unwrap();
    let great = |frame: Frame| Curve::Circle(Circle { frame, radius: 5.0 });
    assert_eq!(rate_on(&globe, &great(meridian)), Some(Vec2::new(0.0, 1.0)));
    let leaning = Frame::from_axis(Vec3::ZERO, Vec3::new(0.0, -1.0, lean).unit().unwrap());
    assert_eq!(rate_on(&globe, &great(leaning.unwrap())), None);
    // The meridian itself, on an edge that runs on over the pole.
    let uv = globe.project(great(meridian).eval(0.3));
    let over = super::write::steady_rate(&globe, &great(meridian), (0.3, uv), 2.0);
    assert_eq!(over, None);

    let ring = Surface::Torus(crate::geom::Torus {
        frame: Frame::WORLD,
        major: 8.0,
        minor: 3.0,
    });
    assert_eq!(
        rate_on(&ring, &circle(up, Vec3::Z, 8.0 + 5f64.sqrt())),
        along
    );
    assert_eq!(rate_on(&ring, &circle(up, tilt, 8.0 + 5f64.sqrt())), None);
    let tube = |origin: Vec3, z: Vec3| {
        let frame = Frame::new(origin, z.unit().unwrap(), Vec3::X).unwrap();
        Curve::Circle(Circle { frame, radius: 3.0 })
    };
    let hub = Vec3::new(8.0, 0.0, 0.0);
    assert_eq!(
        rate_on(&ring, &tube(hub, -Vec3::Y)),
        Some(Vec2::new(0.0, 1.0))
    );
    assert_eq!(rate_on(&ring, &tube(hub, Vec3::new(lean, -1.0, 0.0))), None);
    assert_eq!(
        rate_on(&ring, &tube(hub + Vec3::X * (4.0 * lin), -Vec3::Y)),
        None
    );

    // Two balls whose centres lie two tolerances off one axis meet in a circle that strays from a
    // latitude of each by about that much.
    let mut small = ball(4.0);
    small.transform(&crate::math::Transform::translation(Vec3::new(
        2e-7, 0.0, 6.0,
    )));
    let got = pcurve_kinds(&boolean::union(&ball(5.0), &small).unwrap());
    assert_eq!(got, kinds(&[("B_SPLINE_CURVE_WITH_KNOTS", 2), ("LINE", 4)]));
}

#[test]
fn an_ellipse_on_a_plane_is_written_as_an_ellipse_pcurve() {
    let rod = extrude(&Profile::circle(Vec2::ZERO, 4.0), &Frame::WORLD, 9.0, "rod").unwrap();
    let up = Vec3::new(0.3, 0.1, 1.0).unit().unwrap();
    let (_, below) = boolean::split(&rod, Vec3::new(0.0, 0.0, 5.0), up).unwrap();
    let got = pcurve_kinds(&below);
    assert!(got.contains(&("ELLIPSE".to_string(), 1)), "{got:?}");
}

// Fitted pcurves on a ball and on a rod, each crossed obliquely by another rod: read at three
// points a span, the surface over them strayed 5.33e-8 and 5.35e-8 mm from their edges between
// those points, past the half tolerance the writer promises.
#[test]
fn fitted_pcurves_stay_within_half_the_tolerance_between_their_checks() {
    let rod = |base: [f64; 3], axis: [f64; 3], r: f64| {
        let f = Frame::from_axis(Vec3::from_array(base), Vec3::from_array(axis)).unwrap();
        extrude(&Profile::circle(Vec2::ZERO, 2.0 * r), &f, 40.0, "r").unwrap()
    };
    let shapes = [
        boolean::difference(&ball(6.0), &rod([1.0, -20.0, -15.0], [0.1, 1.0, 0.8], 2.5)),
        boolean::union(
            &rod([0.0, 0.0, -20.0], [0.0, 0.0, 1.0], 4.0),
            &rod([-20.0, 0.5, 0.0], [1.0, 0.0, 0.3], 3.0),
        ),
    ];
    for b in shapes.iter().map(|b| b.as_ref().unwrap()) {
        for c in b.coedge_ids() {
            let edge = b.edge(b.coedge(c).edge);
            let Some(curve) = &edge.curve else { continue };
            let s = &b.face(b.coedge_face(c)).surface;
            let c2 = super::write::pcurve(b, c).unwrap();
            let (t0, t1) = edge.range;
            for i in 0..=1000 {
                let t = t0 + (t1 - t0) * i as f64 / 1000.0;
                let d = on_step_surface(s, &c2, t).dist(curve.eval(t));
                assert!(d <= 0.5 * crate::tol::LINEAR, "{} {d:e}", s.kind());
            }
        }
    }
}

/// Two 2 mm blocks overlapping by 1 mm each way, written as two solids, then made one solid whose
/// one shell holds both blocks' faces: it passes through itself where their faces cross, the
/// faces reaching 1 mm either side of each other.
fn blocks_in_one_shell() -> String {
    let mut b = block("a", [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    b.append(&block("b", [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]));
    let text = write(&b, "board").unwrap();
    let shells: Vec<&str> = text
        .lines()
        .filter(|l| l.contains("=CLOSED_SHELL("))
        .collect();
    assert_eq!(shells.len(), 2);
    let id = |l: &str| l[..l.find('=').unwrap()].to_string();
    let faces = |l: &str| l[l.find("(#").unwrap() + 1..l.rfind("))").unwrap()].to_string();
    let merged = shells[0].replace(
        &faces(shells[0]),
        &format!("{},{}", faces(shells[0]), faces(shells[1])),
    );
    let second = text
        .lines()
        .find(|l| l.contains(&format!("MANIFOLD_SOLID_BREP('',{})", id(shells[1]))))
        .unwrap();
    let one: String = text
        .lines()
        .filter(|&l| l != shells[1] && l != second)
        .map(|l| {
            let l = if l == shells[0] {
                merged.clone()
            } else {
                l.to_string()
            };
            l.replace(&format!(",{}", id(second)), "") + "\n"
        })
        .collect();
    one
}

#[test]
fn a_solid_passing_through_itself_is_kept_with_where_and_how_deep() {
    let imp = read(blocks_in_one_shell().as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert_eq!(imp.brep.solid_ids().count(), 1);
    assert!(check::check_structure(&imp.brep).is_empty());
    assert!(!check::check(&imp.brep).is_empty());
    assert!(!imp.crossed.is_empty());
    let deepest = imp.crossed.iter().map(|x| x.depth).fold(0.0, f64::max);
    assert!((deepest - 1.0).abs() < 1e-9, "{deepest}");
    for x in &imp.crossed {
        assert_eq!((x.index, x.step_id), (0, imp.crossed[0].step_id));
        // Every crossing lies on the overlap's box.
        assert!(
            (0..3).all(|k| x.at[k] >= 1.0 - 1e-9 && x.at[k] <= 2.0 + 1e-9),
            "{:?}",
            x.at
        );
    }
}

#[test]
fn a_tessellated_solid_reads_as_the_polyhedron_of_its_triangles() {
    // The corner of the unit cube cut off by x + y + z = 1, its four triangles counter-clockwise
    // seen from outside: volume 1/6.
    let text = "ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('t','',(''),(''),'','','');
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF'));
ENDSEC;
DATA;
#1=COORDINATES_LIST('',4,((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.)));
#2=TRIANGULATED_FACE('',#1,0,(),$,(),((1,3,2),(1,2,4),(1,4,3),(2,3,4)));
#3=TESSELLATED_SOLID('',(#2),$);
ENDSEC;
END-ISO-10303-21;
";
    let imp = read(text.as_bytes(), "t").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(check::check(&imp.brep).is_empty());
    assert!((measure::volume(&imp.brep) - 1.0 / 6.0).abs() < 1e-12);
    assert_eq!(imp.brep.face_count(), 4);
}

#[test]
fn a_tessellated_solid_of_triangle_strips_and_fans_reads_as_its_triangles() {
    // The unit cube, corners 1 to 8: the bottom, top, front and back as strips, the back's with a
    // triangle of no area at its end, as strips are joined; the left and right as fans. Each quad
    // a, b, c, d counter-clockwise from outside is the strip b, c, a, d and the fan a, b, c, d.
    let text = "ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('t','',(''),(''),'','','');
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF'));
ENDSEC;
DATA;
#1=COORDINATES_LIST('',8,((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),(0.,0.,1.),(1.,0.,1.),\
(1.,1.,1.),(0.,1.,1.)));
#2=COMPLEX_TRIANGULATED_FACE('',#1,0,(),$,(),((4,3,1,2),(6,7,5,8),(2,6,1,5),(4,8,3,7,7)),());
#3=COMPLEX_TRIANGULATED_FACE('',#1,0,(),$,(),(),((1,5,8,4),(2,3,7,6)));
#4=TESSELLATED_SOLID('',(#2,#3),$);
ENDSEC;
END-ISO-10303-21;
";
    let imp = read(text.as_bytes(), "t").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(check::check(&imp.brep).is_empty());
    assert!((measure::volume(&imp.brep) - 1.0).abs() < 1e-12);
}

/// The triangles as loose plane faces, each its own open shell with its own vertices and edges.
fn loose_faces(tris: &[[[f64; 3]; 3]]) -> String {
    let mut lines: Vec<String> = Vec::new();
    let mut add = |s: String| {
        lines.push(s);
        lines.len()
    };
    let num = |x: f64| format!("{x:.17E}");
    let pt = |p: [f64; 3]| format!("({},{},{})", num(p[0]), num(p[1]), num(p[2]));
    for t in tris {
        let vs: Vec<usize> = t
            .iter()
            .map(|&p| {
                let c = add(format!("CARTESIAN_POINT('',{})", pt(p)));
                add(format!("VERTEX_POINT('',#{c})"))
            })
            .collect();
        let mut oes = Vec::new();
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let c = add(format!("CARTESIAN_POINT('',{})", pt(a)));
            let dir = add(format!("DIRECTION('',{})", pt(d.map(|x| x / l))));
            let vec = add(format!("VECTOR('',#{dir},1.)"));
            let line = add(format!("LINE('',#{c},#{vec})"));
            let e = add(format!(
                "EDGE_CURVE('',#{},#{},#{line},.T.)",
                vs[k],
                vs[(k + 1) % 3]
            ));
            oes.push(add(format!("ORIENTED_EDGE('',*,*,#{e},.T.)")));
        }
        let lp = add(format!(
            "EDGE_LOOP('',({}))",
            oes.iter()
                .map(|o| format!("#{o}"))
                .collect::<Vec<_>>()
                .join(",")
        ));
        let bound = add(format!("FACE_OUTER_BOUND('',#{lp},.T.)"));
        let (a, b, c) = (t[0], t[1], t[2]);
        let (u, v) = (
            [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
            [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
        );
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let nl = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        let ul = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        let o = add(format!("CARTESIAN_POINT('',{})", pt(a)));
        let z = add(format!("DIRECTION('',{})", pt(n.map(|x| x / nl))));
        let x = add(format!("DIRECTION('',{})", pt(u.map(|x| x / ul))));
        let ax = add(format!("AXIS2_PLACEMENT_3D('',#{o},#{z},#{x})"));
        let plane = add(format!("PLANE('',#{ax})"));
        let face = add(format!("ADVANCED_FACE('',(#{bound}),#{plane},.T.)"));
        add(format!("OPEN_SHELL('',(#{face}))"));
    }
    let mut out = String::from(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\n\
         FILE_NAME('t','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n",
    );
    for (i, l) in lines.iter().enumerate() {
        out.push_str(&format!("#{}={l};\n", i + 1));
    }
    out.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    out
}

#[test]
fn loose_faces_a_hair_apart_and_some_inside_out_sew_into_their_solid() {
    // The tetrahedron on the unit cube's corner, each face shifted up to 5e-5 along its normal
    // and two of them written inside out: sewn within 1e-3, it is the tetrahedron, volume 1/6 to
    // what the shifts move it.
    let (o, x, y, z) = ([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    let shift = |t: [[f64; 3]; 3], d: f64| {
        let (u, v) = (
            [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]],
            [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]],
        );
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        t.map(|p| [0, 1, 2].map(|k| p[k] + d * n[k] / l))
    };
    let tris = [
        shift([o, y, x], 3e-5),
        [x, o, z],
        shift([o, z, y], -5e-5),
        [x, z, y],
    ];
    let text = loose_faces(&tris);
    let b = super::sew(text.as_bytes(), 1e-3).unwrap();
    assert!(check::check(&b).is_empty());
    assert_eq!(b.face_count(), 4);
    assert!((measure::volume(&b) - 1.0 / 6.0).abs() < 1e-4);
}

/// A cylinder of radius 2 and height 3 whose side is a rational spline surface, straight along its
/// axis and an exact circle round it, so periodic round it: its side is bounded by the two circles
/// alone, a band with no seam, between two plane ends. `top_first` lists the top circle's bound
/// first, so that the seam's first end comes last along the surface's iso line.
fn spline_band_cylinder(top_first: bool) -> String {
    let (r, h, w) = (2.0, 3.0, std::f64::consts::FRAC_1_SQRT_2);
    let round = [
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
    let mut out = String::new();
    let mut ids = [Vec::new(), Vec::new()];
    for (row, z) in [0.0, h].iter().enumerate() {
        for (k, (x, y)) in round.iter().enumerate() {
            let id = 100 + 10 * row + k;
            out += &format!("#{id}=CARTESIAN_POINT('',({},{},{z:?}));\n", x * r, y * r);
            ids[row].push(format!("#{id}"));
        }
    }
    let weights: Vec<String> = (0..9)
        .map(|k| {
            if k % 2 == 0 {
                "1.".into()
            } else {
                format!("{w:?}")
            }
        })
        .collect();
    let row = |i: usize| format!("({})", ids[i].join(","));
    let wrow = format!("({})", weights.join(","));
    out += &format!(
        "#1=(BOUNDED_SURFACE() B_SPLINE_SURFACE(1,2,({},{}),.UNSPECIFIED.,.F.,.T.,.F.) \
         B_SPLINE_SURFACE_WITH_KNOTS((2,2),(3,2,2,2,3),(0.,1.),(0.,0.25,0.5,0.75,1.),.UNSPECIFIED.) \
         GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(({wrow},{wrow})) \
         REPRESENTATION_ITEM('') SURFACE());\n",
        row(0),
        row(1)
    );
    let bounds = if top_first { "#25,#24" } else { "#24,#25" };
    out += &format!(
        "#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=CARTESIAN_POINT('',(0.,0.,0.));
#5=CARTESIAN_POINT('',(0.,0.,{h:?}));
#6=AXIS2_PLACEMENT_3D('',#4,#2,#3);
#7=AXIS2_PLACEMENT_3D('',#5,#2,#3);
#8=CIRCLE('',#6,{r:?});
#9=CIRCLE('',#7,{r:?});
#10=CARTESIAN_POINT('',({r:?},0.,0.));
#11=CARTESIAN_POINT('',({r:?},0.,{h:?}));
#12=VERTEX_POINT('',#10);
#13=VERTEX_POINT('',#11);
#14=EDGE_CURVE('',#12,#12,#8,.T.);
#15=EDGE_CURVE('',#13,#13,#9,.T.);
#20=ORIENTED_EDGE('',*,*,#14,.T.);
#21=EDGE_LOOP('',(#20));
#22=ORIENTED_EDGE('',*,*,#15,.F.);
#23=EDGE_LOOP('',(#22));
#24=FACE_BOUND('',#21,.T.);
#25=FACE_BOUND('',#23,.T.);
#26=ADVANCED_FACE('',({bounds}),#1,.F.);
#30=ORIENTED_EDGE('',*,*,#14,.F.);
#31=EDGE_LOOP('',(#30));
#32=FACE_OUTER_BOUND('',#31,.T.);
#33=PLANE('',#6);
#34=ADVANCED_FACE('',(#32),#33,.F.);
#40=ORIENTED_EDGE('',*,*,#15,.T.);
#41=EDGE_LOOP('',(#40));
#42=FACE_OUTER_BOUND('',#41,.T.);
#43=PLANE('',#7);
#44=ADVANCED_FACE('',(#42),#43,.T.);
#50=CLOSED_SHELL('',(#26,#34,#44));
#51=MANIFOLD_SOLID_BREP('',#50);
"
    );
    format!("ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n{out}ENDSEC;\nEND-ISO-10303-21;\n")
}

#[test]
fn a_band_of_a_spline_surface_closed_round_one_parameter_is_read() {
    // The side's two circles each go once round the surface's period, opposite ways; the reader
    // joins them with a seam along the surface, as it does for a cylinder of its own.
    for top_first in [false, true] {
        let imp = read(spline_band_cylinder(top_first).as_bytes(), "band").unwrap();
        assert!(imp.rejected.is_empty(), "{top_first}: {:?}", imp.rejected);
        // The file is exact: it reads with nothing healed and no tolerance widened.
        assert!(
            imp.healed.is_empty() && imp.tolerant.is_empty(),
            "{top_first}: healed {:?}, tolerant {:?}",
            imp.healed,
            imp.tolerant
        );
        let b = imp.brep;
        assert_eq!(b.solid_ids().count(), 1);
        assert!(check::check(&b).is_empty(), "{:?}", check::check(&b));
        let want = PI * 4.0 * 3.0;
        assert!((measure::volume(&b) - want).abs() < 1e-6 * want);
    }
}

const STRICT: ReadOptions = ReadOptions { strict: true };

#[test]
fn options_are_set_one_by_name() {
    assert_eq!(ReadOptions::default().strict(true), STRICT);
    assert_eq!(
        ReadOptions::default().strict(true).strict(false),
        ReadOptions::default()
    );
    assert!(WriteOptions::default().pcurves);
    assert!(!WriteOptions::default().pcurves(false).pcurves);
}

/// What a strict read of `text` rejects: its one solid, with the reason.
fn strictly_rejected(text: &str) -> String {
    let imp = read_with(text.as_bytes(), "imp", &STRICT).unwrap();
    assert_eq!(imp.brep.solid_ids().count(), 0);
    assert_eq!(imp.rejected.len(), 1, "{:?}", imp.rejected);
    imp.rejected[0].clone()
}

#[test]
fn a_strict_read_keeps_a_file_that_agrees_with_itself_as_it_is() {
    // Another system's block, its geometry meeting exactly: read strictly, it is the plain read's
    // solid, unmoved and with nothing to report.
    let b = block("b", [0.0; 3], [10.0, 6.0, 3.0]);
    let text = write(&b, "b").unwrap().replace("'springrep'", "'other'");
    let imp = read_with(text.as_bytes(), "imp", &STRICT).unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(imp.healed.is_empty() && imp.tolerant.is_empty() && imp.crossed.is_empty());
    assert!(check::check(&imp.brep).is_empty());
    assert!((measure::volume(&imp.brep) - 180.0).abs() < 1e-9);
    let part = read_part_with(text.as_bytes(), "imp", &STRICT).unwrap();
    assert!((measure::volume(&part) - 180.0).abs() < 1e-9);
}

#[test]
fn a_strict_read_rejects_what_the_plain_read_heals() {
    // A vertex 3.7e-5 mm off its corner, and rounds grown 3e-5 mm off the sides they are tangent
    // to: the plain read moves both back.
    let blk = block("b", [0.0; 3], [10.0, 4.0, 3.0]);
    let off = with_vertex_moved(&blk, Vec3::new(2e-5, -1e-5, 3e-5));
    assert_eq!(read(off.as_bytes(), "imp").unwrap().healed.len(), 1);
    let why = strictly_rejected(&off);
    assert!(
        why.contains("ends") && why.contains("from vertex #"),
        "{why}"
    );
    let grown = with_rounds_grown(&rounded_block(), 3e-5);
    assert_eq!(read(grown.as_bytes(), "imp").unwrap().healed.len(), 1);
    let why = strictly_rejected(&grown);
    assert!(why.contains("strict read moves nothing"), "{why}");
}

#[test]
fn a_strict_read_rejects_what_the_plain_read_keeps_with_tolerances_or_moves_between_faces() {
    // Rounds held 0.2 mm off their sides in a file that states 0.4 mm, which the plain read keeps
    // with local tolerances, and the bent edge, which it straightens between its faces.
    let held = with_uncertainty(&with_rounds_grown(&rounded_block(), 0.2), 0.4);
    assert_eq!(read(held.as_bytes(), "imp").unwrap().tolerant.len(), 1);
    let why = strictly_rejected(&held);
    assert!(why.contains("off face #"), "{why}");
    let bent = with_uncertainty(
        &with_edge_bent(&block("b", [0.0; 3], [10.0, 6.0, 3.0]), 0.1),
        0.05,
    );
    let why = strictly_rejected(&bent);
    assert!(why.contains("off face #"), "{why}");
}

#[test]
fn a_strict_read_rejects_what_the_plain_read_heals_loosely() {
    // A 1 mm cube with a corner 0.35 mm off: past what healing may move, so the plain read heals
    // it loosely or rejects it; the strict read rejects it whatever the plain one does.
    let cube = block("b", [0.0; 3], [1.0, 1.0, 1.0]);
    let text = with_vertex_moved(&cube, Vec3::new(0.35, 0.0, 0.0));
    let why = strictly_rejected(&text);
    assert!(why.contains("strict read moves nothing"), "{why}");
}

#[test]
fn a_strict_read_rejects_a_solid_passing_through_itself_and_keeps_one_touching_itself() {
    let why = strictly_rejected(&blocks_in_one_shell());
    assert!(why.contains("passes through itself in"), "{why}");
    // A block bored so the bore touches its top along a line: written as one solid touching
    // itself, read back with the contact joined. It moves nothing, so a strict read keeps it.
    let a = block("a", [0.0, 0.0, 0.0], [12.0, 10.0, 8.0]);
    let frame = Frame::new(Vec3::new(-1.0, 5.0, 5.0), Vec3::X, Vec3::Y).unwrap();
    let bore = extrude(&Profile::circle(Vec2::ZERO, 6.0), &frame, 14.0, "c").unwrap();
    let d = boolean::difference(&a, &bore).unwrap();
    let text = write(&d, "part").unwrap();
    let imp = read_with(text.as_bytes(), "imp", &STRICT).unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(imp.crossed.iter().all(|x| x.depth <= crate::tol::LINEAR));
    let (v0, v1) = (measure::volume(&d), measure::volume(&imp.brep));
    assert!((v1 - v0).abs() <= 1e-9 * v0, "{v1} vs {v0}");
}

#[test]
fn a_face_of_thousands_of_edges_reads_back_on_a_small_stack() {
    // A prism on a 2000-gon of radius 10, 3 tall: each end face's loop has 2000 edges. Read back
    // on a thread with a quarter of a mebibyte of stack, a quarter of what a WebAssembly program
    // has, or on wasm32-wasip1, which cannot spawn a thread, on that program's own stack: tracing a
    // loop's parameter-plane samples must not take a call per edge.
    let n = 2000;
    let pts: Vec<Vec2> = (0..n)
        .map(|i| {
            let a = TAU * i as f64 / n as f64;
            Vec2::new(10.0 * a.cos(), 10.0 * a.sin())
        })
        .collect();
    let p = Profile::polygon(&pts, &["side"]);
    let b = extrude(&p, &Frame::WORLD, 3.0, "gon").unwrap();
    let text = write(&b, "gon").unwrap().replace("'springrep'", "'other'");
    #[cfg(target_os = "wasi")]
    let read_back = read(text.as_bytes(), "imp").unwrap();
    #[cfg(not(target_os = "wasi"))]
    let read_back = std::thread::Builder::new()
        .stack_size(256 * 1024)
        .spawn(move || read(text.as_bytes(), "imp").unwrap())
        .unwrap()
        .join()
        .unwrap();
    assert!(read_back.rejected.is_empty(), "{:?}", read_back.rejected);
    // The polygon's area, n/2 · r² · sin(2π/n), times the height.
    let want = 0.5 * n as f64 * 100.0 * (TAU / n as f64).sin() * 3.0;
    let got = measure::volume(&read_back.brep);
    assert!((got - want).abs() <= 1e-9 * want, "{got} vs {want}");
}

/// The piece a plane cuts from the cone of base radius 2 at z = 0 and apex at z = 2, as three
/// faces: the base, the cut, and the cone between them. `#30` is the cut's plane, `#33` its curve
/// across the cone from `#11` to `#13`, the two ends of the chord x = `x` of the base circle, and
/// `arc` the ends and sense of the base circle's edge. `uses` are the senses of the arc and the
/// chord round the base, of the chord and the curve round the cut, and of the arc and the curve
/// round the cone, whose seam is at `seam` along x.
fn cone_cut(x: f64, cut: &str, arc: &str, seam: f64, uses: [&str; 6]) -> String {
    let y = (4.0 - x * x).sqrt();
    let [ba, bl, cl, ch, ka, kh] = uses;
    step_file(&format!(
        "#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=DIRECTION('',(0.,0.,-1.));
#6=DIRECTION('',({seam:?},0.,0.));
#7=AXIS2_PLACEMENT_3D('',#1,#5,#6);
#8=DIRECTION('',(0.,1.,0.));
#10=CARTESIAN_POINT('',({x:?},{:?},0.));
#11=VERTEX_POINT('',#10);
#12=CARTESIAN_POINT('',({x:?},{y:?},0.));
#13=VERTEX_POINT('',#12);
#20=CIRCLE('',#4,2.);
#21=EDGE_CURVE('',{arc},#20,.T.);
#22=VECTOR('',#8,1.);
#23=LINE('',#10,#22);
#24=EDGE_CURVE('',#11,#13,#23,.T.);
{cut}
#34=EDGE_CURVE('',#11,#13,#33,.T.);
#40=PLANE('',#7);
#41=CONICAL_SURFACE('',#7,2.,{:?});
#50=ORIENTED_EDGE('',*,*,#21,{ba});
#51=ORIENTED_EDGE('',*,*,#24,{bl});
#52=EDGE_LOOP('',(#50,#51));
#53=FACE_OUTER_BOUND('',#52,.T.);
#54=ADVANCED_FACE('',(#53),#40,.T.);
#60=ORIENTED_EDGE('',*,*,#24,{cl});
#61=ORIENTED_EDGE('',*,*,#34,{ch});
#62=EDGE_LOOP('',(#60,#61));
#63=FACE_OUTER_BOUND('',#62,.T.);
#64=ADVANCED_FACE('',(#63),#30,.T.);
#70=ORIENTED_EDGE('',*,*,#21,{ka});
#71=ORIENTED_EDGE('',*,*,#34,{kh});
#72=EDGE_LOOP('',(#70,#71));
#73=FACE_OUTER_BOUND('',#72,.T.);
#74=ADVANCED_FACE('',(#73),#41,.T.);
#80=CLOSED_SHELL('',(#54,#64,#74));
#81=MANIFOLD_SOLID_BREP('',#80);",
        -y,
        PI / 4.0,
    ))
}

/// Reads `text` as one sound solid whose volume and area are `volume` and `area`, and reads it
/// again from what the writer makes of it.
fn assert_reads_as(text: &str, volume: f64, area: f64, what: &str) {
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{what}: {:?}", imp.rejected);
    assert!(check::check(&imp.brep).is_empty(), "{what}");
    let (a, v) = measure::area_and_volume(&imp.brep);
    assert!(
        (v - volume).abs() < 1e-9 * volume,
        "{what}: {v} vs {volume}"
    );
    assert!((a - area).abs() < 1e-9 * area, "{what}: {a} vs {area}");
    let back = read(write(&imp.brep, "p").unwrap().as_bytes(), "back").unwrap();
    assert!(check::check(&back.brep).is_empty(), "{what}");
    let v = measure::volume(&back.brep);
    assert!(
        (v - volume).abs() < 1e-9 * volume,
        "{what}: {v} vs {volume}"
    );
}

#[test]
fn a_hyperbola_edge_is_read_as_its_exact_arc() {
    // The plane x = d cuts the cone r = 2 − z in the hyperbola (2 − z)² − y² = d², and the piece
    // beyond it stands on the base's segment of half angle t = acos(d/2). In polar coordinates the
    // column over (r, θ) is 2 − r high from r = d·sec θ to 2, so the volume is
    // ∫(4/3 − d²sec²θ + d³sec³θ/3)dθ over ±t. The cone's face is √2 times the segment under it,
    // and the cut is ∫z dy = ∫(2 − d·cosh u)·d·cosh u du over ±acosh(2/d).
    // At d = 0.1 the arc spans 7.4 in u and is read in four pieces; at d = 1, in two.
    for d in [1.0_f64, 0.1] {
        let t = (d / 2.0).acos();
        let (sec, tan) = (1.0 / t.cos(), t.tan());
        let volume =
            8.0 / 3.0 * t - 2.0 * d * d * tan + d.powi(3) / 3.0 * (sec * tan + (sec + tan).ln());
        let segment = 2.0 * (2.0 * t - (2.0 * t).sin());
        let u = (2.0 / d).acosh();
        let cut = 4.0 * d * u.sinh() - d * d * (u + (2.0 * u).sinh() / 2.0);
        let area = segment * (1.0 + 2.0_f64.sqrt()) + cut;
        let text = cone_cut(
            d,
            &format!(
                "#25=CARTESIAN_POINT('',({d:?},0.,0.));
#26=DIRECTION('',(-1.,0.,0.));
#27=AXIS2_PLACEMENT_3D('',#25,#26,#8);
#30=PLANE('',#27);
#31=CARTESIAN_POINT('',({d:?},0.,2.));
#32=AXIS2_PLACEMENT_3D('',#31,#3,#5);
#33=HYPERBOLA('',#32,{d:?},{d:?});"
            ),
            "#11,#13",
            -1.0,
            [".F.", ".T.", ".F.", ".T.", ".T.", ".F."],
        );
        assert_reads_as(&text, volume, area, &format!("hyperbola at {d}"));
    }
    // At d = 1 the volume is 8π/9 − 4√3/3 + ln(2 + √3)/3.
    let want = 8.0 * PI / 9.0 - 4.0 * 3.0_f64.sqrt() / 3.0 + (2.0 + 3.0_f64.sqrt()).ln() / 3.0;
    let t = 0.5_f64.acos();
    let got =
        8.0 / 3.0 * t - 2.0 * t.tan() + (t.tan() / t.cos() + (1.0 / t.cos() + t.tan()).ln()) / 3.0;
    assert!((got - want).abs() < 1e-14, "{got} vs {want}");
}

#[test]
fn a_parabola_edge_is_read_as_its_exact_arc() {
    // The plane x + z = 1 is parallel to the cone's side at y = 0 and cuts it in the parabola
    // y² = 1 + 2x, whose focus is on the axis: in polar coordinates r = 1/(1 − cos θ). The piece
    // on the apex's side has volume 8π/9 (the columns are 4/3 − sec²θ/6 for |θ| < π/3 and
    // 1/(6(1 − cos θ)²) beyond, whose integrals are 8π/9 − √3/3 and √3/3), so the piece read here,
    // the rest of the cone's 8π/3, has 16π/9.
    // Its base is the disc less the segment beyond x = 1, 8π/3 + √3. Its cone face stands over
    // the disc less the parabola's inside, 8π/3 − √3, and is √2 times that. Its cut is the
    // parabola's segment of chord 2√3 and height 3/√2: by Archimedes two thirds of their product.
    let r3 = 3.0_f64.sqrt();
    let r2 = 2.0_f64.sqrt();
    let volume = 16.0 * PI / 9.0;
    let area = 8.0 * PI / 3.0 + r3 + r2 * (8.0 * PI / 3.0 - r3) + 2.0 * 6.0_f64.sqrt();
    let h = r2 / 2.0;
    let text = cone_cut(
        1.0,
        &format!(
            "#25=CARTESIAN_POINT('',(1.,0.,0.));
#26=DIRECTION('',({h:?},0.,{h:?}));
#27=AXIS2_PLACEMENT_3D('',#25,#26,#8);
#30=PLANE('',#27);
#31=CARTESIAN_POINT('',(-0.5,0.,1.5));
#28=DIRECTION('',({h:?},0.,{:?}));
#32=AXIS2_PLACEMENT_3D('',#31,#26,#28);
#33=PARABOLA('',#32,{:?});",
            -h,
            r2 / 4.0,
        ),
        "#13,#11",
        1.0,
        [".F.", ".F.", ".T.", ".F.", ".T.", ".T."],
    );
    assert_reads_as(&text, volume, area, "parabola");
}

#[test]
fn an_open_conic_that_closes_on_itself_is_refused() {
    let text = step_file(
        "#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=PARABOLA('',#4,1.);
#6=CARTESIAN_POINT('',(1.,2.,0.));
#7=VERTEX_POINT('',#6);
#8=EDGE_CURVE('',#7,#7,#5,.T.);
#9=ORIENTED_EDGE('',*,*,#8,.T.);
#10=EDGE_LOOP('',(#9));
#11=FACE_OUTER_BOUND('',#10,.T.);
#12=PLANE('',#4);
#13=ADVANCED_FACE('',(#11),#12,.T.);
#14=CLOSED_SHELL('',(#13));
#15=MANIFOLD_SOLID_BREP('',#14);",
    );
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(
        imp.rejected[0].contains("an open conic cannot close on itself"),
        "{:?}",
        imp.rejected
    );
}

/// A hole of radius 2.4 drilled along a radius into a cylinder of radius 56.4, the hole turned
/// `turn` about the cylinder's axis.
fn drilled_cylinder(turn: f64) -> Brep {
    let big = crate::primitive::cylinder(&Frame::WORLD, 56.4, 20.0, "big").unwrap();
    let dir = Vec3::new(turn.cos(), turn.sin(), 0.0);
    let frame = Frame::new(Vec3::new(0.0, 0.0, 10.0) + dir * 40.0, dir, Vec3::Z).unwrap();
    let hole = crate::primitive::cylinder(&frame, 2.4, 30.0, "hole").unwrap();
    boolean::difference(&big, &hole).unwrap()
}

#[test]
fn an_edge_past_the_end_of_a_closed_spline_is_written_on_a_spline_that_starts_at_it() {
    // At no turn the hole's mouth is two edges, split where it crosses the cylinder's seam, and
    // one of them runs past where the mouth's curve ends; turned, it is one edge all the way
    // round, from a vertex a quarter of the way along the curve.
    for (turn, closed) in [(0.0, false), (0.3, true)] {
        let b = drilled_cylinder(turn);
        let dir = Vec3::new(turn.cos(), turn.sin(), 0.0);
        let mut past = 0;
        for e in b.edge_ids() {
            let edge = b.edge(e);
            let Some(curve @ Curve::Intersection(_)) = &edge.curve else {
                continue;
            };
            let Some(turned) = super::write::turn_from(curve, edge.range) else {
                assert!(edge.range.1 <= curve.domain().1 + 1e-9);
                continue;
            };
            past += 1;
            assert_eq!(edge.start == edge.end, closed);
            let (t0, t1) = edge.range;
            assert_eq!(turned.knots[0], t0);
            assert!(turned.knots[turned.knots.len() - 1] >= t1);
            assert!(turned.ctrl[0].dist(b.vertex(edge.start).point) < crate::tol::LINEAR);
            if closed {
                let last = turned.ctrl[turned.ctrl.len() - 1];
                assert!(last.dist(b.vertex(edge.end).point) < crate::tol::LINEAR);
            }
            let spline =
                crate::geom::BSpline::new(turned.degree, turned.knots, turned.ctrl).unwrap();
            for i in 0..=64 {
                let p = spline.eval(t0 + (t1 - t0) * i as f64 / 64.0);
                // On the cylinder, 56.4 from its axis, and on the hole, 2.4 from its own.
                assert!((p.x.hypot(p.y) - 56.4).abs() < crate::tol::LINEAR);
                let off = p - Vec3::new(0.0, 0.0, 10.0);
                let out = off - dir * off.dot(dir);
                assert!((out.norm() - 2.4).abs() < crate::tol::LINEAR);
            }
        }
        assert_eq!(past, 1, "turn {turn}");
        assert_pcurves_follow(&b);
        round_trip(&b);
    }
}

#[test]
fn a_plane_s_edge_past_the_end_of_a_closed_spline_has_that_spline_for_its_pcurve() {
    // A ring of tube radius 2 round a circle of radius 8, tilted, cut out of a slab through its
    // top at z = 10: each mouth it leaves in that plane is a closed curve of the two.
    let tilt = 0.6f64;
    let axis = Vec3::new(tilt.sin(), 0.0, tilt.cos());
    let centre = Vec3::new(0.0, 0.0, 10.0);
    let frame = Frame::new(centre, axis, Vec3::Y).unwrap();
    let ring = crate::primitive::torus(&frame, 8.0, 2.0, "ring").unwrap();
    let slab = block("slab", [-20.0, -20.0, 0.0], [20.0, 20.0, 10.0]);
    let b = boolean::difference(&slab, &ring).unwrap();
    let (mut whole, mut part) = (0, 0);
    for c in b.coedge_ids() {
        let edge = b.edge(b.coedge(c).edge);
        let Some(curve @ Curve::Intersection(_)) = &edge.curve else {
            continue;
        };
        if !matches!(b.face(b.coedge_face(c)).surface, Surface::Plane(_)) {
            continue;
        }
        let Some(turned) = super::write::turn_from(curve, edge.range) else {
            continue;
        };
        if edge.start == edge.end {
            whole += 1;
        } else {
            part += 1;
        }
        let Some(super::write::Curve2::Spline { knots, ctrl, .. }) = super::write::pcurve(&b, c)
        else {
            panic!("a spline on a plane has a spline for its pcurve");
        };
        assert_eq!(knots, turned.knots);
        assert_eq!(knots[0], edge.range.0);
        assert_eq!(ctrl.len(), turned.ctrl.len());
        assert!(turned.ctrl[0].dist(b.vertex(edge.start).point) < crate::tol::LINEAR);
        let spline = crate::geom::BSpline::new(turned.degree, turned.knots, turned.ctrl).unwrap();
        let (t0, t1) = edge.range;
        for i in 0..=64 {
            let p = spline.eval(t0 + (t1 - t0) * i as f64 / 64.0);
            assert!((p.z - 10.0).abs() < crate::tol::LINEAR);
            let off = p - centre;
            let rise = off.dot(axis);
            let reach = (off - axis * rise).norm();
            assert!(((reach - 8.0).hypot(rise) - 2.0).abs() < crate::tol::LINEAR);
        }
    }
    assert_eq!((whole, part), (1, 1));
    assert_pcurves_follow(&b);
    round_trip(&b);
}

#[test]
fn a_rational_circle_is_turned_to_start_where_its_edge_does() {
    // The unit circle as a rational quadratic in four quarter arcs, each a unit of parameter: by
    // symmetry it is at 45° half a unit in, and on an axis at each whole one.
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let square = [
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
    let ctrl: Vec<Vec3> = square.iter().map(|&(x, y)| Vec3::new(x, y, 0.0)).collect();
    let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0];
    let circle = crate::geom::NurbsCurve::new(2, knots, ctrl, Some(weights)).unwrap();
    let curve = Curve::Nurbs(Box::new(circle));
    assert!(super::write::turn_from(&curve, (0.5, 3.0)).is_none());
    let turned = super::write::turn_from(&curve, (0.5, 4.5)).unwrap();
    assert_eq!(turned.knots[0], 0.5);
    assert_eq!(turned.knots[turned.knots.len() - 1], 4.5);
    let back =
        crate::geom::NurbsCurve::new(turned.degree, turned.knots, turned.ctrl, turned.weights)
            .unwrap();
    let at = [
        (0.5, (w, w)),
        (1.0, (0.0, 1.0)),
        (2.0, (-1.0, 0.0)),
        (3.0, (0.0, -1.0)),
        (4.0, (1.0, 0.0)),
        (4.5, (w, w)),
    ];
    for (t, (x, y)) in at {
        assert!(back.eval(t).dist(Vec3::new(x, y, 0.0)) < 1e-14, "at {t}");
    }
    for i in 0..=64 {
        let p = back.eval(0.5 + 4.0 * i as f64 / 64.0);
        assert!((p.x.hypot(p.y) - 1.0).abs() < 1e-14 && p.z == 0.0);
    }
}

/// Everything written lies in the domain of what it is written on, as a reader that wraps only
/// what STEP itself wraps sees it: each edge in its curve's, each pcurve in its surface's and in
/// its own, and each loop's pcurves, its poles' lines among them, meet end to end along every
/// direction the surface has ends in. A seam's two pcurves on its face lie a whole number of
/// periods apart, and not on one another. A face kept as the kernel has it has no pcurves to check.
/// Returns how many faces had such a direction, and how many of those had their surface written
/// from where the face starts.
fn assert_within_step_domains(b: &Brep) -> (usize, usize) {
    let inside = |x: f64, (lo, hi): (f64, f64)| {
        let slack = 1e-9 * (1.0 + lo.abs().max(hi.abs()));
        x >= lo - slack && x <= hi + slack
    };
    for e in b.edge_ids() {
        let edge = b.edge(e);
        let Some(curve) = &edge.curve else { continue };
        if matches!(curve, Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_)) {
            continue;
        }
        let domain = match super::write::turn_from(curve, edge.range) {
            Some(s) => (s.knots[0], s.knots[s.knots.len() - 1]),
            None => curve.domain(),
        };
        let (t0, t1) = edge.range;
        assert!(
            inside(t0, domain) && inside(t1, domain),
            "a {} edge over {t0}..{t1} on a curve written over {domain:?}",
            curve.kind()
        );
    }
    let (mut ended, mut rebased) = (0, 0);
    for f in b.face_ids() {
        let form = super::write::face_form(b, f);
        let (du, dv) = form.ended;
        if du.is_none() && dv.is_none() {
            continue;
        }
        ended += 1;
        rebased += usize::from(form.surface.is_some());
        if form.kept {
            continue;
        }
        let surface = &b.face(f).surface;
        let kind = surface.kind();
        // Each seam's two pcurves on this face, at the middle of its edge: whole periods apart
        // along some direction, and so not one on the other.
        for &lp in &b.face(f).loops {
            for &c in &b.lp(lp).coedges {
                let Some(p) = b.partner(c).filter(|&p| p > c && b.coedge_face(p) == f) else {
                    continue;
                };
                let edge = b.edge(b.coedge(c).edge);
                let (Some(a), Some(z)) = (form.pcurve(b, c), form.pcurve(b, p)) else {
                    continue;
                };
                let t = [0.5 * (edge.range.0 + edge.range.1)];
                let gap = a.track(&t)[0] - z.track(&t)[0];
                let whole = |g: f64, period: Option<f64>| {
                    period.is_some_and(|q| {
                        let n = (g / q).round();
                        n != 0.0 && (g - n * q).abs() <= 1e-6 * q
                    })
                };
                let apart = whole(gap.x, surface.u_period()) || whole(gap.y, surface.v_period());
                assert!(apart, "{kind}: a seam's two pcurves {gap:?} apart");
            }
        }
        for &lp in &b.face(f).loops {
            let mut ends: Vec<(Vec2, Vec2)> = Vec::new();
            // A coedge the writer finds no pcurve for within the tolerance has none written, and
            // its loop's pcurves cannot be followed end to end.
            let mut whole = true;
            for &c in &b.lp(lp).coedges {
                let co = b.coedge(c);
                let edge = b.edge(co.edge);
                if edge.is_degenerate() {
                    // A pole's line across the parameter plane, moved as the writer moves it.
                    let by = form.shift[&c];
                    let (t0, t1) = edge.range;
                    let ts = if co.reversed { [t1, t0] } else { [t0, t1] };
                    let e = b.coedge_uvs(c, &ts);
                    ends.push((e[0] + by, e[1] + by));
                    continue;
                }
                let Some(c2) = form.pcurve(b, c) else {
                    whole = false;
                    continue;
                };
                let (t0, t1) = edge.range;
                if let super::write::Curve2::Spline { knots, .. } = &c2 {
                    let own = (knots[0], knots[knots.len() - 1]);
                    assert!(
                        inside(t0, own) && inside(t1, own),
                        "{kind}: pcurve's own domain"
                    );
                }
                let ts: Vec<f64> = (0..=64).map(|i| t0 + (t1 - t0) * i as f64 / 64.0).collect();
                let uvs = c2.track(&ts);
                // The box the pcurve is known to lie in, as well as points along it.
                let (_, bounds) = super::write::course(b, c, Some(&c2));
                for &uv in uvs.iter().chain(&bounds) {
                    if let Some(d) = du {
                        assert!(inside(uv.x, d), "{kind}: u {} outside {d:?}", uv.x);
                    }
                    if let Some(d) = dv {
                        assert!(inside(uv.y, d), "{kind}: v {} outside {d:?}", uv.y);
                    }
                }
                let (a, z) = (uvs[0], uvs[64]);
                ends.push(if co.reversed { (z, a) } else { (a, z) });
            }
            for k in 0..if whole { ends.len() } else { 0 } {
                let (to, from) = (ends[k].1, ends[(k + 1) % ends.len()].0);
                let gap = to - from;
                for (g, d) in [(gap.x, du), (gap.y, dv)] {
                    if let Some((lo, hi)) = d {
                        assert!(
                            g.abs() <= 1e-6 * (hi - lo),
                            "{kind}: a loop's pcurves {g} apart"
                        );
                    }
                }
            }
        }
    }
    (ended, rebased)
}

/// A cam: a closed spline round `centre`, of radius `r` waved three times, its turn starting at
/// `start` radians round from `+x`.
fn cam(centre: Vec2, r: f64, start: f64) -> Profile {
    let pts: Vec<Vec2> = (0..9)
        .map(|i| {
            let a = start + TAU * i as f64 / 9.0;
            let rr = r * (1.0 + 0.2 * (3.0 * a).cos());
            centre + Vec2::new(rr * a.cos(), rr * a.sin())
        })
        .collect();
    Profile::points(&pts, PointsInterp::PeriodicSpline, "side").unwrap()
}

/// The half of `b` on the side of the plane `x = 0` that `keep` points to.
fn half(b: &Brep, keep: f64) -> Brep {
    let lo = if keep > 0.0 { -100.0 } else { 0.0 };
    boolean::difference(
        b,
        &block("cut", [lo, -100.0, -100.0], [lo + 100.0, 100.0, 100.0]),
    )
    .unwrap()
}

#[test]
fn a_face_across_an_extruded_splines_join_is_written_inside_its_surface() {
    // The cam's turn starts at +x: the half towards +x has its side across the join, the half
    // towards -x has not.
    let tall = extrude(&cam(Vec2::ZERO, 10.0, 0.0), &Frame::WORLD, 4.0, "cam").unwrap();
    for (keep, across) in [(1.0, 1), (-1.0, 0)] {
        let b = half(&tall, keep);
        assert_eq!(assert_within_step_domains(&b), (1, across));
        round_trip(&b);
    }
}

#[test]
fn a_face_across_a_revolved_splines_join_is_written_inside_its_surface() {
    // A ring whose tube is a cam turned about z, the cam's turn starting outermost: the ring's
    // outer half has its face across the join.
    let ring = revolve(
        &cam(Vec2::new(20.0, 0.0), 3.0, 0.0),
        &xz(),
        Vec3::ZERO,
        Vec3::Z,
        TAU,
        "r",
    )
    .unwrap();
    let core = extrude(
        &Profile::circle(Vec2::ZERO, 40.0),
        &Frame::new(Vec3::new(0.0, 0.0, -10.0), Vec3::Z, Vec3::X).unwrap(),
        20.0,
        "core",
    )
    .unwrap();
    for (b, across) in [
        (boolean::difference(&ring, &core).unwrap(), 1),
        (boolean::intersection(&ring, &core).unwrap(), 0),
    ] {
        assert_eq!(assert_within_step_domains(&b), (1, across));
        round_trip(&b);
    }
}

#[test]
fn a_face_across_a_lofted_splines_join_is_written_inside_its_surface() {
    let sections: Vec<crate::loft::Section> = (0..2)
        .map(|i| crate::loft::Section {
            profile: cam(Vec2::ZERO, 8.0 + 2.0 * i as f64, 0.0),
            frame: Frame::new(Vec3::new(0.0, 0.0, 6.0 * i as f64), Vec3::Z, Vec3::X).unwrap(),
        })
        .collect();
    let loft = crate::loft::loft(&sections, false, "l").unwrap();
    for (keep, across) in [(1.0, 1), (-1.0, 0)] {
        let b = half(&loft, keep);
        assert_eq!(assert_within_step_domains(&b), (1, across));
        round_trip(&b);
    }
}

#[test]
fn an_unclamped_closed_nurbs_is_turned_to_start_where_its_edge_does() {
    // A cubic on even knots round eight points, its first three again at the end, as a file may
    // carry a closed curve: it is closed, and has no knot of more than one at its ends.
    let ring: Vec<Vec3> = (0..8)
        .map(|i| {
            let a = TAU * i as f64 / 8.0;
            Vec3::new(5.0 * a.cos(), 3.0 * a.sin(), 0.5 * (2.0 * a).sin())
        })
        .collect();
    let ctrl: Vec<Vec3> = (0..11).map(|i| ring[i % 8]).collect();
    let knots: Vec<f64> = (0..15).map(|i| i as f64).collect();
    let n = crate::geom::NurbsCurve::new(3, knots, ctrl, None).unwrap();
    assert_eq!(n.domain(), (3.0, 11.0));
    assert_eq!(n.period(), Some(8.0));
    let curve = Curve::Nurbs(Box::new(n.clone()));
    let turned = super::write::turn_from(&curve, (5.5, 13.5)).expect("past the end");
    assert_eq!(turned.knots[0], 5.5);
    assert_eq!(turned.knots[turned.knots.len() - 1], 13.5);
    assert!(turned.ctrl[0].dist(n.eval(5.5)) < 1e-12);
    let back = crate::geom::BSpline::new(turned.degree, turned.knots, turned.ctrl).unwrap();
    for i in 0..=64 {
        let t = 5.5 + 8.0 * i as f64 / 64.0;
        assert!(back.eval(t).dist(n.eval(t)) < 1e-12, "at {t}");
    }
}

#[test]
fn a_closed_nurbs_with_unequal_end_weights_is_turned_to_start_where_its_edge_does() {
    // The unit circle in four quarter arcs, the last arc's weights (1, w, 1) taken to
    // (1, w·λ, λ²): the same arc, run at another rate, so the circle ends on a weight of λ².
    let (w, lambda) = (std::f64::consts::FRAC_1_SQRT_2, 2.0);
    let square = [
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
    let ctrl: Vec<Vec3> = square.iter().map(|&(x, y)| Vec3::new(x, y, 0.0)).collect();
    let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w * lambda, lambda * lambda];
    let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0];
    let n = crate::geom::NurbsCurve::new(2, knots, ctrl, Some(weights)).unwrap();
    let curve = Curve::Nurbs(Box::new(n.clone()));
    let turned = super::write::turn_from(&curve, (0.5, 4.5)).expect("past the end");
    let back =
        crate::geom::NurbsCurve::new(turned.degree, turned.knots, turned.ctrl, turned.weights)
            .unwrap();
    for i in 0..=64 {
        let t = 0.5 + 4.0 * i as f64 / 64.0;
        let p = back.eval(t);
        assert!((p.x.hypot(p.y) - 1.0).abs() < 1e-14 && p.z == 0.0, "at {t}");
        assert!(p.dist(n.eval(t)) < 1e-14, "at {t}");
    }
}

#[test]
fn whole_closed_spline_solids_are_written_inside_their_surfaces() {
    // Uncut, each closed side is one face round its seam: in its surface's own turn, and with no
    // surface written anew.
    let tall = extrude(
        &cam(Vec2::new(1.0, -2.0), 10.0, 0.7),
        &Frame::WORLD,
        4.0,
        "cam",
    )
    .unwrap();
    let ring = revolve(
        &cam(Vec2::new(20.0, 0.0), 3.0, 2.0),
        &xz(),
        Vec3::ZERO,
        Vec3::Z,
        TAU,
        "r",
    )
    .unwrap();
    let sections: Vec<crate::loft::Section> = (0..3)
        .map(|i| crate::loft::Section {
            profile: cam(Vec2::ZERO, 8.0 + i as f64, 0.4 * i as f64),
            frame: Frame::new(Vec3::new(0.0, 0.0, 5.0 * i as f64), Vec3::Z, Vec3::X).unwrap(),
        })
        .collect();
    let loft = crate::loft::loft(&sections, false, "l").unwrap();
    for b in [tall, ring, loft] {
        let (ended, rebased) = assert_within_step_domains(&b);
        assert!(ended > 0);
        assert_eq!(rebased, 0);
        round_trip(&b);
    }
}

#[test]
fn a_hole_on_the_join_of_a_face_across_it_is_written_inside_its_surface() {
    // The halves whose sides run across their surfaces' joins, each drilled right on the join: a
    // rod along x through the cam's and the loft's side, and one along z through the ring's tube.
    // The hole's loop lies across the join too, inside the face's own.
    let core = extrude(
        &Profile::circle(Vec2::ZERO, 40.0),
        &Frame::new(Vec3::new(0.0, 0.0, -10.0), Vec3::Z, Vec3::X).unwrap(),
        20.0,
        "core",
    )
    .unwrap();
    let rod = |frame: Frame, d: f64| {
        extrude(&Profile::circle(Vec2::ZERO, d), &frame, 60.0, "rod").unwrap()
    };
    let along_x = Frame::new(Vec3::new(-30.0, 0.0, 2.0), Vec3::X, Vec3::Y).unwrap();
    let tall = extrude(&cam(Vec2::ZERO, 10.0, 0.0), &Frame::WORLD, 4.0, "cam").unwrap();
    let sections: Vec<crate::loft::Section> = (0..2)
        .map(|i| crate::loft::Section {
            profile: cam(Vec2::ZERO, 8.0 + 2.0 * i as f64, 0.0),
            frame: Frame::new(Vec3::new(0.0, 0.0, 4.0 * i as f64), Vec3::Z, Vec3::X).unwrap(),
        })
        .collect();
    let loft = crate::loft::loft(&sections, false, "l").unwrap();
    let ring = revolve(
        &cam(Vec2::new(20.0, 0.0), 3.0, 0.0),
        &xz(),
        Vec3::ZERO,
        Vec3::Z,
        TAU,
        "r",
    )
    .unwrap();
    let along_z = Frame::new(Vec3::new(23.6, 0.0, -30.0), Vec3::Z, Vec3::X).unwrap();
    let outer = boolean::difference(&ring, &core).unwrap();
    for b in [
        boolean::difference(&half(&tall, 1.0), &rod(along_x, 2.0)).unwrap(),
        boolean::difference(&half(&loft, 1.0), &rod(along_x, 2.0)).unwrap(),
        boolean::difference(&outer, &rod(along_z, 1.0)).unwrap(),
    ] {
        assert_eq!(assert_within_step_domains(&b), (1, 1));
        round_trip(&b);
    }
}

#[test]
fn a_piece_across_both_of_a_rings_joins_is_written_inside_its_surface() {
    // A block over the ring's outer side round where both its angle and its tube's turn start,
    // so the piece of the tube it keeps runs across both.
    let ring = revolve(
        &cam(Vec2::new(20.0, 0.0), 3.0, 0.0),
        &xz(),
        Vec3::ZERO,
        Vec3::Z,
        TAU,
        "r",
    )
    .unwrap();
    let b =
        boolean::intersection(&ring, &block("k", [20.5, -6.0, -10.0], [30.0, 6.0, 10.0])).unwrap();
    assert_eq!(assert_within_step_domains(&b), (1, 1));
    round_trip(&b);
}

#[test]
fn an_assembly_writes_each_part_on_its_own_surfaces() {
    // A part with a face across its surface's join, then another: the second's faces are written
    // on their own surfaces, not on the first's turned one, and both come back whole.
    use super::write::{Instance, Part, write_assembly};
    let tall = extrude(&cam(Vec2::ZERO, 10.0, 0.0), &Frame::WORLD, 4.0, "cam").unwrap();
    let cam_half = half(&tall, 1.0);
    let sections: Vec<crate::loft::Section> = (0..2)
        .map(|i| crate::loft::Section {
            profile: cam(Vec2::ZERO, 8.0 + 2.0 * i as f64, 0.0),
            frame: Frame::new(Vec3::new(0.0, 0.0, 6.0 * i as f64), Vec3::Z, Vec3::X).unwrap(),
        })
        .collect();
    let loft_half = half(&crate::loft::loft(&sections, false, "l").unwrap(), 1.0);
    let blk = block("blk", [0.0, 0.0, 0.0], [10.0, 4.0, 3.0]);
    let count = |t: &str, k: &str| t.matches(k).count();
    for (x, y) in [
        (&cam_half, &blk),
        (&blk, &cam_half),
        (&cam_half, &loft_half),
    ] {
        let (alone_x, alone_y) = (write(x, "x").unwrap(), write(y, "y").unwrap());
        let parts = [Part::new("x", x, None), Part::new("y", y, None)];
        let at = crate::math::Transform::translation(Vec3::new(100.0, 0.0, 0.0));
        let instances = [
            Instance::new("i", 0, crate::math::Transform::IDENTITY),
            Instance::new("j", 1, at),
        ];
        let both = write_assembly("a", &parts, &instances).unwrap();
        for k in [
            "PLANE(",
            "SURFACE_OF_LINEAR_EXTRUSION(",
            "PCURVE(",
            "ADVANCED_FACE(",
        ] {
            assert_eq!(
                count(&both, k),
                count(&alone_x, k) + count(&alone_y, k),
                "{k}"
            );
        }
        let back = read(both.as_bytes(), "imp").unwrap();
        assert!(back.rejected.is_empty(), "{:?}", back.rejected);
        assert!(check::check(&back.brep).is_empty());
        let want = measure::volume(x) + measure::volume(y);
        assert!((measure::volume(&back.brep) - want).abs() <= 1e-9 * want);
    }
}

#[test]
fn a_loft_to_a_point_from_a_closed_spline_is_written_across_its_pole() {
    // Its side is closed along `u` and meets the point in a pole, whose line crosses the closed
    // direction in the parameter plane: the loop runs along it and closes.
    let sections = vec![crate::loft::Section {
        profile: cam(Vec2::ZERO, 8.0, 0.0),
        frame: Frame::WORLD,
    }];
    for ruled in [false, true] {
        let apex = Vec3::new(0.0, 0.0, 9.0);
        let b = crate::loft::loft_to_point(&sections, apex, ruled, "l").unwrap();
        assert!(assert_within_step_domains(&b).0 > 0);
        round_trip(&b);
    }
}

/// A file of a solid whose one face is a rational spline surface closed along `u`, with a pole at
/// each end of `v` and one seam: a circle of radius 2 turned along a parabola from (0, 0, −1) to
/// (0, 0, 1) through the circle, so the solid is the revolution of ρ = 1 − z².
fn lemon_text() -> String {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let ring = [
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (-2.0, 2.0),
        (-2.0, 0.0),
        (-2.0, -2.0),
        (0.0, -2.0),
        (2.0, -2.0),
        (2.0, 0.0),
    ];
    let mut t = String::from(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\n\
         FILE_NAME('','',(''),(''),'','','');\n\
         FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n",
    );
    t += "#1=CARTESIAN_POINT('',(0.,0.,-1.));\n#2=CARTESIAN_POINT('',(0.,0.,1.));\n";
    for (i, (x, y)) in ring.iter().enumerate() {
        t += &format!("#{}=CARTESIAN_POINT('',({x:?},{y:?},0.));\n", 10 + i);
    }
    let net: Vec<String> = (0..9).map(|i| format!("(#1,#{},#2)", 10 + i)).collect();
    let ws: Vec<String> = (0..9)
        .map(|i| {
            let x = if i % 2 == 0 { 1.0 } else { w };
            format!("({x:?},{x:?},{x:?})")
        })
        .collect();
    t += &format!(
        "#30=(BOUNDED_SURFACE() B_SPLINE_SURFACE(2,2,({}),.UNSPECIFIED.,.F.,.F.,.F.) \
         B_SPLINE_SURFACE_WITH_KNOTS((3,2,2,2,3),(3,3),(0.,1.,2.,3.,4.),(0.,1.),.UNSPECIFIED.) \
         GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(({})) \
         REPRESENTATION_ITEM('') SURFACE());\n",
        net.join(","),
        ws.join(",")
    );
    t += "#31=B_SPLINE_CURVE_WITH_KNOTS('',2,(#1,#10,#2),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),\
          .UNSPECIFIED.);\n";
    t += "#32=VERTEX_POINT('',#1);\n#33=VERTEX_POINT('',#2);\n";
    t += "#34=EDGE_CURVE('',#32,#33,#31,.T.);\n";
    t += "#35=ORIENTED_EDGE('',*,*,#34,.T.);\n#36=ORIENTED_EDGE('',*,*,#34,.F.);\n";
    t += "#37=EDGE_LOOP('',(#35,#36));\n#38=FACE_OUTER_BOUND('',#37,.T.);\n";
    t += "#39=ADVANCED_FACE('',(#38),#30,.T.);\n#40=CLOSED_SHELL('',(#39));\n";
    t += "#41=MANIFOLD_SOLID_BREP('',#40);\n";
    t += "#50=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n";
    t += "#51=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n";
    t += "#52=(NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT());\n";
    t += "#53=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-7),#50,'distance_accuracy_value',\
          '');\n";
    t += "#54=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#53)) \
          GLOBAL_UNIT_ASSIGNED_CONTEXT((#50,#51,#52)) REPRESENTATION_CONTEXT('','3D'));\n";
    t += "#55=ADVANCED_BREP_SHAPE_REPRESENTATION('lemon',(#41),#54);\n";
    t += "ENDSEC;\nEND-ISO-10303-21;\n";
    t
}

#[test]
fn a_face_between_two_poles_keeps_its_seam_a_period_apart() {
    // Volume π∫(1 − z²)² dz over −1..1, 16π/15; its seam's two pcurves a turn of `u` apart.
    let imp = read(lemon_text().as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    let b = imp.brep;
    let want = 16.0 * PI / 15.0;
    assert!((measure::volume(&b) - want).abs() <= 1e-9 * want);
    assert!(assert_within_step_domains(&b).0 > 0);
    round_trip(&b);
}

#[test]
fn a_face_that_fits_no_one_turn_is_written_as_the_kernel_has_it() {
    // A ring with a cam for its tube, less a slab tilted through it: one face runs a whole turn
    // of the tube and on past its seam, which no one turn of the tube's curve holds. It is kept on
    // the surface as it is, without pcurves, as before, and reads back whole.
    let ring = revolve(
        &cam(Vec2::new(20.0, 0.0), 3.0, 0.0),
        &xz(),
        Vec3::ZERO,
        Vec3::Z,
        TAU,
        "r",
    )
    .unwrap();
    let (tilt, turn) = (0.21f64, 2.2f64);
    let z = Vec3::new(tilt.sin() * turn.cos(), tilt.sin() * turn.sin(), tilt.cos());
    let frame = Frame::new(Vec3::new(0.0, 0.0, 0.3), z, Vec3::Y.cross(z)).unwrap();
    let slab = extrude(
        &Profile::rect(Vec2::ZERO, 300.0, 300.0),
        &frame,
        100.0,
        "cut",
    )
    .unwrap();
    let b = boolean::difference(&ring, &slab).unwrap();
    assert!(check::check(&b).is_empty());
    let kept = b
        .face_ids()
        .filter(|&f| super::write::face_form(&b, f).kept)
        .count();
    assert_eq!(kept, 1);
    round_trip(&b);
}

#[test]
fn a_cut_spline_face_s_trimmed_edges_are_placed_within_their_pcurves() {
    // A lofted cam cut by planes along x and then across z, so its intersection edges start part
    // way along their pcurves: each stretch the writer bounds lies within its pcurve's domain.
    let sections: Vec<crate::loft::Section> = (0..3)
        .map(|i| crate::loft::Section {
            profile: cam(Vec2::ZERO, 8.0 + 2.0 * i as f64, 0.0),
            frame: Frame::new(Vec3::new(0.0, 0.0, 3.0 * i as f64), Vec3::Z, Vec3::X).unwrap(),
        })
        .collect();
    let loft = crate::loft::loft(&sections, false, "l").unwrap();
    let mut written = 0;
    for off in [-4.0, -0.93, 2.4] {
        let first = boolean::difference(
            &loft,
            &block("x", [off, -100.0, -100.0], [off + 100.0, 100.0, 100.0]),
        )
        .unwrap();
        for zc in [1.31, 2.53, 4.97] {
            let b = boolean::difference(
                &first,
                &block("z", [-100.0, -100.0, zc], [100.0, 100.0, 100.0]),
            )
            .unwrap();
            assert!(check::check(&b).is_empty());
            assert_within_step_domains(&b);
            write(&b, "p").unwrap();
            written += 1;
        }
    }
    assert_eq!(written, 9);
}

#[test]
fn a_spline_whose_kind_implies_its_knots_has_them() {
    let knots = |kind: &str, degree: usize, n: usize| {
        super::geometry::implied_knots(0, kind, degree, n).unwrap()
    };
    assert_eq!(
        knots("UNIFORM_CURVE", 2, 5),
        [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    );
    assert_eq!(
        knots("QUASI_UNIFORM_CURVE", 2, 5),
        [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0]
    );
    assert_eq!(
        knots("BEZIER_SURFACE", 2, 5),
        [0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0]
    );
    assert_eq!(knots("QUASI_UNIFORM_CURVE", 1, 2), [0.0, 0.0, 1.0, 1.0]);
    assert!(super::geometry::implied_knots(0, "BEZIER_CURVE", 2, 4).is_err());
    assert!(super::geometry::implied_knots(0, "UNIFORM_CURVE", 3, 3).is_err());
    assert!(super::geometry::implied_knots(0, "B_SPLINE_CURVE", 2, 5).is_err());
}

#[test]
fn quadratic_splines_whose_kind_implies_their_knots_pass_through_their_closed_form_points() {
    let points = "#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(1.,0.,0.));\
                  #3=CARTESIAN_POINT('',(2.,2.,0.));#4=CARTESIAN_POINT('',(4.,1.,0.));\
                  #5=CARTESIAN_POINT('',(5.,3.,1.));#6=CARTESIAN_POINT('',(0.,1.,0.));\
                  #7=CARTESIAN_POINT('',(1.,1.,4.));#8=CARTESIAN_POINT('',(2.,1.,0.));\
                  #9=CARTESIAN_POINT('',(0.,2.,0.));#10=CARTESIAN_POINT('',(2.,0.,0.));\
                  #16=CARTESIAN_POINT('',(1.,2.,0.));";
    let curves = "#11=UNIFORM_CURVE('',2,(#1,#2,#3,#4,#5),.UNSPECIFIED.,.F.,.F.);\
                  #12=QUASI_UNIFORM_CURVE('',2,(#1,#2,#3,#4,#5),.UNSPECIFIED.,.F.,.F.);\
                  #13=BEZIER_CURVE('',2,(#1,#2,#3,#4,#5),.UNSPECIFIED.,.F.,.F.);";
    // A 3 × 3 net on the grid (i, j, 0) with its middle point raised to z = 4.
    let net = "((#1,#6,#9),(#2,#7,#16),(#10,#8,#3))";
    let surfaces = format!(
        "#14=UNIFORM_SURFACE('',2,2,{net},.UNSPECIFIED.,.F.,.F.,.F.);\
         #15=QUASI_UNIFORM_SURFACE('',2,2,{net},.UNSPECIFIED.,.F.,.F.,.F.);"
    );
    let text =
        format!("ISO-10303-21;\nDATA;\n{points}{curves}{surfaces}ENDSEC;\nEND-ISO-10303-21;\n");
    let file = super::part21::parse(text.as_bytes()).unwrap();
    let mut g = super::geometry::Geometry::new(&file);
    let at =
        |c: &Curve, t: f64, p: Vec3| assert!(c.eval(t).dist(p) < 1e-14, "{t}: {:?}", c.eval(t));
    // Uniform: a span starts halfway between its first two control points, and its middle lies at
    // (P0 + 6 P1 + P2) / 8.
    let uniform = g.curve(11).unwrap();
    assert_eq!(uniform.domain(), (0.0, 3.0));
    at(&uniform, 0.0, Vec3::new(0.5, 0.0, 0.0));
    at(&uniform, 0.5, Vec3::new(1.0, 0.25, 0.0));
    at(&uniform, 3.0, Vec3::new(4.5, 2.0, 0.5));
    let quasi = g.curve(12).unwrap();
    assert_eq!(quasi.domain(), (0.0, 3.0));
    at(&quasi, 0.0, Vec3::new(0.0, 0.0, 0.0));
    at(&quasi, 3.0, Vec3::new(5.0, 3.0, 1.0));
    // Bézier: two pieces joined at P2, the first's middle at (P0 + 2 P1 + P2) / 4.
    let bezier = g.curve(13).unwrap();
    assert_eq!(bezier.domain(), (0.0, 2.0));
    at(&bezier, 0.5, Vec3::new(1.0, 0.5, 0.0));
    at(&bezier, 1.0, Vec3::new(2.0, 2.0, 0.0));
    at(&bezier, 2.0, Vec3::new(5.0, 3.0, 1.0));
    // A uniform net of one span starts at the mean of its first 2 × 2 points; a quasi-uniform one
    // of one span is a Bézier patch, through its corner and, in the middle, a quarter of z = 4.
    let uniform = g.surface(14).unwrap();
    assert!(uniform.eval(0.0, 0.0).dist(Vec3::new(0.5, 0.5, 1.0)) < 1e-14);
    let quasi = g.surface(15).unwrap();
    assert!(quasi.eval(0.0, 0.0).dist(Vec3::new(0.0, 0.0, 0.0)) < 1e-14);
    assert!(quasi.eval(0.5, 0.5).dist(Vec3::new(1.0, 1.0, 1.0)) < 1e-14);
}

/// `text`, written without pcurves, with each edge's line and each face's plane given in a form
/// whose kind implies its knots, the kinds taken in turn: a line as the degree-1 spline between its
/// edge's ends, a plane as the bilinear patch over a square of side 200 about its origin.
fn with_implied_knots(text: &str) -> String {
    let file = super::part21::parse(text.as_bytes()).unwrap();
    let point = |r: u64| -> Vec3 {
        let p = file.get(r).unwrap().part("CARTESIAN_POINT").unwrap()[1]
            .as_list()
            .unwrap();
        Vec3::new(
            p[0].as_f64().unwrap(),
            p[1].as_f64().unwrap(),
            p[2].as_f64().unwrap(),
        )
    };
    let dir = |r: u64| -> Vec3 {
        let d = file.get(r).unwrap().part("DIRECTION").unwrap()[1]
            .as_list()
            .unwrap();
        Vec3::new(
            d[0].as_f64().unwrap(),
            d[1].as_f64().unwrap(),
            d[2].as_f64().unwrap(),
        )
    };
    let id_of = |l: &str| l.strip_prefix('#')?.split('=').next()?.parse::<u64>().ok();
    let mut next = text.lines().filter_map(id_of).max().unwrap() + 1;
    let mut added = String::new();
    let mut new_point = |q: Vec3, added: &mut String| {
        let id = next;
        next += 1;
        *added += &format!(
            "#{id}=CARTESIAN_POINT('',({:?},{:?},{:?}));\n",
            q.x, q.y, q.z
        );
        id
    };
    // Each line's edge, by the line's id: its two vertices' points.
    let mut ends: std::collections::HashMap<u64, (u64, u64)> = Default::default();
    for line in text.lines() {
        let Some(id) = id_of(line) else { continue };
        let Some(e) = file.get(id).and_then(|e| e.part("EDGE_CURVE")) else {
            continue;
        };
        let vertex = |k: usize| {
            let v = e[k].as_ref().unwrap();
            file.get(v).unwrap().part("VERTEX_POINT").unwrap()[1]
                .as_ref()
                .unwrap()
        };
        let geometry = e[3].as_ref().unwrap();
        if e[4].as_bool() == Some(true) {
            ends.insert(geometry, (vertex(1), vertex(2)));
        } else {
            ends.insert(geometry, (vertex(2), vertex(1)));
        }
    }
    let curves = [
        "UNIFORM_CURVE",
        "QUASI_UNIFORM_CURVE",
        "BEZIER_CURVE",
        "RATIONAL",
    ];
    let surfaces = ["UNIFORM_SURFACE", "QUASI_UNIFORM_SURFACE", "BEZIER_SURFACE"];
    let (mut nc, mut ns) = (0, 0);
    let mut out = Vec::new();
    for line in text.lines() {
        let id = id_of(line);
        let e = id.and_then(|id| file.get(id));
        if let (Some(id), Some(_)) = (id, e.and_then(|e| e.part("LINE"))) {
            let (a, b) = ends[&id];
            let kind = curves[nc % curves.len()];
            nc += 1;
            out.push(if kind == "RATIONAL" {
                format!(
                    "#{id}=(BOUNDED_CURVE() B_SPLINE_CURVE(1,(#{a},#{b}),.UNSPECIFIED.,.F.,.F.) \
                     CURVE() GEOMETRIC_REPRESENTATION_ITEM() QUASI_UNIFORM_CURVE() \
                     RATIONAL_B_SPLINE_CURVE((2.,2.)) REPRESENTATION_ITEM(''));"
                )
            } else {
                format!("#{id}={kind}('',1,(#{a},#{b}),.UNSPECIFIED.,.F.,.F.);")
            });
            continue;
        }
        if let (Some(id), Some(p)) = (id, e.and_then(|e| e.part("PLANE"))) {
            let place = file.get(p[1].as_ref().unwrap()).unwrap();
            let axes = place.part("AXIS2_PLACEMENT_3D").unwrap();
            let o = point(axes[1].as_ref().unwrap());
            let z = dir(axes[2].as_ref().unwrap());
            let x = dir(axes[3].as_ref().unwrap());
            let y = z.cross(x);
            let mut net = Vec::new();
            for i in [-100.0, 100.0] {
                let row: Vec<String> = [-100.0, 100.0]
                    .iter()
                    .map(|&j| format!("#{}", new_point(o + x * i + y * j, &mut added)))
                    .collect();
                net.push(format!("({})", row.join(",")));
            }
            let kind = surfaces[ns % surfaces.len()];
            ns += 1;
            out.push(format!(
                "#{id}={kind}('',1,1,({}),.UNSPECIFIED.,.F.,.F.,.F.);",
                net.join(",")
            ));
            continue;
        }
        out.push(line.to_string());
    }
    assert!(nc >= curves.len() && ns >= surfaces.len());
    out.join("\n")
        .replace("ENDSEC;\nEND-ISO", &format!("{added}ENDSEC;\nEND-ISO"))
}

#[test]
fn splines_whose_kind_implies_their_knots_are_read() {
    // The block's every edge and face, in turn uniform, quasi-uniform, Bézier and rational: the
    // same block, 10 by 6 by 3, read back.
    let b = block("blk", [0.0, 0.0, 0.0], [10.0, 6.0, 3.0]);
    let text = write_with(&b, "blk", &WriteOptions::default().pcurves(false)).unwrap();
    let implied = with_implied_knots(&text);
    for kind in [
        "UNIFORM_CURVE(",
        "QUASI_UNIFORM_SURFACE(",
        "BEZIER_CURVE(",
        "QUASI_UNIFORM_CURVE()",
    ] {
        assert!(implied.contains(kind), "{kind}");
    }
    let read_back = read(implied.as_bytes(), "imp").unwrap();
    assert!(read_back.rejected.is_empty(), "{:?}", read_back.rejected);
    let c = read_back.brep;
    assert!(check::check(&c).is_empty(), "{:?}", check::check(&c));
    assert_eq!(c.face_count(), 6);
    assert!((measure::volume(&c) - 180.0).abs() <= 1e-9 * 180.0);
    assert!((measure::area(&c) - 2.0 * (60.0 + 30.0 + 18.0)).abs() <= 1e-9 * 216.0);
}

/// `text` with the placement of each `CONICAL_SURFACE` moved onto its other nappe: the same apex,
/// the axis turned over, and the origin as far past the apex as it was before it, so the cone is
/// the same point set and the face on it lies on the nappe the new placement calls far.
fn cones_from_their_far_nappes(text: &str) -> String {
    let file = super::part21::parse(text.as_bytes()).unwrap();
    let ids = text
        .lines()
        .filter_map(|l| l.strip_prefix('#')?.split('=').next()?.parse::<u64>().ok());
    let mut next = ids.max().unwrap() + 1;
    let mut added = String::new();
    let mut out = text.to_string();
    for line in text.lines() {
        let Some((id, body)) = line.strip_prefix('#').and_then(|l| l.split_once('=')) else {
            continue;
        };
        let Some(args) = body.strip_prefix("CONICAL_SURFACE('',#") else {
            continue;
        };
        let place: u64 = args.split(',').next().unwrap().parse().unwrap();
        let cone = file
            .get(id.parse().unwrap())
            .unwrap()
            .part("CONICAL_SURFACE")
            .unwrap();
        let (radius, half_angle) = (cone[2].as_f64().unwrap(), cone[3].as_f64().unwrap());
        let axes = file
            .get(place)
            .unwrap()
            .part("AXIS2_PLACEMENT_3D")
            .unwrap()
            .to_vec();
        let point = |r: u64| -> Vec3 {
            let p = file.get(r).unwrap().part("CARTESIAN_POINT").unwrap()[1]
                .as_list()
                .unwrap();
            Vec3::new(
                p[0].as_f64().unwrap(),
                p[1].as_f64().unwrap(),
                p[2].as_f64().unwrap(),
            )
        };
        let dir = |r: u64| -> Vec3 {
            let d = file.get(r).unwrap().part("DIRECTION").unwrap()[1]
                .as_list()
                .unwrap();
            Vec3::new(
                d[0].as_f64().unwrap(),
                d[1].as_f64().unwrap(),
                d[2].as_f64().unwrap(),
            )
        };
        let (o, z) = (
            point(axes[1].as_ref().unwrap()),
            dir(axes[2].as_ref().unwrap()),
        );
        let x = axes[3].as_ref().unwrap();
        let apex = o - z * (radius / half_angle.tan());
        let o2 = apex * 2.0 - o;
        let (p, d, a) = (next, next + 1, next + 2);
        next += 3;
        added += &format!(
            "#{p}=CARTESIAN_POINT('',({:?},{:?},{:?}));\n#{d}=DIRECTION('',({:?},{:?},{:?}));\n\
             #{a}=AXIS2_PLACEMENT_3D('',#{p},#{d},#{x});\n",
            o2.x, o2.y, o2.z, -z.x, -z.y, -z.z
        );
        out = out.replace(
            &format!("#{id}=CONICAL_SURFACE('',#{place},"),
            &format!("#{id}=CONICAL_SURFACE('',#{a},"),
        );
    }
    out.replace("ENDSEC;\nEND-ISO", &format!("{added}ENDSEC;\nEND-ISO"))
}

#[test]
fn a_face_on_the_far_nappe_of_its_cone_is_read_and_cut() {
    // A frustum of radii 6 and 3 and height 5, its cone written about the other nappe: its volume,
    // 105π, read back, and 100π once bored through by a coaxial rod of radius 1.
    let pts: Vec<Vec2> = [(0.0, 0.0), (6.0, 0.0), (3.0, 5.0), (0.0, 5.0)]
        .iter()
        .map(|&(x, y)| Vec2::new(x, y))
        .collect();
    let p = Profile::polygon(&pts, &["a", "b", "c", "d"]);
    let frustum = revolve(&p, &xz(), Vec3::ZERO, Vec3::Z, TAU, "f").unwrap();
    for pcurves in [true, false] {
        let text = write_with(&frustum, "f", &WriteOptions::default().pcurves(pcurves)).unwrap();
        let far = cones_from_their_far_nappes(&text);
        assert_ne!(far, text);
        let read_back = read(far.as_bytes(), "far").unwrap();
        assert!(read_back.rejected.is_empty(), "{:?}", read_back.rejected);
        // Read whole, and sewn from its faces alone.
        for b in [read_back.brep, super::sew(far.as_bytes(), 1e-6).unwrap()] {
            assert!(check::check(&b).is_empty(), "{:?}", check::check(&b));
            assert!((measure::volume(&b) - 105.0 * PI).abs() <= 1e-9 * 105.0 * PI);
            let rod = extrude(
                &Profile::circle(Vec2::ZERO, 2.0),
                &Frame::new(Vec3::new(0.0, 0.0, -1.0), Vec3::Z, Vec3::X).unwrap(),
                7.0,
                "rod",
            )
            .unwrap();
            let bored = boolean::difference(&b, &rod).unwrap();
            assert!(check::check(&bored).is_empty());
            assert!((measure::volume(&bored) - 100.0 * PI).abs() <= 1e-9 * 100.0 * PI);
            // A coaxial rod of radius 4 cuts the cone face itself, at the height 10/3 where the
            // frustum's radius 6 − 0.6·z is 4: what is left is π times the integral of
            // (6 − 0.6·z)² − 16 up to there, 280π/9.
            let rod = extrude(
                &Profile::circle(Vec2::ZERO, 8.0),
                &Frame::new(Vec3::new(0.0, 0.0, -1.0), Vec3::Z, Vec3::X).unwrap(),
                7.0,
                "rod",
            )
            .unwrap();
            let ring = boolean::difference(&b, &rod).unwrap();
            assert!(check::check(&ring).is_empty());
            let want = 280.0 * PI / 9.0;
            assert!((measure::volume(&ring) - want).abs() <= 1e-9 * want);
            // The plane x = 1 along the axis cuts the cone face in a hyperbola. Each section of
            // what the cut leaves is the disc of radius R = 6 − 0.6·z less its segment past x = 1,
            // R²·acos(1/R) − √(R² − 1), integrated up the height by Simpson's rule.
            let cut = boolean::difference(&b, &block("x", [1.0, -10.0, -1.0], [10.0, 10.0, 6.0]))
                .unwrap();
            assert!(check::check(&cut).is_empty());
            let section = |z: f64| {
                let r = 6.0 - 0.6 * z;
                PI * r * r - (r * r * (1.0 / r).acos() - (r * r - 1.0).sqrt())
            };
            let n = 2000;
            let h = 5.0 / n as f64;
            let simpson: f64 = (0..=n)
                .map(|k| {
                    let w = if k == 0 || k == n {
                        1.0
                    } else if k % 2 == 1 {
                        4.0
                    } else {
                        2.0
                    };
                    w * section(k as f64 * h)
                })
                .sum::<f64>()
                * h
                / 3.0;
            assert!((measure::volume(&cut) - simpson).abs() <= 1e-9 * simpson);
        }
    }
}

#[test]
fn an_edge_of_no_length_between_coincident_vertices_is_dropped() {
    // A 10 mm cube whose corner at the origin is written as two vertices at one point, joined by a
    // line edge of length 0 that the faces x = 0 and y = 0 use, one each way. Joined, the two
    // vertices are the cube's corner and the edge is no edge: the cube is read with its 8
    // vertices, 12 edges and 6 faces, volume 1000 and area 600, as the file's own geometry.
    let text = include_bytes!("fixtures/zero-length-edge.step");
    for options in [ReadOptions::default(), ReadOptions::default().strict(true)] {
        let imp = read_with(text, "cube", &options).unwrap();
        assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
        assert!(imp.healed.is_empty() && imp.tolerant.is_empty());
        let b = &imp.brep;
        assert!(check::check(b).is_empty(), "{:?}", check::check(b));
        assert_eq!(
            (b.vertex_ids().count(), b.edge_ids().count(), b.face_count()),
            (8, 12, 6)
        );
        assert!((measure::volume(b) - 1000.0).abs() <= 1e-9 * 1000.0);
        assert!((measure::area(b) - 600.0).abs() <= 1e-9 * 600.0);
    }
    round_trip(&read(text, "cube").unwrap().brep);
}

#[test]
fn a_band_on_a_surface_of_revolution_gets_a_seam() {
    // Tubes between z = 0 and 10 about the z axis whose outer wall is a band bounded by its two
    // rims and no seam; fixtures/band-fixtures.py writes them and works their volumes.
    //
    // The outer wall the revolution of the cubic Bezier r = 10, 13, 13, 10 at z = 10 t round a
    // cylinder of radius 6: π (∫ r² dz − 6² · 10) = 967π. The rims are written as an outer and an
    // inner bound, as two outer bounds, and with their vertices a quarter turn apart. A window a
    // twelfth of a turn wide through the wall between z = 3 and 7 takes (π / 12) ∫₃⁷ (r² − 36) dz:
    // two off the seam's line leave 66967579π / 75000, one across it, from 160° to 190°, leaves
    // 139492579π / 150000. The two through a tube of radii 6 and 10, its outer wall a cylinder,
    // leave 1792π / 3.
    //
    // The outer wall the extrusion of the ellipse of semi-axes 10 and 8 round a cylinder of radius
    // 5: π (80 − 25) · 10 = 550π; with a window from 135° to 225°, across the seam's line,
    // 575π − 320 atan(5/4). The extrusion of a rational circle of radius 10 whose parameter runs
    // over [2, 3]: π (100 − 25) · 10 = 750π; with a window from 170° to 200°, 725π.
    let cases: [(&str, &[u8], f64); 10] = [
        (
            "band",
            include_bytes!("fixtures/revolution-band.step"),
            967.0 * PI,
        ),
        (
            "two outer",
            include_bytes!("fixtures/revolution-band-two-outer.step"),
            967.0 * PI,
        ),
        (
            "apart",
            include_bytes!("fixtures/revolution-band-apart.step"),
            967.0 * PI,
        ),
        (
            "windows",
            include_bytes!("fixtures/revolution-band-windows.step"),
            66967579.0 * PI / 75000.0,
        ),
        (
            "window on the seam",
            include_bytes!("fixtures/revolution-band-window-on-seam.step"),
            139492579.0 * PI / 150000.0,
        ),
        (
            "cylinder",
            include_bytes!("fixtures/cylinder-band-windows.step"),
            1792.0 * PI / 3.0,
        ),
        (
            "extrusion",
            include_bytes!("fixtures/extrusion-band.step"),
            550.0 * PI,
        ),
        (
            "extrusion window",
            include_bytes!("fixtures/extrusion-band-window.step"),
            575.0 * PI - 320.0 * 1.25f64.atan(),
        ),
        (
            "nurbs extrusion",
            include_bytes!("fixtures/extrusion-nurbs-band.step"),
            750.0 * PI,
        ),
        (
            "nurbs extrusion window",
            include_bytes!("fixtures/extrusion-nurbs-band-window.step"),
            725.0 * PI,
        ),
    ];
    for (name, text, volume) in cases {
        let imp = read(text, "tube").unwrap();
        assert!(imp.rejected.is_empty(), "{name}: {:?}", imp.rejected);
        let b = &imp.brep;
        assert!(check::check(b).is_empty(), "{name}: {:?}", check::check(b));
        let v = measure::volume(b);
        assert!(
            (v - volume).abs() <= 1e-9 * volume,
            "{name}: volume {v}, not {volume}"
        );
        round_trip(b);
        // A seam across a window would make the band's boundary overlap itself, which the mesher
        // refuses. Every mesh point lies within the deflection of the boundary, so the mesh's
        // volume is within the area times the deflection of the solid's.
        let linear = 0.01;
        let m = mesh::tessellate(b, &mesh::MeshOptions::within(linear))
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let rep = mesh::report(b, &m, b.solid_ids().next().unwrap());
        assert!(rep.watertight(), "{name}: {rep:?}");
        assert_eq!(rep.degenerate_triangles, 0, "{name}");
        let area = measure::area(b);
        assert!(
            (rep.volume - volume).abs() <= area * linear,
            "{name}: mesh volume {}",
            rep.volume
        );
    }
}

/// The octant x, y, z ≥ 0 of a ball of radius 10 at the origin: three quarter discs and an
/// eighth of the sphere, bounded by three quarter circles, its sphere framed by `axis` and `ref`.
fn octant(axis: [f64; 3], reference: [f64; 3]) -> String {
    let d = |v: [f64; 3]| format!("({:?},{:?},{:?})", v[0], v[1], v[2]);
    step_file(&format!(
        "#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(1.,0.,0.));
#3=DIRECTION('',(0.,1.,0.));
#4=DIRECTION('',(0.,0.,1.));
#5=DIRECTION('',(-1.,0.,0.));
#6=DIRECTION('',(0.,-1.,0.));
#7=DIRECTION('',(0.,0.,-1.));
#10=VERTEX_POINT('',#1);
#11=CARTESIAN_POINT('',(10.,0.,0.));
#12=VERTEX_POINT('',#11);
#13=CARTESIAN_POINT('',(0.,10.,0.));
#14=VERTEX_POINT('',#13);
#15=CARTESIAN_POINT('',(0.,0.,10.));
#16=VERTEX_POINT('',#15);
#20=LINE('',#1,#25);
#21=LINE('',#1,#26);
#22=LINE('',#1,#27);
#25=VECTOR('',#2,1.);
#26=VECTOR('',#3,1.);
#27=VECTOR('',#4,1.);
#30=EDGE_CURVE('',#10,#12,#20,.T.);
#31=EDGE_CURVE('',#10,#14,#21,.T.);
#32=EDGE_CURVE('',#10,#16,#22,.T.);
#33=AXIS2_PLACEMENT_3D('',#1,#4,#2);
#34=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#35=AXIS2_PLACEMENT_3D('',#1,#3,#4);
#36=CIRCLE('',#33,10.);
#37=CIRCLE('',#34,10.);
#38=CIRCLE('',#35,10.);
#40=EDGE_CURVE('',#12,#14,#36,.T.);
#41=EDGE_CURVE('',#14,#16,#37,.T.);
#42=EDGE_CURVE('',#16,#12,#38,.T.);
#50=ORIENTED_EDGE('',*,*,#31,.T.);
#51=ORIENTED_EDGE('',*,*,#40,.F.);
#52=ORIENTED_EDGE('',*,*,#30,.F.);
#53=EDGE_LOOP('',(#50,#51,#52));
#54=FACE_OUTER_BOUND('',#53,.T.);
#55=AXIS2_PLACEMENT_3D('',#1,#7,#2);
#56=PLANE('',#55);
#57=ADVANCED_FACE('',(#54),#56,.T.);
#60=ORIENTED_EDGE('',*,*,#32,.T.);
#61=ORIENTED_EDGE('',*,*,#41,.F.);
#62=ORIENTED_EDGE('',*,*,#31,.F.);
#63=EDGE_LOOP('',(#60,#61,#62));
#64=FACE_OUTER_BOUND('',#63,.T.);
#65=AXIS2_PLACEMENT_3D('',#1,#5,#3);
#66=PLANE('',#65);
#67=ADVANCED_FACE('',(#64),#66,.T.);
#70=ORIENTED_EDGE('',*,*,#30,.T.);
#71=ORIENTED_EDGE('',*,*,#42,.F.);
#72=ORIENTED_EDGE('',*,*,#32,.F.);
#73=EDGE_LOOP('',(#70,#71,#72));
#74=FACE_OUTER_BOUND('',#73,.T.);
#75=AXIS2_PLACEMENT_3D('',#1,#6,#4);
#76=PLANE('',#75);
#77=ADVANCED_FACE('',(#74),#76,.T.);
#80=ORIENTED_EDGE('',*,*,#40,.T.);
#81=ORIENTED_EDGE('',*,*,#41,.T.);
#82=ORIENTED_EDGE('',*,*,#42,.T.);
#83=EDGE_LOOP('',(#80,#81,#82));
#84=FACE_OUTER_BOUND('',#83,.T.);
#85=DIRECTION('',{});
#86=DIRECTION('',{});
#87=AXIS2_PLACEMENT_3D('',#1,#85,#86);
#88=SPHERICAL_SURFACE('',#87,10.);
#89=ADVANCED_FACE('',(#84),#88,.T.);
#90=CLOSED_SHELL('',(#57,#67,#77,#89));
#91=MANIFOLD_SOLID_BREP('',#90);",
        d(axis),
        d(reference)
    ))
}

#[test]
fn a_sphere_face_whose_edge_passes_through_a_pole_is_split_there() {
    // An eighth of a ball of radius 10: volume π·10³/6, and area three quarter discs and an eighth
    // of the sphere, 3π·10²/4 + π·10²/2 = 125π. The file frames its sphere five ways: its north
    // pole in the middle of the quarter circle in y = 0, its south pole there, its north pole in
    // the middle of the quarter circle in z = 0, at a corner of the face, and inside the face.
    // Six edges with curves bound the octant. One that passes through a pole between its vertices
    // is split there, at a vertex of its own, for the plane face beside it too: seven. A pole at
    // a corner needs nothing, and one inside the face is reached by a seam: seven again.
    let h = std::f64::consts::FRAC_1_SQRT_2;
    let (volume, area) = (1000.0 * PI / 6.0, 125.0 * PI);
    let cases = [
        ([h, 0.0, h], [h, 0.0, -h], 1),
        ([-h, 0.0, -h], [h, 0.0, -h], 1),
        ([h, h, 0.0], [0.0, 0.0, 1.0], 1),
        ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 0),
        ([1.0, 1.0, 1.0], [1.0, -1.0, 0.0], 1),
    ];
    for (axis, reference, more) in cases {
        let imp = read(octant(axis, reference).as_bytes(), "imp").unwrap();
        assert!(imp.rejected.is_empty(), "{axis:?}: {:?}", imp.rejected);
        let b = imp.brep;
        assert!(
            check::check(&b).is_empty(),
            "{axis:?}: {:?}",
            check::check(&b)
        );
        let (v, a) = (measure::volume(&b), measure::area(&b));
        assert!((v - volume).abs() < 1e-9 * volume, "{axis:?}: volume {v}");
        assert!((a - area).abs() < 1e-9 * area, "{axis:?}: area {a}");
        let edges = b.edge_ids().filter(|&e| b.edge(e).curve.is_some()).count();
        assert_eq!(edges, 6 + more, "{axis:?}");
        round_trip(&b);
    }
}

#[test]
fn a_spline_face_that_stands_still_along_an_end_line_is_no_pole() {
    // A prism of height 4 over the region bounded by the segment O A, O = (0, 0) and A = (10, 0), a
    // cubic Bezier from A to B = (5, 10), and the segment B O. Its side on the Bezier is a spline,
    // cubic in u along the Bezier and linear in v up the prism. Where two end control points of the
    // cubic coincide, the side's u derivative vanishes along the whole end line u = 1 (or u = 0),
    // though that line is the segment up from B (or A), not a point. The side is also written
    // with its u, its v, and both running the other way. The volume is 4 times the region's area,
    // ∮x dy over its polynomial boundary, exact: 230 with the end still, 248 with the start
    // still, 251 for four distinct control points. Where the last three control points coincide
    // the cubic is the segment A B, run as B + s³(A − B) back from B, so the side's second
    // derivative across its end line vanishes too: the triangle O A B, 200, also written cubic in
    // v. The four planar faces are the region top and
    // bottom, 10 × 4 on y = 0, and |OB| × 4 = 20√5 through O and B.
    let cases: [(&str, &[u8], f64); 8] = [
        ("end", include_bytes!("fixtures/still-line-end.step"), 230.0),
        (
            "end, u reversed",
            include_bytes!("fixtures/still-line-end-u.step"),
            230.0,
        ),
        (
            "end, v reversed",
            include_bytes!("fixtures/still-line-end-v.step"),
            230.0,
        ),
        (
            "end, both reversed",
            include_bytes!("fixtures/still-line-end-uv.step"),
            230.0,
        ),
        (
            "start",
            include_bytes!("fixtures/still-line-start.step"),
            248.0,
        ),
        (
            "regular",
            include_bytes!("fixtures/still-line-regular.step"),
            251.0,
        ),
        (
            "three",
            include_bytes!("fixtures/still-line-three.step"),
            200.0,
        ),
        (
            "three, cubic in v",
            include_bytes!("fixtures/still-line-three-t.step"),
            200.0,
        ),
    ];
    for (name, text, volume) in cases {
        // The file agrees with itself exactly: a strict read, which moves nothing, keeps it, and
        // the plain read heals nothing.
        let strict = read_with(text, "prism", &ReadOptions::default().strict(true)).unwrap();
        assert!(strict.rejected.is_empty(), "{name}: {:?}", strict.rejected);
        let imp = read(text, "prism").unwrap();
        assert!(imp.rejected.is_empty(), "{name}: {:?}", imp.rejected);
        assert!(imp.healed.is_empty(), "{name}: {:?}", imp.healed);
        let b = &imp.brep;
        assert!(check::check(b).is_empty(), "{name}: {:?}", check::check(b));
        assert_eq!(b.face_count(), 5, "{name}");
        let v = measure::volume(b);
        assert!(
            (v - volume).abs() <= 1e-12 * volume,
            "{name}: volume {v}, not {volume}"
        );
        let mut planar: Vec<f64> = b
            .face_ids()
            .filter(|&f| matches!(b.face(f).surface, Surface::Plane(_)))
            .map(|f| measure::face_area(b, f))
            .collect();
        planar.sort_by(f64::total_cmp);
        let mut expected = [volume / 4.0, volume / 4.0, 40.0, 20.0 * 5f64.sqrt()];
        expected.sort_by(f64::total_cmp);
        assert_eq!(planar.len(), 4, "{name}");
        for (a, e) in planar.iter().zip(expected) {
            assert!(
                (a - e).abs() <= 1e-12 * e,
                "{name}: planar areas {planar:?}, not {expected:?}"
            );
        }
        round_trip(b);
    }
}

/// `text`, a file of solids written one after another, as one solid with voids: the first solid's
/// shell is its outer shell, and the shell of each other one, turned by an `ORIENTED_CLOSED_SHELL`
/// as a void's is, one of its voids. The file's representation then lists that solid alone. A file
/// that writes a compound of bodies this way lists as voids bodies that lie outside the outer one.
fn written_as_voids(text: &str) -> (String, Vec<u64>) {
    let solids: Vec<(String, String)> = text
        .lines()
        .filter_map(|l| {
            let (id, rest) = l.split_once("=MANIFOLD_SOLID_BREP('',")?;
            Some((id.to_string(), rest.trim_end_matches(");").to_string()))
        })
        .collect();
    let (first, voids) = solids.split_first().unwrap();
    let listed: Vec<&str> = voids.iter().map(|v| v.0.as_str()).collect();
    let mut out = String::new();
    for l in text.lines() {
        let at = |id: &str| l.starts_with(&format!("{id}="));
        let line = if at(&first.0) {
            format!(
                "{}=BREP_WITH_VOIDS('',{},({}));",
                first.0,
                first.1,
                listed.join(",")
            )
        } else if let Some(v) = voids.iter().find(|v| at(&v.0)) {
            format!("{}=ORIENTED_CLOSED_SHELL('',*,{},.F.);", v.0, v.1)
        } else if l.contains("=ADVANCED_BREP_SHAPE_REPRESENTATION(") {
            let mut l = l.to_string();
            for v in &listed {
                l = l
                    .replace(&format!(",{v},"), ",")
                    .replace(&format!(",{v})"), ")");
            }
            l
        } else {
            l.to_string()
        };
        out.push_str(&line);
        out.push('\n');
    }
    let ids = listed.iter().map(|v| v[1..].parse().unwrap()).collect();
    (out, ids)
}

/// The volume of each solid of `b`, in its order.
fn solid_volumes(b: &Brep) -> Vec<f64> {
    b.solid_ids()
        .map(|s| measure::volume(&boolean::solids_of(b, &[s])))
        .collect()
}

/// A cube of 10 mm listed as the void of another beside it, touching it along a face or 2 mm
/// apart: neither is a void of the other, so each is read as a solid of its own, of 1000 mm³, and
/// the import says which shells were. A strict read rejects the solid, saying why.
#[test]
fn voids_outside_the_outer_shell_are_read_as_solids_of_their_own() {
    for gap in [0.0, 2.0] {
        let mut b = block("a", [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
        b.append(&block(
            "b",
            [10.0 + gap, 0.0, 0.0],
            [20.0 + gap, 10.0, 10.0],
        ));
        let (text, voids) = written_as_voids(&write(&b, "board").unwrap());
        assert_eq!(instances(&text, "BREP_WITH_VOIDS"), 1);
        assert_eq!(instances(&text, "MANIFOLD_SOLID_BREP"), 0);
        let imp = read(text.as_bytes(), "imp").unwrap();
        assert!(imp.rejected.is_empty(), "gap {gap}: {:?}", imp.rejected);
        assert!(check::check_each(&imp.brep).is_empty());
        let volumes = solid_volumes(&imp.brep);
        assert_eq!(volumes.len(), 2);
        for v in volumes {
            assert!((v - 1000.0).abs() < 1e-9, "gap {gap}: {v}");
        }
        assert_eq!(imp.separated.len(), 1);
        let s = &imp.separated[0];
        assert_eq!((s.index, &s.shells), (0, &voids));
        assert_eq!(imp.bodies.len(), 2);
        assert!(imp.bodies.iter().all(|x| x.step_id == s.step_id));
        // As one part, the two cubes are one box of 2000 mm³ where they touch.
        let part = super::read_part(text.as_bytes(), "imp").unwrap();
        assert!((measure::volume(&part) - 2000.0).abs() < 1e-9);
        let solids = if gap == 0.0 { 1 } else { 2 };
        assert_eq!(part.solid_ids().count(), solids);
        let strict = read_with(text.as_bytes(), "imp", &ReadOptions::default().strict(true));
        let strict = strict.unwrap();
        assert_eq!(strict.brep.solid_ids().count(), 0);
        assert!(strict.separated.is_empty());
        let why = &strict.rejected[0];
        let named = format!("void shell #{} lies outside the outer shell", voids[0]);
        assert!(why.contains(&named), "{why}");
    }
}

/// Of two voids listed, one inside the outer cube and one beside it: the first stays a void, and
/// the second is a solid of its own. The cube of 10 mm less the void of 3 mm holds 973 mm³.
#[test]
fn a_void_inside_stays_a_void_beside_one_read_as_a_solid() {
    let mut b = block("a", [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
    b.append(&block("v", [3.0, 3.0, 3.0], [6.0, 6.0, 6.0]));
    b.append(&block("b", [10.0, 0.0, 0.0], [20.0, 10.0, 10.0]));
    let (text, voids) = written_as_voids(&write(&b, "board").unwrap());
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    assert!(check::check_each(&imp.brep).is_empty());
    let volumes = solid_volumes(&imp.brep);
    assert_eq!(volumes.len(), 2);
    assert!((volumes[0] - 973.0).abs() < 1e-9, "{volumes:?}");
    assert!((volumes[1] - 1000.0).abs() < 1e-9, "{volumes:?}");
    assert_eq!(imp.separated.len(), 1);
    assert_eq!(imp.separated[0].shells, voids[1..]);
    // The void inside alone is a solid with a void, as the file means it.
    let mut c = block("a", [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
    c.append(&block("v", [3.0, 3.0, 3.0], [6.0, 6.0, 6.0]));
    let (text, _) = written_as_voids(&write(&c, "board").unwrap());
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty() && imp.separated.is_empty());
    assert!((measure::volume(&imp.brep) - 973.0).abs() < 1e-9);
    assert_eq!(imp.brep.solid_ids().count(), 1);
}

/// A block listed as the void of a cube it passes through, as one body of a compound may through
/// another: it lies partly outside the cube, so it is no void, and is read as a solid of its own.
/// The cube of 10 mm holds 1000 mm³ and the block of 4 × 3 × 3 mm, 36 mm³.
#[test]
fn a_void_passing_through_the_outer_shell_is_read_as_a_solid_of_its_own() {
    let mut b = block("a", [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
    b.append(&block("v", [8.0, 3.0, 3.0], [12.0, 6.0, 6.0]));
    let (text, voids) = written_as_voids(&write(&b, "board").unwrap());
    let imp = read(text.as_bytes(), "imp").unwrap();
    assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
    // Read as a void, the block would make a solid that passes through itself.
    assert!(imp.crossed.is_empty());
    assert!(check::check_each(&imp.brep).is_empty());
    let volumes = solid_volumes(&imp.brep);
    assert_eq!(volumes.len(), 2);
    assert!((volumes[0] - 1000.0).abs() < 1e-9, "{volumes:?}");
    assert!((volumes[1] - 36.0).abs() < 1e-9, "{volumes:?}");
    assert_eq!(imp.separated.len(), 1);
    assert_eq!(imp.separated[0].shells, voids);
}

/// Two cubes of 10 mm, 2 mm apart, the second listed as a void of the first, with a vertex of the
/// second moved off its corner: 20 nm, which healing moves back, or 30 µm in a file that states
/// 10 µm, which it keeps with a tolerance. The second cube is read as a solid of its own, and the
/// healing or the tolerance is reported for it, the second solid read, not the first.
#[test]
fn healing_and_tolerances_of_a_void_read_as_a_solid_name_that_solid() {
    let mut b = block("a", [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
    b.append(&block("b", [12.0, 0.0, 0.0], [22.0, 10.0, 10.0]));
    let (text, voids) = written_as_voids(&write(&b, "board").unwrap());
    let vp = text
        .lines()
        .rfind(|l| l.contains("=VERTEX_POINT("))
        .unwrap();
    for (d, stated) in [(2e-5, None), (0.03, Some(0.01))] {
        let mut text = with_point_moved(&text, vp, Vec3::new(d, 0.0, 0.0));
        if let Some(u) = stated {
            text = with_uncertainty(&text, u);
        }
        let imp = read(text.as_bytes(), "imp").unwrap();
        assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
        assert_eq!(imp.separated.len(), 1);
        assert_eq!(imp.separated[0].shells, voids);
        assert_eq!(imp.brep.solid_ids().count(), 2);
        let volumes = solid_volumes(&imp.brep);
        assert!((volumes[0] - 1000.0).abs() < 1e-9, "{volumes:?}");
        // A corner moved by d changes the volume by at most d times the area of its three faces.
        assert!((volumes[1] - 1000.0).abs() <= d * 300.0, "{volumes:?}");
        if stated.is_none() {
            assert_eq!(imp.healed.len(), 1, "{:?}", imp.healed);
            let h = &imp.healed[0];
            assert_eq!(h.index, 1, "{h:?}");
            assert!(h.at.x >= 12.0 - 1e-3, "{h:?}");
            assert!(imp.tolerant.is_empty(), "{:?}", imp.tolerant);
        } else {
            assert_eq!(imp.tolerant.len(), 1, "{:?}", imp.tolerant);
            let t = &imp.tolerant[0];
            assert_eq!(t.index, 1, "{t:?}");
            assert!(t.vertices >= 1 && t.largest >= d, "{t:?}");
        }
    }
}

/// Reads `text`, a hand-made solid of torus or sphere faces closed by discs, by default and
/// strictly, and checks it against its closed forms: volume `volume`, area `area`, and the torus
/// or sphere faces whose bound the read turned against the file's flag, by their STEP ids.
fn read_turned(text: &[u8], volume: f64, area: f64, turned: &[u64]) {
    for options in [ReadOptions::default(), ReadOptions::default().strict(true)] {
        let imp = read_with(text, "band", &options).unwrap();
        assert!(imp.rejected.is_empty(), "{:?}", imp.rejected);
        assert!(imp.healed.is_empty() && imp.tolerant.is_empty());
        let b = &imp.brep;
        assert!(check::check(b).is_empty(), "{:?}", check::check(b));
        assert_eq!(b.solid_ids().count(), 1);
        assert!((measure::volume(b) - volume).abs() <= 1e-9 * volume);
        assert!((measure::area(b) - area).abs() <= 1e-9 * area);
        let ids: Vec<u64> = imp.turned_bounds.iter().map(|t| t.face_step_id).collect();
        assert_eq!(ids, turned);
        for t in &imp.turned_bounds {
            assert_eq!((t.index, t.bounds), (0, 1));
            assert!(matches!(
                b.face(t.face).surface,
                Surface::Torus(_) | Surface::Sphere(_)
            ));
        }
        round_trip(b);
    }
}

#[test]
fn a_torus_bound_flagged_against_its_edges_takes_the_sense_its_neighbours_agree_with() {
    // A quarter torus, R = 10 and r = 2 about Z from u = 0 to pi/2, closed by two discs. Its torus
    // face's one loop is listed the way that bounds the quarter band, and flagged .F.: on a closed
    // surface either sense of a loop bounds a region, the band or the rest of the torus, so the
    // flag alone turns the face against both discs across every edge. Read the way the discs
    // agree with, it is the band: volume (pi r^2)(R pi/2) = 20 pi^2, area a quarter of the
    // torus's 4 pi^2 R r plus two discs, 20 pi^2 + 8 pi; the read says it turned face #52's bound.
    let text = include_bytes!("fixtures/torus-bound-against-its-edges.step");
    let (volume, area) = (20.0 * PI * PI, 20.0 * PI * PI + 8.0 * PI);
    read_turned(text, volume, area, &[52]);
    // Flagged .T., or flagged .F. over the loop listed the other way, the file means the same band
    // as it says, and nothing is turned.
    let text = std::str::from_utf8(text).unwrap();
    let agreed = text.replace(
        "#46=FACE_OUTER_BOUND('',#45,.F.)",
        "#46=FACE_OUTER_BOUND('',#45,.T.)",
    );
    read_turned(agreed.as_bytes(), volume, area, &[]);
    let listed_back = text
        .replace(
            "#41=ORIENTED_EDGE('',*,*,#10,.T.)",
            "#41=ORIENTED_EDGE('',*,*,#10,.F.)",
        )
        .replace(
            "#42=ORIENTED_EDGE('',*,*,#22,.T.)",
            "#42=ORIENTED_EDGE('',*,*,#22,.F.)",
        )
        .replace(
            "#43=ORIENTED_EDGE('',*,*,#16,.F.)",
            "#43=ORIENTED_EDGE('',*,*,#16,.T.)",
        )
        .replace(
            "#44=ORIENTED_EDGE('',*,*,#22,.F.)",
            "#44=ORIENTED_EDGE('',*,*,#22,.T.)",
        )
        .replace("(#41,#42,#43,#44)", "(#44,#43,#42,#41)");
    assert_ne!(listed_back, text);
    read_turned(listed_back.as_bytes(), volume, area, &[]);
}

#[test]
fn torus_bounds_flagged_against_their_edges_take_their_sense_from_face_to_face() {
    // A half torus, R = 10 and r = 2 about Z from u = 0 to pi, in three torus faces of pi/3 each,
    // closed by two discs, each torus bound flagged .F. over a loop listed its band's way. The
    // middle face touches only the other two, so its sense follows theirs once the discs have set
    // them. Volume (pi r^2)(R pi) = 40 pi^2, area 40 pi^2 + 8 pi; the read turned all three, faces #80, #87 and #94.
    let text = include_bytes!("fixtures/half-torus-bounds-against-their-edges.step");
    let (volume, area) = (40.0 * PI * PI, 40.0 * PI * PI + 8.0 * PI);
    read_turned(text, volume, area, &[80, 87, 94]);
}

#[test]
fn a_sphere_cap_flagged_against_its_edge_takes_its_pole_from_its_neighbour() {
    // The upper half of a ball of radius r = 5 at the origin, closed by the disc at z = 0. The
    // sphere face's one loop, the equator, is listed anticlockwise about +Z, the way that bounds
    // the upper half, and flagged .F.: on a sphere either sense of a loop bounds a cap, so the flag
    // alone puts the face on the lower half, against the disc. Read the way the disc agrees with,
    // it is the upper half, over the pole at +Z: volume 2 pi r^3 / 3 = 250 pi / 3, area
    // 2 pi r^2 + pi r^2 = 75 pi; the read says it turned face #18's bound.
    let text = include_bytes!("fixtures/hemisphere-bound-against-its-edge.step");
    let (volume, area) = (250.0 * PI / 3.0, 75.0 * PI);
    read_turned(text, volume, area, &[18]);
    let text = std::str::from_utf8(text).unwrap();
    // Flagged .T., the file means the same half as it says, and nothing is turned.
    let agreed = text.replace(
        "#16=FACE_OUTER_BOUND('',#15,.F.)",
        "#16=FACE_OUTER_BOUND('',#15,.T.)",
    );
    read_turned(agreed.as_bytes(), volume, area, &[]);
    // On a sphere whose axis lies in the disc's plane, the equator runs through the sphere's own
    // poles and is read as a latitude about its circle's normal: the same half, and turned.
    let tilted = text.replace(
        "#17=SPHERICAL_SURFACE('',#6,5.0);",
        "#17=SPHERICAL_SURFACE('',#37,5.0);\n#37=AXIS2_PLACEMENT_3D('',#3,#5,#38);\n\
         #38=DIRECTION('',(0.0,1.0,0.0));",
    );
    assert_ne!(tilted, text);
    read_turned(tilted.as_bytes(), volume, area, &[18]);
    // With the disc turned to close the lower half, the flag agrees with it: the lower half, of
    // the same volume and area, read as written.
    let lower = text
        .replace(
            "#9=ORIENTED_EDGE('',*,*,#8,.F.)",
            "#9=ORIENTED_EDGE('',*,*,#8,.T.)",
        )
        .replace(
            "#13=ADVANCED_FACE('',(#11),#12,.F.)",
            "#13=ADVANCED_FACE('',(#11),#12,.T.)",
        );
    read_turned(lower.as_bytes(), volume, area, &[]);
}
