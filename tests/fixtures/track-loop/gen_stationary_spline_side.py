"""Reproducer for a spline face whose parameterisation stands still along one of its boundary
lines: a prism of height H over a plane region bounded by the segments B-O-A and a cubic Bezier
from A to B, its side on the Bezier a B-spline surface, cubic along u and linear along v. Where two
end control points of the cubic coincide, the surface's u derivative vanishes along that whole end
line u = 0 or u = 1, though the line itself is not a point: the parameterisation is stationary
there, not singular.

    python3 -I gen_stationary_spline_side.py OUTDIR

Expected values, from Green's theorem on the polynomial boundary (exact rationals, printed):
volume H * area, area the region's.
"""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stepgen import Step, bspline_curve, bspline_surface  # noqa: E402

H = 4
O, A = (0, 0), (10, 0)


def bezier_area_term(P):
    """The integral of x dy along the cubic Bezier with control points P, exactly."""
    # x(t), y(t) in the Bernstein basis; integrate x(t) y'(t) over [0, 1] by expanding to monomials.
    def mono(c):  # Bernstein coefficients -> monomial coefficients, cubic
        c0, c1, c2, c3 = (F(v) for v in c)
        return [c0, 3 * (c1 - c0), 3 * (c2 - 2 * c1 + c0), c3 - 3 * c2 + 3 * c1 - c0]
    x = mono([p[0] for p in P])
    y = mono([p[1] for p in P])
    dy = [k * y[k] for k in range(1, 4)]  # coefficients of t^0, t^1, t^2
    return sum(x[i] * dy[j] / (i + j + 1) for i in range(4) for j in range(3))


def region_area(P):
    # O -> A along y = 0 contributes nothing to the integral of x dy; A -> B along the Bezier;
    # B -> O along the segment.
    B = P[-1]
    back = -F(B[0]) * F(B[1]) / 2  # integral of x dy from B to O on the segment
    return bezier_area_term(P) + back


def prism(P, name, flip_u=False, flip_v=False):
    """flip_u, flip_v: the spline side's parameters run the other way (u from B to A, v down over
    [-1, 0]), its natural normal turned inward by each, the face's sense written to match."""
    s = Step()
    B = P[-1]
    lift = lambda p, z: (float(p[0]), float(p[1]), float(z))  # noqa: E731
    v = {(k, z): s.vertex(lift(p, z)) for k, p in (("O", O), ("A", A), ("B", B)) for z in (0, H)}
    pts = {"O": O, "A": A, "B": B}

    def line(a, b, za, zb):
        pa, pb = lift(pts[a], za), lift(pts[b], zb)
        d = tuple(pb[i] - pa[i] for i in range(3))
        return s.edge(v[(a, za)], v[(b, zb)], s.line(pa, d))

    bot_OA, top_OA = line("O", "A", 0, 0), line("O", "A", H, H)
    bot_OB, top_OB = line("O", "B", 0, 0), line("O", "B", H, H)
    up_O, up_A, up_B = line("O", "O", 0, H), line("A", "A", 0, H), line("B", "B", 0, H)
    # The Bezier side's bottom and top curves, A to B, which are its surface's v = 0 and v = 1 lines.
    bot_AB = s.edge(v[("A", 0)], v[("B", 0)], bspline_curve(s, [lift(p, 0) for p in P]))
    top_AB = s.edge(v[("A", H)], v[("B", H)], bspline_curve(s, [lift(p, H) for p in P]))
    rows = [[lift(p, 0), lift(p, H)] for p in P]
    if flip_u:
        rows = rows[::-1]
    if flip_v:
        rows = [r[::-1] for r in rows]
    side = bspline_surface(s, rows, (-1.0, 0.0) if flip_v else (0.0, 1.0))
    n_ob = (-float(B[1]), float(B[0]), 0.0)  # outward normal of the side through O and B
    faces = [
        # bottom, seen from below: O -> B -> A
        s.face([(s.loop([(bot_OB, True), (bot_AB, False), (bot_OA, False)]), True)],
               s.plane((0, 0, 0), (0, 0, -1), (1, 0, 0))),
        # top, seen from above: O -> A -> B
        s.face([(s.loop([(top_OA, True), (top_AB, True), (top_OB, False)]), True)],
               s.plane((0, 0, H), (0, 0, 1), (1, 0, 0))),
        # side y = 0, seen from -y: O -> A up A, back along the top, down O
        s.face([(s.loop([(bot_OA, True), (up_A, True), (top_OA, False), (up_O, False)]), True)],
               s.plane((0, 0, 0), (0, -1, 0), (1, 0, 0))),
        # side through O and B, seen from outside: B -> O at the bottom, up O, O -> B on top, down B
        s.face([(s.loop([(bot_OB, False), (up_O, True), (top_OB, True), (up_B, False)]), True)],
               s.plane((0, 0, 0), n_ob, (float(B[0]), float(B[1]), 0.0))),
        # the spline side, its natural normal (u along the Bezier, v up) outward: A -> B along the
        # bottom, up B, B -> A along the top, down A
        s.face([(s.loop([(bot_AB, True), (up_B, True), (top_AB, False), (up_A, False)]), True)], side,
               flip_u == flip_v),
    ]
    return s.text([s.solid(s.shell(faces))], name)


cases = {
    # name: control points, and whether the side's u and v run the other way
    # the last two control points coincide: the surface stands still along u = 1
    "spline_side_stationary_end": [A, (10, 5), (5, 10), (5, 10)],
    # the same, the side's parameters each way round, as one exporter writes it
    "spline_side_stationary_end_u": ([A, (10, 5), (5, 10), (5, 10)], True, False),
    "spline_side_stationary_end_v": ([A, (10, 5), (5, 10), (5, 10)], False, True),
    "spline_side_stationary_end_uv": ([A, (10, 5), (5, 10), (5, 10)], True, True),
    # the first two: along u = 0
    "spline_side_stationary_start": [A, A, (10, 8), (5, 10)],
    # control: four distinct control points
    "spline_side_regular": [A, (10, 5), (8, 9), (5, 10)],
}
out = sys.argv[1]
for name, case in cases.items():
    P, fu, fv = case if isinstance(case, tuple) else (case, False, False)
    open(os.path.join(out, name + ".step"), "w").write(prism(P, name, fu, fv))
    a = region_area(P)
    print(f"{name}: volume {H * a} = {float(H * a)!r}")
