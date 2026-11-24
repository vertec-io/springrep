"""Reproducer for a band on a SURFACE_OF_REVOLUTION: a thick tube between z = 0 and z = H, its outer
wall the revolution of a cubic Bezier profile about the z axis (a barrel), its inner wall a
cylinder of radius RI. The barrel face is bounded by its two rim circles, each a closed edge, and
no seam, as some exporters write it; optionally two windows go through the wall, each bounded by
two planes through the axis and two planes across it, which leave two more loops on the barrel.

    python3 -I gen_revolution_band.py OUTDIR

The profile is z = H t, r(t) the cubic Bezier with control radii R, so every expected value is an
exact rational times pi (printed): the tube, pi (int r^2 dz - RI^2 H); each window, of angle DU,
(DU / 2) int_{z1}^{z2} (r^2 - RI^2) dz.
"""
import math
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stepgen import Step, bspline_curve, face_bounds, revolution  # noqa: E402

H, RI = 10, 6
Z1, Z2 = 3, 7  # the windows' floor and ceiling
WINDOWS = [(math.pi / 6, math.pi / 3), (7 * math.pi / 6, 4 * math.pi / 3)]  # each pi/6 wide


def bez(c, t):
    return sum(math.comb(3, k) * (1 - t) ** (3 - k) * t ** k * c[k] for k in range(4))


def r2_integral(R, z0, z1):
    """The integral of r(z)^2 dz over [z0, z1], exactly, with z = H t."""
    c0, c1, c2, c3 = (F(x) for x in R)
    m = [c0, 3 * (c1 - c0), 3 * (c2 - 2 * c1 + c0), c3 - 3 * c2 + 3 * c1 - c0]
    sq = [sum(m[i] * m[k - i] for i in range(4) if 0 <= k - i < 4) for k in range(7)]
    prim = lambda t: sum(sq[k] * t ** (k + 1) / (k + 1) for k in range(7))  # noqa: E731
    return H * (prim(F(z1, H)) - prim(F(z0, H)))


def expected(R, windows):
    k = r2_integral(R, 0, H) - RI * RI * H
    for _ in range(windows):
        k -= F(1, 12) * (r2_integral(R, Z1, Z2) - RI * RI * (Z2 - Z1))  # (pi/6)/2 = pi/12
    return k  # volume = k pi


def tube(R, name, revolved=True, windows=0, rims=("outer", "outer"), seam_at=(0.0, 0.0)):
    """rims: how the file flags the barrel's bottom and top rim loops. seam_at: the angle of each
    rim circle's one vertex."""
    s = Step()
    rad = lambda z: bez(R, z / H)  # noqa: E731
    at = lambda rho, u, z: (rho * math.cos(u), rho * math.sin(u), float(z))  # noqa: E731

    def circle_edge(rho, z, u):
        vtx = s.vertex(at(rho, u, z))
        c = s.circle((0, 0, z), (0, 0, 1), (math.cos(u), math.sin(u), 0), rho)
        return s.edge(vtx, vtx, c)

    if revolved:
        profile = bspline_curve(s, [(float(R[k]), 0.0, H * k / 3) for k in range(4)])
        outer_surface = revolution(s, profile, (0, 0, 0), (0, 0, 1))
    else:
        assert len(set(R)) == 1
        outer_surface = s.cylinder((0, 0, 0), (0, 0, 1), (1, 0, 0), R[0])
    inner_surface = s.cylinder((0, 0, 0), (0, 0, 1), (1, 0, 0), RI)
    o_bot, o_top = circle_edge(rad(0), 0, seam_at[0]), circle_edge(rad(H), H, seam_at[1])
    i_bot, i_top = circle_edge(RI, 0, 0.0), circle_edge(RI, H, 0.0)
    barrel = [(s.loop([(o_bot, True)]), rims[0]), (s.loop([(o_top, False)]), rims[1])]
    cyl = [(s.loop([(i_bot, False)]), "outer"), (s.loop([(i_top, True)]), "inner")]
    walls = []
    for a, b in WINDOWS[:windows]:
        V = {(rho, u, z): s.vertex(at(r, u, z)) for rho, r in (("i", RI), ("o", None)) for u in (a, b)
             for z in (Z1, Z2) if r is not None}
        for u in (a, b):
            for z in (Z1, Z2):
                V[("o", u, z)] = s.vertex(at(rad(z), u, z))

        def arc(rho, z):  # +u from a to b
            r = RI if rho == "i" else rad(z)
            c = s.circle((0, 0, z), (0, 0, 1), (1, 0, 0), r)
            return s.edge(V[(rho, a, z)], V[(rho, b, z)], c)

        def meridian(u):  # on the barrel, z1 to z2
            if revolved:
                c = bspline_curve(s, [at(float(R[k]), u, H * k / 3) for k in range(4)])
            else:
                c = s.line(at(rad(Z1), u, Z1), (0, 0, 1))
            return s.edge(V[("o", u, Z1)], V[("o", u, Z2)], c)

        def upright(u):  # on the cylinder, z1 to z2
            return s.edge(V[("i", u, Z1)], V[("i", u, Z2)], s.line(at(RI, u, Z1), (0, 0, 1)))

        def radial(u, z):  # inner to outer
            p, q = at(RI, u, z), at(rad(z), u, z)
            return s.edge(V[("i", u, z)], V[("o", u, z)], s.line(p, tuple(q[k] - p[k] for k in range(3))))

        oa1, oa2, ia1, ia2 = arc("o", Z1), arc("o", Z2), arc("i", Z1), arc("i", Z2)
        ma, mb, ua, ub = meridian(a), meridian(b), upright(a), upright(b)
        ra1, ra2, rb1, rb2 = radial(a, Z1), radial(a, Z2), radial(b, Z1), radial(b, Z2)
        barrel.append((s.loop([(oa1, False), (ma, True), (oa2, True), (mb, False)]), "inner"))
        cyl.append((s.loop([(ia1, True), (ub, True), (ia2, False), (ua, False)]), "inner"))
        walls += [
            # floor z = z1, normal +z
            face_bounds(s, [(s.loop([(oa1, True), (rb1, False), (ia1, False), (ra1, True)]), "outer")],
                        s.plane((0, 0, Z1), (0, 0, 1), (1, 0, 0))),
            # ceiling z = z2, normal -z
            face_bounds(s, [(s.loop([(oa2, False), (ra2, False), (ia2, True), (rb2, True)]), "outer")],
                        s.plane((0, 0, Z2), (0, 0, -1), (1, 0, 0))),
            # side at angle a, normal along +u
            face_bounds(s, [(s.loop([(ra1, False), (ua, True), (ra2, True), (ma, False)]), "outer")],
                        s.plane((0, 0, 0), (-math.sin(a), math.cos(a), 0), (math.cos(a), math.sin(a), 0))),
            # side at angle b, normal along -u
            face_bounds(s, [(s.loop([(rb1, True), (mb, True), (rb2, False), (ub, False)]), "outer")],
                        s.plane((0, 0, 0), (math.sin(b), -math.cos(b), 0), (math.cos(b), math.sin(b), 0))),
        ]
    faces = [
        face_bounds(s, barrel, outer_surface),
        face_bounds(s, cyl, inner_surface, same=False),
        # bottom annulus, normal -z; top, normal +z
        face_bounds(s, [(s.loop([(o_bot, False)]), "outer"), (s.loop([(i_bot, True)]), "inner")],
                    s.plane((0, 0, 0), (0, 0, -1), (1, 0, 0))),
        face_bounds(s, [(s.loop([(o_top, True)]), "outer"), (s.loop([(i_top, False)]), "inner")],
                    s.plane((0, 0, H), (0, 0, 1), (1, 0, 0))),
    ] + walls
    return s.text([s.solid(s.shell(faces))], name)


BARREL = (10, 13, 13, 10)
STRAIGHT = (10, 10, 10, 10)
cases = {
    # name: (radii, revolved, windows, rim flags, rim vertex angles)
    "revolution_band": (BARREL, True, 0, ("outer", "inner"), (0.0, 0.0)),
    "revolution_band_two_outer": (BARREL, True, 0, ("outer", "outer"), (0.0, 0.0)),
    "revolution_band_apart": (BARREL, True, 0, ("outer", "inner"), (0.0, math.pi / 2)),
    "revolution_band_windows": (BARREL, True, 2, ("outer", "outer"), (0.0, 0.0)),
    # controls: the same topology on a cylinder, which the reader seams already
    "cylinder_band_windows": (STRAIGHT, False, 2, ("outer", "outer"), (0.0, 0.0)),
}
out = sys.argv[1]
for name, (R, revolved, w, rims, seam_at) in cases.items():
    open(os.path.join(out, name + ".step"), "w").write(tube(R, name, revolved, w, rims, seam_at))
    k = expected(R, w)
    print(f"{name}: volume {k} pi = {float(k) * math.pi!r}")
