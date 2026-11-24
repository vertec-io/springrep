"""Repro: a frustum (z from -10 to -5, radii 10 and 5, apex at the origin, 45 degrees) whose
CONICAL_SURFACE is placed at z = +5 with radius 5 about +Z, so the face lies where
radius + v tan(angle) < 0: on the cone's other nappe, past the apex (ISO 10303-42 gives the cone
for every v). Control: the same frustum on a cone placed at z = -5, radius 5, about -Z.
Expected volume 875 pi / 3. usage: python3 -I r_cone_far_nappe.py OUTDIR"""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stepgen import Step


def frustum(name, cone_origin, cone_axis, same):
    s = Step()
    A, B = (10, 0, -10), (5, 0, -5)
    vA, vB = s.vertex(A), s.vertex(B)
    cb = s.edge(vA, vA, s.circle((0, 0, -10), (0, 0, 1), (1, 0, 0), 10))
    ct = s.edge(vB, vB, s.circle((0, 0, -5), (0, 0, 1), (1, 0, 0), 5))
    seam = s.edge(vA, vB, s.line(A, (-5, 0, 5)))
    bottom = s.face([(s.loop([(cb, False)]), True)], s.plane((0, 0, -10), (0, 0, -1), (1, 0, 0)))
    top = s.face([(s.loop([(ct, True)]), True)], s.plane((0, 0, -5), (0, 0, 1), (1, 0, 0)))
    cone = s.e("CONICAL_SURFACE", "''", s.ax(cone_origin, cone_axis, (1, 0, 0)), 5.0, math.pi / 4)
    side = s.face([(s.loop([(cb, True), (seam, True), (ct, False), (seam, False)]), True)], cone, same)
    return s.text([s.solid(s.shell([bottom, top, side]))], name)


out = sys.argv[1]
open(os.path.join(out, "cone_far_nappe.step"), "w").write(frustum("cone_far_nappe", (0, 0, 5), (0, 0, 1), True))
open(os.path.join(out, "cone_far_nappe_F.step"), "w").write(frustum("cone_far_nappe_F", (0, 0, 5), (0, 0, 1), False))
open(os.path.join(out, "cone_near_nappe.step"), "w").write(frustum("cone_near_nappe", (0, 0, -5), (0, 0, -1), True))
print("expected volume", 875 * math.pi / 3)
