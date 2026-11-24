"""Repro: a cylinder (R 5, height 10) whose side face, bounded by its two circles with no seam edge,
also carries a VERTEX_LOOP (a point bound, as exporters write where another body touches the face)
at a point of the surface's seam, u = 0. Control: the same point bound at u = pi.
Expected volume 250 pi in both. usage: python3 -I r_vertex_loop_on_seam.py OUTDIR"""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stepgen import Step

R, H = 5.0, 10.0


def cyl(point, name):
    s = Step()
    A, B = (R, 0, 0), (R, 0, H)
    vA, vB = s.vertex(A), s.vertex(B)
    cb = s.edge(vA, vA, s.circle((0, 0, 0), (0, 0, 1), (1, 0, 0), R))
    ct = s.edge(vB, vB, s.circle((0, 0, H), (0, 0, 1), (1, 0, 0), R))
    bottom = s.face([(s.loop([(cb, False)]), True)], s.plane((0, 0, 0), (0, 0, -1), (1, 0, 0)))
    top = s.face([(s.loop([(ct, True)]), True)], s.plane((0, 0, H), (0, 0, 1), (1, 0, 0)))
    vl = s.e("VERTEX_LOOP", "''", s.vertex(point))
    side = s.face([(s.loop([(cb, True)]), True), (s.loop([(ct, False)]), False), (vl, False)],
                  s.cylinder((0, 0, 0), (0, 0, 1), (1, 0, 0), R))
    return s.text([s.solid(s.shell([bottom, top, side]))], name)


out = sys.argv[1]
open(os.path.join(out, "vertex_loop_on_seam.step"), "w").write(cyl((R, 0, H / 2), "vertex_loop_on_seam"))
open(os.path.join(out, "vertex_loop_off_seam.step"), "w").write(cyl((-R, 0, H / 2), "vertex_loop_off_seam"))
print("expected volume", math.pi * R * R * H)
