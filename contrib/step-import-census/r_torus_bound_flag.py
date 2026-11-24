"""Repro: a quarter torus (R 10, r 2, axis Z, u from 0 to pi/2) closed by two discs. Its torus face
has one loop: the tube circle at u = 0, the seam arc, the tube circle at u = pi/2, the seam back.
On a closed surface both senses of that loop bound a region (the band, or the rest of the torus),
so the FACE_BOUND flag picks it; some exporters write .F. while listing the loop the band's way.
- torus_bound_ok:          flag .T., loop the band's way (control)
- torus_bound_reversed_ok: flag .F., loop listed the other way (same meaning, valid; control)
- torus_bound_flag_contradicts: flag .F., loop the band's way (the exporter case): the faces then
  disagree across every edge; reading the face as the band gives the solid back.
Expected volume (pi r^2)(R pi/2) = 20 pi^2. usage: python3 -I r_torus_bound_flag.py OUTDIR"""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stepgen import Step, fmt

R, r = 10.0, 2.0


def quarter(name, flag, reverse_loop, c0=True, c1=False):
    s = Step()
    P0, P1 = (R + r, 0, 0), (0, R + r, 0)
    v0, v1 = s.vertex(P0), s.vertex(P1)
    C0 = s.edge(v0, v0, s.circle((R, 0, 0), (0, 1, 0), (1, 0, 0), r))
    C1 = s.edge(v1, v1, s.circle((0, R, 0), (-1, 0, 0), (0, 1, 0), r))
    seam = s.edge(v0, v1, s.circle((0, 0, 0), (0, 0, 1), (1, 0, 0), R + r))
    d0 = s.face([(s.loop([(C0, not c0)]), True)], s.plane((R, 0, 0), (0, -1, 0), (1, 0, 0)))
    d1 = s.face([(s.loop([(C1, not c1)]), True)], s.plane((0, R, 0), (-1, 0, 0), (0, 1, 0)))
    uses = [(C0, c0), (seam, True), (C1, c1), (seam, False)]
    if reverse_loop:
        uses = [(e, not o) for e, o in reversed(uses)]
    oes = [s.e("ORIENTED_EDGE", "''", "*", "*", e, o) for e, o in uses]
    lp = s.e("EDGE_LOOP", "''", tuple(oes))
    bd = s.e("FACE_OUTER_BOUND", "''", lp, flag)
    tor = s.e("TOROIDAL_SURFACE", "''", s.ax((0, 0, 0), (0, 0, 1), (1, 0, 0)), R, r)
    tf = s.e("ADVANCED_FACE", "''", (bd,), tor, True)
    return s.text([s.solid(s.shell([d0, d1, tf]))], name)


out = sys.argv[1]
if len(sys.argv) > 2:  # probe the senses
    for c0 in (True, False):
        for c1 in (True, False):
            open(os.path.join(out, f"probe_{int(c0)}{int(c1)}.step"), "w").write(quarter("p", True, False, c0, c1))
    sys.exit()
open(os.path.join(out, "torus_bound_ok.step"), "w").write(quarter("torus_bound_ok", True, False))
open(os.path.join(out, "torus_bound_reversed_ok.step"), "w").write(quarter("torus_bound_reversed_ok", False, True))
open(os.path.join(out, "torus_bound_flag_contradicts.step"), "w").write(quarter("torus_bound_flag_contradicts", False, False))
print("expected volume", math.pi * r * r * R * math.pi / 2)
