"""Repro: a ball octant whose SPHERICAL_SURFACE axis puts the sphere's pole in the middle of a
boundary arc (not at a vertex). Control files: the axis through a vertex, and through the face's
interior. usage: python3 -I r_sphere_pole_on_edge.py OUTDIR"""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stepgen import Step

R = 10.0


def octant(axis, ref, name):
    s = Step()
    O, X, Y, Z = (0, 0, 0), (R, 0, 0), (0, R, 0), (0, 0, R)
    vO, vX, vY, vZ = (s.vertex(p) for p in (O, X, Y, Z))
    eOX = s.edge(vO, vX, s.line(O, (1, 0, 0)))
    eOY = s.edge(vO, vY, s.line(O, (0, 1, 0)))
    eOZ = s.edge(vO, vZ, s.line(O, (0, 0, 1)))
    aXY = s.edge(vX, vY, s.circle(O, (0, 0, 1), (1, 0, 0), R))
    aYZ = s.edge(vY, vZ, s.circle(O, (1, 0, 0), (0, 1, 0), R))
    aZX = s.edge(vZ, vX, s.circle(O, (0, 1, 0), (0, 0, 1), R))
    f = [
        s.face([(s.loop([(eOY, True), (aXY, False), (eOX, False)]), True)], s.plane(O, (0, 0, -1), (1, 0, 0))),
        s.face([(s.loop([(eOZ, True), (aYZ, False), (eOY, False)]), True)], s.plane(O, (-1, 0, 0), (0, 1, 0))),
        s.face([(s.loop([(eOX, True), (aZX, False), (eOZ, False)]), True)], s.plane(O, (0, -1, 0), (0, 0, 1))),
        s.face([(s.loop([(aXY, True), (aYZ, True), (aZX, True)]), True)], s.sphere(O, axis, ref, R)),
    ]
    return s.text([s.solid(s.shell(f))], name)


out = sys.argv[1]
h = 1 / math.sqrt(2)
cases = {
    # expected volume of every case: pi R^3 / 6
    "sphere_pole_on_edge": ((h, 0, h), (h, 0, -h)),     # pole at (R/sqrt2, 0, R/sqrt2), mid-arc ZX
    "sphere_pole_at_vertex": ((0, 0, 1), (1, 0, 0)),  # control: pole at vertex Z
    "sphere_pole_inside": ((1, 1, 1), (1, -1, 0)),    # control: pole inside the face
}
for n, (ax, ref) in cases.items():
    open(os.path.join(out, n + ".step"), "w").write(octant(ax, ref, n))
print("expected volume", math.pi * R**3 / 6, "area", 3 * math.pi * R * R / 4 + math.pi * R * R / 2)
