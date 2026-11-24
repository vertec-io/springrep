"""Tiny hand-written AP214 writer for the census's minimal reproducers.

Written from the ISO 10303-42 / -21 entity definitions only (no CAD program involved): each repro
script builds a solid entity by entity, so the file holds exactly the construct under test.
"""
import math


def fmt(x):
    if isinstance(x, str):
        return x
    if isinstance(x, bool):
        return ".T." if x else ".F."
    if isinstance(x, float) or isinstance(x, int):
        s = repr(float(x))
        if "e" in s or "E" in s:
            m, e = s.lower().split("e")
            if "." not in m:
                m += "."
            return f"{m}E{int(e)}"
        return s if "." in s else s + "."
    if isinstance(x, (list, tuple)):
        return "(" + ",".join(fmt(v) for v in x) + ")"
    raise TypeError(x)


class Step:
    def __init__(self, uncertainty=1e-7):
        self.lines = []
        self.n = 0
        self.unc = uncertainty

    def add(self, text):
        self.n += 1
        self.lines.append(f"#{self.n}={text};")
        return f"#{self.n}"

    def e(self, name, *args):
        return self.add(f"{name}(" + ",".join(fmt(a) for a in args) + ")")

    # geometry
    def pt(self, p):
        return self.e("CARTESIAN_POINT", "''", tuple(float(c) for c in p))

    def dir(self, d):
        n = math.sqrt(sum(c * c for c in d))
        return self.e("DIRECTION", "''", tuple(c / n for c in d))

    def ax(self, o, z=(0, 0, 1), x=(1, 0, 0)):
        return self.e("AXIS2_PLACEMENT_3D", "''", self.pt(o), self.dir(z), self.dir(x))

    def line(self, p, d):
        n = math.sqrt(sum(c * c for c in d))
        return self.e("LINE", "''", self.pt(p), self.e("VECTOR", "''", self.dir(d), n))

    def circle(self, o, z, x, r):
        return self.e("CIRCLE", "''", self.ax(o, z, x), float(r))

    def plane(self, o, z, x):
        return self.e("PLANE", "''", self.ax(o, z, x))

    def sphere(self, o, z, x, r):
        return self.e("SPHERICAL_SURFACE", "''", self.ax(o, z, x), float(r))

    def cylinder(self, o, z, x, r):
        return self.e("CYLINDRICAL_SURFACE", "''", self.ax(o, z, x), float(r))

    # topology
    def vertex(self, p):
        return self.e("VERTEX_POINT", "''", self.pt(p))

    def edge(self, v1, v2, curve, same=True):
        return self.e("EDGE_CURVE", "''", v1, v2, curve, same)

    def loop(self, oriented):
        """oriented: list of (edge, orientation)."""
        oes = [self.e("ORIENTED_EDGE", "''", "*", "*", e, o) for e, o in oriented]
        return self.e("EDGE_LOOP", "''", tuple(oes))

    def face(self, loops, surface, same=True):
        """loops: [(loop, outer?)]; the first should be the outer one."""
        bs = [self.e("FACE_OUTER_BOUND" if outer else "FACE_BOUND", "''", l, True) for l, outer in loops]
        return self.e("ADVANCED_FACE", "''", tuple(bs), surface, same)

    def shell(self, faces):
        return self.e("CLOSED_SHELL", "''", tuple(faces))

    def solid(self, shell):
        return self.e("MANIFOLD_SOLID_BREP", "''", shell)

    def brep_with_voids(self, outer, voids):
        ov = [self.e("ORIENTED_CLOSED_SHELL", "''", "*", v, False) for v in voids]
        return self.e("BREP_WITH_VOIDS", "''", outer, tuple(ov))

    def text(self, items, name="repro"):
        """items: the solids for one part's ADVANCED_BREP_SHAPE_REPRESENTATION, in millimetres."""
        head = len(self.lines)
        mm = self.add("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) )")
        rad = self.add("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )")
        sr = self.add("( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() )")
        unc = self.e("UNCERTAINTY_MEASURE_WITH_UNIT", f"LENGTH_MEASURE({fmt(self.unc)})", mm, "'distance_accuracy_value'", "''")
        ctx = self.add(f"( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(({unc})) "
                       f"GLOBAL_UNIT_ASSIGNED_CONTEXT(({mm},{rad},{sr})) REPRESENTATION_CONTEXT('',''))")
        ac = self.e("APPLICATION_CONTEXT", "'automotive design'")
        self.e("APPLICATION_PROTOCOL_DEFINITION", "'international standard'", "'automotive_design'", "2000", ac)
        pc = self.e("PRODUCT_CONTEXT", "''", ac, "'mechanical'")
        pr = self.e("PRODUCT", f"'{name}'", f"'{name}'", "''", (pc,))
        pdf = self.e("PRODUCT_DEFINITION_FORMATION", "''", "''", pr)
        pdc = self.e("PRODUCT_DEFINITION_CONTEXT", "'part definition'", ac, "'design'")
        pd = self.e("PRODUCT_DEFINITION", "'design'", "''", pdf, pdc)
        pds = self.e("PRODUCT_DEFINITION_SHAPE", "''", "''", pd)
        origin = self.ax((0, 0, 0))
        rep = self.e("ADVANCED_BREP_SHAPE_REPRESENTATION", f"'{name}'", tuple([origin] + list(items)), ctx)
        self.e("SHAPE_DEFINITION_REPRESENTATION", pds, rep)
        return ("ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('springrep import census reproducer'),'2;1');\n"
                f"FILE_NAME('{name}.step','2026-11-24T00:00:00',(''),(''),'','','');\n"
                "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n"
                + "\n".join(self.lines) + "\nENDSEC;\nEND-ISO-10303-21;\n")


def block(s, x0, y0, z0, dx, dy, dz, zero_edge=False):
    """An axis-aligned block's six faces, outward. With zero_edge, its corner (x0, y0, z0) is
    written as two VERTEX_POINTs at the same point joined by a zero-length LINE edge, which the
    front (y = y0) and left (x = x0) faces both use, one each way."""
    P = [(x0 + dx * (i & 1), y0 + dy * ((i >> 1) & 1), z0 + dz * ((i >> 2) & 1)) for i in range(8)]
    V = [s.vertex(p) for p in P]
    if zero_edge:
        P.append(P[0])
        V.append(s.vertex(P[0]))  # vertex 8: corner 0 again
    E = {}

    def ed(a, b):
        if (a, b) in E:
            return E[(a, b)], True
        if (b, a) in E:
            return E[(b, a)], False
        d = tuple(P[b][i] - P[a][i] for i in range(3))
        if max(abs(c) for c in d) == 0:
            d = (0, 0, 1)
        E[(a, b)] = s.edge(V[a], V[b], s.line(P[a], d))
        return E[(a, b)], True

    faces = [
        ((0, 2, 3, 1), (0, 0, -1), (1, 0, 0)),  # bottom
        ((4, 5, 7, 6), (0, 0, 1), (1, 0, 0)),  # top
        ((0, 1, 5, 4), (0, -1, 0), (1, 0, 0)),  # front y = y0
        ((2, 6, 7, 3), (0, 1, 0), (1, 0, 0)),  # back
        ((0, 4, 6, 2), (-1, 0, 0), (0, 1, 0)),  # left x = x0
        ((1, 3, 7, 5), (1, 0, 0), (0, 1, 0)),  # right
    ]
    if zero_edge:
        faces[2] = ((0, 1, 5, 4, 8), (0, -1, 0), (1, 0, 0))
        faces[4] = ((0, 8, 4, 6, 2), (-1, 0, 0), (0, 1, 0))
    out = []
    for cyc, nz, nx in faces:
        oes = [ed(cyc[i], cyc[(i + 1) % len(cyc)]) for i in range(len(cyc))]
        out.append(s.face([(s.loop(oes), True)], s.plane(P[cyc[0]], nz, nx), True))
    return out
