"""Repros built from blocks. usage: python3 -I r_blocks.py OUTDIR
- zero_length_edge: a 10 mm cube whose corner is two VERTEX_POINTs at one point, joined by a
  zero-length LINE edge used by two faces (as an exporter wrote it). Expected: volume 1000.
- void_outside_touching: BREP_WITH_VOIDS whose 'void' is a second cube outside the outer one,
  sharing a face with it (a compound written as voids). Not valid ISO 10303-42 (a void lies inside
  its outer shell); read as two solids of 1000 each, or refused with that reason.
- void_outside_apart: the same with the second cube 5 mm away.
- void_inside: control, a proper void (2 mm cube inside a 10 mm cube): volume 992.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stepgen import Step, block

out = sys.argv[1]
s = Step(); f = block(s, 0, 0, 0, 10, 10, 10, zero_edge=True)
open(os.path.join(out, "zero_length_edge.step"), "w").write(s.text([s.solid(s.shell(f))], "zero_length_edge"))
for name, x in (("void_outside_touching", 10), ("void_outside_apart", 15)):
    s = Step(); a = s.shell(block(s, 0, 0, 0, 10, 10, 10)); b = s.shell(block(s, x, 0, 0, 10, 10, 10))
    open(os.path.join(out, name + ".step"), "w").write(s.text([s.brep_with_voids(a, [b])], name))
s = Step(); a = s.shell(block(s, 0, 0, 0, 10, 10, 10)); b = s.shell(block(s, 4, 4, 4, 2, 2, 2))
open(os.path.join(out, "void_inside.step"), "w").write(s.text([s.brep_with_voids(a, [b])], "void_inside"))
