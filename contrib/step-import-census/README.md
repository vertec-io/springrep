# STEP import census: hand-written reproducers

Reproducers for the causes listed in the springrep issue "STEP import on real-world files: a census
of 19 failures in 13 causes". Each `.step` file in `fixtures/` was written entity by entity from the
ISO 10303-42 definitions by the Python scripts here (`stepgen.py` plus one `r_*.py` per cause); no
exporter or CAD program was involved. Each failing file comes with a control that reads correctly
and has a closed-form volume.

These files and scripts are our own work and are offered under this repository's licences
(MIT OR Apache-2.0), so they may be copied into springrep's `fixtures/` as tests.

| Cause | Failing file | Control(s) |
|---|---|---|
| 1 cone face past the apex | `cone_far_nappe.step`, `cone_far_nappe_F.step` | `cone_near_nappe.step` |
| 2 zero-length edge | `zero_length_edge.step` | — (closed-form cube, V = 1000) |
| 3 sphere boundary through a pole | `sphere_pole_on_edge.step` | `sphere_pole_at_vertex.step`, `sphere_pole_inside.step` |
| 4 "voids" outside the outer shell | `void_outside_touching.step` | `void_inside.step`, `void_outside_apart.step` |
| 5 torus bound flag contradicts its edges | `torus_bound_flag_contradicts.step` | `torus_bound_ok.step`, `torus_bound_reversed_ok.step` |
| 6 point bound on a periodic seam | `vertex_loop_on_seam.step` | `vertex_loop_off_seam.step` |

Regenerate with `python3 r_<cause>.py` (Python 3, standard library only).
