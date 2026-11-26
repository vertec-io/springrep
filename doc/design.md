# springrep design

How the kernel represents solids, which tolerances it holds them to, how face names follow history,
and which module owns what. Every later stage builds on these rules; change one here first.

## Representation

A solid is exact geometry under an arena of topology, all in one [`Brep`](../src/topo/mod.rs):

```text
Solid ─ Shell ─ Face (surface, reversed, loops, names) ─ Loop ─ Coedge (pcurve)
      ─ Edge (curve) ─ Vertex
```

- **Ids, not pointers.** Each entity lives in a `Vec` and is named by a typed id (`FaceId`,
  `EdgeId`, …). Builders append; nothing is removed in place. An operation that changes topology
  builds a new arena (or appends one with `Brep::append`) rather than patching the old one.
- **Solid.** Its first shell is the outer boundary; any others bound voids.
- **Face.** A bounded piece of one surface. Its outward normal is the surface's natural normal
  `∂S/∂u × ∂S/∂v`, or the opposite when `reversed`. The first loop is the outer one.
- **Loop.** A closed chain of coedges with the face's material on the left, seen from outside the
  solid. In the parameter plane the outer loop therefore runs counter-clockwise when the face is not
  reversed and clockwise when it is; holes run the other way. The checker enforces this.
- **Coedge.** One face's use of an edge, along the edge's parameter or against it (`reversed`). It
  carries a **pcurve**: where it runs in the face's parameter plane, as a function of the edge
  parameter. `PCurve::Affine` (`uv = origin + t·dir`) covers every edge a sweep makes along or
  across a surface's parameter lines; `PCurve::Projected` projects the 3-D point and continues
  periodic parameters from a seed, which is exact on a plane, halving any step between samples that
  turns a periodic parameter by over a quarter period, since beside a pole a short step can turn it
  by nearly half a period either way; at a pole, where every `u` is the same point, `u` is the angle
  the curve's tangent leaves the pole at, so every reader finds the same one however it samples the
  curve. At an end of its edge a projected pcurve takes the parameters the coedge beside it in the
  loop gives the shared vertex, where that one is affine or a spline and reaches the vertex: where
  the map nearly drops rank, as along the row of a swept sheet where the motion runs along its
  curve, projection lands anywhere along the direction the map barely moves in, and the loop would
  not close in parameters. Booleans give every intersection edge a projected pcurve.
  `PCurve::Spline` reads `uv` off a cubic B-spline in `t` (its control points in the plane `z = 0`):
  for an edge whose exact parameters are known at points along it, on a surface whose map nearly
  drops rank there, where projection lands on the wrong foot. It holds over its own edge's range
  alone, so two edges whose coedges carry one are not merged into one edge; an affine change of the
  parameter plane (a period's shift, a mirror, a scale) maps its control points.
- **Edge.** A curve and a range `t0 < t1`, from `start` at `curve(t0)` to `end` at `curve(t1)`. A
  closed curve's edge starts and ends at the same vertex. A **degenerate edge** has no curve: it is
  a pole (a cone's apex, a sphere's pole) that its face's loop spans in the parameter plane. It is
  used by exactly one coedge.
- **Seams.** A face that wraps all the way round a periodic parameter carries a seam edge, used
  twice by that face with pcurves one period apart. So every face's loops close in its parameter
  plane, which tessellation, measurement and point classification rely on, and every closed critical
  curve of a coordinate (a torus's top circle) meets a boundary edge.

### Geometry

Curves ([`geom/curve.rs`](../src/geom/curve.rs)): line (arc-length parameter), circle and ellipse
(angle from the frame's `x`, counter-clockwise about its `z`), non-rational B-spline up to degree 3
([`geom/bspline.rs`](../src/geom/bspline.rs), open or periodic, interpolated at chord-length
parameters: not-a-knot cubic when open, C² when closed; a profile's spline has a knot at every point
and at each end the tangent of the cubic through the four nearest points, the ends other modellers
give a spline through points; the collocation systems are banded, so a fit through thousands of
points costs linear time and memory), and intersection curves
([`geom/intersection.rs`](../src/geom/intersection.rs)): a cubic through points refined onto both
surfaces, with points added until the spline stays within half the linear tolerance of both; a span
beside a split one is split too while it is over three times as long as that one's halves, since
uneven spacing makes the interpolant swing. Where the surfaces cross at a tiny angle θ, a point on
both to within ε is placed across the curve only to within ε/θ (two tori 1e-3 apart cross near their
tops at about 1e-7 rad, so points scatter across a band about 1e-3 wide), and the fit drops any
point whose chords to its neighbours point against each other: a traced branch never reverses
between neighbouring points, and through such scatter the interpolant would fold back on itself.
Rational B-splines of any degree up to 15 ([`geom/nurbs.rs`](../src/geom/nurbs.rs)), as STEP files
carry them, are curves and surfaces too; one whose ends meet within the linear tolerance is periodic
in that parameter, so an edge or loop may run across its join.

Surfaces ([`geom/surface.rs`](../src/geom/surface.rs)): plane, cylinder, cone, sphere and torus in
closed form, and extrusion, revolution and ruled surfaces of any curve. On every surface swept about
an axis, `u` is the angle about it from the frame's `x`; `v` is height (cylinder), arc length along
a generator (cone), latitude (sphere), tube angle (torus), or the profile curve's own parameter
(revolution). Every curve and surface has evaluation, first derivatives, projection (closest point)
and its periods; curves also have second derivatives, breakpoints and coordinate extremes.
Projection onto a B-spline searches only the knot spans whose control points' bounding ball could
hold a nearer point, so its cost does not grow with the span count; an extrusion of one projects
along its direction onto the curve. A B-spline's coordinate extremes are each span's quadratic
roots. A rational curve's are roots found by sampling each span; those of a whole span along x, y
and z are kept with the curve from the first box that reaches the span, so a box over a range works
out only the stretches that end inside a span, and is the box it would be with nothing kept. A moved
curve keeps nothing of its source's. Projection onto a rational surface runs damped Newton steps
from the nearest few of a grid of seed points, spaced about as closely along a long narrow patch (a
twisted loft's side between close corners) as across it, so the nearest seed lies in the nearest
point's basin. The four nearest seeds are the scan's: least in distance, the earlier of two equally
far, one for each point where seeds lie on each other (a seam, a pole). Boxes round rows and runs of
the seeds, built once a surface has been scanned sixteen times, leave out the seeds no nearer than
the fourth so far, which gives the same four in any order wherever none of the seeds that near lies
on another; where one does, the seeds are scanned in turn. Newton takes the full second-order step
where the squared distance's Hessian is positive definite, so it converges quadratically off the
surface too. A parameter at an end of its range, with the nearer point beyond it, is held there, and
the search goes on along that boundary alone: a step solved in both parameters and cut back to the
range stalls short of the boundary's nearest point. At a pole, where a whole row of the surface is
one point and one parameter moves nothing, the line out of the pole that heads for the point most
directly is found first, and the search starts down it. A row is one point where its control
points, blended with their weights, all lie within `LINEAR` of one (the convex hull holds the row);
a parameter that moves nothing along a row that is a curve, as where two end control points
coincide, is a line along which the surface stands still, not a pole. Beside such a boundary line
the squared distance grows only with the fourth power of the step across it and Newton stops short,
so the line is also tried, the search clamped to it, and the nearer foot kept; that foot is a foot
where the squared distance does not bend down into the patch across the line, the first-order
condition saying nothing there. A search that stops where no fraction of its step comes nearer
has found a point and not a foot, and is not taken for one. Projection from a hint
keeps the foot Newton finds when it lies within a twentieth of the seed spacing of the surface: the
global search already assumes no two sheets come that close. A foot farther off than that, a search
that stopped short of a foot, and a hint on a pole all go to the global search, and the nearer of
its point and the hint's is kept. What a rational surface works out of itself on first use (the
seeds, the patches, the span polynomials) is shared by its copies, whichever first asks, and a
surface with other control points starts afresh; a surface moved rigidly takes its source's seeds
and polynomials moved with it. A rational surface also carries small patches, each knot span cut to
about a sixteenth of the surface's longer side, whose control points after knot insertion bound it;
rays and curves are sampled against the surface only where they pass through a patch's box, and a
face's box only searches the patches that poke out of it.

Sweeps choose the analytic surface whenever one exists: an extruded line is a plane and an extruded
circle a cylinder; a revolved line is a cylinder, plane or cone, and a revolved circle a sphere
(centre on the axis) or torus. An arc whose circle's centre lies across the axis from it turns
into the inner sheet of a spindle torus, where the torus's closed-form projection, which takes the
half-plane the point lies in, finds the outer sheet; it is the surface of revolution of the arc's
circle instead, which projects onto the arc. The torus primitive (`primitive::torus`) takes only a
tube clear of its axis, and refuses one that reaches it. The tori whose tube reaches the axis have
calls of their own: `primitive::horn_torus`, the tube as wide as its circle, turned whole and
pinched at the centre; and `primitive::spindle_torus`, a wider tube counted once, its outer arc
turned between the two points it meets the axis at. Only other curves give extrusion or revolution
surfaces. The ellipsoid primitive is the unit ball stretched along its frame's axes by the affine
map (`edit::affine`), so its face is the exact rational image of the sphere's rational form; equal
radii give the sphere itself.

A periodic B-spline of any degree up to 3 is made from its knots and control points
(`BSpline::periodic`, the first `p` points repeated after the rest and the knots a period apart at
both ends) or from a control polygon alone (`BSpline::periodic_uniform`, one unit knot span to each
point). A periodic quadratic on a square is four parabolic arcs, each from one side's middle to the
next; as a profile segment it is a closed loop on its own, and extruded or turned it gives a side
closed round on itself, which a seam cuts as on any periodic surface.

A path of legs (`sweep::along`) is laid out frame by frame first. A mitered corner (`Leg::Miter`)
turns the heading on the spot like a sharp turn, but fills nothing: each straight leg beside it is
extruded on past the corner by the profile's reach times `1 + tan(θ/2)`, θ the turn, and cut by
the plane through the corner square to the sum of the two headings, which bisects the turn. The
section meets that plane in the same outline from either leg, so the two cut faces coincide and
the union removes them. A path whose last leg is a miter that brings it back to its start, heading
and section axes as they started, is closed, and that miter cuts the first leg's start too.

A sweep up to a solid (`sweep::extrude_to_next`, `sweep::revolve_to_next`) runs the profile on
past the solid's box (a turn on to 359°, so its end cap never lies on its start) and cuts that long
sweep by the solid: by difference for a sweep into the solid (`Meet::Into`, a boss joined to it),
by intersection for one out of its material (`Meet::OutOf`, a hole cut from it). The pieces kept
are those that carry a face of the start cap, whose names follow every piece through the boolean:
each runs from the sketch plane to where its part of the profile first crosses the solid's
boundary. A kept piece that also carries the end cap crossed nothing, and is an error.

Lofts ([`loft/`](../src/loft/mod.rs)) match sections edge to edge when every section has as many,
otherwise split at the polar angles of every section's vertices. Edge to edge, each section's loop
starts where its vertices, taken about its centre, best fit the last section's, read either along
world axes or along each section's own axes. The world reading keeps one outline drawn on planes
with different axes, or on a plane facing back down the loft, a straight prism; the frame reading
keeps a helical gear's tooth space, turned past about 4.6 rad over nine sections and moving further
between sections than it is wide, from pairing each vertex with the next one round and coming out
inside out. On a tie the world reading wins. Each column of matched edges is brought to one degree
and knot vector and skinned into a rational B-spline surface: linear across a ruled loft or two
sections, otherwise the lowest-degree Bézier through the end sections that passes within `LINEAR` of
the others. A flat ruled side between two straight edges is a plane. Each skinned surface runs a
tenth of its parameter range on past every edge that is not a seam, as a swept surface runs on past
its profile's ends, so an intersection that ends on a loft edge is traced across it and meets the
edge's vertex.

Every loop of a section is lofted: the sections have as many holes, and each hole is matched and
skinned with the hole at its place in the next section by the same rules as the outer loop, its
sides facing into the hole; the caps carry a loop for each. `loft_to_point` adds the apex as a last
row whose every piece is the one point: each column's rail there has all its control points at the
apex, so a ruled side between two straight edges and the apex is a flat triangle, and any other side
is a surface with a pole along its top, which its face meets in a degenerate edge of its own.

`loft_developable` joins two convex sections on parallel planes by rulings between points where
their tangents run the same way, and so builds their convex hull: the plane through a ruling and
the two parallel tangents at its ends touches the side all along the ruling. Each outline is read
by its tangent's direction, which sweeps once round the circle: a straight edge holds one direction,
a curved edge a range, a corner the range between its two edges. Both outlines are cut at every
direction where either has a corner or the end of an edge. A straight edge pairs with what the
other outline holds at its direction, and each range between cuts pairs with what the other holds
over it, a curved edge being split where its tangent reaches the cut by bisection on its direction.
A corner that pairs with an edge becomes a piece that is one point, as an apex is: a straight edge
gives a flat triangle, a curved one a cone over it with its apex at the corner, a ruled rational
surface whose pole is the corner. Where both outlines hold a corner over a range, the two corners
are joined by a side edge alone.

`loft_guided` lofts along guide curves, one through each vertex of every loop of every section.
Each guide is cut into a rail, the plane of each inner section found on it by sign changes of the
distance over 512 samples and bisection, and the vertex there matched; the guides of a loop are
ordered by the first section's vertices, and every section must hold them in that order. Each guide
is then reparameterized piecewise linearly, its crossing parameters made knots of full multiplicity
and moved to the guides' mean, so that all of a loop's guides cross each section at one parameter.
The side between two neighbouring guides is the Gordon surface of its network (Gordon 1969; Piegl
and Tiller 1997, §10.5): the ruled surface between the guides, plus the section edges interpolated
across the sections (linearly through two rows, by the not-a-knot cubic through more), less the
same interpolation of the corners ruled across; the three are brought to common degrees and knots
in both directions and their control points summed. It then runs along both guides and through
every section edge, and where the edges are straight it is the ruled surface between the guides,
each ruling joining the points of one parameter. The sum is only the summed surface for polynomial
nets, so a guided loft refuses rational curves. A loft to an apex ends every guide there.
`Helix::spiral` gives a guide that is the path of a point screwed about an axis while its distance
from the axis changes linearly with the turn (a helix when it does not change), in degree-7 Hermite
pieces through its points and first three derivatives, spaced as a screw sweep's columns are with
the radius raised to `r + 8|growth|` to bound the eighth derivative.

## Tolerance policy

There is one linear tolerance, [`tol::LINEAR`](../src/tol.rs) = 1e-7 mm, and one angular tolerance,
`tol::ANGULAR` = `LINEAR / MODEL_SIZE` (1e-10 rad for a 1000 mm model), so a direction within it
moves a point at most `LINEAR` across the model. Every geometric decision uses one of these or a
value derived from one where it is used (a parameter tolerance is `LINEAR / speed`). Nothing the
kernel makes carries its own tolerance, and no module defines a private epsilon for geometry. There
are two exceptions. Geometry imported from another system that meets only to that system's tolerance
and cannot be healed to `LINEAR` (below, under STEP) may carry local tolerances on its edges and
vertices (`Brep::edge_tolerance`, `Brep::vertex_tolerance`), as in tolerant modelling (Jackson
1995). And a vertex a boolean places where two solids graze carries the uncertainty of the graze
(below, under booleans). A vertex's sphere holds the tube of every edge ending at it, a tolerance is
the gap measured and never more than a quarter of an edge it holds, two entities coincide within the
sum of their tolerances less the one `LINEAR` two exact entities share (`tol::between`), and
geometry stays exact. Numerical routines (Newton steps, quadrature convergence) have their own
stopping rules, which are about arithmetic, not geometry.

Tangency is one decision too: two surfaces (or a ray and a surface) are tangent where the sine of
the angle between them is at most `tol::TANGENT` = `LINEAR / RESOLVED`, the angle at which they stay
within `LINEAR` of each other for `RESOLVED` (0.1 mm) along their meeting. In closed form a plane
and a cylinder (or two circles) touch when they overlap by at most `LINEAR`; past that they cross,
however thin the sliver between the crossings. A sliver thinner than twice `LINEAR` is cut along
both its edges, so where a test at one of its points falls within `LINEAR` of an edge, the exact
side of that edge (the winding number, or the sign of a ray's hit) decides it.

Two kinds of constant remain that are not tolerances: parameter-plane gaps between images of one
vertex (`arrange.rs`, `sew.rs`), which only need to be far smaller than a period (or twice a
tolerant vertex's tolerance, where its edges end anywhere in its sphere), and orderings and guards
against rounding (a flat triangle in the mesher, a direction `march` cannot take). One such guard is
in sampled curve–surface roots (`csi.rs`), for the pairs with no closed form (see
[Analytic contact](#analytic-contact)): where a curve only touches a surface, rounding in the
distance makes sign changes up to √(noise / curvature) along the curve from the touching point, 3e-8
mm on a 3 mm round, and more where the curve sits a hair inside the surface: a seam 3e-11 mm outside
a cylinder a fillet's tube touches dips into the tube and out again 9.5e-6 mm either side. A sign
change beside an extremum within `LINEAR` of the surface, where the curve leaves the surface on one
side both ways within a sample step, is moved onto the extremum, found from the sign of the slope,
which pins it to rounding; a dip between two samples is found at its extremum too, a touch when
within `LINEAR` and two crossings, one either side, past it. Without it a round's edge that ends
where it touches its pillar gave a second vertex 1.3e-7 mm from the first, and a sliver edge between
them. A boolean takes an edge's meeting with a face to be a vertex already within a cell of it that
lies on both, when the edge stays within `LINEAR` of the face all the way there: two crossings with
no excursion past the tolerance between them are one touching point. A box's edge along a loft's
side, tangent at a ruling, was met 1.3e-7 mm from where the ruling pierces the box, and the side was
never split there. An edge that meets a face on the face's boundary while running tangent to it
(sine at most `TANGENT`) runs along a boundary edge of that face: that is an edge–edge contact,
whose vertices are each edge's ends lying on the other, and the sampled roots along the band are
rounding, so they are not meetings. A circle lofted to a square met by its copy turned 90° about the
axis had its bottom arcs split at 46 such roots, 62 edges where 16 do.

Where two faces nearly coincide (a body and a copy moved a hair off every axis), an edge of one
grazes the other and its meeting with it can lie anywhere along a band LINEAR / sine long, off the
curve the two faces are cut along, which is then never split there. So an edge's meetings with faces
are placed after every face pair's intersection curves are known: a meeting is moved to where an
intersection curve of the edge's own face and the face it meets crosses the edge (Gauss–Newton on
both parameters, then along the curve onto the edge's other face, so the point lies on every surface
of the meeting), the steepest such crossing within the band where the edge stays on the face. A
fitted curve carries its fit error into that point, so where the three surfaces (the edge's two
faces and the face it meets) fix a point within `LINEAR` of it, by Newton steps on their tangent
planes, the meeting goes there: a chamfer's end on a keyway wall otherwise moved 4e-10 mm along the
bore and changed the chamfer's volume by 1e-9 of itself. Inside its own face the edge crosses that
curve at least as steeply as it crosses the face, so the crossing is no worse placed, and it lies on
the curve. It is then taken to be a vertex already there only if that vertex lies on the same curve.
Meetings are placed steepest first, ties broken by position, so where several computations find one
grazing meeting the vertex does not depend on which ran first.

A curve the kernel fits, where no exact form exists, is fitted finer than `LINEAR` by the smallest
angle it will be cut at. A screwed profile's helix crosses a plane square to its axis, or a cone
about it, at its lead angle or less, a hundredth of a radian on a fine thread, and its crossing
slides along it by the fit's error over that angle: [`helix`](../src/sweep/helix.rs) holds each
column of a screwed face within a ten-thousandth of `LINEAR` of its helix. Where a screwed edge is
no helix but an arc about the axis or a line along it, its face is the exact cylinder.

Validity is measured against `LINEAR`: an edge's ends lie within it of its vertices, and every
coedge lies within it of its face's surface. A result that cannot meet that is an error, never a
solid with a widened tolerance. An imported edge or vertex with a local tolerance is measured
against that instead. A boolean cuts such solids as Jackson's does
([Local tolerances in booleans](#local-tolerances-in-booleans)); its own geometry stays exact
wherever it does not touch a tolerant entity.

## Validity

[`check::check`](../src/check.rs) returns the issues that make a solid invalid; an empty list means:

- every non-degenerate edge is used by two coedges, in opposite directions, of faces in the same
  shell (closed, consistently oriented), or, where pieces of the solid touch along it, by an even
  number paired off into partners that run opposite ways in one shell
  ([Results that touch](#results-that-touch)); every degenerate edge by one; and no two faces beside
  each other across an empty wedge round such an edge lie on each other, every interior point of
  one on the other's surface inside it: they would be faces inside the solid (a union that kept
  both of a spline bore and the cylinder it lies on passed without this);
- every loop is connected vertex to vertex and closes in its face's parameter plane;
- each edge's curve ends at its vertices, and each coedge lies on its surface, within `LINEAR`;
- no intersection curve turns back on itself: its direction at neighbouring samples a quarter of a
  knot span apart never differs by more than a right angle (a face bounded by a fold touches itself,
  and no mesher can triangulate it). Imported splines are not tested, since a designer may draw a
  cusp;
- each loop runs the way its face's orientation requires, except a vertex loop, which is a point,
  and a loop that runs only along slits and back, a line contact inside the face, which encloses
  nothing;
- each inner loop lies in the region the face's other loops bound, by the exact winding number of
  those loops about a point on it (a face whose loops wind round a period, as a band round a
  cylinder, has no such region and is not tested), and each void shell lies inside its solid's outer
  shell, by classifying a point inside one of its faces. A sewn blend that kept a hole or void the
  blend took away failed only these two;
- vertices are finite and no two lie within `LINEAR`, and every edge is longer than `LINEAR`;
- each solid's volume is positive;
- no two faces that share no edge meet, away from the vertices they share: the solid neither passes
  through itself nor touches itself anywhere its edges and vertices do not hold the contact, and no
  two solids of the arena pass through or touch each other;
- no solid of the arena lies inside another, found by classifying a point inside one of its faces.
  `check_each`, for an import's placed parts, skips this and pairs only faces of one solid. Pairing
  faces across solids and the nesting probe changed no result on the designs they were run over.

The last is a search over pairs of faces whose boxes meet, sorted along `x`: each edge of one face
is intersected with the other face's surface and every meeting located in that face. An edge is
first halved into stretches, and a stretch is dropped where its box lies clear of the surface (a
plane's side, an analytic surface's distance from the box's centre, or the projection of the
stretch's middle against the box's farthest corner), so an edge that only passes near a face costs a
few distances. An edge that lies in the other face's surface is located at its quarter points; it
can only enter the face elsewhere by crossing the face's boundary, whose edges then meet the first
face. Two faces can also cross along a curve that meets no edge of either away from vertices they
share, as two balls turned apart about a line through both overlap between the points where it
pierces them: where two faces share a vertex but no edge and their tangent planes there cross, a
step along the line they cross in, either way, is brought onto both surfaces by alternating
projections, and if it lands inside both faces, they cross. That step is taken only where the
curves the two surfaces meet in are not all known; where they are, the search below replaces it. A
pair found apart is recorded in both faces' memos. Two faces cut from one input of a boolean are
searched as any others: the input held every contact between its faces in an edge or a vertex, and
the result must hold it still. A block with a hole through it and a rod standing on its front face
held the contact as a slit in that face, from the top edge to the hole and on below the hole to the
rod's end; cut by a plane away from the contact, the face lost the stretch below the hole, the rod's
rim touched it at a point nothing held, and the result passed, since both faces came from one input
(the cause, in the face split, is under [Booleans](#booleans)). It refuses two results the checker
used to pass: a 90° revolve cut by its copy moved 1e-3 or 1e-5 mm along `x` leaves the copy's corner
edge 5e-8 (5e-12) mm from the original's inner wall, a face it shares no edge with, so the result
touches itself there.

The checker settles the sign of a solid's volume without the integral where a ray does, and
reports what the integral alone would. For a solid of one shell whose faces are all joined through
partner coedges, no edge used by more than two, a ray leaves the middle of an edge along the sum of
its two faces' outward normals: a point just along it lies on the outward side of both faces,
whether the edge is convex or concave. The sheet winds round that point once for each face the ray
then leaves the material through, less once for each it enters through, and round a point outside
a closed sheet it winds no times. A count of zero therefore means the normals point out of the
sheet, and a sheet whose faces meet only along the edges and at the vertices they share bounds one
region, so its volume is positive. Up to twelve edges are tried for a ray that meets no face's
surface within the face's box, which needs no face's trimmed domain: the crossings of a plane,
cylinder, sphere or cone come from a closed form, a torus's from its quartic along the line, and a
spline surface's from the patch search of [`csi.rs`](../src/boolean/csi.rs), none of which steps
over a crossing. The first four rays are then followed through the faces they cross, each crossing
located in its face; only a crossing of the four closed forms is followed, since its foot on the
surface is certain. A ray decides nothing where it passes the box of a face on any other surface,
meets a face on its boundary or at a tangent, or starts where a third face lies, and an edge kept
with a local tolerance, or between faces that fold back on each other, starts none.

That the faces meet only where they share an edge or a vertex is searched, not taken, and a ray
settles a solid only where the search leaves nothing open.

Two faces can pass through each other with no edge of either meeting the other: along a closed
curve inside both, as two balls that overlap in a lens, or along a curve from one vertex they share
to another, as two half cylinders whose surfaces touch at both. For two faces on analytic surfaces
whose every curve of intersection a closed form gives, those curves are cut at the vertices the
faces share and the middle of each stretch is located in both faces: a stretch lies inside both
whole or not at all, since it could leave either only across an edge, which the edge search would
show. A curve the surfaces only touch along is left out; analytic surfaces that touch along a curve
stay to one side of each other there. Where the surfaces cut across each other at an end of an edge
the faces share, they meet near it along one curve, the edge's own, and the stretch out of it is
the edge or the edge carried on past its end, which a face whose corner there is convex does not
hold: such a stretch is not located. The search of faces that share no edge makes this test after
its edges, and reports a crossing whatever shells the faces lie in; a pair it clears is recorded in
both faces' memos and not searched again. The curves are worked out once a check for each two
surfaces, within a box that holds every face. Where a pair's curves would have to be traced from
seeds (a cylinder and a cone, two tori) or fitted along a sweeping family (two cylinders askew, a
cylinder and a sphere off its axis: the fit costs what cutting the pair did), or a face lies off the
analytic surfaces, the step out of a shared vertex is taken as before and the pair's solid is left
open. A ball turned one radian about a line through it was built as two hemispheres sewn to a piece
of torus, which pass through each other beyond the axis along an arc between the two points the axis
pierces: the check refuses it, and the sweep is built in parts.

Faces that share an edge are searched before a ray is believed, in the same two ways: no edge of
one that the other does not share may meet the other face, and no stretch of the surfaces' curves
may lie inside both. Two flat faces with straight edges are told apart outright where every vertex
of one off the shared edges stands to one side of the other's plane. A pair found to meet nowhere
else is recorded in both faces' memos with the
places of the shared edges in each face, and is not searched again while both keep their geometry
and share the same edges. An edge against a face off the analytic surfaces costs what the integral
does, so a sheet with such a pair not yet recorded is left to the integral. Two block-shaped lobes
joined corner to corner, one facing out and the larger facing in, whose L-shaped faces share an edge
and pass through each other beyond it, were taken as facing outward by a ray off the first lobe
before edge-sharing faces were searched.

A face must not cross itself either. On an analytic surface a face whose boundary does not cross
itself covers no part of the surface twice: those surfaces fold onto themselves only at a cone's
apex, a sphere's poles and the points where a torus that reaches its axis meets it, where a face
can touch itself and not pass through. Two edges of a face that cross are shown by the search of
edge-sharing faces wherever some other face uses one and not the other. What that leaves is tested
directly: two edges the same other face uses, or the face itself along a seam, are intersected
where both are lines, circles or ellipses, and every edge must be seen not to cross itself, a
spline by its control points (it runs one way along its chord, or turns one way about an axis
through their middle for no more than a turn). A spline face is not searched for crossing itself,
so a sheet with one is left to the integral, as is one where any of this cannot be shown.

So a ray settles a solid whose faces all lie on analytic surfaces, whose boundaries are seen not to
cross themselves, and whose every pair of faces with meeting boxes was searched with every curve of
intersection known. A pair an earlier check found apart by its edges alone, on surfaces whose curves
cannot be asked for, leaves its solid open.

Whenever a check that a ray settled goes on to find any issue, it is made again by the integral
alone, so the issues are those the integral gives first; a solid a ray settled that the search
leaves open takes the integral after all. What no ray settles, a solid of several shells, a sheet
that touches itself along an edge, and a check without the search for crossings (an import's
structure) take the volume integral.

Every operation runs the checker on its result and returns `Error::InvalidResult` instead of a
broken solid. Booleans and blends must keep that contract.

`check::check_each` holds the solids of an arena to the same rules one at a time: vertices of
different solids may lie within `LINEAR` of each other. That is the contract of a STEP import, whose
placed solids are separate parts and may touch (a part seated on another shares its corners); every
solid the reader returns has passed the checker on its own.

## History and face names

Every face carries a set of [`FaceName`](../src/topo/names.rs)s, `feature.face`. A sweep names each
side face after its profile segment and its caps `start` (on the sketch plane) and `end`; a full
revolution has no caps. The rules every operation follows:

- **A split face keeps its names.** Every piece of a face cut in two carries all the names it had.
- **Merged faces carry the union.** When same-domain faces merge, the result carries every name of
  every face merged into it.
- **Blend faces are named by their operation.** A chamfer or fillet face is `<op>.bevel_<i>` or
  `<op>.round_<i>` for the edge pair `i`, and a vertex blend `<op>.corner`.
- A face the result loses simply disappears; nothing is renamed to fill its place.
- **A reshaped face keeps its names.** A face a later operation trims, holes, tilts (a draft) or
  moves along its normal keeps every name it had; an outward shell's inside is the modelled faces
  where they were, so each keeps its names and takes the shell's `inner` too.

Names never depend on geometry or on arena order, and no operation names a face by the surface it
lies on: a face lying on another feature's surface is not that feature's face. So a face with no
names was made by a face its feature left unnamed, and an edit after a feature never moves that
feature's names.

One rule reads surfaces, and it reads them from faces history names: a blend's edge pair stands for
the edges where the surfaces of the faces carrying its two names meet.

## Booleans

[`boolean`](../src/boolean/mod.rs) combines every solid of one arena with every solid of another, in
four steps, each in its own module:

1. **Meeting points.** Every edge of each solid is intersected with every face of the other
   ([`csi.rs`](../src/boolean/csi.rs)), every edge with every edge of the other (which places the
   vertices where coplanar boundaries cross), and every pair of faces' surfaces
   ([`ssi.rs`](../src/boolean/ssi.rs)). Each point goes into one pool, merged within `LINEAR`.
2. **Pieces.** Every input edge and intersection curve is split at every pooled vertex that lies on
   it; pieces outside either face are dropped, and pieces that run along each other become one edge,
   an input edge in preference to an intersection curve. So the two faces of every result edge use
   the same edge by construction, as sweeps do.
3. **Face pieces.** Each face is cut by the intersection pieces inside it and by the boundary of any
   face of the other solid that lies on its surface ([`arrange.rs`](../src/boolean/arrange.rs)). The
   cuts and the boundary form a planar graph in the parameter plane; its regions, traced with the
   material on the left, are the pieces.
4. **Selection and sewing.** Each piece is classified at an interior point, one off the other
   solid's faces where the piece has one (a point on a tangent contact is on both sides): on a
   coincident face of the other solid (same or opposite outward normal), or inside or outside it. A
   piece bounded by a transversal cut is read at the cut's middle where that lies inside the other
   face: it runs into the other solid when the direction into the piece, square to the cut along its
   surface, opposes that face's outward normal (sine at least 1e-3). A piece whose point still lies
   on the other solid's boundary, where no ray has a parity, is decided by that contact: inside one
   face, by the side of it the piece stands on where it stands clear of the face by more than
   rounding (a thousandth of the tolerance), as a sliver cut where two faces cross at a small angle
   does all over, and otherwise by the side the piece's surface bends to (the sign of the gap form
   of [`contact.rs`](../src/boolean/contact.rs), semidefinite along a tangent line); on a seam, as
   inside the face; on an edge the piece's tangent plane holds, by whether that plane runs into the
   wedge of material the edge's two faces bound on both sides of it (along a face it is tangent to,
   by the gap again). A corner, or surfaces that osculate, leave it to the rays. Any other piece
   takes the parity of a ray's crossings of the exact faces
   ([`classify.rs`](../src/boolean/classify.rs)). In a symmetric design the piece's best interior
   point is often on the contact: all 40 of seq seed 2's
   `could not tell whether a face piece is inside` failures were such a point, 18 on a tangent line
   inside a face, 16 on an edge and 6 on a cylinder's seam where it touches a plane. A point is
   located on a face by every image of it in the face's parameter window, since a face that winds
   along a helix can span more than a period, where only one image is inside. The operation keeps
   pieces by that class; [`sew.rs`](../src/boolean/sew.rs) merges kept faces on one surface that
   share an edge only where every use of it lies on that surface (an edge more faces use is a
   contact, and stays), pairs the uses of each edge, groups shells into solids, and runs the
   checker. Where the result touches itself or a piece of it touches another, the contact is kept in
   it ([Results that touch](#results-that-touch)).

### Results that touch

A boolean is the regularized set operation, and its result may touch itself or have pieces that
touch, along a curve or at points: two blocks sharing an edge, a box less a cylinder touching all
four of its sides, two balls touching. The result is one solid, the point set's, and it keeps the
contact:

- **Along a curve**, the contact is an edge `2k` faces use. Round it, faces alternate between wedges
  of material and empty wedges. Each coedge has a _partner_, the coedge across its wedge of
  material, and a coedge _beside_ it across the empty wedge (`Brep::partner`, `Brep::beside`). The
  order round the edge is read a short step into each face, so faces that leave the edge tangent to
  each other (a cylinder lying on a plane) still part. A tangent curve is a cut like any other; a
  cut that dangles inside a face, as a groove's apex edge resting on a wall does, stays in it as a
  slit the face's loop runs along and back.
- **At a point inside a face**, the contact is a _vertex loop_: a loop of one coedge on a degenerate
  edge at the point (`Brep::is_vertex_loop`). The boolean makes one where two faces' surfaces meet
  at an isolated point, or where an edge of the result runs tangent to a face, and the point is a
  vertex both sides share.
- **Sheets and solids.** Faces joined through partners are a sheet; sheets joined through shared
  edges or vertices are one solid, whose first `Solid::outer` shells are outer sheets and the rest
  voids.

[`measure::betti`](../src/measure.rs) gives the point set's Betti numbers. Resolving each contact
into the manifold boundary of the interior gives a surface `M`, so χ(X) = χ(B) − χ(M)/2, `B` the
boundary. An edge `2k` faces use is `k` edges of `M`. A vertex where `f` fans of faces meet is
replaced by its link, the material on a small sphere round it, whose regions the fans' circles part
by turns into material and empty: two balls touching leave two disks, and a cavity touching the skin
from inside an annulus. Which it is follows from whether each fan lies on the other's material side,
read a step from the vertex. `b2` counts the groups of faces joined across empty wedges that enclose
negative volume.

Operations that cannot keep a contact refuse it with `Error::Unsupported`: a blend along a contact
or up to one, and moving a face that runs up to one (a draft, a thicken, a shell). STEP holds a
touching result as the manifold boundary of its interior ([STEP](#step)).

### Solids that are each other

Before any of this, a solid of one input that is a solid of the other, every vertex, edge and face
of one coinciding with one of the other's within their tolerances
([`same.rs`](../src/boolean/same.rs)), is taken from topology: the pair is its own union and
intersection and leaves nothing of a difference, and only the rest is intersected, the facts the
caller gave renumbered onto the copies of the rest as they are where solids pass through untouched.
Faces that lie on each other everywhere are where tracing is least stable; an imported assembly of
605 solids united with itself takes 12 s this way. A distance at exactly the tolerance is read with
a few ulps of the coordinates' size to spare, so a copy moved by `LINEAR` coincides however its
coordinates round: `stress` had 51 `FAIL`s and 21 `REFUSED` among its 81 copies moved 1e-7 mm, and
has none.

Solids a hair apart, farther than the tolerance but close enough that boundaries touching in one lie
within it of the other's over a stretch, are where the face split is least stable. A loop that runs
out along a piece and straight back loses the spike, unless the piece is of an edge the face came
using on both sides: that is a contact the face holds, a slit its loop runs down and back, and a
face that is cut keeps it. Departures along nearly one tangent are
compared by chords a twentieth of the shortest trim compared, and ones equal to the last bit by
where their far ends lie. A loop whose chord area is not clearly larger than its error takes its
sign from the exact integral along its trims. Vertices that split a curve less than the tolerance
apart along it become one. Edges of the two solids running within 1e-5 mm of each other meet where
the offset between them turns round, a crossing at angles no step resolves. Where two surfaces cross
at a tiny angle, an intersection point the Newton steps miss is found by the gap's sign change along
a line across the curve.

Faces of the second input that lie within 1e-6 mm of a face of the first, all along their boundary,
are moved onto that face's surface first ([`boolean/snap.rs`](../src/boolean/snap.rs)), the rest of
the solid rebuilt round them by [`tweak::move_faces`](../src/edit/tweak.rs). Such surfaces cross, if
at all, at an angle the tolerance cannot place: an insert placed by a pose rounded to 0.1 µm sits
0.5 µm off its pocket's axis, and the crescent between them would run within the tolerance of both
cylinders for half a millimetre at each tip. Moved onto one surface, the part sits in its seat, and
the volume moved is at most 1e-6 mm times the face's area. Faces up to 1e-4 mm apart are moved too
where they cross at an angle below `TANGENT`: a cylinder 1 µm off its seat's axis at 24.5 mm radius
crosses it at 4e-8 rad. Two planes are not: they meet along a line found in closed form at any
angle, and the wedge between them can be all the material there is (a cube's top and a plane 1e-6
off it across the top share a wedge 5e-6 mm thick at its end, which moving the plane took away). Where moving all of them at once leaves a corner whose surfaces do not meet,
each is tried alone and kept if the solid still closes; a face tangent to a neighbour cannot move
off it, and is left for the boolean. Planes within the tolerance of each other over both faces are
one plane however small the angle between them, and planes that cross a hair off parallel meet along
a line found without the difference of two numbers `1 − c²` apart, which put it 1e-5 mm off both at
1.5e-5 rad.

Copies moved by one to a few tolerances, 1.1e-7 to 7e-7 mm, failed in 84 of `stress`'s 468 runs of
that band, all for one reason: something of the copy was left a little over the tolerance from the
original, and the boolean cut a piece of edge or face that thin between them. Four rules close it:

- **Vertices free on a curve.** A seam's end lies on a circle, and moving its faces puts it the
  least distance onto them, as far round the circle as the copy was moved. Each vertex of a moved
  face takes the original's vertex within reach where that lies on all its new surfaces
  (`tweak::move_faces_pinned`).
- **Faces already within the tolerance.** One left where it is makes the vertex it shares with a
  nearly parallel neighbour, which does move, slide along their line by the gap over the angle
  between them. Such faces move with the others, and on their own where they lie half the tolerance
  off or more all over, where a distance read with a few ulps to spare and one read without decide
  differently. A face off by that much at only some of its boundary, a part tilted by rounding,
  stays: moving such a cover that way broke its pair with the plate under it.
- **Spline faces.** Extrusions, revolutions, ruled and B-spline faces of one kind snap by the same
  sampled distance; one within the tolerance all over already lies on the other, as a symmetric part
  turned onto itself does, and moves only with the rest.
- **The checker** reads a coedge's distance from its surface with the same few ulps to spare
  (`tol::within`): an edge at exactly the tolerance from a face, as two faces turned 5e-6 rad apart
  leave one, measured 1.0000000028e-7 mm.

Copies moved by 1.1e-7, 2e-7 and 5e-7 mm combine too. A face tangent to a neighbour that does not
move with it still cannot be moved: a plate whose outline is a chain of tangent arcs 1e-6 mm off
another's is the one known pair left. Holding the gap as a local tolerance on its edges, as an
import does, was tried and dropped: it broke 15 `stress` runs that pass when the face is left to the
boolean.

Step 1's pairs are taken from boxes, and where they are many (256 face pairs or 4096 edge pairs
whose boxes overlap) most of them are far apart: a crank and a cam whose boxes overlap had 1,174
face pairs that never meet, each traced from up to 162 seeds at 6 ms a pair. There
[`boolean/near.rs`](../src/boolean/near.rs) meshes the faces of each input whose boxes reach a face
of the other, and drops every face pair, edge pair and edge–face pair whose mesh pieces stay farther
apart than twice both meshes' deflections plus the tolerances and the 1e-5 mm the edge-graze search
looks within: a mesh lies within its deflection of the exact geometry, so where the exact pieces
meet their meshes come that close. A face the mesher leaves bare counts as near everything. The
crank and cam took 11 s and take 0.8 s; an imported assembly of 25,676 faces against a plate took
368 s, mostly in the edge-graze search over 64 samples of every edge pair. Below the threshold
meshing costs more than it saves: always on, it made builds of small parts 40% slower.

A large input touched in one place pays per face for work that is only needed where it is touched.
Three steps skip that work and decide everything as before:

- **Boxes.** [`measure::face_hull_box`](../src/measure.rs) holds a face's loose box and costs a pass
  over its edges: control points for splines, a sphere's or torus's whole extent. Solids, and the
  faces the snap step compares, are ruled out by hull boxes first. The loose box is worked out only
  where hull boxes meet.
- **Splitting curves.** Each curve's vertices are found through a hierarchy over the vertex pool
  rather than a pass over all of it.
- **Classification.** Where the meshes are known, a face clear of the other solid and joined to
  another such face by an edge that is also clear lies on the same side of it, so one ray classifies
  each connected group. A face whose box misses the other solid is outside it with no domain built.

That assembly against a dowel took 0.9 s and takes 0.02 s. Against the plate it takes 3.7 s, where
most of its faces lie in the plate's box: their domains are still built (1.4 s), their boxes worked
out for the snap step (0.9 s) and their meshes made (0.6 s).

A piece of a face narrower than twice the tolerance has no point clear of its boundary, so pieces of
curve are judged by more than their middles: one lies on a face's boundary only if its ends and
middle all lie along one trim, and otherwise the exact side of its middle decides whether it cuts
the face. A piece on a plane that runs into the other solid on one side of its edge only is on the
side of the edge's other face that its point is. Without these, faces between the tolerance and
twice it apart failed, though a copy moved by either value passed.

A boolean's result has its redundant vertices dissolved
([`boolean/simplify.rs`](../src/boolean/simplify.rs)): where only two edges meet, both between the
same two faces and on one curve (one line, one circle, or the same curve running on), an imprint
split an edge the result keeps whole, and the two become one. Two unit cubes united face to face are
a box with 8 vertices and 12 edges. A vertex or edge with a local tolerance is never merged.

### Local tolerances in booleans

An imported solid's tolerant edges and vertices ([STEP](#step), healing step 4) go through a boolean
as tolerant modelling sets out (Jackson 1995, section 4): faces stay exact, and every coincidence
test between two entities uses the sum of their tolerances (`tol::between`).

- **Vertices.** Every pool vertex carries a tolerance: an input vertex its own, a meeting `LINEAR`,
  and a meeting on a tolerant edge (an edge crossing a face, two edges crossing) that edge's. Two
  vertices within the sum of theirs are one; where either is tolerant the one kept grows to hold the
  other's sphere, and a curve split at a vertex its piece ends short of grows it to that distance.
  Two vertices that split one curve within `LINEAR` of each other along it are one, and a piece of
  any other curve between the two, no longer than they lie apart, is dropped: that curve passes
  the vertex once. No tolerance grows past the widest of the inputs', which an import keeps within
  its file's uncertainty, or the reach of a grazing meeting at the vertex; past that, the checker
  refuses the result and names the gap.
- **Edges.** A tolerant edge lies on a face of the other solid where it stays within its own
  tolerance of the face's surface; it is split at every vertex within the sum of the two tolerances;
  and two pieces between the same vertices that run within the sum of theirs are one.
- **Faces.** An input face's boundary is within each trim's own tolerance, so a point on a tolerant
  edge locates on the boundary although it lies off the surface. In the parameter plane, the images
  of a tolerant vertex at the ends of its trims lie up to twice its tolerance apart: they are one
  node when splitting faces and one corner when merging them.
- **Face pairs.** Where a tolerant edge of one face lies on the other face, the two meet along it:
  within twice its tolerance of it their surfaces lie that close (two B-spline faces of one imported
  solid touch rather than cross, which is why healing failed there), so no trace seeds there, and
  where a trace of the pair still fails the edge is taken as their whole contact, as Jackson's
  face–face step leaves contact already known to edge–face coincidence. The checker still reads the
  result for faces that pass through each other.
- **Result.** Each result edge carries its input edge's tolerance, intersection curves `LINEAR`, and
  each vertex its pool tolerance, so the checker holds each to its own.

Every boolean of an imported assembly's three tolerant solids with a box or a rod through them, at
three places each, builds a valid solid (54 of 54; each union, intersection and difference). A cut
through a tolerant edge also splits faces whose edges miss them, so the pieces' volumes add up to
the whole only to within the solid's area times its widest tolerance: the box that cuts one of them
leaves pieces 0.014 mm³ (9e-6) over it.

### Imprint order

A boolean on tolerant models imprints in strict order of dimension (Jackson 1995, section 4.3):
vertex–vertex, vertex–edge, edge–edge, vertex–face, edge–face, face–face, and face–face need not
consider contact already known from edge–edge or edge–face coincidence, since near-tangent curves
are unstable to compute. `Work::meet` and the steps after it run in this order:

| Step               | Meets | How                                                                   |
| ------------------ | ----- | --------------------------------------------------------------------- |
| `Work::new`        | VV    | every input vertex pooled; two within `LINEAR` are one                |
| `pierce`           | EF    | `csi` per edge and face of the other; held, not yet placed            |
| `pierce`           | EF on | an edge lying in a face's surface is recorded as lying on it, not cut |
| `cross_edges`      | EE    | `curve_curve`; each crossing pooled                                   |
| `intersect_faces`  | VF    | per traced face pair, every pooled vertex on both surfaces            |
| `intersect_faces`  | FF    | closed form, or traced from seeds; uses the facts above (below)       |
| `place_hits`       | EF    | each held meeting moved onto the face pair's curve, then pooled       |
| `split_curves`     | VE    | every pooled vertex within `LINEAR` of a curve splits it              |
| `imprint_on_faces` | EF on | an edge lying on a face that its own faces cross there cuts that face |
| `dedup`            | EE on | pieces with the same ends that run along each other become one        |

A grazing meeting is one vertex. Two surfaces each known to `LINEAR` fix where they cross only to
`LINEAR` over the sine between them, so where an edge of each solid crosses the other's face at a
shallow angle, the two crossings found for one meeting lie several tolerances apart along the graze,
and so do the intersection curves of the face pairs round it (two meshing gear flanks at a sine of
0.012 give crossings 3e-7 mm apart). The vertex an edge–edge crossing or an unplaced hit left there
moves to where a later hit is placed on its curve, every curve of the two faces a hit met at is
split at that vertex within the graze's reach (`LINEAR` over the sine, at most a hundred times
`LINEAR`), and the vertex's tolerance grows to hold where each such curve passes it.

Two pieces with the same ends are one where they run within the sum of their tolerances, or, where
one is an intersection curve, within the reach to which its two surfaces fix their crossing:
`LINEAR` plus the other's tolerance over the sine between them, or the tangent band where they are
tangent. Near a point where the two faces of a curve become tangent, the trace and an edge of the
inputs along the same crossing part by more than `LINEAR` (a V-wheel turned a whole turn, built as
the union of two half turns, gave 1e-7 to 4e-6 mm). Left as two pieces between one pair of vertices
they bound a region no wider than that, whose sides the face split orders by rounding, so the result
depended on the platform's floating point.

So the order differs from Jackson's in four places: edge–face meetings are found before edge–edge
ones but placed after face–face (a grazing meeting is placed where the face pair's curve crosses the
edge, which is better conditioned); vertex–edge incidence is read by distance when curves are split,
not before; vertex–face incidence is read only for face pairs that are traced; and coincident edges
are merged last, as pieces. What face–face takes from the steps below it:

- **A vertex on both surfaces where they touch is that touch.** A known vertex (an input's, or one
  an earlier step pooled) that lies on both surfaces, where their normals are parallel within
  `TANGENT` and the gap between them is a definite quadratic form (the difference of the two shape
  operators on the common tangent plane, [`contact.rs`](../src/boolean/contact.rs)), is a point
  where the surfaces touch and share nothing else: no intersection curve passes through it. Within
  `√(2·LINEAR / κ)` of it, `κ` the gap's smaller curvature, the surfaces lie within `LINEAR` of each
  other, so a trace there finds only noise: no seed starts there, and a traced curve that stays
  there is dropped. Where the form is indefinite the branches through the point are traced along its
  null directions as before; where it is semidefinite the surfaces may touch along a curve and
  nothing is decided. A blend's straight round meeting the next corner's round at a tangent section
  touches the corner's wall cylinder at the section's end, a vertex of both; tracing the two
  cylinders from it left a curve 3.3e-5 mm long through it, an edge four faces used.
- **An edge lying on a face is that contact.** Where an edge of one input lies in the surface of a
  face of the other, the contact between the edge's faces and that face along it is the edge. A
  curve of the face pair that runs along it all the way is dropped. Along is within the edge's
  tolerance, within `LINEAR` plus that tolerance over the sine where the two surfaces cross at a
  shallow angle, or within the tangent band where they are tangent. Where the faces are tangent
  along the edge, a trace that runs on past both faces along the edge's curve beyond its ends is
  dropped too: the region it is traced in reaches past the faces, and two tori tangent along a
  meridian (a bull-nose mill's sweeps meeting at an S-bend) stay tangent round the whole circle.
  Kept, such a trace lay a few 1e-8 mm beside the edge and cut off a sliver no point of which could
  be classified. The edge itself cuts the face wherever it lies inside it and one of its own faces crosses the face there
  (`imprint_on_faces`); where its faces only touch the face it is a tangent contact, not a cut, as
  before. Seeds on the edge still start traces: where the two faces nearly coincide (a loft and its
  copy moved 1e-5 mm), an edge of one lies on the other's surface while their intersection crosses
  it, and only a seed there finds that curve.
- **Facts from the operation.** A caller may hand the boolean what its construction makes true
  ([`Facts`](../src/topo/facts.rs)): an edge of one input that lies on a face of the other. The
  boolean checks each at the edge's ends and quarter points against `LINEAR`, and an untrue fact is
  `Error::InvalidInput`, never trusted; a true one is recorded as `pierce` would record it, without
  intersecting the edge with that face. Facts name the inputs' ids; where solids whose boxes meet
  nothing of the other input pass through, the facts about the rest are renumbered into the copies
  that are combined. A blend tool's section sides are swept along the faces its edge joins, so every
  edge of a side face lies on each face of the body on that face's surface; the blend hands those
  over with each tool.

### Tangent contact, in order

Where two faces touch, three sources can say how, and they are consulted strongest first, so a
number never overrides what topology or a closed form already decided (Sugihara and Iri 1989):

1. **Facts.** What the inputs' topology or the operation states: an edge lying on a face (found by
   `pierce` or handed over as a [`Facts`](../src/topo/facts.rs) entry) bounds the face pair's
   contact along it, and a pair's curve that runs along such an edge all the way is dropped,
   whichever layer below produced it.
2. **Exact classification.** A pair with a closed form ([Analytic contact](#analytic-contact))
   decides crossing, touching at a point, tangency along a curve or coincidence from distances in
   closed form. It never asks for seeds, so nothing numeric below runs for it.
3. **Numeric tracing, last.** Only a pair with no closed form is traced, and only then are the known
   vertices where its surfaces touch at a point worked out (`known_touches`): no seed starts within
   such a touch's reach and a trace that stays there is dropped. A trace that stalls where the
   surfaces become tangent ends on the point where they touch, and where a known point of both
   surfaces lies ahead of it in the same band, on that point exactly
   ([`march.rs`](../src/boolean/march.rs)). Two blend sections that start as straight as the face
   they share touch each other there to the fourth order: no search finds such a touching to better
   than a few millionths, but the faces' edges that cross at it do. Where the surfaces are tangent
   at a known vertex and cross there (their gap form is indefinite), two branches of the
   intersection cross at the vertex and a trace passes it only to within `√(2·tol/|λ|)`, `λ` the
   form's larger eigenvalue, where that bends the gap by at least a thousandth per mm (flatter,
   the faces lie along each other over a patch, as a milled floor's do, and a vertex there is on
   no branch in particular): within that reach the vertex still splits the branch, and its
   tolerance grows to cover the gap. A sheet a sweep makes from an edge leaves the face beside it
   at such a point, where the face's grazing curve meets the edge, and comes back into a concave
   face along the second branch. A walk from a known vertex of both surfaces that comes back to
   it with no two points the tolerance apart found no curve and is dropped: where a curve of the
   two surfaces runs out of a bounded one at its corner, the points beside the corner are on both
   only to the tolerance, and the faces meet at the vertex alone. A sheet a sweep makes from one
   edge of a block's top leaves the block's side this way at the top corner, a fraction of a
   degree off the corner's own path.

A stub of this kind, from a fillet's rim meeting a pocket's corner round, came from a pair that both
later layers now cover: the r 1.5 round's cylinder against the r 4 wall cylinder, axes not parallel,
is swept in closed form ([`ssi/sweep.rs`](../src/boolean/ssi/sweep.rs)) and the touching vertex
comes back as a point. The part builds with `known_touches` switched off too; the touch rule guards
the pairs still traced (cones and tori against each other and against cylinders), where the same
stub can arise.

### Face memos

Each face of an arena carries a memo ([`topo/memo.rs`](../src/topo/memo.rs)) of what is derived from
its geometry alone: its rough volume share about its own anchor and its vector area (`∬ n dA`, from
`½∮ p × dp` round its loops), that the checker found its loops on its surface, closed and running
its way, the trims, samples and box a boolean reads, and its share of the loose box. A boolean hands
a result face the memo of the input face it copies when no cut crosses that face, none of its edges
is split or stands in for another, and the arrangement gives back its loops unchanged; so a part
built by many features works each face out once, not once a feature. Adding a loop to a face, A
change to a face and `transform` drop memos; `append` and `clone` share them. Since the vector areas
of a closed shell sum to zero, the volume sign is the sum of the memoed shares moved to one centre.

An arena also keeps what queries derive from the whole of it: the mesh that seeds
[`distance`](../src/distance.rs) and [`rays`](../src/rays.rs), that mesh's triangles under a box
tree, and the mesh pieces and box tree of every solid and of every face. A `clone` shares them. Any
method that adds an entity, or changes one's geometry, tolerance or place in the topology, drops
them (`transform` and `transform_similar` among them), so nothing kept outlives the geometry it was
made from.

[`unify`](../src/boolean/unify.rs) runs step 4's merge and sewing on any valid solid, so faces no
boolean made (an imported file's coplanar triangles) merge as a boolean's do. A solid with nothing
to merge comes back unchanged. A solid kept with local tolerances merges too: faces merge on exact
surfaces, within `LINEAR`, and the edges and vertices that remain keep their tolerances. Planes of
one such solid that lie 1e-4 mm apart across a tolerant edge (in one imported solid) are two faces,
where a kernel merging by the edges' tolerance would make them one.

A solid whose box meets no solid of the other input skips all four steps: it is copied into the
result when the operation keeps it whole, and dropped when it does not. A face whose box meets
nothing of the other solid is outside it without a ray. Where the graph of step 3 has departures
from one node whose tangents differ by less than they bend over a short chord (a fitted curve's
tangent is only as good as its fit), the chords order them; every chord from one node reaches as
far, a twentieth of its shortest curved trim, so they compare how the trims bend rather than how
long they are. The chord's end is found at that distance from the node, not at that fraction of the
trim's length: a chord's angle grows with its length, and two trims that leave tangentially and bend
alike (a cut and an edge 1e-4 mm apart) are told apart by less than a parameter's unevenness changes
it. This holds between any two departures, not only those a micro-radian apart: two trims of the
graph cross only at nodes, so their first crossings of a small sphere about the node come in the
order they leave it. Where a handle's end touches a teapot's body at one point, the cut leaves the
corner tangent to the end's edge; its fitted tangent lay 1.6e-6 rad on the wrong side of the edge
while its chord bent 9e-3 rad into the face, and the turn read from the tangents closed no region.

Surface pairs with a closed form use it, and decide their contact (crossing, touching at points,
tangent along a curve, coincident or apart) from distances in closed form, never from a trace; see
[Analytic contact](#analytic-contact) below. Every other pair is traced
([`march.rs`](../src/boolean/march.rs)) from seed points (the edge–face meetings on both surfaces,
then a grid) along `n_a × n_b`, each step refined onto both surfaces, and fitted as an intersection
curve. A step is taken only where it goes the way the tangent points and the curve passes near its
chord's middle, so refinement that lands behind the last point or on another stretch of the curve
shortens the step instead of being taken. Each point carries its parameters on both surfaces, and
the next step is refined from them where the surfaces cross at a sine of 0.05 or more; where they
nearly coincide, refinement from nearby parameters can settle anywhere along them, so it starts from
a global projection, as it does whenever the hinted refinement fails. The fit likewise measures its
error from the nearest point's parameters, and an inserted point is refined from its neighbour's
first and from scratch only when that lands on the other branch of a lens (below). A walk that the
edge of a spline's parameters stops ends on that edge exactly, where the iso-curve there meets the
other surface: its last point, as near as the shortest step got, can miss the vertex its face's
boundary makes there by more than the tolerance, and the curve's last piece would then be lost.

Where two branches cross, the surfaces touch: two blends each tangent to one face along circles that
cross, or a solid and its copy moved a hair along a shared edge. A walk that stops short of the
region there ends on the point where the normals are parallel, found by Newton steps, when it lies
on both surfaces; the point becomes a vertex, and every branch through it is traced from it along
the two directions in which the surfaces' gap, a quadratic form there, is zero. A seed within the
tolerance of such a point can walk straight into it both ways; that walk is a point, not a branch,
and is dropped, since the branches are traced from the point itself. Where no curve fits through a
branch and such a point, the walk stalled in a band the surfaces share within the tolerance and the
point lies off its curve (a loft's sides against a copy moved 3e-9 mm, which touch along a ruling):
the point and the branches traced from it are dropped, and the branch ends where the walk stopped.
Walks that stall in one such band each find their point anywhere along it, so two crossings within a
step of each other with every point between them on both surfaces are one, unless the surfaces pass
through each other across the line between them (points of one as far out either side as the
crossings are apart lie on opposite sides of the other, beyond the tolerance): two surfaces that
hold one straight line can touch at two points of it and cross along it in between. A branch that
ends on the later of one is moved onto the first: a loft's side and its copy moved 1e-5 mm cross
where they coincide for 2e-7 mm, and the cut down the side otherwise missed the vertex its top edge
was cut at. A closed branch closes only across a gap no wider than its first step, so a curve that
folds back within micrometres (a fillet's tube against a copy moved 1e-5 mm) is followed round the
fold. The fit inserts each new point in the plane square to the curve and keeps it only where
`n_a × n_b` runs the way its neighbours' do: where the surfaces are nearly tangent, the two branches
bounding the thin lens between them run opposite ways, and a point off the curve by little is near
both.

### Analytic contact

Contact between two analytic surfaces is decided once, from distances computed in closed form, and
the result's topology follows from that decision: a touching point comes back as a point, two
branches that cross where the surfaces touch end on that point, and a stretch along which the
surfaces lie within `LINEAR` of each other is a contact, not a sliver to trace. A distance within
`LINEAR` of zero is zero, as a plane and a cylinder that overlap by no more touch along a line.

| Pair                                                                                    | Result                                                                                                                        | Contact decided by                                                        |
| --------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------- |
| plane–plane, plane–cylinder, plane–sphere                                               | line, lines, circle or ellipse; tangent line; point                                                                           | the plane's distance from the axis or centre                              |
| surfaces of revolution about one axis (incl. a sphere's centre)                         | circles from their meridians; a point on the axis                                                                             | the meridians' distances                                                  |
| parallel cylinders; equal cylinders whose axes cross                                    | lines; two ellipses and their two touching points                                                                             | the axes' distance                                                        |
| plane–cone, axis not square to the plane                                                | exact ellipse; hyperbola or parabola along the generators; through the apex two generators, one tangent generator or the apex | the plane's distance from the apex and its angle to the nearest generator |
| plane–torus holding the axis                                                            | the two meridian circles                                                                                                      | the plane's distance from the axis                                        |
| cylinder–cylinder (axes not parallel), cylinder–sphere, plane–torus, sphere–torus       | branches evaluated along a sweeping family, fitted; touching points; nodes                                                    | the family's margins (below)                                              |
| plane–rational spline, its control points on the plane or to one side of it             | nothing; a one-span patch's side, tangent or crossing, or its corners; otherwise traced by `march.rs`                         | the control points' heights over the plane (below)                        |
| cylinder–cone, cone–cone, cylinder–torus, cone–torus, torus–torus, sphere–cone off axis | traced by `march.rs`                                                                                                          | Newton steps where the normals are parallel: sensitive near tangency      |

A rational spline's weights are positive, so each of its points is a mean of its control points and
its height over a plane the same mean of theirs. With every control point clear of the plane on one
side by more than `LINEAR`, the surface does not meet it, whatever its spans. On a patch of one span
with its knots at its ends, a point inside the domain has a positive share of every control point,
and a point of a side between its corners of that side's row alone: where one side's whole row lies
on the plane and every control point off the plane lies to one side, the patch meets the plane in
that side and nowhere else. The contact is the side's own iso-curve, tangent where the row beside
it lies wholly on the plane too and a crossing where none of that row does. Where no side's row
lies on the plane, or only that of a side drawn to a point, the patch meets the plane at its corners
whose control points are on it, and those come back as touching points: a blend's corner patch
touches the face it runs out on at the point one of its sides is drawn to. Control points either
side of the plane, one off it by less than ten times `LINEAR`, a row partly on it beside a side
that is on it, a further corner or side on it, a side between two touching corners whose middle is
within `LINEAR` of the plane and a patch lying in the plane are traced.

The swept pairs ([`ssi/sweep.rs`](../src/boolean/ssi/sweep.rs)) take one surface as a family of
lines or circles (a cylinder's generators, a torus's meridians), each meeting the other surface in
at most two points given by a square root. The family is laid out so that how far its member at
angle `θ` crosses the other surface depends on `c = cos θ` alone, as the least of two functions of
`c` that are each monotone on [-1, 1]: how far the member's nearest point lies inside the other
surface, and how far its farthest point lies outside it, both true distances in mm (a generator's
distance from the other cylinder's axis against its radius; a meridian circle's nearest and farthest
distances from a sphere's centre against its radius; a meridian's reach across a plane). So the
margin has one maximum, found where the two parts cross or at an end, and the contact is decided
from three distances: the margin at that maximum and at `c = ±1` (θ = 0 and π, where the family is
symmetric). A maximum within `LINEAR` of zero is an isolated touching point (or, where the member
itself lies on the other surface, a circle of contact: a ball rolling in a tube); a margin within
half `LINEAR` of zero at θ = 0 or π with the members crossing either side is a node where two
branches cross, which both end on; anywhere else the arcs of `θ` where members cross end where a
member turns tangent, found by bisecting the margin. Each arc's two meetings are then one closed
loop, or two branches between nodes. Where a node's margin is not exactly zero, the squared half
chord there is taken out where the chord is as short, and in proportion to the node's over the
chord's where it is longer, so both branches end on the node; that moves no point off the other
surface by more than the node's margin, and points away from the node by far less. A node's margin
is held to half `LINEAR` because the fit holds half `LINEAR` on top of it: rods of radius 2 and 2 +
1e-7 whose axes cross 1e-7 apart have a margin of 1e-7 at one end, and branches bent onto a node
there could not be fitted, so the result was traced instead and came back open. Past half `LINEAR`
the branches pass each other, as far apart as the chord the margin leaves (6e-4 mm there). For the
same reason two cylinders are one radius whose axes cross (the closed form with two ellipses) only
while the radii's and the axes' gaps sum to `LINEAR`: the ellipses, about the axes' midpoint at the
mean radius, then lie within half `LINEAR` of each. The points come from the closed form and are
fitted with points added until, a quarter, a half and three quarters along every span, the spline
lies within half `LINEAR` of both surfaces and of the closed form itself
(`IntersectionCurve::fit_sampled_seeded`), so nothing is refined numerically and no step can lose a
touching point. A span that strays is parted in two, four or eight by how far it strays, a
cubic's error falling as the fourth power of its spans, and every span is tried again on the
spline returned; a span's three points are the ones it is parted at, so no point of the closed
form is worked out twice. Within the tolerance of both surfaces is not
enough where they cross at a shallow angle: a fillet's tube against a copy's corner ball 1e-3 mm off
its axis crosses it at 1e-4 rad, and a spline 5.8e-5 mm along that band from the loop still lay
within `LINEAR` of both, so the vertices placed on the exact loop fell off it and the loop was
dropped. A plane against a cone ([`ssi/conic.rs`](../src/boolean/ssi/conic.rs)) meets each generator
once: an ellipse within the region's reach of the apex is built exactly from its vertices and the
point a quarter turn round; any other section is followed along the generators only as far as the
region reaches, on the cone's own nappe, the one its placement's origin lies on, which every cone
face does (the STEP reader turns over a file's cone whose face lies past the apex). A
plane all but parallel to a generator cuts an ellipse some 1e16 mm long, which is followed as a
parabola is. That part of a parabola or of a hyperbola's branch is stored as the curve it is, not
fitted: a conic arc with no point at infinity is the image of a circle's arc, so a rational
quadratic `Curve::Nurbs` of Bézier pieces with positive weights (Piegl and Tiller 1997, ch. 7).
The pieces are cut where one would bend through more than a quarter turn or run far out toward an
asymptote, and each is given the parameter in which it starts and ends at unit speed, so the
parameter is close to arc length and the first derivative is continuous across the double knots.
The fitted curve remains only where a middle weight would not be positive.

A spline surface swept by a steady motion keeps, in a private field, the curve it was swept from
and the motion ([`geom/steady.rs`](../src/geom/steady.rs)): a line and how fast a point turns about
it and moves along it, so a slide, a turn and a screw are one form. The surface at `(u, v)` is the
curve's point at `v` after `u` of the motion, within `SCREW_FIT`, all over its own domain. It
keeps them through a copy, through a rigid move, which moves the curve and the line with it, and
through a cut to a rectangle of its parameters, which keeps the parameters; a surface made from
its parts or from new control points starts without them, since it may be any shape. The sheets a
[sweep of a solid](#sweeps-of-solids) builds keep them, and so does each face a profile's segment
screws in [`sweep::helix`](../src/sweep/helix.rs), a thread's flank, where the segment's weights
are all positive; a segment that screws a strip of a cylinder needs neither.

Two such sheets of one motion, `A = M(u)·g(v)` and `B = M(u′)·h(v′)`, meet along the motion's orbits
([`ssi/orbits.rs`](../src/boolean/ssi/orbits.rs)): lines under a slide, circles under a turn,
helices under a screw (a surface the motions of a one-parameter group carry along itself;
Pottmann and Wallner 2001). A common point stays common as `u` and `u′` fall together, until one
reaches the start of its range, so every piece of the intersection is an orbit that starts on
`A`'s first curve where it meets `B`, or on `B`'s first curve where it meets `A`, and runs while
both parameters stay in range: no curve leaves either sheet's rectangle. Those starts are
curve–surface meetings against the stored surfaces, in closed form where the curve is a line. A
start may lie on the other sheet at several of its parameters, where that sheet passes through
itself, and each bounds its own run; they are found on the other sheet's curve, as its points at
the orbit's distance from the axis (under a slide, in a plane along it through the start) that the
motion carries to the start. A meeting the curve–surface search reports on a spline surface run
on past its edge has none, and starts nothing. Runs of one orbit that overlap are joined. A sheet
against a surface the motion carries onto itself whole (a cylinder about a screw's axis; a plane,
cylinder, cone, sphere, torus or revolution about a turn's; a plane, cylinder or extrusion along a
slide) meets it in the orbits from its own first curve, each as long as the sheet. An orbit is
emitted as a line, an arc of a circle, or a helix in the same Hermite pieces as a screw sheet's
columns, so it lies within `SCREW_FIT` of the true orbit and twice that of each sheet. Where the
first curve touches the other surface the orbit is one the surfaces touch along, and is marked
tangent. The pair is traced instead wherever that is not certain to be the whole intersection: the
motions differ by more than `SCREW_FIT` over the region; a first curve lies in the other surface
(two parts of one sheet); a touch has another meeting within 3e-4 mm, which is one touch found
twice; an orbit passes the axis; the other sheet's curve runs along the orbit's cylinder where it
crosses it, so its parameters there are not sharp; the sheet's own nearest point to a start is not
among the parameters found; or a turn comes round to its start.

Curve–surface meetings in closed form: a line against a plane, cylinder, sphere, cone or ring torus;
a circle or ellipse against a plane, a circle against a cylinder whose axis is square to it, and a
circle against a sphere (the sphere's section in the circle's plane). A line against a torus is the
zero set of the quartic `(|p|² + R² − r²)² − 4R²ρ²` along it, whose sign is the distance's: its
stationary points, each the one root of its derivative between the derivative's own stationary
points, split the line into stretches with one crossing at most, and every place the line can touch
the tube is one of them, so a touch is one point whatever the rounding. Other curve–surface pairs
are sampled (`numeric_hits`, below).

A torus or sphere face is boxed by its boundary and the points of its parameter window where a world
coordinate is stationary, not by the whole surface, so faces of one blend ring far apart are not
traced against each other.

A result with no material (a solid minus itself or a copy within tolerance, the common part of
disjoint solids) is `Ok` with an empty arena: no solids, faces or edges. It is valid, and callers
that need a solid check for one; measuring it gives zero. An empty intersection or difference is
refused where a vertex or an edge's middle of one input lies inside the other (outside it, for the
first input of a difference) farther than `LINEAR` from every surface of its faces: material is
there, and an intersection the boolean missed took it away.

Degenerate contact follows two rules:

- **A tangent curve is not a cut.** Where two surfaces touch along a curve, the pieces on its two
  sides lie on the same side of the other solid, so nothing is split there. If the result keeps both
  faces that touch, it would touch itself along that curve, which is not a manifold solid: the
  operation returns `Error::InvalidResult` saying so.
- **An edge used by more than two faces is an error.** Solids that touch along an edge only, such as
  two boxes meeting at a corner line, give a non-manifold edge and an error, never a solid.

A cut that ends on a pole splits the pole's degenerate edge where it arrives in the parameter plane.
Merging faces maps each one's parameters onto the first's (planes in any frame; surfaces of
revolution about one axis, whose angles and heights differ by sign and offset) and moves it by whole
periods to meet its neighbours: a shared edge at one place becomes interior, and one shared a period
apart becomes the merged face's seam.

## Blends

[`blend`](../src/blend.rs) chamfers and fillets edges picked by the names of their two faces: each
pair `[A, B]` takes every edge where a face standing for A meets one standing for B, and all pairs
go in one operation. A name stands for the faces that carry it and every face on the same surface as
one of them: an edge is defined by the surfaces its faces lie on. Each picked edge then takes every
edge that continues it tangentially between faces that continue its own, so a fillet takes the whole
tangent chain; such an edge's blend face carries the names of the pair that reached it. An edge two
pairs pick is blended once, its face named for both. An exact pair takes only the edges where faces
carrying the two names themselves meet, not every face on their surfaces: where a boss's underside
lies on the block's top plane, the pair of the boss's side and the underside picks the convex edge
below the boss and not the concave one where the side stands on the block. The forms that take edges
themselves pick only those: another edge between the same two faces (a bore's other rim) is left.

Where concave and convex picked edges meet at a vertex (other than at a saddle, below), the concave
ones are blended first, and the convex ones are picked again on the result, where they run on round
the concave rounds as tangent chains; each keeps its pair's index, which names its faces.

Each picked edge's blend is a piece: the region between the edge and the bevel line or round arc,
swept along the edge, a prism along a straight edge, a ring or a sector of one about a circular
edge. The section's sides follow the two faces: lines on planes, and on a cylinder along a straight
edge its circle, so a round between a plane and a bore is the circle tangent to both. A chamfer's
setback is the straight distance from the edge in the section, so on a circular side (a bore, or a
sphere about a circular edge's axis) it ends a chord that long from the edge, and the bevel is the
line or cone between the two ends. A chamfer by a distance and an angle sets back the distance on
one face and leaves it at the angle to it, turned back toward the edge, to where the line first
meets the other face's trace (a line or circle); it is the chamfer of two distances whose second is
the straight distance from the edge to that point. Each piece is swept as a tool solid, so its blend face, the contact curves where that face leaves the two faces
the edge joins, and its end sections are exact by construction.

The result is made from the pieces by surgery, the blend sewn in as a face (Jackson 1995, section
5), in [`blend/sew.rs`](../src/blend/sew.rs). Each piece's blend faces are taken whole from its
tool, facing out of the body. The two faces the edge joins give up the edge for the tool's contact
curves; every other edge at the edge's ends is cut back to where its own piece's contact curve ends
on it; a face a piece ends on takes the tool's end section across the corner it loses; an end
section two pieces share (a tangent joint, a mitre) is one edge of both; and a sphere corner is a
face bounded by its three pieces' end arcs. Every curve and point of the result is one the sweeps or
the tools' own trims made, or one already on the body. The sewn solid is checked, testing only face
pairs the blend made or changed for crossings, and blend faces on one surface that share an edge (a
circle split in two) are merged. On one real part's 17 blends the sewing took 0.19 s where the
tools' booleans took 1.44 s (wall time on a shared machine, the sweeps common to both left out).

Three ends need more than the edge's own corner, and each is still read from topology:

- **Past a short edge** ([`blend/run.rs`](../src/blend/run.rs)). Where a contact curve does not meet
  the next edge of its face within that edge's range, the face's loop is walked from the vertex,
  edge by edge, until an edge the contact curve crosses. Every edge before it lies wholly in the
  strip the blend takes off, so it goes with the blend, and so does each vertex where two of them
  meet; the faces across the walked edges bound the blend's end, so the piece is also trimmed by
  each of their surfaces (a plane, cylinder or cone), and their edges at the vanished vertices are
  cut back to where the trimmed tool's end sections meet them. The faces past the end must meet at
  convex edges, so the trims together keep what the body keeps; a concave piece does not walk. A
  chamfer of 1 on a notch's rim runs past a 0.8 chamfer of the flange below it onto the flange's
  cylinder, and the notch's 2 mm round runs past that chamfer's 0.28 mm ellipse edge onto its cone.
  A contact curve that ends exactly on the far vertex of a walked edge (an earlier chamfer's corner
  that a later chamfer's contact line runs through) leaves that vertex on the contact; the trims
  are the same, and the tools cut the piece.
- **Opening the face beyond.** Where the face a piece ends on turns more than a half turn at the
  vertex (a lug's edge meeting the wall it stands on, a chamfer into an inner corner), the end
  section lies past the face's corner, on the face's own surface. The line or circle edge between
  that face and the blend's own face runs on past the vertex to the end section, and the face takes
  the end section as new area.
- **At a crease.** Two convex rounds that each run through the other's wall into air, at a vertex of
  three faces sharing one, cross along a curve no face of the body holds. For sewing, each ends on
  the other's wall instead (the body lies inside both walls, so nothing it keeps is lost) and the
  two tools are joined into one by a union. The union finds the crease where the two blend surfaces
  cross, exactly as any boolean does (a cylinder and a torus along a half disc's rim, two cylinders
  at a tower's rim), and each piece keeps its own blend and side faces by marks set before the
  union; the crease is one edge of both pieces' faces. The joined tool's faces past the crease are
  the end sections where one round crosses the other's wall.
- **Along a face the rounds close.** A full round (every edge of a plate filleted at half its
  thickness, a rib's top at half its width) leaves nothing of the faces between its rounds. The
  balls at the two ends of an edge take the whole edge, so that edge has no piece, and the two
  corners are one ball: its patch is one face, bounded by the end arcs of the other pieces at both
  corners, on a sphere whose poles and seam lie on the great circle square to the edge, which the
  patch meets only between the edge's two faces. The rounds along a face's two long edges then run
  on one axis between the same two balls, so their contact curves on that face are one line: it
  is one edge of both rounds, the face's loop runs out along it and back and is dropped with the
  face, an edge cut back from both ends to one point is dropped, and the two rounds merge into one
  half cylinder. Balls within the tolerance of taking an edge whole would leave a strip, or an
  overlap, narrower than two points the kernel holds to be one. Each ball stands back from its
  corner in proportion to the radius, so the radius at which the balls' centres are one point is
  the one asked for times the edge's length over the two set-backs' sum; the fillet is built at
  that radius, and is the full round exactly. Past the tolerance a narrower radius keeps its strip,
  however thin, and a wider one is refused, since its ball does not fit.

Where the pieces do not close up that way the tools are used instead: each is cut away at a convex
edge and added at a concave one with a boolean, so the result is split, sewn, merged and checked as
every boolean's is, and the side faces hand the boolean the fact that they lie on the body's faces.
The tools are used for an edge at a piece's end that no contact curve or end section meets (a
concave piece running past a short edge), for a body with local tolerances, since sewing would cut
its tolerant edges back exactly while the booleans hold each to its own tolerance
([local tolerances in booleans](#local-tolerances-in-booleans)), and for any sewn solid the checker
refuses. Setting `SPRINGREP_BLEND_VERIFY` makes every blend run both ways, prints a
`springrep-blend:` line for each on stderr naming its feature, and fails a blend whose sewn result
differs from the tools' by more than 1e-9 of its area or volume or in its names. On the designs it
was run over, every sewn blend agrees to 3e-16, except the two whose ends are trimmed past a short
edge, where both ways intersect a cone or cylinder with a traced curve: 1.8e-12 in area and 5.3e-12
in volume (a notch's chamfer) and 4.9e-11 and 1.4e-10 (the notch's round).

Each piece ends where its edge does, in one of ten ways:

- **On a neighbour's section.** Pieces whose edges continue each other tangentially share their end
  section exactly. They are glued into one tool along it without a boolean: the two blend surfaces
  touch tangentially there, which intersection cannot resolve.
- **On a mitre.** Two pieces on one face whose other faces are planes mirrored across the plane
  between them, where that mirror also takes the shared face onto itself: a plane square to it (a
  box's top edges, a pocket's corners), or a sphere centred on it (a square post standing on the
  middle of a ball). The mirror then takes one blend surface onto the other, so they cross on the
  mirror plane, and each piece runs past the corner and is cut there. Two rounds between a plane
  and a sphere are tori about different axes, and their crease is the mirror plane's cut of
  either.
- **Crossing at a corner** ([`blend/strip.rs`](../src/blend/strip.rs)). Two pieces, both convex
  or both concave, sharing one face at a vertex of three faces, across a sharp edge between two
  plane walls that no plane mirrors (a post off the middle of a ball, or longer than it is wide; a
  post on a plate with one wall leaning more than the next; a square hole's rim; a bar's end). The
  two blend surfaces cross on a curve of their own, the crease, and both its ends are known
  outright. The far end is where the two contact curves on the shared face meet: both surfaces
  touch that face there, and so each other. At the near end one piece, the lower, comes to the
  edge between the walls nearer the corner. Where that edge bends the other way from the blends
  (the post, the hole), its contact curve runs on past the edge over its own wall's plane until
  it meets the other piece's surface; where the edge bends the same way (the bar), the higher
  piece has taken the edge that far and the lower's contact stops short of it, on the higher's
  surface. The crease starts at that point, and the higher piece ends on the lower's wall plane
  from the edge to it. Between its ends the crease is held on both surfaces, through points each
  solved on a line across the lower piece's strip laid flat, within the box the two ends are
  corners of: the two surfaces cross a second time through the far end, where they touch, and
  that second curve lies outside the box. Two rounds that come to the edge within 1.5e-7 of each
  other share one corner there, the point between the two, which lies on both surfaces since
  each is tangent to its wall; other blends do so only within the tolerance. On a post 5 by 3
  off the middle of a ball, rounded at 0.5 to 12, the solid holds the volume and area worked by
  hand to 8e-12, each corner's volume being worked two ways that agree to 1e-10.

Where every piece ends square on its own section, on one plane (another face's, or a mitre's
mirror) or at such a corner, no tool is made and nothing is cut: each piece's section is swept a
little past both ends, and its blend face is rebuilt alone, each contact curve stopped where the
end's curve meets it and the end's curves put in place of the sweep's own
([`blend/strip.rs`](../src/blend/strip.rs)). A plane's cut of a cylinder or a plane is the conic;
of a torus or a cone it is held on both surfaces as a crease is. These faces are sewn as the
tools' are. Two pieces whose edges are stretches of one line or circle parted at a vertex, as
a wall's foot is where the ball's seam meets it, are swept as one strip: its contact with the
wall is one edge, and its contact with the ball is parted where the seam ends on it, so no two
faces are left to merge. The strips' own sweeps are not checked, since only their blend faces
are kept and the sewn solid is checked whole. A square post's four mitred foot rounds on a plate
take 0.23 ms this way against 3.8 cut, and on a ball 3.1 ms against 26. About 1.3 ms of the
ball's is the fit of its four creases, each a plane's cut of a torus, which no closed form
gives; 1.5 ms is the sewn solid's check, which settles which way the solid faces by its volume
and searches each two faces that meet at a corner for a crossing, where a plate's faces are
settled by a ray: no closed form gives where two tori cross, so the ray's proof is not to be
had.

An edge of the body that leaves the joint between two pieces along their end section there, and
stops short of the contact curve, lies wholly under the blend: a ball's seam from the foot of a
wall to the pole, where the blend reaches over the pole. It goes with the blend, and so does the
pole. A point just past a pole lies half a turn of longitude from one just short of it, so no
hint places the contact curves' parameters there: each rebuilt loop's projected coedges are
continued one from the next round the loop instead. A blend's edge on a face that would cross
another edge of that face away from a vertex the two share, as a contact reaching over a seam
its own edge does not meet, is refused: the face's boundary would cross itself.
- **Through into air.** Two convex pieces sharing one face, at a vertex of three faces where each
  edge runs on past the vertex through the other's wall into air (a tower's rim where a square side
  meets a leaning cut), each run past the vertex and clear of that wall. Cutting both leaves the
  crease where their blend surfaces meet, as no material lies past either wall to lose; sewing joins
  them at the crease, above.
- **On another face.** A piece whose edge runs into exactly one other face ends square on it when
  that face is a plane square to the edge; otherwise it runs past and is trimmed by that face's
  plane, cylinder or cone, keeping the side its edge lies on. A straight convex piece whose edge
  runs into a corner of several other planes, every edge there convex (an octahedron's apex, or a
  corner earlier chamfers have cut), runs past it and is trimmed by each of those planes: the body
  there is the region inside all of them. The tools cut such a piece; chamfering a solid's edges
  one after another this way mitres each corner.
- **At a sphere.** Three or more fillets on every edge of a corner of planes, all convex or all
  concave, at any angles, stop short of the corner where a sphere closes it: the rolling ball's
  corner, the blend being the union of the balls that fit inside the body (Rossignac and Requicha
  1984). The ball `r` from all the planes on the side it rolls has its centre on each fillet's
  spine, since each spine runs `r` from both its planes; so each fillet stops on the section through
  the centre, and its end is a great circle arc from the ball's foot on one plane to its foot on the
  other. With three planes the centre is where the planes `r` in meet; with more it is their least
  squares point, which must lie `r` from every one (an octahedron's corner). The patch is the spherical polygon of the feet, and
  the material the corner loses (or, concave, gains) is the cell between the end sections, the
  centre and the corner, less the ball: for three planes of solid angle Ω, the cell less Ωr³/3. The
  sphere's poles lie on the great circle square to the direction from the centre to the corner, and
  its seam opposite it, both off the patch, since the corner stands `r` off every plane and so every
  foot lies within a quarter turn of that direction; only a corner so sharp a foot lies on that
  great circle is refused. On a prism whose normals are 68.2° apart the result holds the hand
  volume to 1.5e-16 and the hand area exactly. A concave corner is added with the tools.
- **Along a ridge** ([`blend/ridge.rs`](../src/blend/ridge.rs)). Four convex fillets at a corner
  of four planes that no one ball touches (a pyramid's apex on a rectangle) split it, as the union
  of the balls inside does: the lines `r` in from two opposite planes meet the other two planes `r`
  in at two points, each `r` from three of the planes and clear of the fourth, and a ball rolls
  between them along the line where the two opposite planes `r` in meet. Each fillet stops at the
  section through its own ball's centre, each ball closes its three fillets and the ridge's end as
  a sphere corner, and the ridge is a round between the two opposite planes along the line where
  they meet through the corner, an edge the body does not have, swept from one ball's section to
  the other's. The tools cut it.

- **At a flat corner** ([`blend/flat_corner.rs`](../src/blend/flat_corner.rs)). Three or more
  chamfers on every edge of a convex corner of planes mitre there, each trimmed by the corner's
  other planes as above, and the corner is then cut by the plane through the points where the two
  contact lines on each face cross: the plane a chamfer's own construction gives a vertex, its
  bevel lines' ends. Three points always span it; four or more must lie in it (an octahedron's
  apex), and a corner whose crossings do not (a pyramid's apex on a rectangle) is refused.

- **Over a smaller round** ([`blend/corner.rs`](../src/blend/corner.rs)). A fillet larger than a
  round it runs over, bent the same way (a plate's rim over its rounded corner, a pocket's floor
  past its corner rounds), cannot follow the round. Both pieces beside it share one corner centre,
  the fillet's radius from the plate and both walls, and stop on the planes through it square to
  each wall. The corner between them is a family of quarter ellipses in the half planes through the
  axis square to the plate: a rational biquadratic over a base of circles, exact, and tangent to
  both pieces. A piece along a cylindrical wall is a torus, which ends at the half plane through the
  cylinder's axis, where it and the corner share their tangent plane all along. The material between
  the corner, the walls and the plate is cut (or for a concave rim added) with the piece before it
  as one tool.
- **At a saddle** ([`blend/saddle.rs`](../src/blend/saddle.rs)). Four straight fillets of one radius
  at a vertex of three planes square to each other, alternating concave and convex round it (a boss
  overhanging a block's edge, at each end of its underside's edge), meet where two lines each turn
  from concave to convex, and no ball touches all four. In the vertex's frame each fillet stops a
  radius short of the vertex, square on the cube of side 2r about it. Inside the cube the filleted
  body is nine faces: the cube's six walls, the two convex fillets running on to their ends, and a
  patch bounded by those ends and the concave fillets' end arcs; none of the body's faces reach
  inside. A convex end is a rational quartic on its cylinder that leaves the concave arc square to
  it, along the cylinder's generator, so the patch's corners are right angles. The patch
  ([`blend/saddle_patch.rs`](../src/blend/saddle_patch.rs)) is a rational B-spline of degree 6 whose
  cross derivative along each side is built from the fillet's own homogeneous curve, the direction
  round its cylinder and its axis: whatever the coefficients, it lies in the fillet's tangent plane,
  so the patch meets each fillet G1 to rounding (1e-13 rad), and the corner twists are linear
  conditions on the coefficients. The coefficients and the inner control points fit a Gregory
  surface on the same sides (Gregory 1974). The cube is cut from the body and the nine faces glued
  in along its walls, which match edge for edge; the convex fillets' two pieces then merge.

On the tools' path, pieces that meet at a mitre or run through into air stay separate tools, cut one
after the other: the second cut meets the first's blend surface across it, where a union of the two
pieces would have to resolve their faces lying on the same walls. A run of tangent pieces is glued
into one tool; a piece that will not glue on, as the last piece of a run that closes on itself at a
mitre does (its far end's mitre cap lies on the first piece's), is cut as a tool of its own.

A round about a circular edge as tight as itself (a corner arc of the fillet's own radius) reaches
its axis at one point and turns into a sphere there. A sphere centred on a circular edge's axis
meets each section in a circle about its centre, and a torus about that axis in its tube's circle,
so a round where either meets a plane, cylinder, cone, sphere or torus about that axis (a
hemisphere's equator, a dome fused on a block, a torus section's rim) is the torus tangent to both,
and its side on the sphere or torus sweeps into a zone of that surface itself.

A convex blend about a whole circle whose contact on one face reaches the axis (a chamfer as wide as
a cylinder's radius on its end, a round of the cylinder's own radius) takes that face whole, and the
blend closes to a point on the axis: a cone's apex or a sphere's pole. The section would pinch to a
point there, so it runs out over the face instead, along the face's normal and across to the axis,
where only air lies past a convex edge, and the tool is cut from the body with a boolean.

A fillet along a tangent chain with an edge no constant section sweeps (an ellipse where a bore
meets a sloping face, a circle whose faces do not share its axis) is the surface a ball of the
fillet's radius sweeps rolling along the chain against the faces on both sides (Choi and Ju 1989)
([`blend/rolling.rs`](../src/blend/rolling.rs)). At each point the ball's centre lies in the plane
square to the chain, at the radius from both sides, found by Newton steps whose rows are the two
surface normals and the chain's tangent, from the centre of the ball in the wedge of the faces'
tangent planes, `(na + nb)·r/(1 + na·nb)` off the edge (or from a neighbouring section's centre,
and from the wedge's where that does not settle). A ball larger than the chain's bend (a cylinder
cut aslant and grown past the cut ellipse's least radius of curvature, whose rim the ball rolls
round) has no centre in some of those planes: they cross before they reach the centres' path. The
fillet is then made again with each section in the plane square to the centres' path through the
chain's point, the plane the ball's circle of contact lies in, the two contact normals' cross
product its normal; that row turns with the centre, so it is the plane distance's gradient, by
central differences. Where a side changes face at a joint between faces that
meet tangentially (a slot's flat side into its round end), the contact slides from one onto the
other: the foot is the nearest on that side's faces that lies in its face, by the face's own trimmed
domain (where none does, past an open chain's end, the ball rolls on the link's own surface). The
crossing is where the foot on the face before reaches the curve between the two, found exactly, and
that edge is cut back to it.

The blend is cut into segments at every chain vertex and every such crossing, so each is smooth.
Across, a segment is the exact circular arc between the contacts, a rational quadratic; along, its
rows are interpolated by a cubic, sampled twice as finely until every midway section and both
contact curves lie within a tenth of the linear tolerance of the true ones. The contact curves are
the surface's own edge rows, exact lines in its parameters. Neighbouring segments share their
section arc. A closed chain that crosses a face's seam (a tee's junction over the bar's seam) meets
that face's loop in more than one run of its edges, and each run takes the contact pieces between
the crossings at its two ends. A closed chain of several edges that a face's loop meets in one run
(an elliptic cylinder's rim drawn as two halves) takes them all, starting where the run does, which
need not be the chain's own start.

An open chain ends where one other face meets each end, a plane or any curved face. Past the end the
chain runs on straight along its tangent, so the ball keeps rolling between the same two surfaces (a
spline surface runs on past the edge of its parameters as its end spans' polynomials, the foot found
there by Gauss–Newton steps, so the contact stays smooth instead of stopping at that edge),
and the blend is trimmed by that face's surface: `s*(v)` is where each true section crosses it,
measured along the surface's normal, and the trim curve is fitted like the blend. A trim may cut
across several segments, leaving a crescent of one where it touches a section arc. The end face
loses the corner the blend takes: along its boundary on each side, its edges are dropped up to the
one the contact ends on, which is cut back; where the face turns past a half turn at the vertex, its
first edge (a line, circle or ellipse) runs on past the vertex instead, and the face gains the area.
Forced onto a circular rim, the rolling ball gives the exact torus's volume and area to 1e-9. A
design that picks several such chains has them filleted one after another.

What the region form cannot reach comes back as `Error::Unsupported`, never as a wrong solid:

- a straight edge between faces other than planes and cylinders along it;
- a circular edge between faces other than planes, cylinders, cones and spheres about its axis,
  unless it is filleted as a rolling ball;
- a rolling-ball fillet through a chain vertex another face meets, or whose contact runs past more
  than a run of its end face's edges; and a fillet larger than a convex corner round in its chain
  (the exact path reaches past the corner's axis);
- blends meeting at a vertex other than as a tangent chain, a mirrored mitre, a convex corner they
  run through, a corner they cross at, a sphere corner, a ridge or a flat corner (the error names
  the vertex and the faces), and an edge ending where several other faces meet it other than at a
  convex corner of planes, or where a face other than a plane, cylinder or cone meets it;
- blends that cross at a corner whose faces cannot be sewn (the error names the corner and why):
  beside a piece that must be cut; blends so wide that their contacts on the shared face pass
  without meeting, so the crease has no far end; where a blend reaches over an edge of the face
  it stands on that its own edge does not meet, as a ball's seam beside a post (the error names
  the point); and the largest rounds of a hole's rim;
- a blend wider than a face it sits on, including one whose strip would cross onto the next face
  along its edge, such as where its edge runs on, tangent, as a curve other than a line or circle.

A pair that picks nothing, a size that is not positive, and an edge between tangent faces are
`Error::InvalidInput`.

A result can still fail as an explicit error where a boolean meets tangent contact it cannot
resolve: a round whose surface touches a body face at a single point beside a tangent joint (a
fillet round the root of a boss whose sides are rounded at their corners), or that meets a cone
along the circle where the cone joins the face the round is tangent to.

### Blends past the rolling ball

[`blend/shaped.rs`](../src/blend/shaped.rs) makes the blends a rolling ball does not: a conic or
curvature-continuous section along a straight edge between planes, a fillet between faces that share
no edge, a full round, a ball that grows along its edge, and corners closed by a patch of their own.
A section blend is cut as a prism, the region between the edge and the section, run along the edge;
the curvature-continuous section is a rational quintic, exact. Two corners are not cut but put
together from their faces, each face given by its surface and its bounding edges and the solid
assembled as a file's is ([`step::solid_of`](../src/step/topology.rs)), so no boolean meets the
tangent contact all round them:

- **A set-back corner between curvature-continuous fillets**
  ([`edit/corner.rs`](../src/blend/shaped/corner.rs)). In the distances `u` in from the corner's
  three faces, with `φ` falling from 1 on a face to 0 at the fillet's reach and `B` falling from 1
  to 0 at the setback, the corner is `φ(u₀) + φ(u₁) + φ(u₂) + B(u₀)B(u₁)B(u₂) ≤ 1`: the fillets and
  faces exactly past the setback, and one surface of continuous curvature within it. It is nine
  quintic B-spline faces, cut where `φ` or `B` changes form and fitted to within 1e-8 of the radius;
  the faces along the fillets are rational over the section's denominator, so each ends on its
  fillet's section itself.
- **Three fillets of unequal radii** ([`edit/uneven.rs`](../src/blend/shaped/uneven.rs)). The corner
  is the ball of the largest radius touching all three faces. Along each edge the ball rolls on a
  travel of that radius, shrinking evenly to the edge's own radius, and from there is the fillet's
  ball: its centre runs on a line, so it sweeps a cone, or a cylinder where the radius is the
  largest. Every face is a plane, a cylinder, a cone or a sphere, each tangent to the next. Fillets
  of unequal radii that meet any other way are mitred, the solid less each one's region.

Both first cut the corner away as far as each fillet's start and fillet each edge from inside the
cut, so no fillet has to be cut against another.

## Moving faces

Shells, drafts and thickens move faces and keep the solid's topology where they can
([`tweak.rs`](../src/edit/tweak.rs)): a thicken moves only the faces it names, each by its distance
along its outward normal. Each moved face goes onto a new surface. Each vertex goes where its faces'
new surfaces meet, found by Gauss–Newton steps along their normals from the old point: the least
move where they leave it free (a seam's end lies on only two faces), their one meeting where three
or more pin it. Faces tangent at a vertex (a round and the wall it runs into) stay tangent when both
move one distance, so their equal rows count once. Where only one of two tangent faces moves (a
capsule's side moved in, its dome kept), the steps turn back and forth along their one normal; the
point is then the nearest on the curves the two surfaces meet in, those on a face that stays where
it is first (the moved side's cylinder meets the dome's sphere in two circles, one on the dome and
one below it, equally far). Each edge becomes its two faces' new meeting
between its new vertices: a line where the line through them lies on both, else the branch of their
intersection (closed form where one exists, traced otherwise) through the old edge's middle carried
onto both, trimmed at the vertices and laid the way the old edge ran. Faces moved by a rigid motion
carry the middle by the motion first: settled from where it was, the middle of a round's arc moved
further than its radius lands on the far side of the circle, and the edge would take the rest of it. A seam keeps its parameter
line on the new surface, through its moved ends. Faces whose loops were rebuilt take parameter
curves projected along each edge in turn, so a loop round a period keeps one image of its
parameters. Planar moves keep the older all-planes path, whose edges are lines by construction.

A face's offset is the same kind of surface, exactly: a plane moved along its normal, a cylinder's,
sphere's or torus's radius grown, a cone's radius grown by `d·cos a` as its origin moves `−d·sin a`
along the axis. A shell moves every face by its thickness (the open faces the other way, through the
wall), so a cylinder, sphere or torus is hollowed with its wall offset exactly. A draft tilts a
plane about its line on the neutral plane, and turns a cylinder along the pull into the cone through
its circle there, narrowing a boss along the pull and widening a bore. A draft from a parting line
(`edit::draft_parted`) splits the solid at the parting plane, drafts each side's pieces of the
faces away from it (pulled along the plane's normal on the side it points to, against it on the
other), and joins the sides: both keep the parting line, so the cut faces between them coincide
and the union removes them, and each piece keeps its face's names.

A surface with no exact offset of its own kind gets a certified one: a spline within a tenth of the
linear tolerance of the true offset `p + d·n` at the middle of every span, refined until it is or
refused. An extrusion's offset is the extrusion of its curve's offset and a revolution's the
revolution of its meridian's, so only that curve is fitted: a C¹ cubic through the true offset and
its derivative at nodes that include the curve's own knots (where the offset's curvature jumps, so a
C² interpolant would converge only as h²), with a double knot at each node. Any other surface
(ruled, a spline) is fitted as a bicubic through a grid over the face's parameter box. The box is
read from the face's pcurves, poles among them, so a face that wraps a whole period (a spline ball,
closed on its seam) spans that period, where projecting its edges' points reads only a sliver round
the seam; along a period it does not wrap the box may run past the period's ends. The surface's
own knots cut the box into blocks, each fitted on its own and joined to the next at a triple knot
along the curve both share: across a rational sphere's double knots the offset's curvature jumps,
and a grid through them converges only as h², one cut there as h⁴. Each block's end steps are
halved, where a not-a-knot interpolant errs most. A spline sphere octant offset this way lies on the
grown sphere to 2e-8 wherever it is sampled, and the whole biquadratic ball grown by 2 lies on the
ball of radius 7 to 2e-8.

### Faces that close on the way

Where a move cannot keep the topology because a face closes before the faces get where they go
(walls drafted until they meet over the top, a step raised level with the top so the riser between
them closes, a slot's wall moved onto the other wall, a round offset past its radius), the move is
followed there ([`through.rs`](../src/edit/tweak/through.rs)). Each moved face runs along a path
from its old surface to its new one: a plane through the pencil of planes the two span (so a draft
keeps its hinge and an offset stays parallel), an offset surface by its distance, a cylinder drafted
into a cone by its half angle. The furthest the solid gets with its topology kept is found by
bisection, to a bracket of 1e-15 of the move: a share the move reaches beside one it does not. Where
edges shrink between the last two shares reached, the bracket is cut not at its middle but where the
first of them, run on in a straight line, comes down to the shortest length the checker passes, and
then just to that guess's other side, by a few times what the guesses are still moving by; a guess
only chooses where the bracket is cut, so what it closes on is the same.

Many closings need no bracket. An edge that shrinks in proportion to the move (a riser between two
tops brought level, a slot's floor, a circle about an axis) closes where two shares the move reaches
say it does, so a guess is aimed where the edge is one and a half times the shortest length the
checker passes, or, for an arc, where its radius is 1.2 times the tolerance if that leaves it
longer (an arc is rebuilt down to a radius of about the tolerance), which is just short of where a
bracket ends. The search stops at a guess the move reaches only where
what is made one there is shown to be what is made one at the bracket's end:

- an edge that shrank since the share before is no longer than twice the tolerance there (an arc
  of radius up to one and a half times the tolerance, and never over nine tenths of what is made
  one) and, run on in a straight line, closes before a share the move is known not
  to reach, so the bracket's end lies between the guess and that closing;
- every edge, and every pair of vertices within twice the distance that is made one, is within that
  distance at both ends of that stretch or beyond it at both, by a tenth either way, its length
  run on in a straight line: edges that close together at rates so uneven that the fastest is
  still too long at the guess fail this;
- no vertex lies within that distance, or twice it, of the inside of an edge it is not on, where
  the edge would be split first.

Anything else closes its bracket as before: a guess the move does not reach, a search in which no edge shrinks to its closing (a
torus grown until its hole closes, a ball's cap shrinking through a plane), and an edge already
shorter than what is made one that is not on its way to closing, which stops nothing.

There vertices closer than ten times the tolerance are made one, edges that then run
together one, and edges shorter than that dropped; a face left a strip along a single edge (or an
annulus whose two rims have run together), or lying on another face to face, is taken out, and a
face left running along an edge and back loses it: the edge goes where nothing else uses it, and is
a ridge where two other faces do (an L's top closed between the walls of both its arms, drafted
until they meet). A vertex carried within that distance of an edge it did not touch splits the edge
there first (a T's bar closing, its top edge coming down onto the stem's corners). A closed edge shrunk round
a cone's apex or a sphere's pole becomes that pole on the face whose surface closes there, and a
face left only poles (a ball's cap shrunk onto its pole) is taken out. A sphere shrinking through a
plane reads as touching it while the circle they meet in is still about √(2rε) across, ε the
tolerance, so where nothing is shorter than ten times the tolerance, a closed edge that has shrunk
to under a thousandth of its length since the last event counts as closed. Where every face closes,
the solid is empty, unless a face with more than that distance of its move left faces another
head on, their outward normals opposite: the two would pass through each other (a box's top
pushed through its bottom, a shell's walls thicker than half the part), and the move is refused
as invalid input, as `edit::shell` then refuses its shell. A face that closes onto its own axis or
apex (a frustum's side shrunk past its apex) passes no other face, and leaves the solid empty. The
move goes on from there with every vertex settled again where its faces
meet, a pole put on its moved surface's pole, and faces left on one surface are merged at the end.
Where nothing closes but the move cannot go on, a corner of four faces is being pulled apart (the
T's corners where the bar's ridge meets the stem's walls): it becomes two corners of three faces,
joined by a new edge on the two faces they share, as a straight skeleton's vertex splits
(Aichholzer et al. 1995). Each such corner is tried kept or split each of its four ways (either
pair of opposite faces, the edge laid either way), up to four corners at once, fewest splits first,
and the first that moves on from a little past the event is kept. An event that needs more ends
the move with the error it gave.

### Growing a solid

`edit::offset` grows a solid by a distance `d`: every point within `d` of it. It first moves every
face out and rounds the convex edges at `d`. A hole's cylindrical wall no wider than `d` closes on
the way; over each end of it the rim leaves a dimple, the part of the torus of radius `d` about the
rim inside it, taken from the grown face there. A hole exactly `d` wide (to the tolerance) closes
onto its axis as the move ends, and its dimple comes down to a cusp on the axis. Where that cannot
keep the topology (a face beside a convex edge closes on the way, as a groove's walls cross; a
cone's apex is a convex point; a torus's tube grows past its axis), it builds the set itself
([`near.rs`](../src/edit/near.rs)). A point outside the solid within `d` of it has its foot, its
nearest point, inside a face, on an edge or at a vertex; so the grown solid is the union of the
solid, each face's slab (the points the face sweeps moving out along its normal), a tube about each
convex edge and a ball about each corner (Rossignac and Requicha 1986). A concave edge's foot never
lies off its faces, so it needs no tube; a vertex where the convex edges only run on into each other
tangentially needs no ball. Each piece is exact:

- a plane face's slab is the face extruded along its normal;
- a cylinder's, cone's, sphere's or torus's face, bounded by parallels about its axis, sweeps a
  strip or an annular sector in the meridian half plane, cut into pieces of at most a half turn and
  clipped at the axis, and turned a whole turn. A point the sweep carries past the axis lies nearer
  the same face on the other side of it, so the clip loses nothing. A face bounded by two meridians
  as well is turned from one to the other, the way round its inside lies; there a point carried past
  the axis would lie outside the turn, so such a sweep is refused;
- an extrusion face round a whole closed plane curve, square to its plane, sweeps the region
  between the curve and its parallel at `d`, extruded. Where `d` passes the curve's radius of
  curvature the parallel loops back in a swallowtail; it is cut where it crosses itself, and the
  pieces whose points lie nearer than `d` to the curve are dropped (Farouki and Neff 1990)
  ([`parallel.rs`](../src/edit/parallel.rs)). A shrink past an ellipse's curvature, built as the
  solid less its complement grown, keeps the parallel's arcs between their crossings on the long
  axis, which meet there in a sharp edge. A face on a stretch of an open plane curve, bounded by two
  copies of the stretch and the lines its ends sweep, sweeps the region between the stretch, its
  parallel and the normals at its ends; one that reaches the curve's radius of curvature is
  refused;
- a line edge's tube is a cylinder, and a circle's the disk of radius `d` turned along it, clipped
  at its axis where `d` passes the circle's radius;
- a corner's ball is a sphere.

The booleans resolve where the pieces meet, tangent along the slabs' sides and the tubes' contact
lines. A face of another kind, or an edge other than a line or a circle, is refused with the error
the move gave.

An inward shell whose moved walls cannot be followed (walls that close beside a concave edge, as in
a neck narrower than twice the wall; a bore's moved wall that runs into another wall) is built the
same way: the cavity is the solid less the union of its walls' slabs moved in, a tube about each
concave edge and a ball at each corner of them, and the shell is the solid less the cavity. That is
the sharp shell's cavity where the two agree, so it is used only where every concave edge has a
face beside it that closes (the sharp shell has no corner to meet there, and the points within the
wall of the edge round it) and every open face meets the walls at a right angle or less, so the
walls' slabs reach it as the moved walls would.

### Shrinking a solid, and round shells

`edit::offset` with a negative `d` shrinks a solid: the points of it at least `|d|` from its
outside ([`rounded.rs`](../src/edit/rounded.rs)). It moves every face in by `|d|` and rounds each
concave edge at `|d|`, the round being the cylinder of radius `|d|` about the edge itself; the
convex edges stay sharp. Where a face beside a concave edge closes on the way (a neck narrower than
`2|d|`), or a cone's apex is a concave point, it builds the set: a point inside the solid within
`|d|` of its outside has its foot inside a face, on a concave edge or at a concave corner (a convex
edge's inward normals all lie nearer one of its faces), so the shrunk solid is the solid less the
union of each face's slab moved in, a tube about each concave edge and a ball at each corner of
them (Rossignac and Requicha 1986). Every piece holds only points within `|d|` of the outside, so
a piece more than needed loses nothing. Where the shrink cuts the solid in two, the difference
keeps every piece as a solid of its own.

`edit::shell_rounded` hollows a solid as `edit::shell` does, its wall every point within the
thickness of the faces other than the open ones. An inward wall's cavity is the solid shrunk with
the open faces left where they are; an outward wall's outside is the solid grown that way, each
convex edge between two moved faces rounded and each corner a ball's patch, less the solid. Where a
wall meets an opening, the moved face runs on to the opening's surface, the rim, as in the sharp
shell. An inward wall that cannot be followed is the solid less the union of the walls' slabs, the
tubes about the concave edges beside a wall and the balls at their corners, where each open face
meets the walls at a right angle or less; an outward wall that cannot be followed is refused.

Faces moved by `offset_faces` that run into faces they did not meet (two holes widened until they
merge) cannot be followed either. The moved solid is then the solid with the slabs of the faces
moved out added and those of the faces moved in taken away, where that is the moved solid: every
face beside one that stays meets it square, so the slab's side lies on the other face, and two
faces moved by one distance meet where the slabs leave a sharp corner (moved in at a convex edge,
out at a concave one) or tangentially. The slab's side lies on the face that stays only as far as
that face reaches: where the edge carried along the moved face's normal leaves it (a step raised
past the top its riser ends at), the face has passed the face it meets, and the move is refused
as invalid input.

### Faces taken out

`delete_faces` takes faces out and closes the gap with the faces round it, each run on along its own
surface ([`delete.rs`](../src/edit/delete.rs)). Each neighbour either holds the deleted faces in its
material or has them standing out of it. The general way builds each neighbour's material as a
solid, keeps what is in every neighbour of the first kind and in any of the second, and adds to the
solid, or cuts from it, the pieces of the difference that touch the deleted faces, where those are
bounded by the deleted faces and the neighbours' surfaces alone.

Most gaps close in a way that needs none of that ([`delete/direct.rs`](../src/edit/delete/direct.rs)).
The deleted faces fall into pieces joined edge to edge, and round a piece each neighbour's loop is
cut where it ran along it. Where at most one edge that stays reaches each vertex of the piece's
rim, the loops close by one of:

- whole inner loops of planar faces dropped: a bore's two ends, a boss's foot;
- two rings, each a whole loop or the rim of a face with a seam: on one surface the two faces become
  one (a groove round a shaft), and on two surfaces that meet in a circle each runs on to it (a
  bevelled or rounded rim);
- one chain of two, three or four runs along planar neighbours: two straight edges that face each
  other across the piece become one (a notch, a nick in an edge), three planes meet in a new corner
  (a cut corner), and four either pair off into such edges or meet in two corners joined by a new
  edge (a chamfer, a round between two walls). Where both pairings of four give corners, both are
  tried.

The solid is rebuilt by a splice ([`local.rs`](../src/edit/local.rs)): edges run on along their own
curves to the new vertices, new edges fill what is left, loops are traced through the cuts, faces a
loop passes through become one, and every face whose loops are untouched is copied with its memo.

The result is taken only where it is shown to be the general way's. The region between a piece and
the faces run on is bounded by the piece and by parts of the neighbours' surfaces, so it is a whole
piece of what the general way adds or cuts exactly when it lies on the right side of every
neighbour's surface: in every surface some piece lies within, and, for a piece that stands out of a
neighbour, out of every surface some piece stands out of (a gap filled must lie in one of those,
where there are any). Each side is decided in closed form from the region's rims: its edges, the
stretches edges ran on by, and the new edges. Against a plane the distance is linear, so over a
planar, cylindrical or conical face it is greatest and least on the face's edges, and over a sphere
or torus turned about the plane's normal on its meridian; against a cylinder, cone, sphere or torus
every face of the region must be turned about that surface's axis, and the distance is then a
function on the half plane through the axis, where each face is a curve the rims hold. A rim that
touches a surface it does not lie on, away from its ends, fails the test. The region must also lie
in the box the general way works in. Last, the checker does not pair faces that share an edge, nor
find a shell inside the region: a face that keeps its loops and whose box meets the region must
belong to a changed face's shell and have the region strictly to one side of its surface. A result
that fails any of this, or the checker, or a gap of any other shape, goes the general way, and so
does a choice between two closings that both pass.

### Walls and slabs

A shell's wall and the slab thickened faces sweep are each a solid less the same solid with some of
its faces moved: `shell` takes the cavity (every face but the open ones moved in) from the solid,
or the solid from its outside moved out, and `thicken` the solid from itself with the named faces
moved out, or itself with them moved in from the solid. The two solids share every face that did
not move, but for the strip such a face gained or lost beside a moved one, so their difference is
built from their own faces, with no intersection ([`between.rs`](../src/edit/between.rs)):

- each moved face of the larger solid, as it is;
- each moved face of the smaller, turned to face the other way, its loops run backwards;
- for each face that did not move, the strips the larger has of it and the smaller has not. Where
  a whole loop runs along moved faces (a shell's opening, a cylinder's end beside its thickened
  side), the strip lies between the larger's loop and the smaller's run backwards: the larger's is
  the outer one where it was the face's outer loop, and the hole where it was a hole, the strip
  then lying round it. Elsewhere each run of edges along moved faces bounds a strip of its own,
  closed at both ends by the stretch of the edge that leads up to the run, between the smaller's
  corner on it and the larger's. That stretch is new, on the larger's own curve, and one edge of
  the strips of the two faces it lies between.

The faces must be paired first. Where the move kept the solid's topology, the two solids have the
same faces, loops and coedges in the same order, and every face that did not move is paired with
itself: a corner no moved face is at must be the same point in both. Where a wall of a shell closed
on the way in (rounds thinner than the wall), an open face is paired with the one face of the
cavity on its surface, facing its way and carrying its names, and only where each has a single loop
that runs along moved faces all the way round; an open face the cavity has no face on is kept
whole (a frustum's top above the apex its moved side closes on).

The result is taken only where it is the difference. The faces laid are the two boundaries with
what they share left out, the smaller's turned over, so a solid they bound is the points in exactly
one of the two solids; that is the larger less the smaller where the smaller lies in the larger.
Each strip is laid as a piece the larger has: the smaller's corner must lie on the larger's edge
short of its corner by more than ten times the tolerance, and a strip or a shell laid where the
smaller reached out of the larger is turned inside out, which the checker refuses (a loop that runs
the wrong way, an inner loop outside its face, a void outside its solid, a volume that is not
positive), as it does two faces that meet away from their edges. The faces must hang together as
one shell, or, where no face is shared (a ball shelled, a torus's whole face thickened), as the
larger's faces with the smaller's a void in them. A pairing that fails, a result of more pieces, or
one the checker refuses goes to the boolean difference, as it always did: so does an opening beside
a concave edge, where the wall moves in under itself and the cavity's face on the opening's plane
reaches past the opening as the end of that wall.

## Sheets

A sheet is faces with no solid behind them: one or more shells, each in no solid
([`sheet`](../src/sheet.rs)). `sheet::from_faces` copies faces of a solid, each keeping its surface,
its outward normal and its names; faces that share an edge in the solid share it in the sheet, and
an edge only one of them uses is on the sheet's boundary. Each piece of faces joined edge to edge is
one shell. A sheet's front is the side its faces' normals point to; `sheet::flipped` turns it over,
reversing each face and each loop.

`check::check_sheet` holds a sheet to the rules of [Validity](#validity) with two changes: every
non-degenerate edge is used once, on the boundary, or twice, by coedges that run opposite ways in
faces of one shell, never more; and every shell is in no solid. The volume, inside-out and nesting
tests have no solid to test.

`sheet::thicken(sheet, front, back, corners, feature)` builds the solid the sheet sweeps moving
from `back` behind it to `front` in front of it along its normals. Each side is built alone, and
two sides are united along the sheet.

- **Rounded.** One side thickened by `d` is the points within `d` of the sheet on that side, the
  union of a piece for each face, crease and corner (Rossignac and Requicha 1986)
  ([`sheet/near.rs`](../src/sheet/near.rs)): each face's slab ([Growing a solid](#growing-a-solid));
  at each crease the side turns convex at, the sector of the tube about the crease between its two
  faces' normals (a straight crease's extruded, a circle's turned, each while the normals keep
  their angle to the crease along it); and at each corner where two or more such creases meet other
  than running on tangentially, the part of the ball about it inside the cone its faces' normals
  span, the ball cut by each plane through two normals that has the rest on one side. Where the side
  turns concave the slabs overlap, and two faces that meet tangentially share a slab's side. A
  corner on the boundary or beside a crease turned the other way is refused.
- **Sharp.** The moved faces extended until they meet, as `edit::thicken` moves a solid's faces:
  the rounded thickening of the sheet's back by `d` is a solid whose faces beside the sheet's
  boundary lie along the sheet's normals there, so the sheet's faces in it moved `d` out
  ([Moving faces](#moving-faces)), less it, is the sheet thickened to its front with every crease
  sharp.

The faces on the sheet and those moved from it keep the sheet's faces' names; the walls at its
boundary are `<feature>.side` and the rounds `<feature>.round`.

## Moving whole solids

`Brep::transform` moves a solid rigidly and leaves every parameter alone. `transform_similar`
([`topo/similar.rs`](../src/topo/similar.rs)) takes any similarity: a uniform scale times a rotation
or a mirror, then a shift. A stretch or shear is refused with `Error::Unsupported`, because a
cylinder or sphere it moves is no longer one.

- **Mirror.** Each curve and surface is rebuilt as its mirror image in a right-handed frame. A
  plane, circle or ellipse keeps its parameters, with z flipped. A cylinder, cone, sphere, torus or
  surface of revolution runs its angle backwards, so u becomes period − u. A face whose parameters
  are kept turns its normal over. Every loop is reversed, and the pcurves map the old uv through the
  same change.
- **Scale.** Each curve's parameter scales by the factor its geometry needs: a line's by s (its
  direction stays a unit vector), a conic's and a spline's not at all. Each surface's parameters
  scale the same way: a plane's both, a cylinder's or cone's height, an extrusion's sweep. Edge
  ranges, affine and spline pcurves and the local tolerances scale with them.

## Sweeps of solids

[`sweep::solid`](../src/sweep/solid.rs) builds everything a solid passes through as it moves
rigidly: `solid(body: &Brep, motion: &Motion, feature: &str) -> Result<Brep>`, the union over
`τ ∈ [0, 1]` of `M(τ)·body`. The body is given where it is drawn, and the motion places it:

- `Motion::Slide(d)`: carried straight, `M(τ)` the shift by `τ·d`.
- `Motion::Turn { origin, axis, from, to }`: turned right-handed about the line, from angle `from`
  to angle `to`.
- `Motion::Screw { screw, from, to }`: carried along a [`Helix`](../src/sweep/helix.rs) from turn
  `from` to turn `to`.
- `Motion::Along { path, range, up }`: carried by a frame riding a curve, its origin on the curve,
  its x along the tangent and its z as near `up` as square to the tangent allows; the body is drawn
  in the frame at the range's start.
- `Motion::Then(a, b)`: two motions over one parameter, `M(τ) = B(τ)·A(τ)`: `a` places the body,
  and `b` carries the whole of it.

### The envelope

In the body's own frame each point `p` moves at `V(τ, p) = M(τ)⁻¹·M′(τ)·p`. Where `V·n > 0`,
with `n` the outward normal, the body's surface moves out into space it did not hold (an ingress
point in the terms of Blackmore, Leu and Wang 1997); where `V·n < 0` it moves in; where `V·n = 0`
it grazes. A point the moving body covers at some `τ` either lies in the body at `τ = 0` or the
body reached it through the region `W(τ)` of its surface where `V·n > 0`. So the swept solid is the
body at the start together with what `W` sweeps, and where nothing overlaps its boundary is: the
body's faces outside `W(0)` at the start, the faces of `W(1)` at the end, and the surface that each
curve of the boundary of `W` sweeps between them (Wang and Wang 1986; Martin and Stephenson 1990;
Weld and Leu 1990). The boundary of `W` has two kinds of curve:

- **Characteristic curves**, where `V·n = 0` on a smooth face: what they sweep is the envelope,
  tangent to the moving body all along them. It is the solution set the Jacobian of the map
  `(τ, u, v) ↦ M(τ)·S(u, v)` drops rank on (Abdel-Malek and Yeh 1997).
- **Sharp edges** between a face piece in `W` and one outside it: what they sweep is a ruled or
  swept sheet of the edge itself.

For a slide, a turn or a screw, `V` is the same field at every `τ` (each is a one-parameter group
of motions), so `W` and its curves are fixed on the body and the solid is a product: the boundary
of `W` times `[0, 1]`, capped by the two pieces of the body. For the other motions `W(τ)` slides
over the body as `τ` runs.

### Finding the curves

On each face `g(u, v) = n(u, v)·V(S(u, v))`, the normal taken outward. A face on which `g`
vanishes at every sample (a plane along a slide, a cylinder turned about its own axis) is still:
it moves within its own surface and sweeps nothing new; it is not in `W`. Otherwise the roots of
`g` along each edge of the face, bracketed by samples and pinned by bisection, are where its curves
meet the face's boundary; a root on a seam is two places on the boundary, a period apart.

A rigid velocity is `ω × p + v`, so on most surfaces `g` is simple in one parameter
([`contour.rs`](../src/sweep/solid/contour.rs)): affine in both on a plane; affine in `v` on a
cylinder, cone or extrusion, whose normal is the same all along `v`, so the curve is
`v = −n(u)/d(u)`; and `a(u) cos v + b(u) sin v` on a sphere or torus, so it is the two angles
`atan2(b, a) ± π/2` (on a sphere one of them, a great circle). Where both coefficients vanish at
one `u`, the whole line `u = const` is a curve (a ball's meridian under a turn about its polar
axis, a generator of a cylinder under a slide). Each solution's points are exact. A branch is
clipped to the face by the face's winding number along it, each run inside the face ends at the
root the branch passes through there, and runs that cross inside the face (a whole line across
another branch) end at a vertex there. The winding number is read at even steps along the branch
and, closing in from both sides, beside each root the branch passes through, so a run is found
however little of the branch's parameter it spans: far from a screw's axis a leaning cylinder
grazes nearly along two rulings, each curve crossing the face within a thousandth of a radian. A
root that no curve ends at is refused: the face would be taken whole for one side. A run stops
where the face's boundary is first within the tolerance; a branch that runs into an apex beside a
seam is that near the seam well short of the apex, so an end with no root, where the run stopped
away from every pole, is first looked for at the face's poles, along the branch a little past
where the run stopped, and is the pole where the branch passes through it. The run's points are fitted with a cubic B-spline at
chord-length parameters, points added in each span that strays from them by more than a tenth of
`LINEAR`; points on one line or one circle within rounding are that line or circle exactly.

On other surfaces (revolutions of other curves, B-splines, ruled surfaces) each curve is traced
across the face from a root, a predictor step along the curve and Newton steps on `g` back onto it,
until it reaches another root, and the traced points are fitted the same way. A curve that comes
round to its start through a seam is a periodic spline.

Under a turn, a point of an edge on the axis is held still. Each edge is cut where it touches the
axis, so each such point is a vertex. Where the axis pierces a face, `V` vanishes and changes sign
across the point, so `g` vanishes there and a characteristic curve passes through it; the curve is
cut there, and its two halves sweep apart from a pole at the point.

A turn of a whole turn or more ends where it starts, so its two ends would be one face; it is built
as the union of the sweeps of its two halves. Turns about one line compose, so the second half's
sweep is the first's turned half a turn, and only the first is swept.

The faces are split along the curves by the booleans' arrangement
([`arrange.rs`](../src/boolean/arrange.rs)), each piece is in `W` or not by the sign of `g` at a
point inside it, and the pieces are joined by [`sew`](../src/boolean/sew.rs) into the result, which
merges neighbours on one surface and runs the checker.

### The swept surfaces

What a curve of the boundary of `W` sweeps is a face bounded by the curve at the start, the curve
at the end, and the paths of its two ends. It is exact where the motion has a closed form:

- a slide sweeps an extrusion of the curve: a plane where the curve is a line or lies in a plane
  that holds the slide (a drum's rim pushed across its axis, whose end face slides within its own
  plane: an extrusion there would turn singular where the rim runs along the slide), a cylinder
  where it is a circle square to the slide;
- a turn sweeps a plane where the curve lies square to the axis; a cylinder, cone, sphere or torus
  where the curve is a line or a circle in a plane through the axis; the rational B-spline whose
  columns follow its control points' circles, within `SCREW_FIT`, where it is a line askew to the
  axis; else the revolution of the curve;
- a screw sweeps the rational B-spline whose columns follow the helices of the curve's control
  points, as [`helix`](../src/sweep/helix.rs) builds a screwed face, within its `SCREW_FIT`; a
  line along the screw's axis sweeps the exact cylinder.

A B-spline sheet keeps the curve and the motion that swept it, so that two sheets of one sweep, and
a sheet and a cone or cylinder of the same sweep, meet in closed form along the motion's orbits
(see [Analytic contact](#analytic-contact)).

Under a turn or a screw, a curve that stays at one distance from the axis sweeps a strip of the
cylinder about the axis at that distance. The path of a vertex under a screw is a helix, cut at each
whole turn, so that no edge of the result runs more than once round the axis.

Every piece of one curve (an edge cut at several roots) sweeps a face on the one surface the whole
curve sweeps, so faces that two parts of the sweep share are one surface to the booleans that unite
the parts. A swept face folds where the motion runs along its curve, the velocity crossed with
the curve's direction turning over; past the fold the envelope lies inside the solid, and such a
sweep's faces are trimmed (see [Folds](#folds)).

A characteristic curve is fitted only where no closed form exists, so a case whose curves are lines
and circles under a slide, turn or screw is exact throughout.

The faces are joined by construction: each curve of the boundary of `W` has its piece of the cap
at the start on one side and its piece at the end on the other, and the paths of its ends are
shared with the swept faces of its neighbours. The direction each face runs is fixed by those caps,
so no face's side has to be guessed. A vertex the motion holds still, on a turn's axis, has no
path: it is a pole of the swept faces beside it.

### Overlaps

The product is a solid only where the sweep does not run into itself. It does where a turn carries
a plane past the line it grazes along: the plane where it starts and where it ends cross beyond
that line, and the envelope between them lies inside the solid. It does where a concave edge of `W`
sweeps a sheet into the body, and where a long turn or screw brings the body back over where it
was. The checker finds each of these as faces that meet without sharing an edge. The whole
product's refusal is asked only whether it is one, so its faces are sewn by
`sew::assemble_or_refuse`, whose check (`check::check_to_first`) ends the search for such faces
at the first pair found that does not merely touch; a product with no such pair is searched to the
end, as any result is. Then the sweep is
built in parts: the body at the start, its faces cut along their curves, united one by one with
the solid each piece of face in `W` sweeps on its own (its two ends and the faces its boundary
sweeps). Each part starts on a face of the body and grazes it along its curves there; the boolean
is told so (each such edge lies on the other's face), as a blend tells it of a round's tangent
edges.

The checker pairs only faces that share no edge, and never a face with itself, so it cannot see a
swept face that passes through itself or through a swept face beside it: a body screwed about a
line close by turns back through what it swept. Those are found before the faces are sewn
([`orbit.rs`](../src/sweep/solid/orbit.rs)). A steady motion carries each point along its orbit, a
line, a circle or a helix, and an orbit is named by two lengths the motion keeps: the place in the
plane square to a slide; the distance from the axis and the height along it for a turn; the
distance and `z − c·θ` for a screw of rise `c`, read again a whole turn on. Two points of the swept
curves on one orbit, less than the motion's span apart along it, are one point of two sheets, or
of one sheet twice. So the curves are drawn in the plane of the two lengths, and the sheets cross
where two drawings cross, or one crosses itself, at points that close along their orbit. The
drawings are sampled, and each crossing of their chords is confirmed on the curves by Newton's
method: two curves that leave a shared end along one line, as a grazing curve and the edge it ends
on do, have chords that cross beside it, and the steps run to the end.

A swept face can also pass through an end of the sweep: the body's faces left where it starts, or
those carried outward where it ends. It shares an edge with the ends beside its curve, and may
cross any of them along a curve that meets no edge at all. So the paths of points along each swept
curve are intersected with each end's surface, and a meeting inside the face is a crossing: a drum
lying nearly across a screw's axis carries a sliver of its far end outward and then turns through
where the sliver went.

A product with either kind of crossing is invalid, and is built in parts or refused as any other
overlap is; without these tests it was a solid that bounded the overlap twice.

A body that is itself a union, as a ball on a stem, may instead be swept a piece at a time,
since the sweep of a union is the union of the sweeps.

The unions look for where each edge meets each face's surface by the signed distance to the
surface along the edge ([`csi::Tracker`](../src/boolean/csi.rs)). A turn's sheets are surfaces of
revolution, and the point of one nearest a point lies in the half-plane through the point and the
axis, among the curve's samples there and Newton's steps from each sampled least distance. The
samples and their places in that half-plane (`geom::Meridian`) are the same for every point, so
they are worked out once for an edge, and `Surface::project_on` projects each point from them: a
seed turned to the point's angle lies as far from the point as its place in the half-plane lies
from the point's, so only the seeds within rounding of the nearest there are placed on the surface
and compared, and the seed chosen is the one `Surface::project` chooses. A projection is not
continued from the last point's: a curve that lies in no half-plane through the axis has an image
there that turns back on itself, the surface's two sheets meet along that fold, and no step is
short enough to say the last foot's sheet is still the nearer.

### Folds

Where a swept face folds, the faces built as above bound more than the solid: past the fold, and
where the folded sheet runs through its neighbours and the caps, they lie inside it. They are
trimmed to what the moving body never reaches into, each in its own parameter plane
([`trim`](../src/sweep/solid/trim.rs); Peternell, Pottmann, Steiner and Zhao 2005; Seong, Elber and
Kim 2006). A point `p` lies inside the swept solid where the body at some `τ` holds it, which is
read along its inverted trajectory `M(τ)⁻¹·p` through the body's own frame (Ilies and Shapiro 1999;
Erdim and Ilies 2008): the deepest the body holds it over the motion, `max_τ d(M(τ)⁻¹·p)`, with `d`
the depth inside the body. A point of a face of the sweep is reached where that is more than a
hair (`1e-10` of the body's size): on the true boundary the body touches the point and never holds
it. For a body of revolution (every face a plane square to one axis, or a cylinder, cone, sphere or
torus about it) `d` is read in its meridian section, the distance from the outline of segments and
arcs, its sign from the side of the nearest piece or, at a corner, from the outline's winding. The
depth along the path is sampled finely enough that no point moves more than a fiftieth of the
body's size between samples, and each peak that could top the others is climbed by Brent's search
(Brent 1973), which steps to the top of a parabola where the depth is smooth and by golden section
where it has a corner. The depth is a distance, so it changes no faster than the point moves, and
under a slide, a turn or a screw the point's speed through the body's frame is the same all along
its path: a sample is read only where that bound lets it be the highest or beside a peak that
could be.

- Each edge is cut where the body begins to reach into it, by bisection between samples.
- Each face's plane is searched on a grid as the edges' planes are (marching squares), each grid
  line whose ends are reached differently cut by bisection; each curve found is cut where it leaves
  the face, at the vertex on the face's boundary there, and a closed one gets two vertices.
- The curves are made exact where the grid is coarse, each point moved across the curve onto it by
  Brent's root finder on the depth; each curve is cut at every vertex of the others that lies on it.
- Each curve is where two faces cross, so it has a partner between the same two vertices on another
  face, or on the same face where a sheet runs through itself. The two share one piece, fitted
  through one curve's exact points; each face's pcurve is the spline through its own points in its
  plane.
- Each region of each plane between its boundary and its curves that the body does not reach is a
  face, and the faces are sewn once. A face cut from a spline sheet, or from a swept sheet whose
  curve is a spline, lies on the sheet cut down to the face's box of parameters, the parameters
  kept (knot insertion, Piegl and Tiller 1997, ch. 5): on the whole sheet a point of it has a
  second foot across the fold or the self-crossing, which a projection may land on.

Where a curve's face folds, the sheet folds along the sweep of the place on the curve where the
motion runs along it (the velocity crossed with the curve's direction turns over, pinned by
bisection), and a curve across the cap there reaches the cap's edge at that place, tangent to the
fold. The reach is poorly conditioned there: on the cap side its depth rises only with a power of
the distance from the curve. So where the curve lies on the body to rounding at the fold, the cut
of the edge is moved to the fold, and a curve on a sheet whose partner is a cap is found as the
sign change of the body's own depth at that end of the motion, `d(M(τₑ)⁻¹·p)`, which is smooth; the
cap's pcurve is then the curve projected onto the cap.

Only the sheets of slides, turns and screws are trimmed so far, and only bodies of revolution have
a depth.

### Motions whose `W` moves

Under a frame along a curve or two motions at once ([`moving`](../src/sweep/solid/moving.rs)), `W(τ)`
slides over the body, its curves appear, vanish and pass through vertices as `τ` runs. Rather than
follow them from event to event, the boundary is found in planes of two parameters, where every
event is a crossing of curves:

- **An edge's plane** is `(τ, t)`, `t` the edge's parameter, carried to `M(τ)·e(t)`. In it, the
  curves where a face beside the edge grazes, `g(τ, e(t)) = 0`, are found on a grid as marching
  squares finds a contour (Lorensen and Cline 1987, in two dimensions), each crossing of a grid
  line pinned by bisection and a saddle on the zero set made a vertex. The plane's sides are the
  edge cut where the span starts and ends, and the paths of its two vertices, cut wherever a curve
  reaches them. Where just one face beside the edge moves outward, the region is a swept face of
  the result.
- **A face's plane** is `(τ, x)`, `x` the place along its grazing curve: on a plane face the line
  `A x + B y + C = 0` by the place from its point nearest the face's middle, on a cylinder, cone,
  extrusion or torus the curve's `u`. Where a cylinder's or cone's curves run nearly along its
  rulings it is their `v`, and at each `v` the rate is `a cos u + b sin u + c`, zero at `u =
  atan2(b, a) ± acos(−c / √(a² + b²))`: two curves. Its regions are bounded by the curves of the edges' planes carried
  over, and the face's grazing curves where the span starts and ends; where the curve lies on the
  face, the region is a swept face. It faces the way the face does, or the other way where the
  curve's point drifts toward the side the face moves outward on (the outer side of a turn).

A curve shared by two planes is one piece of the result, so the faces join by construction. The
planes' regions are traced as the faces of a planar graph, each a loop with its inside on its
left. Each plane's surface is the one map fitted once: a plane or an extrusion where the points lie
on one, else the quintic B-spline interpolating a grid of them at the plane's own parameters
(de Boor 1978),
refined until it lies within a hundredth of `LINEAR` of the points between, so the plane's
straight sides are exact parameter lines of it. Every other side is a curve whose points in the
plane are known exactly, and its pcurve is the spline through them at the parameters of the edge's
curve nearest their images, with points added where the surface over it strays from the curve by
more than a quarter of `LINEAR`, each moved from the spline's by the least change that meets the
curve (damped Gauss–Newton steps, Nocedal and Wright 2006, ch. 10), and kept on the plane's line
where all its points lie on one. Projected instead, the curve's points land anywhere along a
direction in which the map nearly drops rank, as where a grazing curve moves a thousandth as fast
as the body. The caps are the body cut at the span's ends as for a steady motion.

Where a face moves within its own surface for an instant (the sides of a key carried round an
inflexion of its path), the side it carries outward turns over at once; the motion is cut there
and each span's end reads the rate of change of the velocity instead. Where the boundary runs into
itself, the span is the union of the body at its start with what each face sweeps where it moves
outward (a point the body reaches is in it at the start, or was crossed by such a face), and
where one of those runs into itself too, the union of the sweeps over its halves. A swept face
whose map turns over inside it, past a cusp of the envelope, sends the span straight to that
union. A torus grazes along two curves of `u`, `v = atan2(−a, b)` and half a turn on; at a place
where both coefficients vanish it grazes along a whole meridian and its curves' `v` turns once
round about the place, which is refused, as a sphere is. Two faces meeting smoothly along an edge
graze along one curve there, which is one piece of both their planes. A region still far thinner
in its plane than it is long is refused rather than fitted across.

## Sewing

[`sew`](../src/sew.rs) builds a solid from given surfaces. A `Piece` is a surface and the loops of
curve spans that bound it, the outer loop first; `sew` joins pieces along the edges they share,
the way a STEP file's open shells are sewn ([`step::sew`](../src/step/mod.rs)), and the result
passes the checker or is an error.

- **Pieces.** `Piece::whole` takes a B-spline surface bounded by its four boundary parameter
  lines, each the exact rational curve of its row of the net (`NurbsSurface::iso_curve`). A line
  whose control points are all one point is a pole and bounds nothing: the loop crosses it as a
  gap in the parameter plane, which becomes a degenerate edge, as at a cone's apex or a sphere's
  pole. Two boundary lines that lie on each other, on a surface closed round a parameter, become
  one seam used twice. `Piece::planar` takes loops of spans that lie in one plane, the plane found
  by Newell's normal of the outer loop's points; its front faces the side the outer loop runs
  counter-clockwise about. `Piece::bounded` takes any surface and loops; which way the outer loop
  runs in the surface's parameters sets the piece's front.
- **Analytic edges.** A span on a whole B-spline whose control points lie in order along a line is
  that line, and a quadratic one whose every knot span has five points on one circle is that
  circle's arc: a conic is fixed by five points. Sweeps choose the analytic form for the same
  reason, and a blend rolled along a circle edge needs it.
- **Joining.** Span ends within the tolerance become one vertex. Two edges between one pair of
  vertices become one where the points a quarter, a half and three quarters along each lie on the
  other: two arcs between the same ends can both pass nearest their chord's middle at an end.
  Pieces are then turned to agree across their edges, healed onto their surfaces as a read solid
  is, and the shell is turned to face out of the volume it closes.
- **Coons patches.** `sew::coons` fills a loop of four polynomial spans with the bilinearly
  blended Coons patch (Coons 1967; Farin 2002):
  `(1 − v)·c₀(u) + v·c₁(u) + (1 − u)·d₀(v) + u·d₁(v)` less the bilinear blend of the corners. Each
  pair of opposite spans is brought to one degree and knot vector, and since a B-spline whose
  coefficients are a linear function at the Greville abscissae is that function (de Boor 1978),
  each linear blend is exact on the other direction's knots, so the patch is one B-spline surface
  through all four spans. A rational span is refused: the blend of rational curves is no B-spline.

## Slabs

[`slab`](../src/primitive/slab.rs) builds the solid between a frame's xy plane and a clamped
B-spline height field over a rectangle. Its coefficients sit at the knots' Greville abscissae, so x
and y run linearly with the parameters and the top is a graph. Each side's top edge is that field's
boundary row as a B-spline curve, and its pcurve on the top is a parameter line. Every coefficient
must be above the floor. By the convex hull property that keeps the whole field above it, so the
solid never pinches.

## STEP

[`step`](../src/step/mod.rs) reads and writes AP214 with its own ISO 10303-21 parser
([`part21.rs`](../src/step/part21.rs)).

- **Writing** puts every solid of the arena in one part's advanced B-rep, in millimetres and
  radians. Each face's STEP name is its names joined by `;`, and the header names the kernel as the
  originating system, so a read gives those names back. Every curve and surface type has an exact
  STEP form: a cone narrowing along its axis is written about the reversed axis, an intersection
  curve as its spline, and a ruled surface as the B-spline surface between its rails (two polynomial
  rails, or two conics over equal angles; any other pair is an explicit error). Voids are written
  turned inside out under `ORIENTED_CLOSED_SHELL(.F.)`, as STEP stores them. STEP's solids are
  manifold, so a solid whose pieces touch is written as the manifold boundary of its interior: each
  outer sheet a solid of its own with the voids inside it, each partner pair of an edge's coedges
  its own edge, each fan of faces round a vertex its own vertex, and neither a slit nor a vertex
  loop, which lie inside one piece's face. Each edge is a `SURFACE_CURVE` (a seam a `SEAM_CURVE`)
  with its pcurve on every face that uses it, in the surface's STEP parameters (a cone's `v` is its
  height along the STEP axis, turned with it). An affine pcurve is a line, a spline one its own
  spline, and on a plane a line or spline edge maps onto it exactly, and a circle or an ellipse that
  turns the way the plane's parameters do is that conic in them (one turning the other way has no
  conic there with the edge's own parameter). An edge that is a parameter line of its surface by its
  geometry is a line from the coedge's own parameters where it starts: a circle with the axis of a
  cylinder, a cone, a sphere, a torus or a surface of revolution (its normal along the axis, its
  centre on it) runs along `u` at one radian a radian; a meridian circle of a sphere or a torus
  (its plane through the axis, its centre the sphere's or the tube's, its radius theirs) runs along
  `v`, on a sphere only short of the poles; a line along the rulings of a cylinder, a cone or an
  extrusion runs along `v`. Each test, the conic on a plane's too, bounds how far the edge can lie
  from the pcurve's image: the tilt times the radius or the length it acts across, the offsets of
  the centre and the radius, and the gap at the coedge's start, added together, must come within
  the 0.45 of the tolerance a fitted pcurve keeps to, so the bound holds at any size of model. A
  curve that only lies near a parameter line is fitted. Any other pcurve is a cubic through the
  coedge's parameter-plane points, with points added where the surface over it strays more than 0.45 of the tolerance from
  the edge at any of seven points a span, which keeps it within half the tolerance between them (a
  part's file grows by up to 7%). A reader given no pcurves must project each traced edge onto its
  face itself. A ruled surface between conics, whose rational rails do not run at the angle's rate,
  gets no pcurves. [`WriteOptions`](../src/step/write.rs) with `pcurves` off
  ([`write_with`](../src/step/write.rs), [`write_assembly_with`](../src/step/write.rs)) writes no
  pcurve at all: every edge, a seam's too, is its curve in space alone, which the standard allows,
  and the file is smaller and written faster for a reader that works out each edge's course on its
  faces itself. The reader here reads either form as the same solid.
- **Closed splines have ends in STEP.** The kernel wraps a closed spline's parameter, and the
  parameter of an extrusion or a revolution of one along it and of a closed NURBS surface; STEP
  wraps only a conic's parameter and the angles of its elementary surfaces and of a surface of
  revolution, and gives every spline two ends. So everything is written inside the ends of what it
  lies on. An edge that runs past the end of its closed spline's turn is written on the spline of
  one turn from the edge's start. A face's coedges' pcurves are moved by whole periods so each
  loop meets end to end within one turn of its surface, a pole's line across the parameter plane
  among them, placed by a box each is known to lie in (a line's ends, a spline's control points
  over eight stretches of its edge), not by points along it, so a seam's two pcurves stay a turn
  apart. A face whose box runs past the end of the surface's turn is written on its surface over the
  turn with the face in the middle, clear of both ends: the extrusion or the revolution of its
  curve's turn from there, or the NURBS surface's net over that turn. A turn from anywhere is cut
  by knot insertion out of two turns of the spline, clamped at the ends of its own turn (a file's
  closed spline need not be) and with the second turn's homogeneous points scaled to meet the
  first's, which leaves the curve as it was however its end weights differ; so the turn keeps the
  kernel's parameters there, and the pcurves keep theirs. A NURBS surface's rows each take their
  own scale, which moves the surface by at most the ratio of the largest to the smallest less one,
  times half the widest box of one patch's control points: a surface the kernel finds closed meets
  at one scale but for rounding, which moves it by nothing measurable. A face that fits no one turn
  is written as before, on its surface as the kernel has it and without pcurves, leaving a reader
  to find its course: a loop that runs round a closed direction and on across its own seam, as a
  slab tilted through a ring's tube can leave, or a face across the join of a surface whose scales
  would move it by more than a tenth of the tolerance. An edge round its closed curve more than
  once, which the kernel never makes, is an explicit error.
- **Writing an assembly** ([`write_assembly`](../src/step/write.rs)) puts several parts in one file.
  Each part is a product with its own advanced B-rep, written once however often it is used, and
  each instance is a use of its part's product under the assembly's, placed by its motion. A part
  given a colour has its solids styled with it. Each part's representation states its own
  uncertainty.
- **What a read says of each body** is in the import's `bodies`, one for each solid: the id of the
  file's solid it came from, the name of the product whose shape holds it, the uses that place it
  from the root down, and its surface colour where a presentation style names the solid. A face
  styled on its own is in `face_colours`. Only a surface style counts: a curve's colour is not a
  body's.
- **Reading** walks the product structure (`NEXT_ASSEMBLY_USAGE_OCCURRENCE` with its placing
  transformation, shape representations joined by plain relationships, mapped items) to every
  solid's placements, converts units, builds each distinct solid once (in parallel, those with the
  most faces and spline control points first; the output does not depend on which finishes first)
  and copies it to each placement. A transformation that scales uniformly builds its solid in scaled
  units. A left-handed placement uses `transform_similar` to reflect the geometry and turn the
  loops, so the reflected solid still faces out.
  Faces of another system's file are named `<feature>.face_<i>` in placement order.
- **Entities.** A B-spline curve or surface whose kind implies its knots (`UNIFORM_CURVE`,
  `QUASI_UNIFORM_CURVE`, `BEZIER_CURVE` and their surfaces, simple or rational) is read with the
  knots its kind gives, whole numbers: a uniform spline's each once, from minus its degree up; a
  quasi-uniform spline's from 0, each end one more time than its degree and each other once; a
  Bézier spline's from 0, each end one more time than its degree and each join as many times as
  its degree. An `OFFSET_SURFACE` of a plane, cylinder, cone, sphere or torus is read as the
  surface of the same kind it is; an offset of any other surface is an explicit error, never an
  approximation. A `VERTEX_LOOP` on a cone or sphere, the apex the face closes on, becomes a seam (a
  generator or a meridian) from the face's other bound up to the apex and back, as the kernel's own
  cones and spheres are bounded; a sphere bounded by a pole's alone, as some files write a whole
  ball, takes a meridian from the other pole. Anywhere else on its face's surface it becomes a
  vertex loop; off the surface it is an explicit error. `POLY_LOOP` faces and `FACETED_BREP` solids
  are read as planes bounded by straight edges, each side shared with the face across it.
- **Contacts and voids.** Within one solid, vertices that coincide and then edges that coincide over
  their whole range are joined: a valid solid has no two vertices that close, so they are the two
  sides of a contact the file wrote apart, and the edge four faces use comes back. An edge whose two
  ends are then one vertex, on a curve that does not close, and no longer than `LINEAR` between them
  is no edge, and each loop goes on without it; one that is all of a loop is kept. A bound that runs
  along an edge and back is a slit and is kept, as a point bound is: each marks where something
  touched the face. A loop tries a seam's second use a period away before taking it as running back
  over the first. A face closed round one parameter whose two rims each go once round it, the ways
  opposite, and which has no seam gets one, so that its loop winds round the period zero times and
  closes in the parameter plane: on a cylinder, cone, sphere or torus a line through a vertex of
  each rim; on any other surface closed in one parameter, a spline, a revolution or an extrusion of
  a closed curve, the surface's own line of the other parameter (an iso curve, the profile turned,
  a line along the direction), where it is farthest round from the rims' vertices and clear of the
  face's holes. A shell is the way the volume it encloses says, outer positive and void negative:
  some files write a void's faces facing into it and turn its shell too, which read as written is a
  void inside out. A void is a hole in the material its outer shell bounds, but some files write a
  compound of bodies as one solid whose "voids" are the other bodies, beside the outer shell and
  often touching it, where joining their contacts gives edges four faces use that pair across no
  material. Every void is placed against the outer shell by its vertices and the middles of its
  edges, those within the tolerance of the shell left out, before the solid is taken as read: the
  outer shell of the solid read whole, or, where that read is rejected, the outer shell read
  alone. A void with a point outside lies beside the outer shell or passes through it, as one body
  of a compound may through another, and is read as a solid of its own, however the file turns
  it, and `Import::separated` names it; its healing and tolerances are reported for that solid. A
  strict read rejects the solid instead, naming those voids. On a surface closed both ways (a torus,
  a sphere) either sense of a loop bounds a region, so a bound's flag alone picks the face, and some
  files flag a loop against the way they list it. A bound of such a face whose every edge shared
  with a settled face runs the same way as that face's use is turned; faces on other surfaces are
  settled as written, and a bound settles once a settled face lies across one of its edges, so the
  sense passes from face to face. A sphere cap the read gave a seam to the pole on its loop's left
  has, turned, the other cap on its left, so its seam is drawn to the antipode instead. The import
  names each face turned (`Import::turned_bounds`); turning moves none of the file's geometry, only
  a seam the read made, so a strict read does it too.
- **Pcurves and poles.** Each loop is followed through its face's parameter plane by projection,
  each point continuing from the last. Where a loop leaves a pole any angle is possible: the one its
  next point implies is tried in each nearby period, a seam's second use must run a period from its
  first, and the loop must close and enclose area the way the face's orientation requires. A gap at
  a pole becomes a degenerate edge. A pole is a line of the parameter plane that is one point, known
  in closed form for each surface, and not only a point where a derivative vanishes. An edge that
  passes through a pole of a sphere between its vertices (a meridian over the pole) is split there
  before its loops are followed: the loop leaves the pole at another value of the free parameter
  than it arrived with, and only at a vertex can it jump along the pole's line. An edge's parameter
  at a vertex is an end of its curve where the two points lie within `LINEAR`, judged on the points:
  a curve that stands still at its end moves only with the square of the parameter there. The free
  parameter where a use reaches or leaves a pole is the one every reader of its pcurve finds
  (`topo::project_on`): where `u` is free, the `u` of a point just along the edge, or, where the
  edge leaves along the surface's axis as every curve through a horn torus's pole does, of a point
  the way the edge bends; where `v` is free (a B-spline blend's pole along `u = const`), the `v`
  whose `u`-line leaves the pole along the edge, from the surface's own derivative, since
  projecting a point beside the pole of a thin patch lands at any `v`.
- **A cone face lies on its cone's own nappe.** A file's cone is the whole double cone, and a face
  may lie past its apex from where the placement sits. The reader gives such a cone the placement
  about the axis turned over, from the point as far past the apex as the old origin is before it:
  the same point set, with the face on its own nappe, as every cone face the kernel makes is.
- **Healing** ([`heal.rs`](../src/step/heal.rs)). Exporters hold vertices, curves and surfaces
  together only to their own tolerance (an assembly's file may declare 0.01 mm). Wherever it can,
  the imported geometry moves instead, until each solid closes at `LINEAR`:
  1. Where an edge misses its faces or vertices by more than `LINEAR`, each vertex involved is moved
     onto every surface that meets there (and onto the curves that already fit), and each edge that
     still misses is refitted through its vertices onto both its faces' surfaces. A refit follows
     the file's curve bent to the moved vertices, so where two faces touch it stays on the line of
     the vertices rather than a second line beside it; a refit that stays on both faces everywhere
     is kept as its spline, otherwise as an intersection curve.
  2. Where faces miss each other, so that step fails (three surfaces with no common point at a
     vertex, two faces that touch rather than cross, a B-spline face that stops short of its
     corner), the surfaces of the faces involved move too ([`adjust.rs`](../src/step/adjust.rs)).
     Each gets a few unknowns: a rigid motion and its radii, or its B-spline control points
     (coincident ones together, so poles and seams stay closed). Every vertex and sample point of an
     edge on those faces is free as well, and damped Gauss–Newton steps of least weighted change
     move them until each point lies on each of its faces; moving a surface costs a hundred times
     what moving a point does, so surfaces move only where points cannot close a gap. Where an
     edge's sample slides along its edge, or a vertex whose faces all touch slides across their
     common normal, it costs as much as a surface: neither direction is fixed to first order, and a
     free point runs far along it. Each step is tried undamped and then damped, since touching faces
     leave the normal equations nearly singular. Where a B-spline face's own control points cannot
     follow the other surface closely enough (a cubic blend boundary along a curve on a cylinder),
     the step starts again with knots inserted until each parameter has 16 spans; the points it
     placed along an edge that would not fit seed that edge's refit, since two touching faces meet
     anywhere across a band. Step 1 then runs again on the moved surfaces, up to three times.
  3. Every solid records how far its geometry moved: the largest distance any vertex, edge or
     surface moved, what moved and where, and the largest surface movement on its own, which bounds
     how far a face moved. `Import::healed` carries it for each placed solid, for the caller to
     report. A solid that would move more than `MAX_HEAL` (0.05 mm) is rejected with an error naming
     the distance and what moved.
  4. Where healing fails, or moves geometry further than the file's stated uncertainty (its
     `distance_accuracy_value`), the solid is kept as the file has it instead, if it can be: each
     edge gets a tolerance of 1.5 times the widest gap sampled between its curve and its faces'
     surfaces, each vertex the widest of its edges' and of its own distance from its edges' ends and
     its faces. The tolerance is the gap measured, whatever the file states: exporters write gaps
     wider than the uncertainty they state, and a file that states none can still be kept. What
     bounds it is the model: a tolerance past a quarter of an edge it holds leaves that edge's ends
     indistinct, and is an error; where there is one and healing succeeded, the healed solid stands.
     `Import::tolerant` lists each such solid, with the uncertainty the file states and how many of
     its tolerances are past it. The writer gives a file holding local tolerances the widest of them
     as its uncertainty, which is what other kernels read; reading such a file of its own, the
     kernel keeps its geometry so rather than healing it within that uncertainty, so a solid held
     with local tolerances reads back as it was written. The conjugate-gradient solves of step 2
     stop once the residual has not halved in 200 iterations: on faces that touch the system is
     nearly singular and the residual stalls far above its goal.

  The read returns the solids that check valid and a reason for each it rejected. A solid is
  rejected for its structure, the geometry each entity lies on, and its loops
  (`check::check_structure`), not for passing through itself: a solid whose faces meet away from
  their shared edges is kept, and each place is reported in `Import::crossed` with the two faces,
  the point, and how far one passes through the other (`check::self_crossings`). One imported mesh
  of 23442 triangles folds through itself at four corners, up to 0.014 mm deep. The data is the part
  as its file gives it, and refusing it refused the whole part. A boolean that cuts through such a
  fold may still fail.

## Measurement and tessellation

[`measure`](../src/measure.rs) integrates the exact geometry. A face integral over its trimmed
parameter domain becomes a boundary integral by Green's theorem; the inner integral runs along a
parameter in which the integrand is a low-order polynomial or trigonometric polynomial (so a fixed
Gauss rule is exact to rounding; on a B-spline surface the held parameter's basis is summed into one
iso-line per boundary point), and the boundary integral is adaptive. Along a plane and round the
axis of a cylinder, cone, sphere or torus the inner integral is taken in closed form: each of
the four is a circle turned about its axis, so the point and `su × sv` are trigonometric
polynomials of degree one in the angle. Its tolerance halves with each
halving of a piece, so where the two rules differ by noise in the integrand rather than by their
error (the inner rule's own error, near 1e-8 of `|su × sv|` on a fitted sheet, jumps each time the
boundary point crosses a knot), they miss it by about as many times at every depth; once a halving
no longer cuts the miss to three quarters of the last and the miss is within a millionth of the
integrand's size, the halves are taken (the round-off test of Piessens et al. 1983), where halving
to the depth limit took minutes. Volume is `⅓∮ (p − c)·n dA`
over the faces. Area and volume agree with closed forms to about 1e-15 relative on analytic solids.
The same integral carries the first and second moments, `¼∬ xᵢ(x·n)` and `⅕∬ xᵢxⱼ(x·n) dA`, for the
centre and the inertia tensor about it, each settled to a tolerance in the size of the loose box
round the geometry, since a moment that symmetry makes zero has no size of its own to settle to. The
bounding box is tight: edge extremes in closed form or by root finding, plus interior stationary
points of surfaces whose critical sets need not meet a boundary.

A plane section (`measure::section`) is the cap the boolean of the solid and a half space leaves,
so it needs a solid and costs a boolean. A slice (`measure::slice`, `measure::slices` for parallel
planes) cuts face by face instead, so it takes a sheet as well: each edge is met with the plane in
closed form where one exists (`csi::curve_surface`), each face's surface with the plane by the
booleans' surface intersection (`ssi::intersect`), each curve split at the face's boundary meetings
and a piece kept where its middle lies inside the face (`boolean::domain`). Where the cut has to be
traced, it starts from the face's boundary meetings and from where the plane crosses the
iso-parametric lines through each critical point of the height along the normal, and the seams. A
closed loop of the cut that touches no edge bounds a disk of the parameter plane, over which the
height has an interior extremum, a critical point (Morse theory), whose iso-parametric lines run out
of the disk across the loop; a loop that winds round a periodic parameter crosses its seam. The
critical points are in closed form on the sphere and the torus; on a surface of revolution they are
the profile's extremes along two directions (`Curve::extremes`, sampled for a `Curve::Nurbs`
profile); a ruled surface has none that a loop could surround, its height being affine along each
ruling. On a spline surface they are the common zeros of two Bernstein polynomials, sought by de
Casteljau subdivision of each Bézier patch with the convex hull of the control values as the
bound, and a piece whose Jacobian's interval determinant keeps one sign solved by Newton's method;
a spline line's meetings with the plane are isolated the same way. The search looks at no more than
65,536 pieces of a patch, seeding the middles of any left, and does not resolve critical points
that are not isolated below cells 2⁻¹² of the patch on a side. Pieces run with the cut on their left seen
from the normal's side, so a solid's outer outline turns counter-clockwise about the normal, and
they join at the plane's meetings with edges, within each edge's and vertex's own tolerance and at
least ten times `LINEAR`. Where more than two pieces meet, as at a vertex where faces pinch, each
arriving piece continues on the first leaving clockwise from it round the point in the plane, so
each wedge of the cut keeps its own corner. A face lying in the plane gives its loops, turned to run with the face on their left seen from the
normal; an edge lying
in it, of no such face, is a piece of the slice; a curve the plane only touches cuts nothing,
nor does a traced loop no wider than ten times the join distance. The
faces' domains, boxes and critical points are worked out once for every plane of `slices`.

[`distance`](../src/distance.rs) finds the closest points of two selections of faces, edges,
vertices and solids. Both sides' meshes give the distance to within their deflections, each between
a ten-thousandth and a thousandth of its part's size, so a part micrometres across has triangles
that stand near its faces and one kilometres across is not meshed to 0.01 mm; every pair of mesh
pieces that close, nearest first, seeds a Newton search on the exact surfaces, curves and points,
damped to keep each step downhill, and a result counts only where both points lie on their face or
edge. The true minimum lies in one of those basins, so the least result is exact, even where a whole
ring of points is equally near. A surface whose nearest point lies off its face is covered by the
face's edges, and a solid wholly inside another, by its winding number in the other's mesh, is at
zero. Before either side is meshed, a vertex of one side or the middle of one of its edges that a
solid of the other holds, farther than the tolerance from the surface of every face whose box
reaches it, is a point of both: the distance is zero, and no mesh is made.

[`rays`](../src/rays.rs) meshes an arena once, at the same scale, and casts rays at it. The mesh
says where to look and no more: a ray's candidate faces are those with a mesh triangle within twice
the deflection of it, so no crossing is missed, and each is searched over the stretches of the ray
that pass its triangles. What is found there is every meeting of the line with the face's exact
surface, by the intersection the booleans use (`csi::curve_surface`), which never steps along the
line looking for a change of sign: two crossings closer than any step would be stepped over, as at a
gear's root. A plane, cylinder, cone, sphere and torus have closed forms, a torus whose tube reaches
its axis in both its sheets. An extrusion, a surface of revolution and a ruled surface each turn on
one parameter of their profile: the profile's distance from the plane through the line, the
difference between a ring's radius and the line's distance from the axis at the ring's height, and
the triple product that vanishes where the line and a ruling lie in one plane. Each is a polynomial
over a piece of a spline, times the spline's denominator on a rational one, and a quadratic in the
tangent of the half angle on a circle or an ellipse; its stationary points split the piece into
stretches with one crossing at most. A spline surface has no such parameter: it is cut into Bézier
patches, each read as its distances from two planes through the line and lying among its control
points, as each half of a patch cut in two does. A patch whose control points all stand to one side
of either plane holds no crossing. A patch whose slopes, bounded by the steps between its control
points, make the map from its parameters to the two distances one to one holds one crossing at most,
and is solved outright. Any other is halved, down to patches the size of the tolerance, which a
touch or two crossings a hair apart comes to. A run of more than 64 such patches along the line is a
line that lies in the surface, as a face's straight edge lies along the blend beside it, and has no
crossings to find. Only a profile that is an intersection curve, and a line that lies in a spline
surface, are still sampled. A crossing counts inside its face and range, by the face's exact domain
(`boolean::domain`), the test the booleans and the checker use. A crossing on an edge counts once. A
grazing ray stays within the tolerance of the surface along a stretch; the touch is where the
surface's normal turns square to the ray, found by bisection.

Whether a solid holds a point (`Rays::holds`) is the parity of a ray's crossings. A ray that passes
through an edge, grazes a face or starts on one does not say, and the next of a fixed set of
directions is tried; a point no ray decides is on the boundary. There is one such query: the
booleans' own classifier (`boolean::classify`) follows the same rules on the same intersection,
without the mesh, for the few points a boolean asks about.

The normal and curvature at a point of the boundary ask whether a point of a face's surface lies on
the face through `measure::Outline`. The face's loops, drawn as polygons in its parameters, decide,
except within their chord sag of the boundary, where a polygon cuts inside a curved edge: there the
nearest edge decides by its own side, and a point on an edge is on both its faces. The edges are
asked only where the point lies that near the chords along them, since projecting onto a spline edge
costs far more than walking the chords.

[`mesh`](../src/mesh/mod.rs) discretises each edge once; both faces of an edge use the same points,
so a valid solid's mesh is watertight by construction. An edge is sampled until its chords stay
within the linear and angular deflection of the curve, then thinned to the fewest of those samples
that still do, and split again wherever the surface over a chord in a face's parameter plane strays
from the curve or its normal turns along the chord, or between either end and its middle, by more
than the angular deflection (a straight edge along a twisted loft's side lies on the surface
everywhere, but its triangles must follow the turn, and the normal can turn away and back between
two ends that face the same way). A chord no longer than the linear deflection is never split for
its turn alone, so a hair-wide circle keeps its samples well apart. A knot within 1e-9 of the edge's
parameter range of either end is not a sample point, so rounding cannot crowd samples onto a vertex.
Edges that join one pair of vertices with one chord each get a sample at their middle: four faces
would otherwise use that one segment (two short arcs where two holes' circles cross on a face).
Chords of one face that cross in its parameter plane, as two edges leaving a point at a hair's angle
can, are halved until they do not (no further than the edge's tolerance, below which their
parameter-plane points are rounding), and so is a curved chord with another of the face's boundary
points between it and its curve, which a coarse deflection would leave on the wrong side (a disc's
rim drawn as a square past its bolt holes). A loop of only two chords gets a point between them.
Deflections must be positive and finite. Where a face's loops pass one vertex twice (a crescent
pinched where two circles touch), images of it that only rounding parts, far closer than a period or
a pole's span, become one point, so the face is meshed through that point once rather than twice a
hair apart. Edges and faces are meshed on as many threads as the machine offers, each result in its
place, so the mesh is the same on any number; a solid of flat faces and at most 128 straight edges
is meshed on the calling thread alone, being less work than starting the others.

[`display`](../src/mesh/display.rs) makes what a viewer draws from that mesh: each face's triangles
with nodes of its own, so shading creases where faces meet, each node's normal the exact surface's
there, faces in the order the STEP writer writes them, each with the id of the face it draws; and
the edges between faces as polylines, each with the faces it joins, the least and greatest angle
between their outward normals along it (positive where the edge is convex, 0 where the faces meet
smoothly), and whether the two lie on one surface, as where a face was cut in two. It reports what
it measured and leaves to the viewer which edges count as creases. Seams, where a face meets itself,
are listed apart, and poles are left out. It packs the faces as one binary glTF mesh (`POSITION`,
`NORMAL`, a `_FACE` index per node, optionally `COLOR_0`), so a viewer draws a part without meshing
it again; `glb_with_edges` adds the edges as lines, each node with its edge's index as `_EDGE` into
a table of faces, turns and splits in the mesh's `extras`, and the seams and split edges as lines of
their own.
`display_within` takes the mesh's fineness, coarser for a part that is only context.

Each face is a constrained Delaunay triangulation of its boundary samples in its parameter plane
([`cdt.rs`](../src/mesh/cdt.rs), with an exact orientation test in
[`predicates.rs`](../src/mesh/predicates.rs)), scaled per face so that one unit along each parameter
is the longest step the deflection limits allow along it, at most the face's size
([`face.rs`](../src/mesh/face.rs)). A cylinder's triangles therefore run the length of its straight
lines. A triangle is refined while the surface at its centroid or an interior edge's middle strays
from it by more than the linear deflection (the centroid may stray by its boundary segments' own sag
on top), or the normals at its corners turn by more than the angular one, and, on a curved face,
while it is skinnier than about 20° and wider than one unit. However small, a triangle across which
the normal turns more than a right angle is refined: a ball narrower than the deflection, sampled at
its poles and seam alone, would otherwise be triangles from pole to pole that fold flat. The point
goes in at the triangle's circumcentre, or at its worst sample when the circumcentre lies outside
the face. Boundary segments are never split. A point within 1e-3 of its triangle's longest side of
one of the triangle's corners, or in a triangle that flat, is not inserted: a triangle whose error
lies along a boundary segment would otherwise be refined ever nearer that segment, leaving slivers.
Nor is a point within `LINEAR` in space of the line through two corners of its triangle, and once
refined, each triangle with a corner within `LINEAR` of its longest side in space has that side
flipped where both new triangles are taller: the scaled plane stretches a face far narrower one way
than the other (a pocket's round of radius 1e-3 mm, 3.8 mm long), so a triangle well shaped there
can be flat in space. An imported assembly of 605 solids meshed with 22 such triangles, and now with
none. Points outside the triangulation's box are never located, since the orientation test overflows
there. The parameter-plane scale reads curvature from one-sided differences, keeping the smaller, so
a singularity a hair outside the face (a spindle torus's apex) cannot shrink it. At a pole, where
one point spans a stretch of the parameter plane, each triangle is sampled with its pole corner
moved along the pole to face the triangle's other corners, and two triangles on the same three
points facing opposite ways, a point joined to the pole on both sides of one chord, are dropped:
they cover nothing. A face whose refinement passes a cap on its points is an error, never a hang.

Faces mesh in parallel. Afterwards an interior edge between two boundary points that some other
triangle of the mesh also uses is split at its middle, so no mesh edge is used more than twice.

## Module map

| Module        | Owns                                                                      |
| ------------- | ------------------------------------------------------------------------- |
| `math`, `tol` | Points, vectors, frames, rigid transforms, boxes; the two tolerances      |
| `numeric`     | Gauss–Legendre rules, adaptive quadrature, dense and banded solves, roots |
| `geom/`       | Curves and surfaces                                                       |
| `topo/`       | The arena, builders, names, pcurves, facts, tolerances, similarities      |
| `profile`     | Planar sketch profiles, the readings of point lists, regions combined     |
| `primitive`   | Block, cylinder, sphere, ellipsoid, torus, cone; a height-field block     |
| `sweep`       | Profiles extruded, turned, along legs, a helix or a curve; moving solids  |
| `loft`        | Solids skinned through planar sections                                    |
| `boolean/`    | Intersection (closed form, swept, traced), contact, splitting, sewing     |
| `blend`       | Chamfers and fillets: one size, per edge, shaped sections, corner patches |
| `edit`        | Faces moved, tilted or deleted; shells, offsets, thickens, affine maps    |
| `pick`        | Faces and edges found by the names faces carry, and convex edges          |
| `step/`       | STEP AP214 import and export, product structure, healing                  |
| `mesh/`       | Watertight tessellation (constrained Delaunay), the mesh report, display  |
| `measure`     | Area, volume, mass properties, boxes, extents, plane sections, slices     |
| `distance`    | Closest distance between selections, and from a point to one              |
| `overlap`     | Whether placed solids share material, and how much                        |
| `draw`        | Hidden-line drawings                                                      |
| `affine`      | A solid stretched, sheared or mirrored exactly                            |
| `codec`       | The binary encoding of built geometry kept on disk                        |
| `rays`        | Where rays cross faces, touches included                                  |
| `sew`         | Solids sewn from bounded pieces of given surfaces; Coons patches          |
| `sheet`       | Sheets: faces in no solid, taken from a solid, turned over, thickened     |
| `check`       | The validity checker                                                      |
| `probe`       | Thin places in a valid solid: narrow faces, slivers and gaps, for hygiene |

Heavy tests stay out of this crate: its own tests are fast unit tests.

## References

The papers this document and the source cite by author and year.

- **Abdel-Malek and Yeh 1997.** K. Abdel-Malek and H.-J. Yeh, "Geometric representation of the
  swept volume using Jacobian rank-deficiency conditions", _Computer-Aided Design_ 29(6), 1997,
  pp. 457–468.
- **Aichholzer et al. 1995.** O. Aichholzer, F. Aurenhammer, D. Alberts and B. Gärtner, "A novel
  type of skeleton for polygons", _Journal of Universal Computer Science_ 1(12), 1995,
  pp. 752–761.
- **Blackmore, Leu and Wang 1997.** D. Blackmore, M. C. Leu and L. P. Wang, "The sweep-envelope
  differential equation algorithm and its application to NC machining verification",
  _Computer-Aided Design_ 29(9), 1997, pp. 629–637.
- **Brent 1973.** R. P. Brent, _Algorithms for Minimization without Derivatives_, Prentice-Hall,
  1973.
- **Choi and Ju 1989.** B. K. Choi and S. Y. Ju, "Constant-radius blending in surface modelling",
  _Computer-Aided Design_ 21(4), 1989, pp. 213–220.
- **Coons 1967.** S. A. Coons, _Surfaces for Computer-Aided Design of Space Forms_, MIT Project MAC
  report MAC-TR-41, 1967.
- **de Boor 1978.** C. de Boor, _A Practical Guide to Splines_, Springer, 1978.
- **Erdim and Ilies 2008.** H. Erdim and H. T. Ilieş, "Classifying points for sweeping solids",
  _Computer-Aided Design_ 40(9), 2008, pp. 987–998.
- **Farouki and Neff 1990.** R. T. Farouki and C. A. Neff, "Analytic properties of plane offset
  curves" and "Algebraic properties of plane offset curves", _Computer Aided Geometric Design_ 7,
  1990, pp. 83–99 and 101–127.
- **Farin 2002.** G. Farin, _Curves and Surfaces for CAGD: A Practical Guide_, 5th ed., Morgan
  Kaufmann, 2002.
- **Gordon 1969.** W. J. Gordon, "Spline-blended surface interpolation through curve networks",
  _Journal of Mathematics and Mechanics_ 18(10), 1969, pp. 931–952.
- **Gregory 1974.** J. A. Gregory, "Smooth interpolation without twist constraints", in R. E.
  Barnhill and R. F. Riesenfeld (eds.), _Computer Aided Geometric Design_, Academic Press, 1974.
- **Hartmann 2001.** E. Hartmann, "Parametric Gⁿ blending of curves and surfaces", _The Visual
  Computer_ 17, 2001, pp. 1–13. The source of the word "thumbweight" for a blend section's shape
  parameter.
- **Ilies and Shapiro 1999.** H. T. Ilieş and V. Shapiro, "The dual of sweep", _Computer-Aided
  Design_ 31(3), 1999, pp. 185–201.
- **Jackson 1995.** D. J. Jackson, "Boundary representation modelling with local tolerances",
  _Proceedings of the Third Symposium on Solid Modeling and Applications_ (Solid Modeling '95), ACM,
  1995, pp. 247–253.
- **Kronrod 1965.** A. S. Kronrod, _Nodes and Weights of Quadrature Formulas_, Consultants
  Bureau, 1965.
- **Lorensen and Cline 1987.** W. E. Lorensen and H. E. Cline, "Marching cubes: a high resolution 3D
  surface construction algorithm", _Computer Graphics_ 21(4) (SIGGRAPH '87), 1987, pp. 163–169.
- **Martin and Stephenson 1990.** R. R. Martin and P. C. Stephenson, "Sweeping of three-dimensional
  objects", _Computer-Aided Design_ 22(4), 1990, pp. 223–234.
- **Nocedal and Wright 2006.** J. Nocedal and S. J. Wright, _Numerical Optimization_, 2nd ed.,
  Springer, 2006.
- **Peternell et al. 2005.** M. Peternell, H. Pottmann, T. Steiner and H. Zhao, "Swept volumes",
  _Computer-Aided Design and Applications_ 2(5), 2005, pp. 599–608.
- **Piegl and Tiller 1997.** L. Piegl and W. Tiller, _The NURBS Book_, 2nd ed., Springer, 1997.
- **Piessens et al. 1983.** R. Piessens, E. de Doncker-Kapenga, C. W. Überhuber and D. K. Kahaner,
  _QUADPACK: A Subroutine Package for Automatic Integration_, Springer, 1983.
- **Pottmann and Wallner 2001.** H. Pottmann and J. Wallner, _Computational Line Geometry_,
  Springer, 2001.
- **Rossignac and Requicha 1984.** J. R. Rossignac and A. A. G. Requicha, "Constant-radius blending
  in solid modelling", _Computers in Mechanical Engineering_ 3(1), 1984, pp. 65–73.
- **Rossignac and Requicha 1986.** J. R. Rossignac and A. A. G. Requicha, "Offsetting operations
  in solid modelling", _Computer Aided Geometric Design_ 3(2), 1986, pp. 129–148.
- **Seong, Elber and Kim 2006.** J.-K. Seong, G. Elber and M.-S. Kim, "Trimming local and global
  self-intersections in offset curves/surfaces using distance maps", _Computer-Aided Design_ 38(3),
  2006, pp. 183–193.
- **Shewchuk 1997.** J. R. Shewchuk, "Adaptive precision floating-point arithmetic and fast robust
  geometric predicates", _Discrete & Computational Geometry_ 18, 1997, pp. 305–363.
- **Sugihara and Iri 1989.** K. Sugihara and M. Iri, "A solid modelling system free from topological
  inconsistency", _Journal of Information Processing_ 12(4), 1989.
- **Wang and Wang 1986.** W. P. Wang and K. K. Wang, "Geometric modeling for swept volume of moving
  solids", _IEEE Computer Graphics and Applications_ 6(12), 1986, pp. 8–17.
- **Wang et al. 2008.** W. Wang, B. Jüttler, D. Zheng and Y. Liu, "Computation of rotation
  minimizing frames", _ACM Transactions on Graphics_ 27(1), 2008.
- **Weld and Leu 1990.** J. D. Weld and M. C. Leu, "Geometric representation of swept volumes with
  application to polyhedral objects", _The International Journal of Robotics Research_ 9(5), 1990,
  pp. 105–117.
