# opencadkernel

A geometry kernel for CAD work: 2D curve algebra with a B-rep solid layer
built on top of it.

Early. The structure and the invariants are settled; most of the operations
are not written yet.

## Layering

```
acis    lift a SatDocument into brep, lower it back        [feature = "acis"]
  │
brep    owned mutable topology, SSI, boolean, blends       [feature = "brep"]
  │
geom2d  curves, intersection, containment, tessellation    [feature = "geom2d"]
  └ offset  parallel offset of a polyline                  [feature = "offset"]
```

`geom2d` is the default feature, compiles alone, and pulls in nothing.

The direction is the point. Splitting a face during a boolean means
projecting the intersection curve into that face's `(u, v)` space and running
a loop boolean there — a 2D problem. The quality of `geom2d` therefore caps
the quality of `brep`, which is why 2D comes first rather than as an
afterthought.

## Where this sits

```
OpenCADStudio ──► opencadcodec      DWG/DXF/ACIS
      ├─────────► opencadkernel
      └─────────► acadifc       IFC ↔ CAD conversion

acadifc ──► opencadcodec
       └──► opencadkernel

opencadkernel ──► opencadcodec          under `acis` only
```

Consumers depend on this crate directly rather than through each other, so
an editor that never touches IFC does not build the IFC subsystem.

`acadifc` is a sibling consumer, not a layer above. IFC conversion is
geometry work — `IfcAdvancedBrep` carries NURBS surfaces, `IfcFacetedBrep`
is polygonal, `IfcExtrudedAreaSolid` is a sweep, and profiles are 2D curves —
so it reaches into `brep` and `geom2d` the same way an editor does.

**This crate must not depend on `acadifc`.** That edge would close a cycle,
and it is the one direction the graph above forbids.

Conversion into CAD is also why `brep::Provenance` has a `Synthesized`
variant: geometry arriving from IFC has no source record to fall back on, so
every node it produces is emitted in full.

## Two invariants

**Local frames.** Every operation lifts its input into a translation-only
local frame before doing arithmetic and shifts results back on the way out.
Drawings live at survey coordinates; an origin at 1.2e6 spends seven of
`f64`'s ~15 significant digits before anything has happened. What follows
does not look like a precision bug — arc fits quietly degrade to polylines,
stepping along a curve overflows to infinity, point-in-polygon flips sign at
an edge and a fill bleeds. Patching each site as it appears does not
converge. Lifting at the boundary does. See `geom2d::frame`.

**Provenance.** Every B-rep node remembers the record it was lifted from and
whether it has been edited since. Lowering writes clean nodes back as their
original bytes, so attributes, parameter-space curves and surface types the
kernel does not model survive an edit that did not touch them. Only dirtied
and newly built nodes are re-emitted. Without this, one boolean rewrites
every body in the file and quietly drops whatever the kernel cannot express.
See `brep::Provenance`.

## Scope

Curves: line, ray, infinite line, circle, arc, elliptical arc, polyline with
bulges, and NURBS. Intersection is closed-form wherever one exists, a polyline
is taken apart into pieces that have one, and what remains is subdivided and
refined until it converges to the tolerance asked for. Nothing is sampled at a
fixed density and no result is left approximate.

Nearest points and containment come off the same dispatch: the point on a
curve nearest a cursor, and whether a point lies inside a boundary — the
latter by casting a ray and counting crossings, retried from another angle
when the ray grazes a corner rather than passing through it.

In: 2D curve intersection, offset, containment, chaining, area and length;
B-rep topology, surface–surface intersection, booleans; exact ACIS lift and
lower.

Out, for now: general fillet, chamfer and shelling. Constant-radius blends
between analytic faces are tractable and planned; corner blends where three
or more meet, variable radius, and blends on spline faces are not, and
pretending otherwise would be a poor use of anyone's time.

Analytic surfaces come first throughout. Plane, cone, cylinder, sphere and
torus intersect each other in closed form, and that covers most of the solids
that turn up in real drawings. NURBS-against-NURBS intersection is a separate
and much larger problem; a body this kernel cannot edit exactly should say so
rather than approximate quietly.

## Dependencies

One so far, and only under `offset`:

| Crate | For |
| --- | --- |
| `cavalier_contours` | polyline offset — the surviving-pieces problem, not the moving-segments one |

Offsetting is behind its own feature because of it: the default build stays
dependency-free. The version is pinned, because the offset path reaches into
the crate's internals rather than only its public API.

The rest are intended but not yet added; each arrives with the code that needs
it:

| Crate | For |
| --- | --- |
| `robust` | exact orientation predicates — classification is where naive kernels break |
| `rstar` | spatial index for intersection candidate pruning and snapping |
| `iTriangle` | triangulation for fills |
| `curvo` | NURBS evaluation, derivatives, knot insertion |

Vector maths is this crate's own — `geom2d::Vec2`, with `Vec3` and a transform
to follow when the B-rep layer starts. Not `glam` or `nalgebra`: what is needed
is add, subtract, scale, dot, cross and normalise, and taking a dependency for
that would pin every consumer to a version on this crate's account. Public
signatures stay `[f64; 2]`, which converts for free and imports nothing. If a
decomposition is ever genuinely needed — fitting a surface through a mesh, say
— that module can take one behind its own feature.

One decision is still open: whether the ACIS bridge depends on the codec crate
directly, or inverts so the codec stays unaware of the kernel.

## Licence

MPL-2.0.
