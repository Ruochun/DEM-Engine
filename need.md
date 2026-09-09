Below is a writeup intended to be handed directly to another engineering/coding agent. I have phrased it as an implementation design rather than as a paper summary, and I distinguish the Feng-derived part from the proposed DEME hybrid where appropriate.

---

# Design Proposal: Feng-Style Intersection-Boundary Contact Geometry for DEM-Engine v3 Mesh–Mesh Contacts

## 1. Objective

Implement an alternative mesh–mesh contact-geometry derivation in the `Mesh_Particles` branch of `Ruochun/DEM-Engine`.

The proposed method should:

1. **retain DEM-Engine v3's existing GPU broad-phase/contact-candidate infrastructure;**
2. **retain its mesh-patch/body-pair contact aggregation architecture;**
3. replace, for mesh–mesh contacts, the current triangle-projection-based derivation and subsequent voting of:

   * contact area,
   * contact normal,
   * contact point,

   with a **Feng-style additive surface-intersection-segment formulation**;
4. optionally retain DEME's existing triangle-projection machinery, or a cheaper derived variant of it, **only for estimating a representative penetration depth** \(\delta\);
5. continue presenting the normal DEME force model with an effective contact record approximately of the form

$$
(A,\mathbf n,\delta,\mathbf x_c),
$$

so that existing user-defined linear, Hertzian, damping, friction, cohesion, etc. force models remain usable.

The motivation is to obtain **more globally geometric, tessellation-insensitive contact area/normal/contact-point estimates**, while preserving the existing flexibility and GPU architecture of DEM-Engine.

---

# 2. Background

## 2.1 Current DEM-Engine v3 approach

In the `Mesh_Particles` branch, mesh surfaces are represented by triangles and the existing GPU infrastructure already performs:

* spatial binning / broad phase;
* candidate triangle-pair generation;
* triangle-prism or triangle-level contact tests;
* mesh-patch grouping;
* computation/reduction of geometric contact quantities;
* patch-based force calculation.

The current triangle–triangle geometric resolution uses a **bidirectional projection method**. One triangle is projected onto the other's plane, the projected polygon is clipped against the reference triangle via Sutherland–Hodgman clipping, and the reverse projection is also evaluated. The code then selects an appropriate projection, currently favoring the direction with smaller projection distance, to determine quantities such as area, normal, depth, and centroid/contact point.

Downstream patch-based force calculation already expects aggregated quantities including:

```cpp
finalAreas
finalNormals
finalPenetrations
finalContactPoints
```

and then passes them to the user-customizable DEME force model.

Therefore the proposed work should **not redesign the entire mesh-particle system**. It should primarily replace the geometric reduction used to produce those final patch-level values.

---

# 3. Feng contact theory relevant to this implementation

The source paper is:

**Y.T. Feng, "An energy-conserving contact theory for discrete element modelling of arbitrarily shaped particles: Contact volume based model and computational issues", CMAME 373 (2021), 113493.** 

Its central observation for triangulated solid surfaces is that the geometric quantities defining a contact can be computed from the **boundary of the penetration region**, rather than from the penetration volume or contact surfaces themselves.

For two intersecting closed bodies,

$$
\Gamma = S_1 \cap S_2
$$

is the intersection curve between their surfaces. For triangulated bodies, \(\Gamma\) becomes a collection of line segments produced by intersecting triangle pairs. 

Critically, those segments **do not have to be explicitly stitched into an ordered closed polyline** if only the additive contact quantities are needed. Each segment's contribution can be calculated independently and summed. 

That property is the main reason the approach is attractive for DEME/GPU implementation.

---

# 4. Contact-area vector

Suppose one intersecting triangle pair produces an oriented intersection segment

$$
\Delta\Gamma_i:
\qquad
\mathbf x_i \rightarrow \mathbf x_{i+1}.
$$

Feng derives the contribution

$$
\boxed{
\mathbf S_i
=
\frac12\,
\mathbf x_i\times\mathbf x_{i+1}
}
$$

and therefore

$$
\boxed{
\mathbf S
=
\frac12
\sum_i
\mathbf x_i\times\mathbf x_{i+1}
}
$$

for the whole contact boundary. 

Here \(\mathbf S\) is a **vector area**.

From it,

$$
A = \|\mathbf S\|
$$

and the contact normal can be obtained from its direction.

Care must be taken with the DEME convention for whether the final normal points from body B toward A or vice versa. Feng's orientation is tied to the ordered intersection segment

$$
\boldsymbol\tau = \mathbf n_1 \times \mathbf n_2.
$$

The implementation must explicitly map this convention onto DEME's current `B2A` convention rather than assuming the signs coincide. Feng defines the oriented intersection direction from the two outward surface normals. 

This sign issue should be covered by unit tests.

---

# 5. Contact geometric moment and contact point

Feng also defines

$$
\mathbf G
=
\int_{S_1}
\mathbf r\times d\mathbf S.
$$

For the segmented intersection curve this can be evaluated additively.

For segment

$$
\mathbf x_i\rightarrow\mathbf x_{i+1},
\qquad
\Delta\mathbf x_i
=
\mathbf x_{i+1}-\mathbf x_i,
$$

the paper gives

$$
\boxed{
\mathbf G_i
=
-\frac13
\left[
\mathbf x_i\cdot\mathbf x_{i+1}
+
\frac13
\Delta\mathbf x_i\cdot\Delta\mathbf x_i
\right]
\Delta\mathbf x_i
}
$$

and

$$
\boxed{
\mathbf G=\sum_i\mathbf G_i.
}
$$

The paper notes that this form is somewhat cheaper than its equivalent surface-triangle formulation. 

The normal contact line is then

$$
\mathbf x_c(\lambda)
=
\frac{\mathbf n\times\mathbf G}{A}
+
\lambda\mathbf n.
$$

Thus the geometry alone determines a **contact line**, rather than one unique point. 

For DEME we probably do **not initially need to implement Feng's more involved minimum-surface criterion for selecting \(\lambda\)**.

A simpler first implementation is proposed below.

---

# 6. Recommended DEME hybrid interpretation

The proposed DEME implementation is **not identical to Feng's original contact constitutive model**.

Feng's recommended linear energy model uses

$$
w(V_c)=k_nV_c,
$$

giving

$$
F_n=k_n A.
$$

In that special case the overlap volume itself disappears and only the surface-intersection boundary is needed. 

However, DEME intentionally supports arbitrary user-written force laws that often depend on a scalar penetration depth:

$$
F_n=f(\delta,\dot\delta,R,E,\ldots).
$$

Therefore the proposed hybrid is:

$$
\boxed{
\text{Feng geometry}
+
\text{DEME penetration measure}
+
\text{DEME arbitrary force law}.
}
$$

In particular:

### Obtain from Feng-style intersection reduction

$$
A,
\qquad
\mathbf n,
\qquad
\mathbf x_c.
$$

### Obtain separately

$$
\delta.
$$

### Then give DEME's existing force model

$$
(A,\mathbf n,\delta,\mathbf x_c).
$$

This preserves DEME's existing contact-model abstraction while improving the geometry supplied to it.

Note that once an arbitrary \(\delta\)-based force law is used, **Feng's proof of exact energy conservation no longer automatically applies**. We are borrowing the geometric construction, not claiming the complete Feng constitutive theory.

---

# 7. Proposed GPU algorithm

The desired GPU pipeline is approximately:

```text
Existing DEME broad phase
        |
        v
candidate mesh triangle pairs
        |
        v
existing cheap candidate / prism test
        |
        v
exact triangle-surface intersection
        |
        +--> intersection endpoints x0, x1
        |
        +--> optional penetration estimate delta_i
        |
        v
per-triangle-pair contributions
        |
        +--> S_i
        +--> G_i
        +--> penetration statistics
        |
        v
group/reduce by DEME patch-pair
        |
        +--> S = sum(S_i)
        +--> G = sum(G_i)
        +--> delta_patch
        |
        v
derive
        A      = |S|
        n      = normalize(S)
        x_c    = function(S,G)
        delta  = retained DEME measure
        |
        v
existing finalAreas/finalNormals/
finalPenetrations/finalContactPoints
        |
        v
existing patch-based force kernel
```

This should reuse as much of the current infrastructure as possible.

---

# 8. Narrow-phase triangle–triangle intersection

## 8.1 Desired output

The new narrow-phase geometry routine should conceptually have an interface such as:

```cpp
template <typename Vec, typename Scalar>
__device__ bool calcTriTriIntersectionSegment(
    const Vec& A0,
    const Vec& A1,
    const Vec& A2,
    const Vec& B0,
    const Vec& B1,
    const Vec& B2,
    Vec& segStart,
    Vec& segEnd
);
```

A successful result means the two *surface triangles themselves* intersect along a line segment.

The output must be consistently oriented.

---

## 8.2 Segment orientation

Let

$$
\mathbf n_A
=
\frac{
(B_A-A_A)\times(C_A-A_A)
}{
\|(B_A-A_A)\times(C_A-A_A)\|
},
$$

and similarly \(\mathbf n_B\).

The segment tangent required by Feng is

$$
\mathbf t
=
\mathbf n_A\times\mathbf n_B.
$$

If the initially generated segment has

$$
(\mathbf x_1-\mathbf x_0)\cdot\mathbf t<0,
$$

swap the endpoints.

That ensures all local segment contributions carry the proper orientation and therefore sum with the correct cancellation properties.

This step is essential.

---

# 9. Local contribution kernel

For every accepted surface-intersection segment:

```cpp
Vec dx = x1 - x0;

Vec S_i = 0.5 * cross(x0, x1);

Scalar q =
    dot(x0, x1)
    + (1.0 / 3.0) * dot(dx, dx);

Vec G_i =
    -(1.0 / 3.0) * q * dx;
```

Store or immediately reduce:

```text
S_i.x
S_i.y
S_i.z

G_i.x
G_i.y
G_i.z
```

and any quantities needed for \(\delta\).

These operations are cheap compared with the present clipping procedure.

One important theoretical benefit of Feng's construction is precisely that the segment contributions are independent and additive; global reconstruction of the contact boundary is unnecessary for the basic geometry. 

The paper specifically notes that the local evaluations can be performed independently and only require summation afterward, making the procedure highly parallelizable. 

---

# 10. Reference-origin issue

The formula

$$
\mathbf S_i=\frac12\mathbf x_0\times\mathbf x_1
$$

is mathematically independent of translation after summing a complete closed oriented contact boundary.

Nevertheless, computing cross products using very large world-space coordinates can unnecessarily amplify floating-point cancellation.

DEME should therefore preferably evaluate the segment formula relative to a convenient local origin.

For example,

$$
\mathbf y_0=\mathbf x_0-\mathbf x_{\rm ref},
$$

$$
\mathbf y_1=\mathbf x_1-\mathbf x_{\rm ref}.
$$

Then calculate

$$
\mathbf S_i
=
\frac12\mathbf y_0\times\mathbf y_1.
$$

A good candidate for \(\mathbf x_{\rm ref}\) is:

* midpoint between the two owner centers;
* one owner center;
* patch reference location.

Feng explicitly notes that a constant coordinate shift does not change the final vector area. 

For \(\mathbf G\), however, the translation transformation must be handled consistently. Do not casually apply the same shift without deriving/updating the corresponding \(G\) transformation.

For a first implementation, one conservative strategy is:

* compute \(\mathbf S\) in locally shifted coordinates;
* compute \(\mathbf G\) in `double3` global/world coordinates;
* optimize this later if profiling shows it matters.

---

# 11. GPU reduction by patch pair

The existing patch architecture should be preserved.

Each triangle-pair segment is associated with the same patch-pair key currently used to aggregate triangle-level contacts.

Conceptually define temporary arrays such as

```cpp
double3 segmentS;
double3 segmentG;
double  segmentDepthValue;
double  segmentDepthWeight;
PatchPairID pairID;
```

Then group/sort/reduce by `pairID`.

Final reduction:

$$
\mathbf S_{\rm patch}
=
\sum_i\mathbf S_i
$$

$$
\mathbf G_{\rm patch}
=
\sum_i\mathbf G_i.
$$

Existing CUB-based sort/reduce infrastructure should be reused where practical.

Do **not** perform atomic additions directly into body-pair records unless benchmarks establish that contention is small. DEME already has a data-parallel reduction philosophy; keeping local outputs contiguous and reducing afterward will probably scale better.

---

# 12. Computing final contact area and normal

After reduction,

$$
A=\|\mathbf S\|.
$$

If

$$
A < A_{\rm tolerance},
$$

the patch should be treated as degenerate or routed through a fallback path.

Otherwise,

$$
\mathbf n_F
=
\frac{\mathbf S}{\|\mathbf S\|}.
$$

Then convert this to DEME's expected normal orientation.

The current patch force kernel comments indicate that its effective normal is stored as:

```cpp
float3 B2A = finalNormals[relativeIndex];
```

that is, the normal felt by A and pointing from B to A.

Therefore explicit body-orientation validation is required before assigning `finalNormals`.

---

# 13. Contact-point calculation

## 13.1 First implementation

Use the zero-\(\lambda\) point on Feng's contact line:

$$
\boxed{
\mathbf x_{c,0}
=
\frac{
\mathbf n\times\mathbf G
}{
A
}
}
$$

assuming that \(\mathbf G\) has been evaluated in the same world coordinate convention as Feng's formula.

This should be implemented and validated first.

---

## 13.2 Optional improved point

Feng notes that any point along

$$
\mathbf x_c
=
\mathbf x_{c,0}
+
\lambda\mathbf n
$$

produces the same normal-force moment.

The paper subsequently derives a minimum-squared-area criterion for selecting a unique point. 

That can be implemented later if needed.

For DEME, a simpler physically intuitive alternative may actually be preferable:

1. calculate Feng's line;
2. calculate a weighted average of triangle-pair contact positions;
3. project that average onto the Feng contact line.

If \(\bar{\mathbf x}\) is the current weighted geometric center,

$$
\lambda
=
(\bar{\mathbf x}-\mathbf x_{c,0})\cdot\mathbf n,
$$

$$
\boxed{
\mathbf x_c
=
\mathbf x_{c,0}
+
\lambda\mathbf n.
}
$$

This preserves the force-equivalent Feng line while selecting the point closest to DEME's intuitive patch location.

This is **a proposed DEME hybrid**, not something stated in Feng's paper.

I recommend this as an optional second-stage improvement rather than part of the initial implementation.

---

# 14. Retaining penetration depth

This is the largest design choice.

Feng's intersection-boundary formulation intentionally does not require conventional penetration depth in its preferred linear contact model. DEME does.

Therefore retain an independent penetration estimate.

Three implementation options should be considered.

---

## Option A — retain existing projection depth

Continue to evaluate the existing bidirectional projection calculation but only extract

$$
\delta_i.
$$

Do not use its:

* projected area;
* projected normal;
* projected centroid.

Advantages:

* smallest behavior change;
* existing force laws retain familiar depth semantics;
* easiest validation against current DEME3.

Disadvantage:

* much of the current projection/clipping cost may remain unless a cheaper depth-only path is created.

**Recommended for the first working implementation.**

---

## Option B — implement a cheaper depth-only routine

Once an intersecting triangle pair is known, calculate candidate vertex-to-opposite-plane penetrations and choose a symmetric shallowest translation-like measure.

This avoids polygon clipping entirely.

One possible local quantity is based on the two directional plane penetrations:

$$
d_{A\rightarrow B},
\qquad
d_{B\rightarrow A}.
$$

Then select an appropriate symmetric value, potentially analogous to DEME's current shorter-projection criterion.

The exact definition must be decided carefully because triangle–triangle penetration depth is not uniquely defined for arbitrary nonparallel triangles.

This should be treated as an optimization after Option A establishes correctness.

---

## Option C — patch-level penetration from existing pair values

Compute per-pair \(\delta_i\) and reduce with a weighted statistic such as

$$
\delta_{\rm patch}
=
\frac{
\sum_i w_i\delta_i
}{
\sum_i w_i
}.
$$

Candidate weights include:

$$
w_i=A_i,
$$

or

$$
w_i=A_i\delta_i.
$$

However, the Feng segment itself does not naturally possess a scalar "area \(A_i\)" comparable to DEME's old projected polygon. Therefore if this weighting is retained, it may require current projection information.

A simpler initial choice may be:

* maximum valid penetration;
* average directional penetration;
* current DEME patch depth behavior.

The first implementation should aim to **preserve existing DEME depth semantics as closely as practical**, rather than simultaneously redesigning the force model.

---

# 15. Important difference: surface intersection vs prism/margin contact

DEME's mesh implementation uses finite triangle "sandwiches"/prisms and margins for contact detection and persistence. The current code therefore has contact states in which the effective collision primitives overlap even if the mathematical zero-thickness surface triangles do not yet intersect.

The Feng geometry only exists when the actual surface triangles intersect.

This means the new method cannot simply replace every geometry calculation unconditionally.

The implementation should distinguish at least:

```text
Case A:
actual surface triangles intersect
    -> use Feng segment geometry

Case B:
only margin/prism overlap exists
    -> no true surface-intersection segment
```

For Case B, several policies are possible.

### Recommended initial policy

Use the current DEME contact geometry as a fallback for margin-only contacts.

Thus:

```cpp
if (surfaceIntersectionExists) {
    use Feng segment contribution;
} else if (marginContactExists) {
    use old DEME geometry/fallback;
}
```

This is particularly important for:

* `familyExtraMarginSize`;
* contact-history persistence;
* friction models;
* approaching contacts;
* contacts that have temporarily separated by a very small amount.

Trying to force the Feng formulation to represent a non-intersecting surface pair would destroy one of the clean theoretical properties we are trying to gain.

---

# 16. Coplanar and nearly coplanar triangles

This is a major edge case.

Feng's presentation states that ordinary intersecting non-coplanar triangle pairs produce two intersection points; coplanar triangles are excluded from the ordinary contact-segment construction. 

Real meshes can nevertheless generate:

* nearly coplanar neighboring surfaces;
* flat-on-flat contacts;
* identical or nearly parallel facets.

If

$$
\|\mathbf n_A\times\mathbf n_B\|
$$

is very small, segment orientation becomes ill-conditioned.

Recommended behavior:

```cpp
if (norm(cross(nA,nB)) < angularTolerance) {
    route to fallback projection geometry;
}
```

Do not try to force an intersection-segment formulation into the degenerate case in the first implementation.

Later, a specialized coplanar overlap algorithm could be added.

---

# 17. Degenerate mesh conditions

The implementation should explicitly handle:

* zero-area triangles;
* reversed winding;
* inconsistent mesh orientation;
* non-manifold surfaces;
* open surfaces;
* duplicated triangles;
* touching exactly at one point;
* intersection only along a shared edge;
* extremely short intersection segments.

Feng's global cancellation properties assume an oriented surface representing a meaningful body boundary.

Therefore the new method will work best for **watertight consistently oriented meshes**.

DEME does not necessarily need to reject other meshes outright, but should have fallback behavior.

Suggested policy:

```text
valid well-oriented nondegenerate pair
    -> Feng

degenerate / ambiguous pair
    -> existing projection method
```

This allows incremental deployment.

---

# 18. Patch structure should remain

One might ask whether Feng's additivity makes DEME's patch system unnecessary.

For now, **do not remove the patch system**.

Reasons:

1. DEME already uses patches for contact-history and force organization.
2. Multiple disconnected contact regions can exist between the same two concave bodies.
3. Friction histories should not necessarily be merged across unrelated regions.
4. The patch abstraction is valuable even if the underlying geometry becomes more globally additive.

Thus the natural new hierarchy is:

```text
triangle intersection segments
        ↓
existing/contact-island patch assignment
        ↓
Feng reduction within each patch
```

rather than:

```text
all segments between two complete bodies
        ↓
one global contact
```

The latter could be investigated later as a separate mode.

---

# 19. Proposed data model

A useful temporary record might be:

```cpp
struct TriTriFengContribution {
    PatchPairID patchPair;

    double3 S;
    double3 G;

    double penetration;

    // optional:
    double penetrationWeight;
    double3 oldContactPoint;
};
```

In practice DEME probably should retain Structure-of-Arrays layout rather than literally instantiate this struct:

```cpp
patchPairID[]
Sx[]
Sy[]
Sz[]
Gx[]
Gy[]
Gz[]
depth[]
weight[]
```

because this will better match existing GPU/CUB workflows.

---

# 20. Suggested new utility functions

Possible API-level device helpers:

```cpp
template <typename T1, typename T2>
__device__ bool triTriSurfaceIntersectionSegment(
    const T1& A0,
    const T1& A1,
    const T1& A2,
    const T1& B0,
    const T1& B1,
    const T1& B2,
    T1& x0,
    T1& x1
);
```

```cpp
template <typename T1, typename T2>
__device__ bool calcFengTriTriContribution(
    ...,
    T1& S,
    T1& G,
    T2& penetration
);
```

and a finalization routine:

```cpp
__device__ bool finalizeFengPatchGeometry(
    const double3& S,
    const double3& G,
    double penetration,
    float3& normal,
    double& area,
    double3& contactPoint
);
```

The names are only suggestions.

---

# 21. Suggested integration point

Do **not** initially modify the downstream force kernel.

The best migration target is to continue producing exactly:

```cpp
finalAreas
finalNormals
finalPenetrations
finalContactPoints
```

because `DEMCalcForceKernels_PatchBased.cu` already consumes those quantities cleanly.

Therefore the implementation should be localized mostly to:

* triangle–triangle geometric contact routines;
* intermediate contact arrays;
* patch voting/reduction kernels;
* corresponding host-side allocations;
* validation tests.

This reduces risk substantially.

---

# 22. Numerical precision

Recommended first implementation:

### Triangle intersection

Use FP32 where safe, consistent with current DEME GPU narrow phase, but include mixed/double fallback near numerical ambiguity if needed.

### Feng accumulation

Use:

```cpp
double3 S;
double3 G;
```

for patch reductions.

The individual triangle operations are cheap enough that accumulation accuracy is probably more important than saving a few arithmetic cycles.

This matters especially because vector-area cancellation between many segments is an essential feature of the theory.

A configuration involving a large contact patch may contain many contributions whose transverse components cancel, leaving a significantly smaller resultant.

FP32 reduction could degrade precisely the global geometric advantage being sought.

Optimization to FP32/FP64 mixed accumulation can follow profiling.

---

# 23. Expected GPU-performance effects

The new geometry potentially removes a substantial amount of local work.

The current projection routine can involve:

* signed-distance calculations;
* projection;
* construction of temporary polygons;
* Sutherland–Hodgman clipping;
* duplicate-vertex checking;
* polygon reordering;
* centroid computation;
* polygon-area calculation;
* evaluation in both projection directions.

The Feng narrow phase needs roughly:

```text
triangle–triangle surface intersection
orientation
cross product
dot products
small fixed arithmetic
```

plus reduction.

So its local arithmetic and register footprint may be considerably smaller.

However, **do not assume an overall performance improvement until benchmarked**.

The actual cost may be controlled by:

* candidate generation;
* number of triangle pairs;
* patch sorting;
* CUB reductions;
* memory traffic;
* contact-history processing.

The paper itself identifies finding intersecting triangle pairs as the dominant cost in its surface-mesh formulation. 

---

# 24. Development plan

I recommend implementation in four stages.

## Stage 1 — Add Feng geometry in diagnostic mode

Do not affect forces.

For each current DEME mesh–mesh contact:

1. run existing DEME geometry;
2. additionally compute Feng segment contributions;
3. reduce them by patch;
4. write diagnostic values:

   * `fengArea`,
   * `fengNormal`,
   * `fengContactPoint`;
5. compare with existing values.

This makes debugging much easier.

---

## Stage 2 — Replace area and normal

Use:

$$
A_F,
\qquad
\mathbf n_F
$$

while retaining:

$$
\delta_{\rm old},
\qquad
\mathbf x_{c,\rm old}.
$$

This isolates errors in the two easiest Feng quantities.

---

## Stage 3 — Replace contact point

Switch to

$$
\mathbf x_{c,F}
$$

or the projected-nearest-on-Feng-contact-line hybrid described above.

Now:

$$
(A,\mathbf n,\mathbf x_c)
$$

all come from Feng geometry.

Penetration remains legacy DEME.

---

## Stage 4 — Optimize penetration path

Once behavior is validated, remove unnecessary clipping work and build a depth-only routine.

At this point the expensive old projection polygon can potentially disappear entirely for ordinary nondegenerate contacts.

---

# 25. Required validation tests

A good implementation should not be judged only by "simulation does not crash."

We need geometric and mechanical tests.

## Test 1: two intersecting planar triangular surfaces

Known intersection segment.

Validate:

* segment endpoints;
* orientation;
* \(\mathbf S_i\);
* \(\mathbf G_i\).

---

## Test 2: two closed cubes

Use a simple offset/intersection.

Compare resulting:

$$
A,\mathbf n,\mathbf x_c
$$

against analytically expected overlap-boundary geometry.

---

## Test 3: triangulation invariance

Represent the same cube face with:

* 2 triangles;
* 8 triangles;
* 32 triangles.

Keep the physical contact state identical.

The Feng-derived

$$
A,\mathbf n,\mathbf x_c
$$

should remain approximately invariant.

This is one of the most important tests.

---

## Test 4: mesh-diagonal flip

Triangulate the same quad with opposite diagonals.

Current triangle-local approaches can change local contact decomposition.

The new method should give essentially identical patch geometry.

---

## Test 5: refined sphere meshes

Contact two approximated spheres using successively refined surface meshes.

Check convergence and sensitivity of:

* area;
* normal;
* force;
* torque.

---

## Test 6: strongly oblique contact

Use two surfaces approaching nearly perpendicular orientation.

The current implementation specifically contains logic intended to stabilize this situation.

Ensure the Feng method remains well behaved.

---

## Test 7: flat-on-flat contact

Nearly coplanar triangles.

Confirm fallback behavior.

---

## Test 8: two disconnected contact patches

A concave object contacting another body in two separate places.

Verify that patch identities remain separate and Feng reductions do not merge them accidentally.

---

## Test 9: energy / elastic impact

Even though arbitrary DEME force laws invalidate Feng's strict conservation proof, compare:

* old DEME3;
* hybrid geometry;

for a frictionless elastic impact.

Check total translational + rotational energy drift.

Feng's paper validates the complete native model on complex triangulated impacts, including large penetrations and irregular concave shapes. 

---

## Test 10: performance

Benchmark:

```text
number candidate tri pairs
number true intersections
narrow-phase time
patch aggregation time
force time
total timestep time
register usage
occupancy
GPU memory traffic
```

Compare old vs new for:

* coarse particles;
* high-resolution particles;
* large particle counts;
* dense packing.

---

# 26. Acceptance criteria

The implementation should be considered successful if:

1. existing non-mesh DEME contacts remain unchanged;
2. mesh–mesh contacts still produce the same downstream DEME contact interface;
3. simple analytic tests give correct vector area and normal;
4. mesh refinement / diagonal changes produce substantially smaller geometric variation than the current implementation;
5. force and torque remain stable over contact-patch transitions;
6. margin-only and degenerate cases fall back safely;
7. no material regression occurs in large GPU simulations;
8. ideally, triangle narrow-phase cost decreases.

---

# 27. What should *not* be done in the first implementation

Avoid expanding scope unnecessarily.

Do not initially:

* replace DEME's spatial binning;
* replace its patch/contact-island architecture;
* redesign force models;
* impose Feng's \(F=kA\) law on all contacts;
* remove penetration depth;
* implement Feng's nonlinear volume-energy version;
* explicitly reconstruct closed intersection loops;
* attempt a general exact coplanar polygon formulation;
* remove old geometry before the new path has side-by-side validation.

The goal is a **controlled replacement of mesh–mesh geometry derivation**, not a rewrite of DEM-Engine.

---

# 28. Core conceptual summary for the implementer

The current DEME3 concept is approximately:

$$
\boxed{
\text{triangle pair}
\rightarrow
(A_i,\mathbf n_i,\delta_i,\mathbf x_i)
\rightarrow
\text{patch voting}
}
$$

where much of the final geometry is constructed from triangle-local projected contacts.

The proposed concept is:

$$
\boxed{
\text{triangle pair}
\rightarrow
\text{oriented intersection segment}
\rightarrow
(\mathbf S_i,\mathbf G_i)
}
$$

followed by

$$
\boxed{
\mathbf S=\sum_i\mathbf S_i,
\qquad
\mathbf G=\sum_i\mathbf G_i
}
$$

and then

$$
\boxed{
A=\|\mathbf S\|,
\qquad
\mathbf n=\frac{\mathbf S}{\|\mathbf S\|},
\qquad
\mathbf x_c
\leftarrow
(\mathbf S,\mathbf G)
}
$$

while obtaining

$$
\boxed{\delta}
$$

from DEME's existing penetration machinery.

Thus the intended final architecture is:

$$
\boxed{
\text{DEME broad phase}
+
\text{Feng surface-intersection geometry}
+
\text{DEME penetration}
+
\text{DEME patches}
+
\text{DEME user force law}.
}
$$

The fundamental implementation insight from Feng is that **the intersection segments need not be assembled into a closed polyline**. Their vector-area and geometric-moment contributions are additive and can be computed independently, which is exactly the property that makes this formulation attractive for GPU DEM. 

I would have the implementation agent start specifically with **Stage 1 diagnostic mode**. That gives us a very clean way to compare Feng's \(A,n,x_c\) against DEME3's existing patch outputs on exactly the same candidate/contact topology before allowing the new geometry to influence forces.
