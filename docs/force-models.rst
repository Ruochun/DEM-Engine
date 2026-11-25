Contact force models
====================

This page specifies the built-in frictional and frictionless Hertzian models,
including their numerical approximations. The executable definitions are
``src/kernel/DEMCustomizablePolicies/FullHertzianForceModel.cu`` and
``FrictionlessHertzianForceModel.cu``; the shared rotational law is
``contactResistanceTorque`` in ``src/kernel/DEMHelperKernels.cuh``.
User scripts for cohesion, electrostatics, fracture and mooring add their own
constitutive terms and are examples, not additional built-in Hertz models.

Frames and force-model outputs
------------------------------

Let :math:`\mathbf n` be ``B2A``, the unit contact normal pointing from B to A,
and :math:`\mathbf r_A,\mathbf r_B` the global vectors from each owner's center
of mass to the contact point. ``locCPA`` and ``locCPB`` supply these vectors in
owner-local coordinates. ``ARotVel`` and ``BRotVel`` are also owner-local;
rotate them with ``AOriQ`` and ``BOriQ`` before taking their difference.

All five writable output vectors start at zero for every contact evaluation.
``force`` is the conventional equal-and-opposite contact force. The four
additional outputs are additive and expressed in global coordinates:

.. math::

   \mathbf F_A &= \mathbf{force}+\mathbf{A\_only\_force},\\
   \mathbf F_B &= -\mathbf{force}+\mathbf{B\_only\_force},\\
   \boldsymbol\tau_A &= \mathbf r_A\times\mathbf F_A+\mathbf{A\_only\_torque},\\
   \boldsymbol\tau_B &= \mathbf r_B\times\mathbf F_B+\mathbf{B\_only\_torque}.

Forces have units N and direct torques have units N m. Independent forces act
at the contact point, not at the center of mass. Custom models are responsible
for conservation when using independent outputs. Built-in rotational resistance
always uses equal-and-opposite direct torques.

Force-model ingredient reference
--------------------------------

A custom model is a CUDA/C++ code fragment supplied through
``DefineCustomModel`` or ``ReadCustomModelFile``. It executes once per accepted
contact pair: a sphere pair, a sphere--analytical pair, or an aggregated mesh
patch pair. Mesh primitive pairs supply geometry to aggregation; they do not
independently execute the patch force law. A and B identify the two contact
sides, not necessarily two particles (one or both can be boundaries).

The tables below describe the named inputs available at model evaluation.
Treat these as inputs rather than overwriting them. Only the five output vectors
and explicitly registered wildcards are intended for model updates. Units below
assume SI; another consistent unit system can be used throughout the simulation.

Geometry and contact inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Ingredient
     - Type / units
     - Meaning and use
   * - ``ContactType``
     - ``deme::contact_t``
     - Contact category: ``deme::SPHERE_SPHERE_CONTACT``, ``SPHERE_TRIANGLE_CONTACT``,
       ``SPHERE_ANALYTICAL_CONTACT``, ``TRIANGLE_TRIANGLE_CONTACT``, or
       ``TRIANGLE_ANALYTICAL_CONTACT`` (all in namespace ``deme``). Mesh categories
       describe the final patch contact. Rejected ``NOT_A_CONTACT`` pairs skip the model.
   * - ``B2A``
     - ``float3`` / dimensionless
     - Global unit contact normal pointing from B toward A. Positive normal force
       on A acts along this vector.
   * - ``overlapDepth``
     - ``double`` / m
     - Signed penetration, positive for overlap. Extra contact margins can admit
       zero or negative penetration; the model must decide whether to apply force.
       Built-in Hertz models require positive penetration.
   * - ``overlapArea``
     - ``double`` / m²
     - Geometric contact area; mesh contacts use the projected patch aggregate.
       The built-in contact radius is ``sqrt(overlapArea / deme::PI)``.
   * - ``contactPnt``
     - ``double3`` / m
     - Global contact point, shared by both sides. For mesh patches this is the
       aggregated point, not necessarily an individual triangle's contact point.
   * - ``locCPA``, ``locCPB``
     - ``float3`` / m
     - Vectors from each owner's center to the contact point, expressed in that
       owner's local frame. Rotate with ``AOriQ`` / ``BOriQ`` for global moment arms.
   * - ``ARadius``, ``BRadius``
     - ``float`` / m
     - Actual component radius for spheres; ``DEME_HUGE_FLOAT`` placeholders for
       meshes and analytical objects. These are **not** the mesh intrinsic radii
       configured with ``SetContactRadius``. Use ``rollingResistanceLength`` for RR.
   * - ``rollingResistanceLength``
     - ``const float`` / m
     - Producer-resolved rolling lever arm, including enabled intrinsic/absolute
       bounds. Finite and nonnegative when resolved; negative when requested bounds
       lack a usable scale. A model must reject a negative value when it actually
       applies RR. Nonfinite values are rejected before model evaluation. See the
       rolling section below for the activation conditions and bound controls.

Owner state and identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Ingredient
     - Type / units
     - Meaning and use
   * - ``AOwnerPos``, ``BOwnerPos``
     - ``double3`` / m
     - Global owner positions (the reference centers used for contact moment arms).
   * - ``bodyAPos``, ``bodyBPos``
     - ``double3`` / m
     - Global component reference locations: sphere centers, mesh patch centers,
       or analytical-entity reference locations. These are not generally the
       contact point or the owner center; a mesh patch center need not lie on its surface.
   * - ``AOriQ``, ``BOriQ``
     - ``float4`` / dimensionless
     - Owner orientation quaternions, stored as ``(x, y, z, w)``. Applying
       ``applyOriQToVector3(vector, quaternion)`` rotates a vector from owner-local
       to global coordinates in place.
   * - ``AOwnerMass``, ``BOwnerMass``
     - ``float`` / kg
     - Owner masses, not individual sphere-component masses. Boundary values follow
       the solver's configured mass conventions; do not infer mobility from mass alone.
   * - ``ALinVel``, ``BLinVel``
     - ``float3`` / m/s
     - Global owner translational velocities. These are not contact-point velocities;
       add the angular-velocity cross moment-arm contribution.
   * - ``ARotVel``, ``BRotVel``
     - ``float3`` / rad/s
     - Owner angular velocities in their respective **local** frames. Rotate both
       to global before subtracting them or taking a cross product with a global vector.
   * - ``AOwnerMOI``, ``BOwnerMOI``
     - ``float3`` / kg m²
     - Diagonal owner moments of inertia in the owner's local principal frame.
   * - ``AOwner``, ``BOwner``
     - ``deme::bodyID_t``
     - Owner indices, suitable for indexing owner wildcard arrays. Several sphere
       components or mesh patches can share an owner.
   * - ``AGeo``, ``BGeo``
     - ``deme::bodyID_t``
     - Geometry indices: sphere component ID, mesh **patch** ID, or analytical
       component ID. These are not interchangeable with owner IDs or mesh triangle IDs.
   * - ``AOwnerFamily``, ``BOwnerFamily``
     - ``deme::family_t``
     - Family numbers of the two owners, usable for family-dependent model rules.
   * - ``bodyAMatType``, ``bodyBMatType``
     - ``deme::materialsOffset_t``
     - Material-table indices for the contacting components/patches.
   * - ``ts``
     - ``float`` / s
     - Current dynamics timestep, used to integrate contact history.
   * - ``time``
     - ``float`` / s
     - Current elapsed simulation time at force evaluation.

Velocity, MOI, owner/geometry IDs, ``ts`` and ``time`` are acquired on demand by
scanning the model text for their names. Use the exact identifiers in the model
body; do not rely on references hidden only inside a prerequisite function or a
macro expansion to request an optional ingredient. Prerequisite device functions
should take per-contact inputs as arguments. ``granData`` and ``simParams`` are
kernel implementation details, not needed to consume the named interface.

Writable outputs
~~~~~~~~~~~~~~~~

All five are ``float3`` vectors in the global frame, initialized to zero at each
evaluation. They are additive as defined by the equations above:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Output
     - Units
     - Effect
   * - ``force``
     - N
     - Adds this force to A and its negative to B, at the shared contact point.
   * - ``A_only_force``, ``B_only_force``
     - N
     - Adds force only to the named side, at the contact point. Its moment about
       that owner's center is included automatically.
   * - ``A_only_torque``, ``B_only_torque``
     - N m
     - Adds a direct torque only to the named owner; no cross product or additional
       force is generated. Assign opposite torques for a pure resisting couple.

Do not include the moment of ``force`` or the side-only forces again in the
side-only torque outputs. ``torque_only_force`` is unsupported and produces a
setup-time deprecation error; it is not an alias for any of these outputs.

Material tables and persistent history
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Register material names with ``SetMustHaveMatProp`` and identify pairwise names
with ``SetMustPairwiseMatProp``. Registered properties become named float tables:
``E[bodyAMatType]`` reads a per-material value, while
``Crr[bodyAMatType][bodyBMatType]`` reads a pairwise value. Names such as ``E``,
``nu``, ``CoR``, ``mu``, ``Crr`` and ``Cspin`` belong to the chosen model, not to
a universally predeclared set of scalar inputs. Derived locals such as
``E_cnt``, ``Cspin_cnt``, ``cnt_rad`` and ``mass_eff`` are computed by the force
script itself; they are not supplied ingredients.

``SetPerContactWildcards`` registers named ``float`` history scalars. Each name
is loaded for the current contact, can be modified directly, and is written back
after evaluation. New contact history starts at zero; retained pairs carry their
history across contact detection updates. Rejected contacts have their history
cleared. For margin-admitted contacts with nonpositive penetration, clearing
history is the model's responsibility (the built-in Hertz model does so).

``SetPerOwnerWildcards`` instead provides ``float*`` arrays. For a registered
name ``state``, ``state``, ``state_A`` and ``state_B`` all alias the same array;
use ``state[AOwner]`` and ``state[BOwner]`` to access the owners. Multiple contacts
can update one owner concurrently, so use appropriate atomic operations or avoid
concurrent writes. Do not assume a global ordering among contact evaluations.
Geometry-wildcard acquisition is not enabled in the current force-model path.

Geometry and effective material properties
------------------------------------------

Contact area uses geometric sphere-intersection formulas and mesh polygon projection/reduction.
For positive penetration :math:`\delta`, the equivalent-area contact radius and effective mass are

.. math::

   a=\sqrt{A/\pi},\qquad m^*=\frac{m_A m_B}{m_A+m_B}.

Intrinsic curvature is :math:`\kappa_i=1/R_i` for a finite radius and zero for a wall.
The scale :math:`R^*=1/(\kappa_A+\kappa_B)` bounds rolling resistance; it does not replace
contact area. Sphere components supply their actual radii and analytical boundaries supply
zero curvature. Mesh volume-equivalent radius is a global size proxy, not local curvature.
Owner masses retain the existing prescribed/fixed-owner conventions.

Mesh radius lifecycle
~~~~~~~~~~~~~~~~~~~~~

``mesh->SetContactRadius(R)`` explicitly overrides the intrinsic radius;
``R=0`` explicitly selects a wall. The implementation stores zero curvature,
not a huge radius in reciprocal arithmetic. Negative and nonfinite inputs are
rejected. ``ClearContactRadius()`` restores automatic selection at the next
``Initialize()`` or full ``Update()``. ``GetContactRadius()`` returns the last
resolved value, with zero denoting a wall/unresolved radius.

With ``SetMeshUniversalContact(true)``, initialization and full updates infer
:math:`R=(3V/(4\pi))^{1/3}` from watertight, nonshell solids unless overridden.
Open meshes, shells, and invalid/nonpositive volumes retain zero curvature.
With universal contact disabled, meshes default to walls and no automatic
volume-radius calculation is performed. Explicit overrides win in both cases.
A closed container should therefore explicitly use ``SetContactRadius(0)``
when universal contact is enabled.

Template overrides are copied into later instances; existing instances keep
their own settings. Explicit radii are fixed length settings: scale the override
explicitly if you scale an overridden mesh. ``ComputeMassProperties`` remains a const query with no
radius side effects. Lightweight host/device mesh deformation leaves the
resolved radius unchanged. A synchronized full ``Update()`` reads the latest
geometry and refreshes automatic radii, while preserving explicit overrides.
Rolling bounds also take effect at this full boundary;
``UpdateSimParams()`` alone does not refresh them.

Patch geometry and diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mesh patches use area-weighted normal voting, maximum projected penetration, and
contact-point weights proportional to projected area times penetration. These are geometric
approximations, not a continuum contact-manifold solver. In flat mesh contact, projected area
may remain finite as penetration tends to zero. The resulting :math:`A/(\pi\delta)` can be
unbounded. RR length bounds limit that rolling lever arm, but do not cap contact area,
normal/tangential stiffness, damping, or spinning torque. Area or normal discontinuities
can therefore still require a smaller timestep or improved geometry.

Per-triangle ``P`` output is an allocated normal **load**, not load divided by triangle area;
``PV`` multiplies allocated load by slip speed. Allocation uses projected geometric weights.
Ordinary contact output uses the force-model area where a custom model records ``overlapArea``.


.. math::

   \frac1{E^*}&=\frac{1-\nu_A^2}{E_A}+\frac{1-\nu_B^2}{E_B},\\
   \frac1{G^*}&=\frac{2(2-\nu_A)(1+\nu_A)}{E_A}
               +\frac{2(2-\nu_B)(1+\nu_B)}{E_B}.

``E`` and ``nu`` belong to each material. ``CoR``, ``mu``, ``Crr`` and ``Cspin``
are pairwise properties. An explicitly supplied pair value overrides the
arithmetic mean of the two material values. Missing ``Cspin`` values default
to zero with the generic missing-material-property warning. The built-in law expects ``Crr`` and ``Cspin`` to be finite
and nonnegative. Both are dimensionless, but they use different length scales.

Normal force (both built-in models)
-----------------------------------

The relative contact-point velocity is

.. math::

   \mathbf v=\mathbf v_A+\boldsymbol\omega_A\times\mathbf r_A
              -\mathbf v_B-\boldsymbol\omega_B\times\mathbf r_B,
   \quad v_n=\mathbf v\cdot\mathbf n,
   \quad\mathbf v_t=\mathbf v-v_n\mathbf n.

With :math:`e=\max(\mathrm{CoR},\mathrm{DEME\_TINY\_FLOAT})`, the implementation uses

.. math::

   \beta&=\frac{\ln e}{\sqrt{(\ln e)^2+\pi^2}},\quad S_n=2E^*a,\\
   k_n&=\frac23 S_n,\quad
   \gamma_n=2\sqrt{5/6}\,\beta\sqrt{S_nm^*},\\
   f_n&=k_n\delta+\gamma_n v_n,\quad\mathbf F_n=f_n\mathbf n.

For the usual :math:`0<e\le1`, :math:`\gamma_n\le0`: damping opposes relative
normal motion. The existing normal law does **not** clamp :math:`f_n` on
unloading; a negative value is possible. Rotational resistance uses only the
compressive load :math:`N=\max(f_n,0)`. No built-in contact force is applied
when :math:`\delta\le0`.

Tangential sliding (frictional model only)
------------------------------------------

The stored displacement :math:`\boldsymbol\xi` comprises ``delta_tan_x/y/z``.
At each step of length :math:`h`, update and project it into the tangent plane:

.. math::

   \boldsymbol\xi'&=\boldsymbol\xi+h\mathbf v_t,\\
   \boldsymbol\xi&=\boldsymbol\xi'-(\boldsymbol\xi'\cdot\mathbf n)\mathbf n,\\
   k_t&=8G^*a,\quad g_t=-2\sqrt{5/3}\,\beta\sqrt{m^*k_t},\\
   \mathbf F_t^{\rm trial}&=-k_t\boldsymbol\xi-g_t\mathbf v_t.

For ``mu > 0``, limit the trial force magnitude to :math:`\mu|f_n|`. If clipped,
back-calculate :math:`\boldsymbol\xi=-(\mathbf F_t+g_t\mathbf v_t)/k_t`.
The implementation zeros a trial force whose magnitude is at most
``DEME_TINY_FLOAT``; back-calculation occurs only in the clipping branch.
The resultant is :math:`\mathbf{force}=\mathbf F_n+\mathbf F_t`.
For ``mu <= 0``, no tangential force is added, although the displacement update
still occurs. History is cleared for nonpenetrating contacts or lost pairs.
This existing projection does not explicitly preserve displacement magnitude
when the normal changes. Its friction cap uses :math:`|f_n|`, not the clamped
load used for rotational resistance.

Rolling and spinning (frictional model only)
--------------------------------------------

Resolve global relative angular velocity into normal and tangential components:

.. math::

   \Delta\boldsymbol\omega&=\boldsymbol\omega_A-\boldsymbol\omega_B,\\
   \boldsymbol\omega_s&=(\Delta\boldsymbol\omega\cdot\mathbf n)\mathbf n,\\
   \boldsymbol\omega_r&=\Delta\boldsymbol\omega-\boldsymbol\omega_s.

The kinetic moment scales are

.. math::

   M_r=C_{rr}N L_r,\qquad
   M_s=\frac{3\pi}{16}C_{\rm spin}N a.

The raw rolling length is :math:`A/(\pi\delta)`.
``SetRollingResistanceLengthLimit(multiplier, absolute_max)`` limits it to
``multiplier * R*`` when an intrinsic scale exists, and additionally to a positive
``absolute_max`` when supplied. The default multiplier is **1**; zero (default)
disables only the absolute cap. Parameters must be finite, with positive multiplier
and nonnegative absolute maximum. With bounds enabled, an active RR contact between
two zero-curvature objects requires an absolute cap or a finite intrinsic radius.

``SetRollingResistanceLengthBoundsEnabled(false)`` disables both bounds and exposes
the raw geometry-derived length, including very large finite values. This is not a
claim that the raw length is physically accurate. Non-finite results still error.
Apply changes with a full ``Update()`` after initialization.

The producer exposes ``rollingResistanceLength`` as a finite nonnegative length,
or a negative value when requested bounds cannot be resolved. It rejects NaN/infinity
before force-model evaluation. A force model using RR must reject a negative ingredient
when positive overlap, positive compressive load and positive ``Crr`` activate RR.
Spin-only contacts need no rolling length. The bundled force models implement this
check and pass the explicit length as the required last argument of
``contactResistanceTorque<float>``. The helper performs no implicit length resolution.

The first is an empirical rolling-resistance moment. The second follows from
integrating fully mobilized circumferential Coulomb traction over a circular
Hertz pressure distribution:

.. math::

   p(r)=p_0\sqrt{1-r^2/a^2},\quad
   N=\frac{2\pi}{3}p_0a^2,\quad
   M_s=2\pi C_{\rm spin}\int_0^a p(r)r^2\,dr.

``Cspin`` is an independently calibrated friction coefficient; it need not
match ``mu``. Using this pressure-profile factor with an equivalent-area
geometric radius is an approximation, especially for noncircular mesh patches.

With :math:`\omega_0=10^{-3}` rad/s, the implemented torque on A is

.. math::

   \mathbf T_A=-M_r\frac{\boldsymbol\omega_r}
                               {\max(|\boldsymbol\omega_r|,\omega_0)}
                -M_s\frac{\boldsymbol\omega_s}
                               {\max(|\boldsymbol\omega_s|,\omega_0)},
   \qquad\mathbf T_B=-\mathbf T_A.

These are added to ``A_only_torque`` and ``B_only_torque``. Below the
regularization speed the torque decreases linearly to zero. There is no stored
rotational spring, static moment, contact-age gate, or rotational history.
The law is active immediately for positive overlap and compressive load.
The frictionless model supplies neither rolling nor spinning resistance.

The added torque pair conserves angular momentum and has nonpositive
continuous-time power:

.. math::

   P=\mathbf T_A\cdot\Delta\boldsymbol\omega
     =-M_r\frac{|\boldsymbol\omega_r|^2}{\max(|\boldsymbol\omega_r|,\omega_0)}
      -M_s\frac{|\boldsymbol\omega_s|^2}{\max(|\boldsymbol\omega_s|,\omega_0)}\le0.

It vanishes for common rigid-body rotation. This energy property does not
ensure stability for an arbitrary integration timestep. Near rest, the
regularization acts like angular damping; strong resistance or small inertia
can require a shorter timestep. There is no per-contact stopping cap, and
multiple contacts contribute simultaneously.

Sliding and spinning are evaluated independently. Their independent limits do
not enforce a shared contact-patch traction budget. Partial torsional slip,
static rotational resistance, and coupling between sliding and spin are not
modeled. The rolling coefficient is empirical and must be recalibrated when
migrating from the former rolling law.

Migration and reporting
-----------------------

``torque_only_force`` is removed from the force-model interface. An identifier
reference in model code or prerequisites raises a setup-time deprecation error;
mentions in comments and string literals are allowed. The former quantity had
units of force and generated a moment through a cross product. Assigning that
value directly to a new torque output is dimensionally wrong.

Bundled custom scripts retain their other constitutive terms and now use the
same direct rolling couple, with spinning disabled unless explicitly added.
The built-in Hertzian model no longer allocates ``delta_time`` history.

Contact records contain total ``ForceA`` and ``ForceB``, plus direct
extra ``TorqueA`` and ``TorqueB``, all global. Torques exclude
:math:`\mathbf r\times\mathbf F`; tracker total-wrench queries include it.
``GetForceA()``, ``GetForceB()``, ``GetTorqueA()`` and ``GetTorqueB()`` are
available in C++ and Python. For legacy compatibility, ``GetForce()`` and
``GetTorque()`` return A's arrays; generic ``Force`` and ``Torque`` keys alias
the same A storage rather than maintaining duplicate arrays.
Legacy tracker ``...ForAll`` contact-list queries still select one side per
contact (A when both owners are selected); use per-owner queries or the detailed
contact container to inspect both independent sides. CSV output
uses ``f_a_x/y/z``, ``f_b_x/y/z``, ``torque_a_x/y/z`` and ``torque_b_x/y/z``
when the corresponding output flags are enabled. These replace the former
``f_x/y/z``, ``torque_x/y/z``, ``ForceB_X/Y/Z`` and ``TorqueB_X/Y/Z`` headers;
CSV readers must migrate their column names; legacy CSV headers are no longer reserved.
Consumers must use the header rather than fixed column offsets.
Asynchronous writers own these fields as part of their completed snapshot.
Contact filtering considers both owners, including contacts with only a direct
torque. The legacy single threshold is applied separately to force and torque
magnitudes (in their respective units).

Tests
-----

``DEMTest_ContactArea`` checks all five projected contact types, both force
reduction paths, mixed-type scratch lifetimes, sphere--plane normalization,
radius lifecycle and configurable rolling bounds.
``--host-only`` checks radius policies and rolling bounds without CUDA initialization.

``DEMTest_SpinningCube`` checks spin-down on analytical and meshed planes against
zero-spinning-friction controls. Its ``--rr-error`` mode intentionally triggers
the unresolved rolling-bound GPU diagnostic. Run deliberate GPU-trap tests last,
especially on WSL where context teardown can block later GPU work.

``DEMTest_ContactResistance`` checks the constitutive limits and deprecation,
then exercises independent outputs, pure rolling/spinning, common rotation,
rotated owner frames, separation and re-contact history reset, both reduction
paths, mesh pressure from a B-only force, and no-record mode.
``--host-only`` runs the constitutive and deprecation checks without GPU setup.
``DEMTest_DeviceDataRetrieval`` checks asymmetric contact-wrench retrieval;
``DEMTest_ContactOutputSnapshot`` checks asynchronous contact output.

For a focused local run (including refreshing JIT source assets)::

   cmake --build build --target sync_kernels
   cmake --build build --target DEMTest_ContactResistance DEMTest_DeviceDataRetrieval DEMTest_ContactOutputSnapshot -- -j2
   build/bin/DEMTest_ContactResistance
   build/bin/DEMTest_DeviceDataRetrieval
   build/bin/DEMTest_ContactOutputSnapshot
