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

Geometry and effective material properties
------------------------------------------

For positive penetration :math:`\delta` (``overlapDepth``) and geometric overlap
area :math:`A` (``overlapArea``), define

.. math::

   a=\sqrt{A/\pi},\qquad R_{\rm eff}=a^2/\delta,\qquad
   m^*=\frac{m_A m_B}{m_A+m_B}.

The radius :math:`a` is an equivalent-area radius. In particular, using actual
geometric overlap area is not identical to solving the elastic Hertz contact
patch. For meshes and other noncircular contacts these are geometric proxies,
not exact local-curvature or pressure-distribution calculations. Masses are
owner masses, including the existing prescribed/fixed-owner mass conventions.

.. math::

   \frac1{E^*}&=\frac{1-\nu_A^2}{E_A}+\frac{1-\nu_B^2}{E_B},\\
   \frac1{G^*}&=\frac{2(2-\nu_A)(1+\nu_A)}{E_A}
               +\frac{2(2-\nu_B)(1+\nu_B)}{E_B}.

``E`` and ``nu`` belong to each material. ``CoR``, ``mu``, ``Crr`` and ``Cspin``
are pairwise properties. An explicitly supplied pair value overrides the
arithmetic mean of the two material values. Missing ``Cspin`` values default
to zero without a warning. Built-in ``Crr`` and ``Cspin`` values must be finite
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

   M_r=C_{rr}N R_{\rm eff},\qquad
   M_s=\frac{3\pi}{16}C_{\rm spin}N a.

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

Contact records contain total ``Force`` on A and ``ForceB`` on B, plus direct
extra ``Torque`` on A and ``TorqueB`` on B, all global. Torques exclude
:math:`\mathbf r\times\mathbf F`; tracker total-wrench queries include it.
``GetForceB()`` and ``GetTorqueB()`` are available in C++ and Python.
Legacy tracker ``...ForAll`` contact-list queries still select one side per
contact (A when both owners are selected); use per-owner queries or the detailed
contact container to inspect both independent sides. CSV output
adds ``ForceB_X/Y/Z`` and ``TorqueB_X/Y/Z`` when the corresponding output flags
are enabled. Consumers must use the header rather than fixed column offsets.
Asynchronous writers own these fields as part of their completed snapshot.
Contact filtering considers both owners, including contacts with only a direct
torque. The legacy single threshold is applied separately to force and torque
magnitudes (in their respective units).

Tests
-----

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
