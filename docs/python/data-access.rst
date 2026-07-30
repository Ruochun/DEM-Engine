Retrieving simulation data
==========================

Trackers expose both ordinary Python-returning methods and zero-host-copy-style
methods that fill caller-owned CUDA memory.

Host access
-----------

Use host methods for inspection, logging, and small result sets:

.. code-block:: python

   solver.DoDynamicsThenSync(0.01)
   position = tracker.Pos()          # first tracked owner
   positions = tracker.Positions()  # every owner in this tracker
   velocities = tracker.Velocities()

Most single-object methods accept an optional owner offset. Bulk methods return
all owners covered by the tracker. Quaternion values use public ordering
``(x, y, z, w)``.

Direct CUDA access
------------------

The ``*ToDevice`` methods synchronously fill a raw CUDA pointer represented by
a Python integer. This avoids first materializing the result as a Python list.
The allocation remains entirely owned by the caller.

This CuPy example retrieves all tracked positions:

.. literalinclude:: examples/cupy_device_retrieval.py
   :language: python
   :linenos:

The direct interface is intentionally low level. Every call must satisfy all
of these requirements:

* ``pointer`` is a valid, writable CUDA-accessible pointer on ``device``;
* the allocation remains alive until the synchronous call returns;
* ``capacity`` is a number of result elements, not bytes;
* dtype and row layout match the method;
* the current simulation state has been synchronized before it is read.

Supplying the wrong pointer, device, dtype, layout, or capacity can cause a
native CUDA error or memory corruption. Prefer host access unless eliminating
the host transfer materially helps the application.

Tracked-owner methods
---------------------

The following table gives the destination representation. A ``float3`` result
maps to a C-contiguous ``(capacity, 3)`` array of ``float32``; ``float4`` maps
to ``(capacity, 4)``.

.. list-table::
   :header-rows: 1

   * - Method
     - Destination element
     - Meaning
   * - ``PositionsToDevice``
     - ``float3``
     - Global position
   * - ``VelocitiesToDevice``
     - ``float3``
     - Global linear velocity
   * - ``AngularVelocitiesLocalToDevice``
     - ``float3``
     - Local angular velocity
   * - ``AngularVelocitiesGlobalToDevice``
     - ``float3``
     - Global angular velocity
   * - ``OrientationQuaternionsToDevice``
     - ``float4``
     - Quaternion ``(x, y, z, w)``
   * - ``FamiliesToDevice``
     - ``uint32``
     - Family number
   * - ``MassesToDevice``
     - ``float32``
     - Mass
   * - ``MOIsToDevice``
     - ``float3``
     - Principal moments of inertia
   * - ``ContactAccelerationsToDevice``
     - ``float3``
     - Global contact acceleration
   * - ``ContactAngularAccelerationsLocalToDevice``
     - ``float3``
     - Local contact angular acceleration
   * - ``ContactAngularAccelerationsGlobalToDevice``
     - ``float3``
     - Global contact angular acceleration
   * - ``OwnerWildcardValuesToDevice``
     - ``float32``
     - Named owner wildcard

The required owner capacity is the number of owners covered by the tracker.
In the common case where a tracker was created from a clump batch, this is the
number of clumps added in that batch.

Contact-pair methods
--------------------

Python exposes three all-owner compacting methods:

* ``GetContactForcesForAllToDevice(points, forces, capacity, device)``
* ``GetContactForcesAndLocalTorqueForAllToDevice(points, forces, torques, capacity, device)``
* ``GetContactForcesAndGlobalTorqueForAllToDevice(points, forces, torques, capacity, device)``

Every destination is a C-contiguous ``(capacity, 3)`` ``float32`` allocation.
The methods return the number of valid compacted contact pairs; only rows
before that count contain results. Capacity must cover the simulation's total
recorded contact count, not merely the final number involving this tracker.
Contact retrieval also requires force recording to remain enabled.

Destination device
------------------

``device`` names the logical CUDA device that owns the destination allocation;
it does not have to be inferred from the pointer. With CuPy, enter
``with cp.cuda.Device(device):`` before allocating, pass ``array.data.ptr``,
and keep the array in scope through the call. DEME handles retrieval when the
destination differs from a worker device, subject to the CUDA capabilities of
the system.
