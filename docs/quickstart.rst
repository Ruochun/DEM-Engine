Quickstart
==========

Python
------

The package follows the same setup flow as the C++ API:

.. code-block:: python

   import deme

   solver = deme.DEMSolver()
   solver.InstructBoxDomainDimension(1.0, 1.0, 1.0)
   solver.SetGravitationalAcceleration([0.0, 0.0, -9.81])

Material, particle, boundary, timestep, and output choices must be supplied
before ``Initialize()``. Once initialized, advance the simulation with
``DoDynamics()`` or the corresponding duration-based helper.

The complete runnable Python examples are still being brought into this branch.
Until then, use the C++ demos as the authoritative examples for simulation
setup and consult the :doc:`python/index` for Python-specific conventions.

C++
---

Include the public API and construct a solver:

.. code-block:: cpp

   #include "DEM/API.h"

   deme::DEMSolver solver;
   solver.InstructBoxDomainDimension(1.0f, 1.0f, 1.0f);
   solver.SetGravitationalAcceleration(make_float3(0.0f, 0.0f, -9.81f));

See :doc:`cpp-api/index` for the generated API reference.
