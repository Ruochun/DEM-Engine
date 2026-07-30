Quickstart
==========

Python
------

The :doc:`python/quickstart` contains a complete runnable sphere simulation.
The Python guide also covers the setup/runtime boundary, multi-GPU device
selection, tracker access, and direct retrieval into CUDA arrays.

C++
---

Include the public API and construct a solver:

.. code-block:: cpp

   #include "DEM/API.h"

   deme::DEMSolver solver;
   solver.InstructBoxDomainDimension(1.0f, 1.0f, 1.0f);
   solver.SetGravitationalAcceleration(make_float3(0.0f, 0.0f, -9.81f));

See :doc:`cpp-api/index` for the generated API reference.
