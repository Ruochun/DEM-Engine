DEM-Engine documentation
========================

DEM-Engine (DEME) is a GPU-accelerated discrete element method solver with a
Chrono-like C++ API and preview Python bindings distributed as ``deme3`` with the import name ``deme``.

This documentation covers installation, the simulation model, the Python and
C++ interfaces, and the constraints that matter when running JIT-compiled CUDA
kernels.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   quickstart
   concepts
   visualization
   troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: Interfaces

   python/index
   cpp-api/index

.. toctree::
   :maxdepth: 2
   :caption: Development

   developer/architecture
   developer/documentation
