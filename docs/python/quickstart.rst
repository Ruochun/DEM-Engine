Python quickstart
=================

Install DEME in an environment with a compatible NVIDIA driver and CUDA
runtime, then verify the package and version:

.. code-block:: console

   python -c "import deme; print(deme.__version__)"

The following complete example creates a material, a bounded domain, and one
spherical clump, then advances the simulation:

.. literalinclude:: examples/sphere_drop.py
   :language: python
   :linenos:

Run it with:

.. code-block:: console

   python docs/python/examples/sphere_drop.py

The solver constructor initializes CUDA worker resources, so even importing
successfully is not sufficient to run a simulation without a visible,
supported NVIDIA GPU. See :doc:`device-selection` when the process can see more
than one GPU.

The example uses ``DoDynamicsThenSync`` because the position is read
immediately afterward. For longer simulations, asynchronous ``DoDynamics``
calls can overlap host work; synchronize before reading results or exiting.

Where to go next
----------------

* :doc:`solver-lifecycle` explains which operations belong before and after
  ``Initialize()``.
* :doc:`data-access` covers host-returning tracker methods and direct retrieval
  into CUDA arrays.
* :doc:`api-overview` maps common tasks to the Python objects and methods that
  implement them.
