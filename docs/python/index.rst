Python interface
================

Import and compatibility
------------------------

Use the lowercase package in new code:

.. code-block:: python

   import deme

``import DEME`` re-exports the same API for backward compatibility. The
capitalized spelling should not be used in new examples.

Version
-------

The runtime version is available as:

.. code-block:: python

   assert deme.__version__ == "3.0.0"

API conventions
---------------

The current bindings intentionally retain the C++ method names, such as
``SetGravitationalAcceleration`` and ``DoDynamics``. Pythonic snake-case aliases
and type stubs are planned, but are not yet part of the stable interface.

CUDA vector types exposed by the API are generally accepted as three- or
four-element Python sequences and returned as tuples. Quaternion ordering is
``(x, y, z, w)``.

Objects returned by solver setup methods frequently own or reference native
simulation state. Keep the solver alive for as long as trackers, inspectors, or
other solver-owned handles are in use.

Python API reference status
---------------------------

The extension contains pybind11 docstrings, but a complete generated Python API
reference and type-stub set are still pending. The generated
:doc:`../cpp-api/index` is currently the detailed method reference; Python-only
differences should be documented on this page as they are audited.
