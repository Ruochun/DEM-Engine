Installation
============

System requirements
-------------------

DEME requires:

* a 64-bit Linux system for the currently supported Python package;
* an NVIDIA GPU;
* an NVIDIA driver compatible with the selected CUDA Toolkit;
* CUDA Toolkit 11 or newer, including NVRTC and CUDA headers;
* CMake 3.18 or newer and a CUDA-compatible C++ compiler when building from
  source.

The exact Python, CUDA, compiler, driver, and GPU architecture matrix is being
validated for DEM-Engine 3. A wheel should not be assumed portable across
CUDA major versions until that matrix is published.

Python package
--------------

Install a released wheel with:

.. code-block:: console

   python -m pip install deme

The canonical import is:

.. code-block:: python

   import deme

The historical ``import DEME`` spelling remains available as a compatibility
alias.

Build a wheel from a checkout
-----------------------------

Initialize the bundled dependencies and use an isolated build environment:

.. code-block:: console

   git submodule update --init --recursive
   python3 -m venv .venv-wheel
   . .venv-wheel/bin/activate
   python -m pip install --upgrade pip build
   python -m build --wheel

Install the wheel in a separate environment and test it from outside the source
tree:

.. code-block:: console

   python3 -m venv /tmp/deme-wheel-test
   . /tmp/deme-wheel-test/bin/activate
   python -m pip install dist/deme-3.0.0-*.whl
   cd /tmp
   python -c "import deme; print(deme.__version__, deme.__file__)"

C++ build
---------

.. code-block:: console

   git submodule update --init --recursive
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build --parallel

Use a focused demo or modular-test target first when validating a change.
