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

Build prerequisites
~~~~~~~~~~~~~~~~~~~

The wheel contains a native CUDA/C++ extension and is compiled on the machine
that creates it. Before building, verify that the intended Python interpreter,
CMake, CUDA compiler, and NVIDIA driver are available:

.. code-block:: console

   python3 --version
   cmake --version
   nvcc --version
   nvidia-smi

Initialize the bundled Git dependencies:

.. code-block:: console

   git submodule update --init --recursive

Run the following commands from the repository root. Start with an empty
``dist/`` directory, or move artifacts from earlier builds elsewhere, so that
validation and installation cannot accidentally select an older wheel.

Create and validate the wheel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use a dedicated virtual environment for packaging:

.. code-block:: console

   python3 -m venv .venv-wheel
   source .venv-wheel/bin/activate

   python -m pip install --upgrade pip
   python -m pip install build twine

   python -m build --wheel
   python -m twine check dist/*

``python -m build --wheel`` invokes the ``scikit-build-core`` backend from
``pyproject.toml``. That backend configures CMake with
``DEME_BUILD_PYTHON=ON``, compiles the native ``deme._deme`` extension in
Release mode, and places the resulting wheel under ``dist/``.

``twine check`` validates the wheel metadata and the rendering of its package
description. A successful build should produce a platform-specific file whose
name resembles:

.. code-block:: text

   dist/deme-3.0.0-<python-tag>-<abi-tag>-linux_<architecture>.whl

This is not a pure-Python or universal wheel. Its filename tags determine which
Python interpreter and operating-system ABI pip will accept, while CUDA and
driver compatibility must also be validated separately.

Test the wheel in a clean environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Leave the packaging environment, create a separate test environment, and
install the wheel there:

.. code-block:: console

   deactivate

   python3 -m venv /tmp/deme-wheel-test
   source /tmp/deme-wheel-test/bin/activate

   python -m pip install --upgrade pip
   python -m pip install dist/deme-3.0.0-*.whl
   python -m pip check

Run import checks from outside the source tree. Otherwise, files in the
checkout could hide missing wheel contents:

.. code-block:: console

   cd /tmp
   python -c "import deme; print(deme.__version__, deme.__file__)"
   python -c "import DEME; print(DEME.__version__)"

The first command should report version ``3.0.0`` and a module path inside
``/tmp/deme-wheel-test``. The second verifies the compatibility import; new
applications should continue to use lowercase ``import deme``.

On a host with a supported visible NVIDIA GPU, also construct a solver to
exercise CUDA initialization and confirm the selected logical device:

.. code-block:: console

   python -c "import deme; s = deme.DEMSolver([0]); print(s.GetGPUDeviceIDs())"

Expected output is ``[0, 0]``. Import checks alone do not exercise solver
construction, CUDA device selection, or worker allocation.

Return to the checkout when testing is complete:

.. code-block:: console

   deactivate
   cd /path/to/DEM-Engine

Optional publication checks
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Publishing changes external package state and is only for authorized
maintainers. TestPyPI can be used before a production upload:

.. code-block:: console

   source .venv-wheel/bin/activate
   python -m twine upload --repository testpypi dist/deme-3.0.0-*.whl

After validating the TestPyPI artifact, an authorized production upload uses:

.. code-block:: console

   python -m twine upload dist/deme-3.0.0-*.whl

Use API-token authentication or a configured trusted-publishing workflow.
Uploading is not part of the normal local build or test procedure.

Wheel portability
~~~~~~~~~~~~~~~~~

Before distributing a wheel, record and test at least:

* the Python and ABI tag in the wheel filename;
* the Linux distribution and minimum compatible ``glibc`` baseline;
* the CUDA Toolkit used for compilation;
* the minimum NVIDIA driver version;
* the GPU architectures included by the CUDA build; and
* imported shared-library dependencies.

Until the supported compatibility matrix is published, build and validate
wheels on the oldest intended deployment platform and test them on each
supported Python, CUDA/driver, and GPU configuration.

C++ build
---------

.. code-block:: console

   git submodule update --init --recursive
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build --parallel

Use a focused demo or modular-test target first when validating a change.
