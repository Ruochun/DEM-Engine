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

The release-wheel policy currently covers 64-bit x86 Linux with glibc 2.28 or
newer, CPython 3.9 through 3.14, and CUDA 12.9. The installed machine must
provide an NVIDIA driver compatible with CUDA 12.9. Source builds can continue
to use other supported CUDA Toolkit versions, but those builds are outside the
binary-wheel compatibility policy.

Python package
--------------

Install a released wheel with:

.. code-block:: console

   python -m pip install deme

The canonical import is:

.. code-block:: python

   import deme

The historical ``import DEME`` spelling remains available as a compatibility
alias. New code should use the canonical lowercase ``deme`` distribution and
import name.

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

Run the environment-creation commands from the repository root. Start with an
empty ``dist/`` directory, or move artifacts from earlier builds elsewhere, so
that validation and installation cannot accidentally select an older wheel.

Create and validate the wheel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use a dedicated virtual environment for packaging:

.. code-block:: console

   python3 -m venv .venv-wheel
   source .venv-wheel/bin/activate

   python -m pip install --upgrade pip
   python -m pip install build twine

   cd ..
   python -m build --wheel --outdir DEM-Engine/dist DEM-Engine
   cd DEM-Engine
   python -m twine check dist/*

The build command deliberately runs from the checkout's parent directory. A
pre-existing local ``build/`` directory in the repository root can otherwise
shadow the PyPA package named ``build`` and cause
``No module named build.__main__``. Replace ``DEM-Engine`` with the checkout
directory name if it differs.

``python -m build --wheel`` invokes the ``scikit-build-core`` backend from
``pyproject.toml``. That backend configures CMake with
``DEME_BUILD_PYTHON=ON``, compiles the native ``deme._deme`` extension in
Release mode, and places the resulting wheel under ``dist/``.

``twine check`` validates the wheel metadata and the rendering of its package
description. A successful build should produce a platform-specific file whose
name resembles:

.. code-block:: text

   dist/deme-3.0.1-<python-tag>-<abi-tag>-linux_<architecture>.whl

This is not a pure-Python or universal wheel. Its filename tags determine which
Python interpreter and operating-system ABI pip will accept, while CUDA and
driver compatibility must also be validated separately.

Build the complete supported Python matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One native wheel must be built by each targeted CPython interpreter. Building
with Python 3.13, for example, creates only the ``cp313`` wheel; it does not
also create wheels for the other supported Python versions. DEM-Engine 3
currently targets CPython 3.9 through 3.14.

The following reproducible Conda workflow creates a separate build environment
for every target. Run the environment-creation commands from the repository
root:

.. code-block:: console

   conda create --yes --name deme-wheel-py39 python=3.9 pip
   conda create --yes --name deme-wheel-py310 python=3.10 pip
   conda create --yes --name deme-wheel-py311 python=3.11 pip
   conda create --yes --name deme-wheel-py312 python=3.12 pip
   conda create --yes --name deme-wheel-py313 python=3.13 pip
   conda create --yes --name deme-wheel-py314 python=3.14 pip

   for env in deme-wheel-py39 deme-wheel-py310 deme-wheel-py311 deme-wheel-py312 deme-wheel-py313 deme-wheel-py314; do
       conda run --name "$env" python -m pip install --upgrade pip build twine auditwheel
   done

Build once with each interpreter. As in the single-version workflow, run the
builder from the checkout's parent directory to prevent the local ``build/``
directory from shadowing the PyPA ``build`` package:

.. code-block:: console

   cd ..
   for env in deme-wheel-py39 deme-wheel-py310 deme-wheel-py311 deme-wheel-py312 deme-wheel-py313 deme-wheel-py314; do
       conda run --name "$env" python -m build --wheel --outdir DEM-Engine/dist DEM-Engine
   done
   cd DEM-Engine

   conda run --name deme-wheel-py313 python -m twine check dist/*.whl

The resulting directory should contain six distinct wheels with ``cp39``,
``cp310``, ``cp311``, ``cp312``, ``cp313``, and ``cp314`` tags.
Confirm that explicitly:

.. code-block:: console

   ls -1 dist/deme-3.0.1-cp3*-linux_*.whl
   for wheel in dist/*.whl; do
       conda run --name deme-wheel-py313 python -m auditwheel show "$wheel"
   done

Generating all six files is only the build step. Each wheel must still be
installed and exercised with its matching Python version before that version
is considered validated. CUDA, Linux ABI, and GPU compatibility also require
separate testing; ``auditwheel show`` reports the native shared-library and
``glibc`` requirements but does not prove runtime compatibility.

Build release wheels with cibuildwheel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Conda commands above are useful for native development builds. Release
wheels use ``cibuildwheel`` and PyPA's CUDA-enabled
``manylinux_2_28_x86_64_cuda12_9`` container so the result does not inherit the
Linux ABI of the maintainer's workstation. The configuration is stored in
``pyproject.toml``.

With Docker available, build the same complete matrix locally from the parent
of the checkout:

.. code-block:: console

   python3 -m venv .venv-cibuildwheel
   source .venv-cibuildwheel/bin/activate
   python -m pip install --upgrade pip
   python -m pip install "cibuildwheel==4.1.1" twine auditwheel

   cd ..
   python -m cibuildwheel --platform linux --output-dir DEM-Engine/wheelhouse DEM-Engine
   cd DEM-Engine

   python -m twine check wheelhouse/*.whl
   for wheel in wheelhouse/*.whl; do
       python -m auditwheel show "$wheel"
   done

This release process uses ``auditwheel repair`` to copy ordinary redistributable
native dependencies into each wheel and assign the
``manylinux_2_28_x86_64`` tag. It explicitly excludes ``libcuda.so.1``,
``libcudart.so.12``, and ``libnvrtc.so.12``. DEME runtime-compiles CUDA kernels,
so a compatible CUDA 12.9 toolkit (including NVRTC, its builtins, and headers)
and NVIDIA driver must be installed on the deployment host. Bundling a driver
stub is incorrect, while bundling NVRTC without all of its dynamically loaded
resources produces an incomplete runtime. Before publishing, inspect the
repaired wheel and ``auditwheel show`` output to confirm that CUDA is the only
non-system external dependency.

Automated wheel builds
~~~~~~~~~~~~~~~~~~~~~~

``.github/workflows/python-wheels.yml`` runs the same policy as six parallel
jobs, one for each CPython ABI. It runs on relevant pull requests, release tags,
or manual dispatch. Every job:

* checks out Git submodules recursively;
* builds in the CUDA 12.9 manylinux 2.28 container;
* repairs the wheel with ``auditwheel``;
* checks package metadata and the expected Python/platform filename tags; and
* uploads the wheel as a workflow artifact for testing or release assembly.

The hosted build runners do not provide a usable NVIDIA GPU. Consequently this
workflow validates compilation, repair, metadata, and tags but deliberately
does not claim GPU runtime validation. Install each artifact on a compatible
GPU host and run the tests below before publishing it.

Test the wheel in a clean environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Leave the packaging environment, create a separate test environment, and
install the wheel there:

.. code-block:: console

   deactivate

   python3 -m venv /tmp/deme-wheel-test
   source /tmp/deme-wheel-test/bin/activate

   python -m pip install --upgrade pip
   python -m pip install dist/deme-3.0.1-*.whl
   python -m pip check

Run import checks from outside the source tree. Otherwise, files in the
checkout could hide missing wheel contents:

.. code-block:: console

   cd /tmp
   python -c "import deme; print(deme.__version__, deme.__file__)"
   python -c "import DEME; print(DEME.__version__)"

The first command should report version ``3.0.1`` and a module path inside
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

Publish the deme distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Publishing changes external package state and is only for authorized
maintainers. ``python-wheels.yml`` uses PyPI Trusted Publishing, so it does not
store a long-lived PyPI token in GitHub.

Before the first upload, create a GitHub environment named ``pypi`` under
``Settings`` then ``Environments``. Configure required reviewers so that the
publication job always pauses for approval. Then sign in to PyPI, open the
account-level ``Publishing`` page, and add a pending GitHub publisher with:

* PyPI project name: ``deme``;
* GitHub owner: ``Ruochun``;
* repository: ``DEM-Engine``;
* workflow filename: ``python-wheels.yml``; and
* environment: ``pypi``.

The pending publisher creates ``deme`` on the first successful upload. It
does not reserve the name before that upload. The project name must exactly
match ``name = "deme"`` in ``pyproject.toml``.

To build without publishing, open the repository's ``Actions`` tab, select
``Build Python wheels``, choose ``Run workflow``, leave
``publish_to_pypi`` disabled, and run it from the intended commit or branch.
Download and test all six artifacts after the jobs succeed.

To publish the already-reviewed source commit, dispatch the same workflow
again with ``publish_to_pypi`` enabled. The six build jobs run again; only if
all succeed does the ``Publish deme wheels to PyPI`` job enter the protected
``pypi`` environment. Approve that deployment after checking the commit and
wheel jobs. The publishing job downloads the six artifacts and uploads them
with a short-lived PyPI OIDC credential.

PyPI does not allow replacing a file or reusing an existing release version.
If any ``deme`` version ``3.0.1`` file has already been uploaded, increment the
project version and rebuild the complete wheel set rather than retrying with
different bytes under the same version.

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
