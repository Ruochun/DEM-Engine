Troubleshooting
===============

CUDA or PTX mismatch
--------------------

Errors such as ``CUDA_ERROR_UNSUPPORTED_PTX_VERSION`` usually indicate that the
installed NVIDIA driver cannot execute code produced by the selected CUDA
Toolkit. Verify the driver/toolkit compatibility and rebuild with a supported
toolkit.

Runtime JIT compilation cannot find headers
-------------------------------------------

DEME compiles kernels at runtime. Confirm that the CUDA Toolkit headers and the
``share/DEME/kernel`` and ``include`` resources installed with DEME are
available. The CUDA header major and minor version must match the loaded NVRTC
library. If several CUDA Toolkits are installed, set ``CUDA_HOME`` to the one
providing that NVRTC version, for example:

.. code-block:: console

   export CUDA_HOME=/usr/local/cuda-12.8

Errors in CUDA or CURAND headers involving undefined internal identifiers often
mean that an unversioned ``/usr/local/cuda`` link selected headers from a newer
Toolkit than the loaded ``libnvrtc``. Current DEME releases reject that mismatch
and report both the required NVRTC version and how to select matching headers.

Import works but initialization fails
-------------------------------------

``import deme`` verifies only that the extension and its immediate shared
libraries can load. A meaningful installation test must construct a solver and
run a small simulation through ``Initialize()`` so NVRTC, kernel resources, and
the GPU driver are exercised.

Conda ``GLIBCXX`` errors
-------------------------

Build with compilers compatible with the target Conda environment. A wheel
built against a newer system ``libstdc++`` may import on the build host but fail
inside another environment.

Stale runtime kernels
---------------------

After changing kernel or force-model sources, use a clean build or ensure the
runtime kernel assets have been refreshed. Old copied text sources can make the
runtime behavior disagree with the compiled host code.
