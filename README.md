# MinOS

MinOS (**Min**imal **O**ptimal-control **S**olver) is a C++ solver for optimal control problem (OCP) using a direct method together with CasADi to generate the C code of the OCP functions and a NLP solver (such as IPOPT) to solve the subsequent large-scale nonlinear programming problem (NLP). Interfaces are available for MATLAB, Python, and C or C++.

This README is developer-oriented. User-oriented README is provided with the release in the documentation.

Build and tested with:

* *Windows 10* with MSVC >=2019 and MATLAB >=R2022b
* *Linux* (Debian 12) with g++ and MATLAB >= R2022b

## Building requirements

* [CMake](https://cmake.org/) build system to build the software.
* A C++ compiler. Preferred option is [MSVC](https://learn.microsoft.com/en-us/cpp) for *Windows* and [GCC](https://gcc.gnu.org/) for *Linux*. Alternativelly, for *Windows* one can employ [MinGW64-w64](https://www.mingw-w64.org/) (not tested). Intel compilers may be also employed.
* A Windows distribution of [GCC](https://gcc.gnu.org/). Make sure that the GCC root directory (i.e. where `bin/gcc.exe` is) is in the `PATH` environment variable. This is needed in *Windows* to compile the problem library in the MATLAB and Python interfaces.
* A NLP solver to solve the large-scale NLP. Possibilities are:
  * [IPOPT](https://github.com/coin-or/Ipopt) version 3.14 with linear solver [MUMPS](https://github.com/coin-or-tools/ThirdParty-Mumps), which are freely available. For *Windows*, dynamic and static libraries should be placed in `./bin` and `./lib`. Alternativelly, one could set a `IPOPT_DIR` environment variable pointing to the IPOPT installation path, which should contain the folders `bin`, `lib`, and `include`. IPOPT headers provided with this repository in `./include/coin-or` may be also employed (current IPOPT version is 3.14). For *Linux*, one should install IPOPT in `/usr/lib` and `/usr/include` (or `/usr/local/lib` and `/usr/local/include`). Again, alternativelly one could use a local installation of IPOPT with and set `$IPOPT_DIR` pointing to the IPOPT installation path that contains `lib` and `include` folders.
  * [KNITRO](https://www.artelys.com/solvers/knitro/), which is a commercial software. Free trial license is freely available for academics. A `KNITRODIR` environment variable pointing to the KNITRO installation path should be automatically set when installing KNITRO. If not, one needs to set `KNITRODIR` pointing to the KNITRO installation path. KNITRO is compatible only with MSVC for *Windows*.
  * [WORHP](https://worhp.de/), whose license is freely available for academics. For *Windows*, one should set a `WORHPDIR` environment variable pointing to the KNITRO installation path.
  * [SNOPT](https://ccom.ucsd.edu/~optimizers/solvers/snopt/), which is a commercial software. Free trial license is freely available for academics. A `SNOPTDIR` environment variable pointing to the SNOPT installation path needs to be set.
* [MATLAB](https://www.mathworks.com/products/matlab.html) to run the MATLAB examples.
* [Python](https://www.python.org/) to run the Python examples and generate the C code for the C++ examples; use `interfaces/python/requirements.txt` to install the Python requirements.
* [CasADi](https://web.casadi.org/) to generate the C code for the examples (with MATLAB and Python interfaces).
* [Doxygen](https://doxygen.nl/) and a latex distribution to generate the documentation.
* A CMake-CPack generator compatible with your system to generate the release.

Make sure that the necessary runtime libraries are found by the system. Basically, you need to add the binary directory of the NLP solver into `PATH` for *Windows* and `LD_LIBRARY_PATH` for *Linux*.

## Installing IPOPT with MUMPS

* *Windows*: pre-compiled binaries compiled with MSVC are available in [IPOPT-releases](https://github.com/coin-or/Ipopt/releases); use *-md* version (non-debug). For MinGW-w64, you need to build IPOPT along with its dependencies youself, as the IPOPT binaries must be compatible with the selected compiler; see also [Installing IPOPT](https://coin-or.github.io/Ipopt/INSTALL.html).
* *Linux*: you should build IPOPT along with its dependencies youself. For Debian-based distro, follow the steps in [Installing IPOPT](https://coin-or.github.io/Ipopt/INSTALL.html). For other distro, just replace the commands with you package manager and the corresponding packages. In some distros pre-built IPOPT packages may be employed, as they meet the supported version (e.g. [IPOPT-aur](https://aur.archlinux.org/packages/coin-or-ipopt) for ArchLinux). Note that MATLAB supports only 64-bit integer version of IPOPT (together with its dependencies), since LAPACK library used by MATLAB employs 64-bit integers; see [Building for 64-bit integers](https://coin-or.github.io/Ipopt/INSTALL.html#INT64_BUILD) for details on building IPOPT with 64-bit integers. Basically, you need to obtain liblapack64-dev and build METIS with 64-bit integers; MUMPS and IPOPT can be built with 64-bit integers using the configure flag `--with-intsize=64`. Note that using a 64-bit integer version of IPOPT has a number of limitations, such such limited linear solvers that can be used.

## Building

Building is performed using CMake. Building should be performed in release mode in *Windows*.

```command
mkdir build
cd build
cmake ..
cmake --build . [--config Release]
```

To employ MinGW-w64 compilers in *Windows*, you need to add `path/to/mingw-w64/bin` to the `PATH` and configure the project using

```command
set PATH=path/to/mingw-w64/bin;%PATH%
cmake -G "MinGW Makefiles"
```

The `CMakeLists.txt` file contains a number of options to change the building settings

* `WITH_IPOPT`: compile IPOPT interface
* `WITH_WORHP`: compile WORHP interface
* `WITH_KNITRO`: compile KNITRO interface
* `WITH_SNOPT`: compile SNOPT interface
* `WITH_MEX_INTERFACE`: compile MEX interface for MATLAB
* `WITH_PYTHON_INTERFACE`: compile Python interface
* `WITH_GCC`: install GCC distribution (for Windows only)
* `WITH_DOCS`: generate docs using Doxygen (require Doxygen and latex installed)
* `WITH_CPACK`: generate CPack target
* `WITH_CPP_EXAMPLES`: compile C++ examples
* `SKIP_HESSIAN`: use approximate Hessian calculation in C++ example
* `NLPSOLVER_EXAMPLE`: NLP solver for C++ examples (ipopt, knitro, worhp, snopt)

## Test

To test the build (run all C++ examples if `WITH_CPP_EXAMPLES` is `ON`)

```command
ctest [-C Release]
```

To run MATLAB and Python examples, you need to build and install the project. For example, install the project in a local directory

```command
cmake --install . --prefix ./../install
```

and run the examples in `<installdir>/examples/matlab` and  `<installdir>/example/python`. See also  `<installdir>/requirements.txt` for the Python requirements.

## Package

To package the software into an installer

```command
cpack
```

You can add the option:

* `-g NSIS` to employ NSIS for *Windows* (default in *Windows*)
* `-g TGZ` to generate a TAR.GZ file (default in *Linux*)
* `-g ZIP` to generate a ZIP file (both *Windows* and *Linux*)

## Versioning

The version is identified by a triplet `X.Y.Z` where:

* `X` is the major version. Change in this number means changes in the OCP formulation or severe changes in the API that for sure break the backward-compatibility.
* `Y` is the minor version. Change in this number means new functionalities added, but that should not break the backward-compatibility.
* `Z` is the release number. Change in this number means minor fixing. which are done in a backward-compatible way.
