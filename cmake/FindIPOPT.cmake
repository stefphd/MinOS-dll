# Find IPOPT
# Try to locate the IPOPT library using the IPOPT_DIR environment variable.
# If IPOPT_DIR not exist, just try to find IPOPT in the ${CMAKE_SOURCE_DIR} directory.
# If not found, search is extended to default system paths.
# IPOPT is found by searching for the library and header files.
#
# Create the following variables:
#
# IPOPT_DIR         - Root directory of IPOPT
# IPOPT_LIBRARY     - Default library to link against to use IPOPT
# IPOPT_LIBRARY_DIR - Directory to the IPOPT library
# IPOPT_INCLUDE_DIR - Directory to include to use IPOPT
# IPOPT_RUNTIME     - IPOPT runtime libraries
# IPOPT_RUNTIME_DIR - Directory to the IPOPT runtime libraries
# IPOPT_FOUND       - If false, IPOPT not found

# Test the IPOPT_DIR environment variable and create cmake IPOPT_DIR variable accordingly
set(IPOPT_DIR_TEST $ENV{IPOPT_DIR})
if(IPOPT_DIR_TEST)
  set(IPOPT_DIR $ENV{IPOPT_DIR}     CACHE PATH "Path to IPOPT install directory")
else()
  set(IPOPT_DIR ${CMAKE_SOURCE_DIR} CACHE PATH "Path to IPOPT install directory")
endif()

# Find library
# first search in ${IPOPT_DIR} with no additional (default) system paths
find_library(IPOPT_LIBRARY 
             NAMES ipopt ipopt.dll libipopt.dll
             PATHS ${IPOPT_DIR}/ ${IPOPT_DIR}/lib ${IPOPT_DIR}/bin
             PATH_SUFFIXES coin coin-or
             NO_DEFAULT_PATH)
# if not found, search in default paths
find_library(IPOPT_LIBRARY 
             NAMES ipopt ipopt.dll libipopt.dll
             PATH_SUFFIXES coin coin-or)
get_filename_component(IPOPT_LIBRARY_DIR "${IPOPT_LIBRARY}" DIRECTORY)

# Find include path
# first search in ${IPOPT_DIR} with no additional (default) system paths
find_path(IPOPT_INCLUDE_DIR 
          NAMES IpIpoptApplication.hpp
          PATHS ${IPOPT_DIR}/ ${IPOPT_DIR}/include
          PATH_SUFFIXES coin coin-or
          NO_DEFAULT_PATH)
# if not found, search in default paths
find_path(IPOPT_INCLUDE_DIR 
          NAMES IpIpoptApplication.hpp
          PATH_SUFFIXES coin coin-or)
          
# Find DLL path
if (WIN32)
  find_file(IPOPT_RUNTIME
            NAMES ipopt-3.dll ipopt.dll libipopt-3.dll libipopt.dll
            PATHS ${IPOPT_DIR}/bin
            PATHS ${IPOPT_DIR}/dll
            PATH_SUFFIXES coin coin-or
            NO_DEFAULT_PATH)
  get_filename_component(IPOPT_RUNTIME_DIR "${IPOPT_RUNTIME}" DIRECTORY)
  file(GLOB_RECURSE IPOPT_RUNTIME 
    "${IPOPT_RUNTIME_DIR}/*.dll"
    )
else()
  set(IPOPT_RUNTIME ${IPOPT_LIBRARY}) 
  get_filename_component(IPOPT_RUNTIME_DIR "${IPOPT_RUNTIME}" DIRECTORY)
  file(GLOB_RECURSE IPOPT_RUNTIME 
    "${IPOPT_RUNTIME_DIR}/lib*ipopt*.so*"
    "${IPOPT_RUNTIME_DIR}/lib*mumps*.so*"
    )
  set(IPOPT_RUNTIME_DIR ${IPOPT_RUNTIME_DIR}/)
endif()

# Handle IPOPT_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPOPT DEFAULT_MSG IPOPT_LIBRARY IPOPT_INCLUDE_DIR IPOPT_RUNTIME)
