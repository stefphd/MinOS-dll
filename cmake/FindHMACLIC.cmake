# Find HMACLIC
# Try to locate the HMACLIC library using the HMACLICDIR environment variable
# and searching in the system paths.
# HMACLIC is found by searching for the static library and header file.
#
# Create the following variables:
#
# HMACLIC_DIR         - Root directory of HMACLIC
# HMACLIC_LIBRARY     - Default library to link against to use HMACLIC
# HMACLIC_LIBRARY_DIR - Directory to the HMACLIC library
# HMACLIC_INCLUDE_DIR - Directory to include to use HMACLIC
# HMACLIC_RUNTIME     - HMACLIC runtime libraries
# HMACLIC_RUNTIME_DIR - Directory to the HMACLIC runtime libraries
# HMACLIC_FOUND       - If false, HMACLIC not found

# Test the HMACLICDIR environment variable and create cmake HMACLIC_DIR variable accordingly
set(HMACLIC_DIR ${CMAKE_SOURCE_DIR} CACHE PATH "Path to HMACLIC install directory")
set(HMACLIC_DIR $ENV{HMACLICDIR})

# Find library
find_library(HMACLIC_LIBRARY 
             NAMES hmaclic hmaclic.dll libhmaclic.dll
             PATHS ${HMACLIC_DIR}/ ${HMACLIC_DIR}/lib ${HMACLIC_DIR}/bin)
get_filename_component(HMACLIC_LIBRARY_DIR "${HMACLIC_LIBRARY}" DIRECTORY)

# Find include path
find_path(HMACLIC_INCLUDE_DIR 
          NAMES hmaclic.h
          PATHS ${HMACLIC_DIR}/ ${HMACLIC_DIR}/include)
          
# Find DLL path
if (WIN32)
  find_file(HMACLIC_RUNTIME
            NAMES hmaclic.dll libhmaclic.dll
            PATHS ${HMACLIC_DIR}/bin)
  get_filename_component(HMACLIC_RUNTIME_DIR "${HMACLIC_RUNTIME}" DIRECTORY)
else()
  set(HMACLIC_RUNTIME ${HMACLIC_LIBRARY}) 
  get_filename_component(HMACLIC_RUNTIME_DIR "${HMACLIC_RUNTIME}" DIRECTORY)
endif()

# Handle HMACLIC_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HMACLIC DEFAULT_MSG HMACLIC_LIBRARY HMACLIC_INCLUDE_DIR HMACLIC_RUNTIME)
