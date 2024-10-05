# Find WORHP
# Try to locate the WORHP C library using the WORHPDIR environment variable
# WORHP is found by searching for the static library and header file.
#
# Create the following variables:
#
# WORHP_DIR         - Root directory of WORHP
# WORHP_LIBRARY     - Default library to link against to use WORHP
# WORHP_LIBRARY_DIR - Directory to the WORHP library
# WORHP_INCLUDE_DIR - Directory to include to use WORHP
# WORHP_FOUND       - If false, WORHP not found

# Test the WORHPDIR environment variable and create cmake WORHP_DIR variable accordingly
set(WORHP_DIR_TEST $ENV{WORHPDIR})
if(WORHP_DIR_TEST)
  set(WORHP_DIR $ENV{WORHPDIR}      CACHE PATH "Path to WORHP install directory")
else()
  set(WORHP_DIR ${CMAKE_SOURCE_DIR} CACHE PATH "Path to WORHP install directory")
endif()

# Find library
find_library(WORHP_LIBRARY 
             NAMES worhp
             PATHS ${WORHP_DIR}/ ${WORHP_DIR}/lib ${WORHP_DIR}/bin
             PATH_SUFFIXES worhp)
get_filename_component(WORHP_LIBRARY_DIR "${WORHP_LIBRARY}" DIRECTORY)

# Find include path
find_path(WORHP_INCLUDE_DIR
        NAMES worhp.h
        PATHS ${WORHP_DIR}/ ${WORHP_DIR}/include
        PATH_SUFFIXES worhp)

# Handle WORHP_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(WORHP DEFAULT_MSG WORHP_LIBRARY WORHP_DIR WORHP_INCLUDE_DIR)