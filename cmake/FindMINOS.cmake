# Find MINOS
# Try to locate the MINOS library using the MINOSDIR environment variable
# or the MINOS_DIR cmake variable.
# MINOS is found by searching for the static library and header file.
#
# Create the following variables:
#
# MINOS_DIR         - Root directory of MINOS
# MINOS_LIBRARY     - Default library to link against to use MINOS
# MINOS_LIBRARY_DIR - Directory to the MINOS library
# MINOS_INCLUDE_DIR - Directory to include to use MINOS
# MINOS_FOUND       - If false, MINOS not found

# Test the MINOSDIR environment variable and create cmake MINOSDIR_DIR variable accordingly
set(MINOS_DIR_TEST $ENV{MINOSDIR})
if(MINOS_DIR_TEST)
  set(MINOS_DIR $ENV{MINOSDIR}      CACHE PATH "Path to MINOS install directory")
endif()

find_library(MINOS_LIBRARY 
             REQUIRED
             NAMES minos libminos
             PATHS ${MINOS_DIR} ${MINOS_DIR}/lib ${MINOS_DIR}/bin
             PATH_SUFFIXES minos)
get_filename_component(MINOS_LIBRARY_DIR "${MINOS_LIBRARY}" DIRECTORY)

# Find include path
find_path(MINOS_INCLUDE_DIR
            NAMES minos.h 
            PATHS ${MINOS_DIR} ${MINOS_DIR}/include
            PATH_SUFFIXES minos)

# Handle MINOS_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MINOS DEFAULT_MSG MINOS_LIBRARY MINOS_DIR MINOS_INCLUDE_DIR)