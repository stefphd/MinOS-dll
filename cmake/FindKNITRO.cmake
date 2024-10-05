# Find KNITRO
# Try to locate the KNITRO C library using the KNITRODIR environment variable
# KNITRO is found by searching for the static library and header file.
#
# Create the following variables:
#
# KNITRO_DIR         - Root directory of KNITRO
# KNITRO_VERSION     - Version of KNITRO
# KNITRO_LIBRARY     - Default library to link against to use KNITRO
# KNITRO_LIBRARY_DIR - Directory to the KNITRO library
# KNITRO_INCLUDE_DIR - Directory to include to use KNITRO
# KNITRO_FOUND       - If false, KNITRO not found

# Test the KNITRODIR environment variable and create cmake KNITRO_DIR variable accordingly
set(KNITRO_DIR $ENV{KNITRODIR} CACHE PATH "Path to KNITRO install directory")

# Extract the version number from KNITRODIR
set(KNITRO_VERSION "")
if (WIN32)
string(REGEX MATCH "([0-9]+)\\.([0-9]+)\\.([0-9]+)" _ KNITRO_VERSION_MATCH ${KNITRO_DIR})
set(KNITRO_VERSION ${CMAKE_MATCH_1}${CMAKE_MATCH_2}${CMAKE_MATCH_3})
endif()

# Find library
find_library(KNITRO_LIBRARY 
             NAMES knitro "knitro${KNITRO_VERSION}"
             PATHS ${KNITRO_DIR}/ ${KNITRO_DIR}/lib
             PATH_SUFFIXES knitro)
get_filename_component(KNITRO_LIBRARY_DIR "${KNITRO_LIBRARY}" DIRECTORY)

# Find include path
find_path(KNITRO_INCLUDE_DIR 
        NAMES knitro.h
        PATHS ${KNITRO_DIR}/ ${KNITRO_DIR}/include
        PATH_SUFFIXES knitro)

# Handle KNITRO_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(KNITRO DEFAULT_MSG KNITRO_LIBRARY KNITRO_DIR KNITRO_INCLUDE_DIR)