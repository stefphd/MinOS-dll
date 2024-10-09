# Find SNOPT
# Try to locate the SNOPT7 C library using the SNOPTDIR environment variable
# and searching in the system paths.
# SNOPT is found by searching for the static library and header file.
#
# Create the following variables:
#
# SNOPT_DIR         - Root directory of SNOPT
# SNOPT_LIBRARY     - Default library to link against to use SNOPT
# SNOPT_LIBRARY_DIR - Directory to the SNOPT library
# SNOPT_INCLUDE_DIR - Directory to include to use SNOPT
# SNOPT_FOUND       - If false, SNOPT not found

# Test the SNOPTDIR environment variable and create cmake SNOPT_DIR variable accordingly
set(SNOPT_DIR ${CMAKE_SOURCE_DIR} CACHE PATH "Path to SNOPT install directory")
set(SNOPT_DIR $ENV{SNOPTDIR})

# Set SNOPT version
set(SNOPT_VERSION 7)

# Find library
find_library(SNOPT_LIBRARY 
             NAMES "snopt${SNOPT_VERSION}"
             PATHS ${SNOPT_DIR}/ ${SNOPT_DIR}/lib ${SNOPT_DIR}/bin
             PATH_SUFFIXES snopt7)
get_filename_component(SNOPT_LIBRARY_DIR "${SNOPT_LIBRARY}" DIRECTORY)

# Find include path
find_path(SNOPT_INCLUDE_DIR
        NAMES snopt.h 
        PATHS ${SNOPT_DIR}/ ${SNOPT_DIR}/include
        PATH_SUFFIXES snopt "snopt${SNOPT_VERSION}")

# Handle SNOPT_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SNOPT DEFAULT_MSG SNOPT_LIBRARY SNOPT_DIR SNOPT_INCLUDE_DIR)