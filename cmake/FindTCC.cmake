# Find Tiny C Compiler (TCC)
# Try to locate the Tiny C compiler searching in system paths.
#
# Create the following variables
#
# TCC_DIR           - Root directory of TCC
# TCC_EXECUTABLE    - TCC executable
# TCC_INCLUDE_DIR   - TCC include directory for TCC library
# TCC_FOUND         - If false, TCC not found

# Find executable
find_program(TCC_EXECUTABLE
    NAMES tcc tcc.exe
    PATH_SUFFIXES tcc)
# Extract root directory
get_filename_component(TCC_DIR "${TCC_EXECUTABLE}" DIRECTORY)
# Find include path
find_path(TCC_INCLUDE_DIR
    NAMES libtcc.h
    PATHS ${TCC_DIR}
    PATH_SUFFIXES libtcc)
# Handle TCC_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TCC DEFAULT_MSG TCC_EXECUTABLE TCC_DIR)