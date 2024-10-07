# Find Windows GCC compiler
# Try to locate the Tiny C compiler searching in system paths.
#
# Create the following variables
#
# GCC_DIR           - Root directory of GCC
# GCC_EXECUTABLE    - GCC executable

# Find executable
find_program(GCC_EXECUTABLE
    NAMES bin/gcc.exe)
# Find root directory
find_path(GCC_DIR
    NAMES bin/gcc.exe)
# Handle GCC_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GCC DEFAULT_MSG GCC_EXECUTABLE GCC_DIR)