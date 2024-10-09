# Find Windows GCC compiler
# Try to locate the GCC using the GCCDIR environemt variable and
# searching in the system paths.
#
# Create the following variables
#
# GCC_DIR           - Root directory of GCC
# GCC_EXECUTABLE    - GCC executable

# Test the GCCDIR environment variable and create cmake GCC_DIR variable accordingly
set(GCC_DIR $ENV{GCCDIR})

# Find executable
find_program(GCC_EXECUTABLE
    NAMES bin/gcc.exe
    PATHS ${GCC_DIR})

# Extract root directory
find_path(GCC_DIR
    NAMES bin/gcc.exe)

# Handle GCC_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GCC DEFAULT_MSG GCC_EXECUTABLE GCC_DIR)