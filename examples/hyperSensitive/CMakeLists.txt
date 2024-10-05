# Example name
set(EXAMPLE_NAME hyperSensitive)

# Generate C code with Python
execute_process(
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/generate.py" "${EXAMPLE_NAME}.c" "${CMAKE_CURRENT_BINARY_DIR}/" ${SKIP_HESSIAN_FLAG}
    RESULT_VARIABLE PYTHON_RESULT)

# Create targets if C code generated
if (NOT PYTHON_RESULT EQUAL 0)
    message(WARNING "Unable to generate ${EXAMPLE_NAME}.c. Skipping example ${EXAMPLE_NAME}.")
else()
    # Build problem dynamic library
    add_library(${EXAMPLE_NAME}-lib SHARED
                "${CMAKE_CURRENT_BINARY_DIR}/${EXAMPLE_NAME}.c"
                )
    set_target_properties(${EXAMPLE_NAME}-lib PROPERTIES
                            OUTPUT_NAME ${EXAMPLE_NAME})
    # Create the executable
    add_executable(${EXAMPLE_NAME}-exe
                  main.cpp ) 
    set_target_properties(${EXAMPLE_NAME}-exe PROPERTIES
                        OUTPUT_NAME ${EXAMPLE_NAME}-main)
    target_link_libraries(${EXAMPLE_NAME}-exe
                          PRIVATE minoslib)
    # Add test
    add_test(NAME ${EXAMPLE_NAME}-exe 
            COMMAND ${EXAMPLE_NAME}-main ${NLPSOLVER_EXAMPLE}
            WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
    # Add IPOPT path to PATH
    if (WIN32)
    set_property(TEST ${EXAMPLE_NAME}-exe PROPERTY ENVIRONMENT_MODIFICATION  
        PATH=path_list_append:${BINARY_DIR})
    else()
    set_property(TEST ${EXAMPLE_NAME}-exe PROPERTY ENVIRONMENT_MODIFICATION
        LD_LIBRARY_PATH=path_list_append:${BINARY_DIR})
    endif()
endif()