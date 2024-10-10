# Set binary directory to add to PATH
set(BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
# Check Hessian
if (SKIP_HESSIAN)
    set(MINOS_BASELIB ${MINOS_BASELIB}-nh)
    # Add --skip_hessian flag if SKIP_HESSIAN
    set(SKIP_HESSIAN_FLAG --skip_hessian)
endif()
# Find python interpreter
find_package(Python3 
            COMPONENTS Interpreter REQUIRED)
# Check CasADi
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import casadi"
    RESULT_VARIABLE CASADI_RESULT)

if (NOT CASADI_RESULT EQUAL 0)
    message(FATAL_ERROR "CasADi not found for C++ examples.")
endif()
message(STATUS "CasADi found.")   
# Logical target to build all example with --target examples
add_custom_target(examples)
add_dependencies(examples minos-${NLPSOLVER_EXAMPLE})
# Add example subdirectories
add_subdirectory(brachistochrone)
add_subdirectory(brachistochrone_auxdata)
add_subdirectory(brachistochrone_c)
add_subdirectory(brysonDenham)
add_subdirectory(hyperSensitive)
add_subdirectory(goddardRocket)
add_subdirectory(directCollocation)
# Configure CMakeLists.txt for release
configure_file(${CMAKE_SOURCE_DIR}/cmake/CMakeLists.in
			   ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
			   @ONLY)
# Install rules
install(DIRECTORY ./
        COMPONENT cpp_examples
        DESTINATION ./examples/cpp
        PATTERN ./CMakeLists.txt EXCLUDE)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
        COMPONENT cpp_examples
        DESTINATION ./examples/cpp)
install(FILES ${CMAKE_SOURCE_DIR}/cmake/FindMINOS.cmake
        COMPONENT cpp_examples
        DESTINATION ./examples/cpp/cmake)  
# CPack
if (WITH_CPACK)
    cpack_add_component(cpp_examples
                        DISPLAY_NAME "C++ examples"
                        DESCRIPTION "MinOS examples for C++"
                        DEPENDS libraries
                        INSTALL_TYPES full static)
endif()