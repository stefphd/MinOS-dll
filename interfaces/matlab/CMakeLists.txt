# Find the MATLAB package
find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM MEX_COMPILER)
# MinOS MEX
get_filename_component(Matlab_MX_LIBRARY_PATH ${Matlab_MX_LIBRARY} DIRECTORY )
find_library(MXLIBUT NAMES libut ut 
             PATHS ${Matlab_MX_LIBRARY_PATH} REQUIRED)
matlab_add_mex(NAME minos-mex
                SRC ${CMAKE_CURRENT_SOURCE_DIR}/minosMex.cpp
                LINK_TO ${MXLIBUT} minos)
set_target_properties(minos-mex PROPERTIES OUTPUT_NAME minosMex)
# P-code buildOCP
set(BUILDOCP_MFILE ${CMAKE_CURRENT_SOURCE_DIR}/buildOCP.m)
set(BUILDOCP_PFILE ${CMAKE_CURRENT_BINARY_DIR}/buildOCP.p)
set(MATLAB_PCODE_CMD "pcode('${BUILDOCP_MFILE}')")
add_custom_command(TARGET minos-mex POST_BUILD
    DEPENDS ${BUILDOCP_MFILE}
    COMMAND ${Matlab_MAIN_PROGRAM} -batch ${MATLAB_PCODE_CMD}                               
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}          
    COMMENT "Running MATLAB to convert ${BUILDOCP_MFILE} to ${BUILDOCP_PFILE}"
    VERBATIM)
# Configure files
configure_file(${CMAKE_SOURCE_DIR}/cmake/buildOCP.in
        ${CMAKE_CURRENT_BINARY_DIR}/buildOCP.m
        @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/cmake/minosMex.in
        ${CMAKE_CURRENT_BINARY_DIR}/minosMex.m
        @ONLY)
# Install files
set(MINOSMEX_DESTDIR lib)
if (WIN32)
        set(MINOSMEX_DESTDIR bin)
endif()
install(TARGETS minos-mex
        RUNTIME 
        DESTINATION ${MINOSMEX_DESTDIR}
        COMPONENT matlab)
install(FILES ${BUILDOCP_PFILE} ${CMAKE_CURRENT_BINARY_DIR}/buildOCP.m ${CMAKE_CURRENT_BINARY_DIR}/minosMex.m
        COMPONENT matlab
        DESTINATION ${MINOSMEX_DESTDIR})
install(DIRECTORY ./examples/
        COMPONENT matlab
        DESTINATION ./examples/matlab
        PATTERN *.cpp EXCLUDE
        PATTERN CMakeLists.txt EXCLUDE
        PATTERN examples/*.txt EXCLUDE
        PATTERN examples/*.cpp EXCLUDE
        PATTERN examples/*.c EXCLUDE
        PATTERN examples/*.asv EXCLUDE
        PATTERN examples/*.mex* EXCLUDE
        PATTERN examples/*.log EXCLUDE)
# CPack
if (WITH_CPACK)
        cpack_add_install_type(matlab 
                        DISPLAY_NAME "Matlab interface")
        cpack_add_component(matlab 
                        DISPLAY_NAME "Matlab interface"
                        DESCRIPTION "MinOS MEX inteface"
                        DEPENDS libraries
                        INSTALL_TYPES full matlab)
endif()