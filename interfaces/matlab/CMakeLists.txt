# Find the MATLAB package
find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM MEX_COMPILER)
# MinOS MEX
get_filename_component( Matlab_MX_LIBRARY_PATH ${Matlab_MX_LIBRARY} DIRECTORY )
find_library(MXLIBUT NAMES libut ut 
             PATHS ${Matlab_MX_LIBRARY_PATH} REQUIRED)
matlab_add_mex(NAME minosMex
                SRC ${CMAKE_SOURCE_DIR}/interfaces/matlab/minosMex.cpp
                LINK_TO ${MXLIBUT} minoslib)
# P-code buildOCP
set(BUILDOCP_MFILE ${CMAKE_CURRENT_SOURCE_DIR}/buildOCP.m)
set(BUILDOCP_PFILE ${CMAKE_CURRENT_BINARY_DIR}/buildOCP.p)
set(MATLAB_PCODE_CMD ${Matlab_MAIN_PROGRAM} -batch "pcode('${BUILDOCP_MFILE}')")
add_custom_command(
    OUTPUT  ${BUILDOCP_PFILE}
    DEPENDS ${BUILDOCP_MFILE}
    COMMAND ${MATLAB_PCODE_CMD}                               
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}          
    COMMENT "Running MATLAB to convert ${BUILDOCP_MFILE} to ${BUILDOCP_PFILE}"
    VERBATIM                                                 
)
add_custom_target(minosMexPcode ALL
    DEPENDS ${BUILDOCP_PFILE}
)

# Install files
set(MINOSMEX_DESTDIR lib)
if (WIN32)
        set(MINOSMEX_DESTDIR bin)
endif()
install(TARGETS minosMex
        RUNTIME 
        DESTINATION ${MINOSMEX_DESTDIR}
        COMPONENT matlab)
install(FILES ${BUILDOCP_PFILE}
        COMPONENT matlab
        DESTINATION ${MINOSMEX_DESTDIR})
install(FILES buildOCP_help.m
        COMPONENT matlab
        RENAME buildOCP.m
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
cpack_add_install_type(matlab 
                DISPLAY_NAME "Matlab interface")
cpack_add_component(matlab 
                DISPLAY_NAME "Matlab interface"
                DESCRIPTION "MinOS MEX inteface"
                DEPENDS libraries
                INSTALL_TYPES full matlab)