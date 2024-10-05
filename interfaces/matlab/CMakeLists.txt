# MinOS MEX
get_filename_component( Matlab_MX_LIBRARY_PATH ${Matlab_MX_LIBRARY} DIRECTORY )
find_library(MXLIBUT NAMES libut ut 
             PATHS ${Matlab_MX_LIBRARY_PATH} REQUIRED)
matlab_add_mex(NAME minosMex
                SRC ${CMAKE_SOURCE_DIR}/interfaces/matlab/minosMex.cpp
                LINK_TO ${MXLIBUT} minoslib)
set(MINOSMEX_DESTDIR lib)
if (WIN32)
        set(MINOSMEX_DESTDIR bin)
endif()
install(TARGETS minosMex
        RUNTIME 
        DESTINATION ${MINOSMEX_DESTDIR}
	COMPONENT matlab)
# Install files
install(FILES buildOCP.m
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
cpack_add_install_type(matlab 
                DISPLAY_NAME "Matlab interface")
cpack_add_component(matlab 
                DISPLAY_NAME "Matlab interface"
                DESCRIPTION "MinOS MEX inteface"
                DEPENDS libraries
                INSTALL_TYPES full matlab)