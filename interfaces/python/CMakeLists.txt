find_package(Python3
                COMPONENTS Development.Module
                REQUIRED)
find_package(Cython MODULE REQUIRED)
include(UseCython)
cython_transpile(${CMAKE_CURRENT_SOURCE_DIR}/minosPy.pyx
                LANGUAGE CXX 
                OUTPUT_VARIABLE MINOSPY)
python3_add_library(minos-py MODULE ${MINOSPY})
set_target_properties(minos-py PROPERTIES OUTPUT_NAME minosPy)
target_link_libraries(minos-py PRIVATE minos)
# Configure files
configure_file(${CMAKE_SOURCE_DIR}/cmake/minosPy.in
        ${CMAKE_CURRENT_BINARY_DIR}/minosPy.pyi
        @ONLY)
# Install rules
set(MINOSPY_DESTDIR lib)
if (WIN32)
        set(MINOSPY_DESTDIR bin)
endif()
install(TARGETS minos-py
        DESTINATION ${MINOSPY_DESTDIR}
	COMPONENT python)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/minosPy.pyi
        COMPONENT python
        DESTINATION ${MINOSPY_DESTDIR})
install(FILES requirements.txt
        COMPONENT python
        DESTINATION ./)
install(DIRECTORY ./examples/
        COMPONENT python
        DESTINATION ./examples/python
        PATTERN *pxd EXCLUDE
        PATTERN *pyx EXCLUDE
        PATTERN CMakeLists.txt EXCLUDE
        PATTERN examples/*.c EXCLUDE
        PATTERN examples/*.cpp EXCLUDE
        PATTERN examples/*.pyd EXCLUDE
        PATTERN examples/*.pyi EXCLUDE
        PATTERN examples/*.so EXCLUDE
        PATTERN examples/*.log EXCLUDE
        PATTERN examples/*.txt EXCLUDE
        PATTERN old-minos.py EXCLUDE
        PATTERN __pycache__ EXCLUDE
        PATTERN examples/build EXCLUDE)

# CPack
cpack_add_install_type(python 
                    DISPLAY_NAME "Python interface")
cpack_add_component(python
                    DISPLAY_NAME "Python interface"
                    DESCRIPTION "MinOS Python interface"
                    DEPENDS libraries
                    INSTALL_TYPES full python)