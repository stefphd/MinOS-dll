find_package(Python
                COMPONENTS Development.Module
                REQUIRED)
find_package(Cython MODULE REQUIRED)
include(UseCython)
cython_transpile(${CMAKE_SOURCE_DIR}/interfaces/python/minosPy.pyx
                LANGUAGE CXX 
                OUTPUT_VARIABLE MINOSPYIPOPT)
python_add_library(minosPy MODULE "${MINOSPYIPOPT}")
target_link_libraries(minosPy PRIVATE minoslib)
install(TARGETS minosPy
        DESTINATION bin
	COMPONENT python)

# Install rules
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