# Find LUA
find_package(LUA REQUIRED)
# Create interface library for LUA
add_library(lua INTERFACE)
target_include_directories(lua INTERFACE ${LUA_INCLUDE_DIR})
target_link_libraries(lua INTERFACE ${LUA_LIBRARIES})
# Command-line interface target
add_executable(minos-cli ${CMAKE_CURRENT_SOURCE_DIR}/minos_cli.cpp)
target_include_directories(minos-cli PRIVATE include)
target_link_libraries(minos-cli PRIVATE minos lua)
if (MSVC)
        # This is to fix warning "LINK : warning LNK4098: defaultlib 'LIBCMT' conflicts with use of other libs" 
        # in MSVC when linking static LUA lib
        target_link_options(minos-cli PRIVATE "/NODEFAULTLIB:LIBCMT")
endif()
# Install target
install(TARGETS minos-cli
        RUNTIME
        COMPONENT cli)
install(DIRECTORY ./examples/
        COMPONENT cli
        DESTINATION ./examples/cli
        PATTERN *.dll EXCLUDE
        PATTERN *.so EXCLUDE
        PATTERN *.c EXCLUDE)
# CPack
if (WITH_CPACK)
    cpack_add_install_type(cli 
                    DISPLAY_NAME "Command-line interface")
    cpack_add_component(cli 
                    DISPLAY_NAME "Command-line interface"
                    DESCRIPTION "Include the command-line inteface"
                    DEPENDS libraries
                    INSTALL_TYPES full cli)
endif()