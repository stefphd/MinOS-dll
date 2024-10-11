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
# Install target
install(TARGETS minos-cli
        RUNTIME
        COMPONENT cli)
# CPack
if (WITH_CPACK)
    cpack_add_install_type(cli 
                    DISPLAY_NAME "Command-line interface")
    cpack_add_component(cli 
                    DISPLAY_NAME "Command-line interface"
                    DESCRIPTION "MinOS command-line inteface"
                    DEPENDS libraries
                    INSTALL_TYPES full cli)
endif()