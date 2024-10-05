# Build MinOS

cmake_minimum_required(VERSION 3.25)

# Set the project name etc
project(MinOS
        LANGUAGES C CXX
        VERSION 1.0.0
        DESCRIPTION "MINimal Optimal-control Solver"
        )

# User options
option(WITH_IPOPT "Compile IPOPT interface" ON)
option(WITH_WORHP "Compile WORHP interface" ON)
option(WITH_KNITRO "Compile KNITRO interface" ON)
option(WITH_SNOPT "Compile SNOPT interface" ON)
option(WITH_MEX_INTERFACE "Compile MEX interface for MATLAB" ON)
option(WITH_PYTHON_INTERFACE "Compile Python interface" ON)
option(WITH_DOCS "Generate docs using Doxygen (require Doxygen and latex)" ON)
option(WITH_CPACK "Generate CPack target" ON)
option(WITH_CPP_EXAMPLES "Compile with C++ examples" ON)
option(SKIP_HESSIAN "Use approximate Hessian calculation in C++ examples" OFF)
set(NLPSOLVER_EXAMPLE "ipopt" CACHE STRING "NLP solver for C++ examples (ipopt, knitro, worhp, snopt)")

# CPP_EXAMPLES_NLPSOL choices
set_property(CACHE NLPSOLVER_EXAMPLE PROPERTY STRINGS
             "ipopt" "knitro" "worhp" "snopt")

# Set output paths
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")

# Set copyright info
set(COPYRIGHT_YEAR "2024")
set(COPYRIGHT_NAME "Stefano Lovato")
set(COPYRIGHT_CONTACT "stefano.lovato@unipd.it")

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

# Default build type is Release
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (default Release)" FORCE)
endif()

# Find doxygen
if (WITH_DOCS)
    find_package(Doxygen REQUIRED)
    find_package(LATEX REQUIRED)
endif()

# Set the C++ standard to C++ 11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set C/C++ flags
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find NLP solver packages
if (WITH_IPOPT)
    # Find IPOPT
    find_package(IPOPT)
    if(NOT IPOPT_FOUND)
        message(WARNING "IPOPT not found, disabling WITH_IPOPT.")
        set(WITH_IPOPT OFF)
    endif()
endif()
if (WITH_WORHP)
    # Find IPOPT
    find_package(WORHP)
    if(NOT WORHP_FOUND)
        message(WARNING "WORHP not found, disabling WITH_IPOPT.")
        set(WITH_WORHP OFF)
    endif()
endif()
if (WITH_KNITRO)
    # Find IPOPT
    find_package(KNITRO)
    if(NOT KNITRO_FOUND)
        message(WARNING "KNITRO not found, disabling WITH_IPOPT.")
        set(WITH_KNITRO OFF)
    endif()
endif()
if (WITH_SNOPT)
    # Find IPOPT
    find_package(SNOPT)
    if(NOT SNOPT_FOUND)
        message(WARNING "SNOPT not found, disabling WITH_IPOPT.")
        set(WITH_SNOPT OFF)
    endif()
endif()

# MinOS sources with no solver source
set(MINOS_SOURCES ${CMAKE_SOURCE_DIR}/src/minos_base.cpp
                  ${CMAKE_SOURCE_DIR}/src/minos_nlp.cpp
                  ${CMAKE_SOURCE_DIR}/src/minos_print.cpp
                  ${CMAKE_SOURCE_DIR}/src/minos_utils.cpp
                  ${CMAKE_SOURCE_DIR}/src/minos_dlopen.cpp)

# Main MinOS base library 
add_library(minoslib SHARED
			${MINOS_SOURCES}
			)
target_compile_definitions(minoslib
							PRIVATE MINOS_VERSION=\"${CMAKE_PROJECT_VERSION}\")
target_include_directories(minoslib
							PUBLIC include
							)
set_target_properties(minoslib PROPERTIES
					OUTPUT_NAME minos)
set_target_properties(minoslib 
					PROPERTIES 
					PUBLIC_HEADER include/minos.h)
install(TARGETS minoslib
		COMPONENT libraries
		PUBLIC_HEADER
		COMPONENT libraries)

# IPOPT interface
if (WITH_IPOPT)
	# Create an interface library for IPOPT
	add_library(ipopt INTERFACE)
	target_include_directories(ipopt 
							INTERFACE ${IPOPT_INCLUDE_DIR})
	target_link_libraries(ipopt
						  INTERFACE ${IPOPT_LIBRARY})
	# Create IPOPT interface library
	add_library(minos-ipopt SHARED src/ipoptsolver.cpp)
	target_include_directories(minos-ipopt PRIVATE include)
	target_link_libraries(minos-ipopt PRIVATE ipopt minoslib)
	# Install target
	install(TARGETS minos-ipopt
			RUNTIME
			COMPONENT libraries)
	if (WIN32)
		install(FILES ${IPOPT_RUNTIME} 
				COMPONENT libraries
				DESTINATION bin)
	endif(WIN32)
endif()

# WORHP interface
if (WITH_WORHP)
	# Create an interface library for WORHP
	add_library(worhp INTERFACE)
	target_include_directories(worhp
							INTERFACE ${WORHP_INCLUDE_DIR})
	target_link_libraries(worhp 
						INTERFACE ${WORHP_LIBRARY})
	# Create WORHP interface library
	add_library(minos-worhp SHARED src/worhpsolver.cpp)
	target_include_directories(minos-worhp PRIVATE include)
	target_link_libraries(minos-worhp PRIVATE worhp minoslib)
	# Install target
	install(TARGETS minos-worhp
			RUNTIME
			COMPONENT libraries)
endif()

# KNITRO interface
if (WITH_KNITRO)
	# Create an interface library for KNITRO
	add_library(knitro INTERFACE)
	target_include_directories(knitro
							INTERFACE ${KNITRO_INCLUDE_DIR})
	target_link_libraries(knitro 
						INTERFACE ${KNITRO_LIBRARY})
	# Create KNITRO interface library
	add_library(minos-knitro SHARED src/knitrosolver.cpp)
	target_include_directories(minos-knitro PRIVATE include)
	target_link_libraries(minos-knitro PRIVATE knitro minoslib)
	# Install target
	install(TARGETS minos-knitro
			RUNTIME
			COMPONENT libraries)
endif()

# SNOPT interface
if (WITH_SNOPT)
	# Create an interface library for SNOPT
	add_library(snopt INTERFACE)
	target_include_directories(snopt
							INTERFACE ${SNOPT_INCLUDE_DIR})
	target_link_libraries(snopt 
						INTERFACE ${SNOPT_LIBRARY})
	# Create SNOPT interface library
	add_library(minos-snopt SHARED src/snoptsolver.cpp)
	target_include_directories(minos-snopt PRIVATE include)
	target_link_libraries(minos-snopt PRIVATE snopt minoslib)
	# Install target
	install(TARGETS minos-snopt
			RUNTIME
			COMPONENT libraries)
endif()

# Other installed files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Copyright.in
			   ${CMAKE_CURRENT_BINARY_DIR}/Copyright.txt
			   @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Copyright.txt
		COMPONENT libraries
		DESTINATION ./)

# Enable testing functionality
enable_testing()

# CPack
if (WITH_CPACK)

	# Set vars
	set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}")
	set(CPACK_PACKAGE_DESCRIPTION "${CMAKE_PROJECT_DESCRIPTION}")
	set(CPACK_PACKAGE_VENDOR "unipd")
	set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_BINARY_DIR}/Copyright.txt")
	set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/cmake\\\\minosLogo.bmp")
	if (WIN32)
		# NSIS
		set(CPACK_NSIS_MUI_ICON "${CMAKE_SOURCE_DIR}/cmake/minosLogo.ico")
		set(CPACK_NSIS_MUI_UNIICON "${CMAKE_SOURCE_DIR}/cmake/minosLogo.ico")
		set(CPACK_NSIS_MODIFY_PATH ON)
		# INNOSETUP
		set(CPACK_INNOSETUP_ICON_FILE "${CMAKE_SOURCE_DIR}/cmake/minosLogo.ico")
		# Wix
		set(CPACK_WIX_UPGRADE_GUID "FBF85747-FC84-4321-937A-A2686D63809A")
		set(CPACK_WIX_PRODUCT_ICON "${CMAKE_SOURCE_DIR}/cmake/minosLogo.ico")
	endif()

	# Include CPack
	include(CPack)
	
	# Installation types
	cpack_add_install_type(full
						DISPLAY_NAME "Full")
	cpack_add_install_type(dynamic 
						DISPLAY_NAME "C++ Libraries")
	
	# Add main application component
	cpack_add_component(libraries 
						DISPLAY_NAME "Application"
						DESCRIPTION "MinOS library"
						INSTALL_TYPES full dynamic matlab python
						REQUIRED)

endif()

# C++ examples
if (WITH_CPP_EXAMPLES)
	add_subdirectory(examples)
endif()

# MATLAB interface
if (WITH_MEX_INTERFACE)
	# Find the MATLAB package
	find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER)
	add_subdirectory(interfaces/matlab)
endif()

# Python interface
if (WITH_PYTHON_INTERFACE) 
	add_subdirectory(interfaces/python)
endif()

# Docs
if (WITH_DOCS) 
    # Configure files
	if(WIN32) # Windows readme
		set(README_IN "${CMAKE_CURRENT_SOURCE_DIR}/cmake/README.win.in")
	else() # Unix assumed otherwise
		set(README_IN "${CMAKE_CURRENT_SOURCE_DIR}/cmake/README.unix.in")
	endif()
    set(README_OUT "${CMAKE_CURRENT_BINARY_DIR}/README.md")
    configure_file(${README_IN} ${README_OUT} @ONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Doxyfile.in 
                    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile 
                    @ONLY)
    # Custom target to generate docs
    add_custom_target( doc_doxygen ALL
                        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                        COMMENT "Generating API documentation with Doxygen"
                        VERBATIM )
    # Install rules
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs 
			COMPONENT docs
			DESTINATION ./)

	# CPack
	cpack_add_component(docs
						DISPLAY_NAME "Documentation"
						DESCRIPTION "HTML documentation"
						INSTALL_TYPES full dynamic matlab python)	
endif()
