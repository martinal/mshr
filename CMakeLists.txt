project( MSHR )

########################## Initialize cmake ############################
cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

# Borrow some cmake modules from cgal
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/3rdparty/CGAL-4.4/cmake/modules")


# Use C++11
SET(CMAKE_CXX_FLAGS "-std=c++0x -Wall")

# Boost
find_package( Boost REQUIRED system filesystem program_options )

if ( NOT Boost_FOUND )
  message(FATAL_ERROR "This project requires the Boost library, and will not be compiled.")
endif()

# Temporary: attempt to find linking problem
find_package(MPI REQUIRED)

#  GMP_INCLUDE_DIR       - the GMP include directory
#  GMP_LIBRARIES_DIR     - directory where the GMP libraries are located
#  GMP_LIBRARIES         - Link these to use GMP
find_package(GMP REQUIRED)

# Try to find the MPFR libraries
# MPFR_FOUND - system has MPFR lib
# MPFR_INCLUDE_DIR - the MPFR include directory
# MPFR_LIBRARIES_DIR - Directory where the MPFR libraries are located
# MPFR_LIBRARIES - the MPFR libraries
# MPFR_IN_CGAL_AUXILIARY - TRUE if the MPFR found is the one distributed with CGAL in the auxiliary folder
find_package(MPFR REQUIRED)

add_subdirectory(3rdparty)
#include_directories(BEFORE ${CGAL_INCLUDE_DIR})
include_directories(BEFORE ${EXTERNAL_INCLUDE_DIRS})
add_definitions("${EXTERNAL_DEFINITIONS}")

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# Set installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")


if(WIN32 AND NOT CYGWIN)
  set(INSTALL_CMAKE_DIR CMake)
else()
  set(INSTALL_CMAKE_DIR lib/CMake/mshr)
endif()

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()


# include for local header directory
include_directories( BEFORE include )

file( GLOB_RECURSE SOURCES src/*.cpp )

add_library(mshr SHARED ${SOURCES})
add_dependencies( mshr CGAL-4.4 )

# Link the executable to (static) CGAL libraries
target_link_libraries( mshr ${Boost_LIBRARIES}
                            ${GMP_LIBRARIES}
                            ${MPFR_LIBRARIES}
                            ${CGAL_LIB_DIR}/libCGAL.a
                            ${CGAL_LIB_DIR}/libCGAL_Core.a
			    tet
                            ${DOLFIN_LIBRARIES}
#                            ${DOLFIN_3RD_PARTY_LIBRARIES}
                     )

# Build the python bindings
add_subdirectory(swig)

# Fix CGAL's random generator (usefull to reproduce bugs)
option(INIT_CGAL_RANDOM "Seed for CGAL's random generator" OFF)
if (INIT_CGAL_RANDOM)
  add_definitions( -DINIT_RANDOM_GENERATOR=0 )
endif()

export(PACKAGE mshr)

add_executable( mshrable app/mshrable.cpp )
target_link_libraries( mshrable mshr 
                                ${Boost_LIBRARIES} 
				${MPI_CXX_LIBRARIES} )

# install library and app
install(TARGETS mshr mshrable
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# install header files
install(DIRECTORY include/ DESTINATION include)

message(STATUS "INSTALL DIR: ${INSTALL_CMAKE_DIR}")
message(STATUS "INCLUDE DIR: ${INSTALL_INCLUDE_DIR}")
message(STATUS "PROJECT SOURCE: ${PROJECT_SOURCE_DIR}")
message(STATUS "INSTALL CMAKE: ${INSTALL_CMAKE_DIR}")
message(STATUS "INSTALL INCLUDE: ${INSTALL_INCLUDE_DIR}")

file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
message(STATUS "REL INCLUDE: ${REL_INCLUDE_DIR}")

### add cmake configuration file ###

# common for both config files
set(CONF_EXTERNAL_INCLUDE_DIRS "${DOLFIN_INCLUDE_DIRS};${DOLFIN_3RD_PARTY_INCLUDE_DIRS};${DOLFIN_PYTHON_INCLUDE_DIRS}")
set(CONF_EXTERNAL_LIBRARIES "${DOLFIN_LIBRARIES};${DOLFIN_3RD_PARTY_LIBRARIES};${DOLFIN_PYTHON_LIBRARIES}")
set(CONF_CXX_DEFINITIONS "${DOLFIN_CXX_DEFINITIONS}")
string(REPLACE "\"" "\\\"" CONF_CXX_DEFINITIONS "${CONF_CXX_DEFINITIONS}")
set(CONF_CXX_FLAGS "${DOLFIN_CXX_FLAGS}")

# for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include")
set(CONF_LIBRARIES_DIRS "${PROJECT_BINARY_DIR}/lib")

configure_file(mshr-config.cmake.in "${PROJECT_BINARY_DIR}/mshr-config.cmake" @ONLY)
file(COPY use-mshr.cmake DESTINATION "${PROJECT_BINARY_DIR}/")



# ... and for the install tree
set(CONF_INCLUDE_DIRS "${INSTALL_INCLUDE_DIR}")
set(CONF_LIBRARIES_DIRS "${INSTALL_LIB_DIR}")
configure_file(mshr-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mshr-config.cmake" @ONLY)

# Install the FooBarConfig.cmake and FooBarConfigVersion.cmake
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mshr-config.cmake" DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
install(FILES "use-mshr.cmake" DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)


# install python module


option(ENABLE_TESTS "Enable testing" OFF)

if(ENABLE_TESTS)
  enable_testing()
  add_subdirectory(test)
else()
  message(STATUS "Testing disabled")
endif()