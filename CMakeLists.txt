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
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/3rdparty/CGAL-4.3/cmake/modules")



# Use C++11
SET(CMAKE_CXX_FLAGS "-std=c++0x -Wall")

# Boost
find_package( Boost REQUIRED system thread)

if ( NOT Boost_FOUND )
  message(FATAL_ERROR "This project requires the Boost library, and will not be compiled.")
endif()

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
include_directories(BEFORE ${CGAL_INCLUDE_DIR})

message(STATUS "Done with 3rd party libs")


find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# include for local header directory
include_directories( BEFORE include )

FILE( GLOB_RECURSE SOURCES src/*.cpp )

add_library(mshr SHARED ${SOURCES})
add_dependencies( mshr CGAL-4.3 )

# Link the executable to (static) CGAL libraries
target_link_libraries( mshr ${Boost_LIBRARIES}
                            ${GMP_LIBRARIES}
                            ${MPFR_LIBRARIES}
                            ${CGAL_LIB_DIR}/libCGAL.a
                            ${CGAL_LIB_DIR}/libCGAL_Core.a
                            ${DOLFIN_LIBRARIES}
                            ${DOLFIN_3RD_PARTY_LIBRARIES}
                     )


add_subdirectory(swig)

# Compile a simple program for testing
add_executable( testprogram2D test/test2.cpp )
target_link_libraries( testprogram2D ${Boost_LIBRARIES}
                                   ${GMP_LIBRARIES}
                                   ${MPFR_LIBRARIES}
				   ${DOLFIN_LIBRARIES}
				   ${DOLFIN_3RD_PARTY_LIBRARIES}
				   mshr
                     )

add_executable( testprogram3D test/test3D.cpp )
target_link_libraries( testprogram3D ${Boost_LIBRARIES}
                                   ${GMP_LIBRARIES}
                                   ${MPFR_LIBRARIES}
				   ${DOLFIN_LIBRARIES}
				   ${DOLFIN_3RD_PARTY_LIBRARIES}
				   mshr
                     )
