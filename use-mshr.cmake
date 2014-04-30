#
# This file sets up include directories, link directories, and
# compiler settings for a project to use mshr. It should not be
# included directly, but rather through the DOLFIN_USE_FILE setting
# obtained from mshr-config.cmake.
#

if (NOT MSHR_USE_FILE_INCLUDED)
  set(MSHR_USE_FILE_INCLUDED 1)

  # Add compiler definitions needed to use mshr
  add_definitions(${MSHR_CXX_DEFINITIONS})

  # Add compiler flags needed to use MSHR
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MSHR_CXX_FLAGS}")

  # Add include directories needed to use MSHR
  include_directories(${MSHR_INCLUDE_DIRS})
  include_directories(SYSTEM ${MSHR_EXTERNAL_INCLUDE_DIRS})

  # Add link directories needed to use MSHR
  message(STATUS "Setting link directories")
  link_directories(${MSHR_LIBRARIES_DIRS})
endif()
