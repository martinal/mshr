

# Make sure cmake picks up the same swig as Dolfin used
if (NOT DEFINED DOLFIN_SWIG_EXECUTABLE )
  message(WARNING "DOLFIN_SWIG_EXECUTABLE not defined. This should be in DOLFINConfig.cmake")
else()
  set(SWIG_EXECUTABLE ${DOLFIN_SWIG_EXECUTABLE})
endif()

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# FIXME: Instead use DOLFIN_PYTHON_EXECUTABLE when
# this has been added to DOLFINConfig.cmake
find_package(PythonInterp 2)

# Don't use these. Instad get the values from Dolfin
# find_package(PythonLibs)
# include_directories(${PYTHON_INCLUDE_PATH})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_SWIG_FLAGS 
  -module mshr
  -shadow
  -modern
  -modernargs
  -fastdispatch
  -fvirtual
  -nosafecstrings
  -noproxydel
  -fastproxy
  -fastinit
  -fastunpack
  -fastquery
  -nobuildnone
  -Iinclude/swig
  -Wextra
  ${DOLFIN_CXX_DEFINITIONS}
  ${DOLFIN_PYTHON_DEFINITIONS}
)

add_definitions( ${DOLFIN_PYTHON_DEFINITIONS} )

set_source_files_properties(mshr.i PROPERTIES CPLUSPLUS ON)
swig_add_module(mshr python mshr.i)
swig_link_libraries(mshr ${DOLFIN_PYTHON_LIBRARIES} mshr )

# Python command string to discover module install location
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(PYTHON_LIB_DISCOVER_STR "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=True))")
else()
  set(PYTHON_LIB_DISCOVER_STR "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))")
endif()

# Probe Python interpreter
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "${PYTHON_LIB_DISCOVER_STR}"
  OUTPUT_VARIABLE INSTALL_PYTHON_MODULE_DIR
  )


set(INSTALL_PYTHON_MODULE_DIR ${INSTALL_PYTHON_MODULE_DIR}
  CACHE PATH "Python extension module installation directory.")

# Don't remove runtime link paths
set_target_properties(${SWIG_MODULE_mshr_REAL_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

# install compiled swig generated module
install(TARGETS ${SWIG_MODULE_mshr_REAL_NAME} DESTINATION ${INSTALL_PYTHON_MODULE_DIR})

# install pure python swig generated module
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mshr.py
  DESTINATION ${INSTALL_PYTHON_MODULE_DIR})