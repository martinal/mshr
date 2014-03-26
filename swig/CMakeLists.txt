find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})

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
swig_link_libraries(mshr ${PYTHON_LIBRARIES} mshr )

