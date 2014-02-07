FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS 
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

ADD_DEFINITIONS( ${DOLFIN_PYTHON_DEFINITIONS} )

SET_SOURCE_FILES_PROPERTIES(mshr.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(mshr python mshr.i)
SWIG_LINK_LIBRARIES(mshr ${PYTHON_LIBRARIES} mshr )