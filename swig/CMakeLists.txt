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
  ${DOLFIN_CXX_DEFINITIONS}
)



set(CMAKE_SWIG_FLAGS "-Wextra;${CMAKE_SWIG_FLAGS} ${DOLFIN_PYTHON_DEFINITIONS}")

message(STATUS "Dolfin cxx defs: ${DOLFIN_CXX_DEFINITIONS}")
message(STATUS "Pythond defs: ${DOLFIN_PYTHON_DEFINITIONS}")
#add_definitions("${DOLFIN_PYTHON_DEFINITIONS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DOLFIN_PYTHON_DEFINITIONS}")


SET_SOURCE_FILES_PROPERTIES(mshr.i PROPERTIES CPLUSPLUS ON)
#SET_SOURCE_FILES_PROPERTIES(mshr.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE(mshr python mshr.i)
SWIG_LINK_LIBRARIES(mshr ${PYTHON_LIBRARIES} mshr )