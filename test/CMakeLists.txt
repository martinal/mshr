# Copyright (C) 2014 Benjamin Kehlet
#
# This file is part of mshr.
#
# mshr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# mshr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with mshr.  If not, see <http://www.gnu.org/licenses/>.

# This file defines simple regression tests. 


# test with a number of input files
message(STATUS "Preparing tests")

set(CGAL_DATA_DIR "${CMAKE_SOURCE_DIR}/3rdparty/CGAL-4.3/demo/Polyhedron/data")

# Files that should pass
set(TESTFILES_VALID "anchor.off;couplingdown.off;cross.off;cube.off;dragknob.off;elephant.off;ellipsoid.off;handle.off;icosahedron.off;joint.off;knot1.off;knot2.off;pinion.off;pipe.off;pyramid.off;rotor.off;sphere.off;spool.off;star.off;translated-cube.off;tripod.off")
set(TESTFILES_SELFINTERSECTING "bones.off;cow.off;man.off;oblong.ogg")
set(TESTFILES_NOTCLOSED "cube-ouvert.off;mushroom.off")

foreach(CURRENT_DATA_FILE IN LISTS TESTFILES_VALID)
  message(STATUS "Adding test: ${CURRENT_DATA_FILE}")
  #add_test("FilePass-${CURRENT_DATA_FILE}" "../mshr-surfacefile" "${CGAL_DATA_DIR}/${CURRENT_DATA_FILE}")
endforeach(CURRENT_DATA_FILE)

# Not sure this is the best way. Is saves the environment at "cmake
# time" and applies it when running the test. Changing PYTHONPATH
# later will not have effect. Not sure how to add tp PYTHONPATH at
# runtime.
set(PYTHON_ENVIR "${CMAKE_BINARY_DIR}/swig:$ENV{PYTHONPATH}")

add_test("Python-DummyImport" "python" "${CMAKE_SOURCE_DIR}/test/dummy.py")
set_property(TEST "Python-DummyImport" PROPERTY
    ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR} "
  )

file(GLOB PYTHONDEMOS "${CMAKE_SOURCE_DIR}/demo/python/*.py")
foreach(CURRENT_DEMO IN LISTS PYTHONDEMOS)
  add_test("PythonDemo-${CURRENT_DEMO}" "python" "${CURRENT_DEMO}")
  set_property(TEST "PythonDemo-${CURRENT_DEMO}" PROPERTY
    ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
  )
endforeach()
message(STATUS "Python demos: ${PYTHONDEMOS}")
