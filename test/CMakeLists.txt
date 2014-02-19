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
# For now we just test


# test with a number of input files
message(STATUS "Preparing tests")

string(COMPARE NOTEQUAL "${TEST_DATA_DIR}" "" RUN_TESTDATA_TESTS)

if (RUN_TESTDATA_TESTS)
  message(STATUS "Enabling test with data in ${TEST_DATA_DIR}")

  file(GLOB TEST_DATA_FILES "${TEST_DATA_DIR}/*")
  message(STATUS "Data files: ${TEST_DATA_FILES}")

  foreach(CURRENT_DATA_FILE IN LISTS TEST_DATA_FILES)
    message(STATUS "Adding test: ${CURRENT_DATA_FILE}")
    add_test("TestFileInput_${CURRENT_DATA_FILE}" "mshr-surfacefile" ${CURRENT_DATA_FILE})
  endforeach(CURRENT_DATA_FILE)

else()
  message(STATUS "No path to test given. Skipping tests.")
endif()
