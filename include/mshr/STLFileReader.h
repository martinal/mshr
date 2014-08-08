// Copyright (C) 2014 Benjamin Kehlet
//
// This file is part of mshr.
//
// mshr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// mshr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mshr.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __STL_FILE_READER_H
#define __STL_FILE_READER_H

#include <string>
#include <vector>
#include <array>

#include <boost/tuple/tuple.hpp>

namespace mshr
{

class STLFileReader
{
public:

  /// @brief Parse STL file.
  /// @param filename The file to be read
  /// @param vertices array to return vertices in
  /// @param facets array to return facets in as indices to the vertex array
  static void read(const std::string filename, 
                   std::vector<std::array<double, 3> >& vertices,
                   std::vector<std::vector<std::size_t> >& facets);
};

}
#endif
