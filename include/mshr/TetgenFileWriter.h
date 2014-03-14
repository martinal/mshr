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

#ifndef __MSHR_TETGENOUTPUT_H_
#define __MSHR_TETGENOUTPUT_H_

#include <dolfin/geometry/Point.h>

#include <string>
#include <vector>
#include <array>


namespace mshr
{

// Forward declaration
class CSGCGALDomain3D;


/// This class can output a geometry in .poly format that Tetgen understands.
class TetgenFileWriter
{
 public:
  static void write(const CSGCGALDomain3D &domain, std::string filename);
  static void write(const std::vector<dolfin::Point> &vertices,
                    const std::vector<std::array<std::size_t, 3> > &facets,
                    std::string filename);
};

}

#endif

