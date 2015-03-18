// Copyright (C) 2015 Benjamin Kehlet
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

#ifndef HOLE_TRIANGULATE_H__
#define HOLE_TRIANGULATE_H__

#include <dolfin/geometry/Point.h>
#include <vector>


namespace mshr
{
class HoleTriangulate
{
 public:
  static void triangulate_polygon_3D(const std::vector<dolfin::Point> polygon);

};
}
#endif
