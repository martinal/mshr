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
// along with mshr. If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __TRIANGULATION_REFINEMENT_H
#define __TRIANGULATION_REFINEMENT_H

#include <dolfin/mesh/MeshEditor.h>

namespace 
{
  // Simple convenience class for converting from barycentric coordinates wrt to 
  // a reference triangle.
  class RefTriangle
  {
   public:
    RefTriangle(dolfin::Point v1, dolfin::Point v2, dolfin::Point v3)
      : v1(v1), v2(v2), v3(v3) {}

    dolfin::Point barycentric2Point(double l1, double l2, double l3)
    {
      const double x = l1*v1.x() + l2*v2.x() + l3*v3.x();
      const double y = l1*v1.y() + l2*v2.y() + l3*v3.y();
      const double z = l1*v1.z() + l2*v2.z() + l3*v3.z();

      return dolfin::Point(x, y, z);
    }

   private:
    dolfin::Point v1, v2 , v3;
  };
  //-----------------------------------------------------------------------------
  inline dolfin::Point get_edge_point(dolfin::Point a, dolfin::Point b, double f)
  {
    // print "Get point {} between {} and {}".format(f, a, b)
    const dolfin::Point e = b-a;
    return a + e*f;
  }
  //-----------------------------------------------------------------------------
  inline std::size_t get_edge_vertex(std::map<std::tuple<std::size_t, std::size_t, std::size_t>, 
                                       std::size_t> edge_vertices, 
                              std::size_t a, std::size_t b, std::size_t i, std::size_t N)
  {
    std::size_t v = edge_vertices[std::make_tuple(std::min(a, b), std::max(a, b), a < b ? i : N-i)];
    // print "getting edge vertex ({}, {}, {}) = {}".format(a, b, i, v)
    return v;
  }
  //-----------------------------------------------------------------------------
  inline void add_cell(dolfin::MeshEditor& editor, std::size_t cell_no, std::size_t v0, std::size_t v1, std::size_t v2)
  {
    // print "    Adding cell {}: {}, {}, {}".format(cell_count, v0, v1, v2)
    editor.add_cell(cell_no, v0, v1, v2);
  }
}
//-----------------------------------------------------------------------------
void refine_triangulation(const std::vector<dolfin::Point> initial_vertices,
                          const std::vector<std::tuple<std::size_t, std::size_t, std::size_t> > initial_triangulation,
                          std::size_t N,
                          std::vector<dolfin::Point>& vertices,
                          std::vector<std::tuple<std::size_t, std::size_t, std::size_t> >& triangles)
{
  const std::size_t num_vertices = initial_vertices.size() + 
    (N-1)*initial_triangulation.size() + (N+1)*N/2.0;
  const std::size_t num_triangles = N*N*initial_triangulation.size();

  vertices.clear();
  triangles.clear();
  vertices.reserve(num_vertices);
  triangles.reserve(num_triangles);
}

#endif
