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
  inline std::size_t get_edge_vertex(const std::map<std::array<std::size_t,3 >, 
                                     std::size_t>& edge_vertices, 
                                     std::size_t a, 
                                     std::size_t b, 
                                     std::size_t i, 
                                     std::size_t N)
  {
    const std::array<std::size_t, 3> key{std::min(a, b), std::max(a, b), a < b ? i : N-i};
    std::size_t v = edge_vertices.at(key);
    // print "getting edge vertex ({}, {}, {}) = {}".format(a, b, i, v)
    return v;
  }
  //-----------------------------------------------------------------------------
  inline void add_cell(std::vector<std::array<std::size_t, 3> >& triangles, std::array<std::size_t, 3> v)
  {
    std::cout << "    Adding cell " << v[0] << ", " << v[1] << ", " << v[2] << std::endl;
    triangles.push_back(v);
  }
}
//-----------------------------------------------------------------------------
// TODO: Use outputiterator to return triangulation instead of vectors to avoid 
// copying
void refine_triangulation(const std::vector<dolfin::Point> initial_vertices,
                          const std::vector<std::array<std::size_t, 3> > initial_triangulation,
                          std::size_t N,
                          std::vector<dolfin::Point>& vertices,
                          std::vector<std::array<std::size_t, 3> >& triangles)
{
  const std::size_t num_vertices = initial_vertices.size() + 
    (N-1)*initial_triangulation.size() + (N+1)*N/2.0;
  const std::size_t num_triangles = N*N*initial_triangulation.size();

  vertices.clear();
  triangles.clear();
  vertices.reserve(num_vertices);
  triangles.reserve(num_triangles);

  std::size_t vertex_count = 0;

  // Add the corner vertices
  for (const dolfin::Point& p : initial_vertices)
  {
    vertices.push_back(p/p.norm());
    vertex_count += 1;
  }

  std::map<std::array<std::size_t, 3>, std::size_t> edge_vertices;

  // Add the "inner" vertices along the edges
  for (std::size_t i=1; i < N; i++)
  {
    for (const std::array<std::size_t, 3>& t : initial_triangulation)
    {
      for (std::size_t j = 0; j < 3; j++)
      {
        if (t[j] < t[(j+1)%3])
        {
          std::cout << "Adding vertex: " << t[j] << ", " << t[(j+1)%3] << std::endl;
          dolfin::Point v = get_edge_point(initial_vertices[t[j]], initial_vertices[t[(j+1)%3]], float(i)/N);
          vertices.push_back(v/v.norm());
          // TODO: Use vertices.size() ?
          edge_vertices[{t[j], t[(j+1)%3], i}] = vertex_count;
          vertex_count += 1;
        }
      }
    }
  }


  /***************************** Add the triangles *********************************/
  std::size_t cell_count = 0;
  for (const std::array<std::size_t, 3>& triangle : initial_triangulation) 
  {
    std::cout << "Processing triangle (" << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << ")" << std::endl;
    RefTriangle ref_triangle(initial_vertices[triangle[0]], initial_vertices[triangle[1]], initial_vertices[triangle[2]]);
    const std::size_t vertex_start = vertex_count;

    for (std::size_t i = 1; i < N; i++)
    {
      std::cout << "i = " << i << std::endl;

      const double l1 = static_cast<double>(i)/N;
      for (std::size_t j=1; j < N-i; j++)
      {
        std::cout << "  j = " << j << std::endl;

        // Don't edge along initial vertices
        if (i+j == N)
          continue;

        const double l2 = static_cast<double>(j)/N;
        const double l3 = 1.0 - l1 - l2;
            
       dolfin::Point p = ref_triangle.barycentric2Point(l1, l2, l3);
        p /= p.norm();
        vertices.push_back(p);
        vertex_count += 1;
      }
    }

    //  add the "corner" facets
    std::cout << "  Adding corner facets" << std::endl;
    add_cell(triangles, { triangle[0],
                          get_edge_vertex(edge_vertices, 
                                          triangle[0],
                                          triangle[1], 
                                          1, N),
                          get_edge_vertex(edge_vertices,
                                          triangle[0], triangle[2],
                                          1, N)});
    cell_count += 1;

    add_cell(triangles, {triangle[1],
                         get_edge_vertex(edge_vertices,
                                         triangle[1], 
                                         triangle[2],
                                         1, N),
                         get_edge_vertex(edge_vertices, 
                                         triangle[1],
                                         triangle[0],
                                         1, N)});
    cell_count += 1;

    add_cell(triangles, {triangle[2],
                         get_edge_vertex(edge_vertices,
                                         triangle[2], 
                                         triangle[0], 
                                         1, N),
                         get_edge_vertex(edge_vertices,
                                         triangle[2], 
                                         triangle[1],
                                         1, N)});
    cell_count += 1;


    if (N == 2)
    {
      std::cout << "  Refinement level 2: Adding interior cell" << std::endl;
      add_cell(triangles, {get_edge_vertex(edge_vertices, triangle[0], triangle[1], 1, N),
                           get_edge_vertex(edge_vertices, triangle[1], triangle[2], 1, N),
                           get_edge_vertex(edge_vertices, triangle[2], triangle[0], 1, N)});
      cell_count += 1;
    }
    else
    {
      std::cout << "  Add facets incident to original edges" << std::endl;
      add_cell(triangles, {vertex_start,
            get_edge_vertex(edge_vertices, triangle[2], triangle[1], 1, N),
            get_edge_vertex(edge_vertices, triangle[2], triangle[0], 1, N)});
      cell_count += 1;

      add_cell(triangles, {vertex_start+N-3,
            get_edge_vertex(edge_vertices, triangle[1], triangle[0], 1, N),
            get_edge_vertex(edge_vertices, triangle[1], triangle[2], 1, N)});
      cell_count += 1;

      add_cell(triangles, {vertex_count - 1,
            get_edge_vertex(edge_vertices,triangle[0], triangle[2], 1, N),
            get_edge_vertex(edge_vertices,triangle[0], triangle[1], 1, N)});
      cell_count += 1;


      std::cout << "  Add the facets along the original edges" << std::endl;

      // along edge(triangle[1] <--> triangle[2])
      add_cell(triangles, {vertex_start,
            get_edge_vertex(edge_vertices,triangle[2], triangle[1], 2, N),
            get_edge_vertex(edge_vertices,triangle[2], triangle[1], 1, N)});
      cell_count += 1;

      for (std::size_t i = 3; i < N; i++)
      {
        add_cell(triangles, {vertex_start+i-3,
              vertex_start+i-2,
              get_edge_vertex(edge_vertices,triangle[2], triangle[1], i-1, N)});
        cell_count += 1;

        add_cell(triangles, {vertex_start+i-2,
              get_edge_vertex(edge_vertices,triangle[2], triangle[1], i, N),
              get_edge_vertex(edge_vertices,triangle[2], triangle[1], i-1, N)});
        cell_count += 1;
      }
            
      // along edge(triangle[2], triangle[0])
      add_cell(triangles, {vertex_start,
            get_edge_vertex(edge_vertices,triangle[2], triangle[0], 1, N),
            get_edge_vertex(edge_vertices,triangle[2], triangle[0], 2, N)});
      cell_count += 1;

      std::size_t current_inner_vertex = 0;
      for (std::size_t i = 3; i < N; i++)
      {
        const std::size_t step = N-2-i+3;
        add_cell(triangles, {vertex_start+current_inner_vertex+step,
              vertex_start+current_inner_vertex,
              get_edge_vertex(edge_vertices,triangle[2], triangle[0], i-1, N)});
        cell_count += 1;

        add_cell(triangles, {get_edge_vertex(edge_vertices,triangle[2], triangle[0], i-1, N),
              get_edge_vertex(edge_vertices,triangle[2], triangle[0], i, N),
              vertex_start+current_inner_vertex+step});
        cell_count += 1;
        
        current_inner_vertex += step;
      }

      // along edge(triangle[1], triangle[0])
      add_cell(triangles, {vertex_start+N-3,
            get_edge_vertex(edge_vertices,triangle[1], triangle[0], 2, N),
            get_edge_vertex(edge_vertices,triangle[1], triangle[0], 1, N)});
      cell_count += 1;

      current_inner_vertex = N-3;
      for (std::size_t i = 3; i < N; i++)
      {
        const std::size_t step = N-3-i+3;
        add_cell(triangles, {vertex_start+current_inner_vertex,
              vertex_start+current_inner_vertex+step,
              get_edge_vertex(edge_vertices,triangle[1], triangle[0], i-3+2, N)});
        cell_count += 1;

        add_cell(triangles, {vertex_start+current_inner_vertex+step,
              get_edge_vertex(edge_vertices,triangle[1], triangle[0], i-3+3, N),
              get_edge_vertex(edge_vertices,triangle[1], triangle[0], i-3+2, N)});
        cell_count += 1;
            
        current_inner_vertex += step;
      }
    }
    
    std::cout << "  Add the inner facets that don't touch initial edges" << std::endl;
    
    std::size_t row_offset = 0;
    for (std::size_t i = 1; i < N-2; i++)
    {
      const std::size_t row_length = N-1-i;
      add_cell(triangles, {vertex_start+row_offset,
            vertex_start+row_offset+row_length,
            vertex_start+row_offset+1});
      cell_count += 1;

      for (std::size_t j = 0; j < N-3-i; j++)
      {
        add_cell(triangles, {vertex_start+row_offset+row_length+j,
              vertex_start+row_offset+row_length+j+1,
              vertex_start+row_offset+j+1});
        cell_count += 1;

        add_cell(triangles, {vertex_start+row_offset+row_length+j+1,
              vertex_start+row_offset+j+2,
              vertex_start+row_offset+j+1});
        cell_count += 1;

        row_offset += row_length;
      }
    }
  }
}


#endif
