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

// OBS! Experimental code

#include <mshr/SurfaceConsistency.h>
#include "Point3FuzzyStrictlyLess.h"

#include <dolfin/log/log.h>
#include <vector>
#include <map>
#include <deque>
#include <iostream>
#include <limits>

// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef Kernel::Point_3 Point_3;
// typedef Kernel::Triangle_3 Triangle_3;

// namespace
// {
// inline bool triangles_intersect(const std::vector<std::array<double, 3> >& vertices,
//                                 const std::vector<std::array<std::size_t, 3> >& facets,
//                                 std::size_t a,
//                                 std::size_t b)
// {
//   const std::array<std::size_t, 3>& ta = facets[a];
//   const std::array<std::size_t, 3>& tb = facets[b];

//   std::cout << "  Intersection test: " << std::endl;
//   std::cout << "    Facet " << a << ": (" << ta[0] << ", " << ta[1] << ", " << ta[2] << ")" << std::endl;
//   std::cout << "    Facet " << b << ": (" << tb[0] << ", " << tb[1] << ", " << tb[2] << ")" << std::endl;

//   for (std::size_t i = 0; i < 3; i++)
//   {
//     for (std::size_t j = 0; j < 3; j++)
//     {
//       if (ta[i] == tb[j] || (ta[i] == tb[(j+1)%3] && ta[(i+1)%3] == tb[j]))
//       {
//         std::cout << "    Neighbor" << std::endl;
//         return false;
//       }
//     }
//   }

//   Triangle_3 t1(Point_3(vertices[ta[0]][0], vertices[ta[0]][1], vertices[ta[0]][2]),
//                 Point_3(vertices[ta[1]][0], vertices[ta[1]][1], vertices[ta[1]][2]),
//                 Point_3(vertices[ta[2]][0], vertices[ta[2]][1], vertices[ta[2]][2]));

//   Triangle_3 t2(Point_3(vertices[tb[0]][0], vertices[tb[0]][1], vertices[tb[0]][2]),
//                 Point_3(vertices[tb[1]][0], vertices[tb[1]][1], vertices[tb[1]][2]),
//                 Point_3(vertices[tb[2]][0], vertices[tb[2]][1], vertices[tb[2]][2]));

//   const bool i = CGAL::do_intersect(t1, t2);
//   std::cout << "    Result: " << (i ? "True" : "False") << std::endl;
//   return i;
// }
// }

namespace mshr
{

void SurfaceConsistency::checkConnectivity(std::vector<std::array<std::size_t, 3> >& facets,
                                           std::set<std::size_t>& duplicating,
                                           bool error)
{
  log(dolfin::TRACE, "Checking connectivity");

  // Store halfedges
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> halfedges;

  for (std::size_t facet_no = 0; facet_no < facets.size(); facet_no++)
  {
    std::array<std::size_t, 3>& f = facets[facet_no];
    // Check for (topologically) degenerate facets
    if ( f[0] == f[1] || f[0] == f[2] || f[1] == f[2] )
      dolfin::dolfin_error("SurfaceConsistency.cpp",
                           "confirm surface connectivity",
                           "Facet %d is topologically degenerate", facet_no);

    if (halfedges.count(std::make_pair(f[0], f[1])) > 0 ||
        halfedges.count(std::make_pair(f[1], f[2])) > 0 ||
        halfedges.count(std::make_pair(f[2], f[0])) > 0)
    {
      duplicating.insert(facet_no);
    }
    else
    {
      halfedges[std::make_pair(f[0], f[1])] = facet_no;
      halfedges[std::make_pair(f[1], f[2])] = facet_no;
      halfedges[std::make_pair(f[2], f[0])] = facet_no;
    }
  }
}
//-----------------------------------------------------------------------------
void SurfaceConsistency::filterFacets(const std::vector<std::array<std::size_t, 3> >& facets,
                                      const std::vector<std::array<double, 3> >& vertices,
                                      std::size_t start, std::set<std::size_t>& skip)
{
  // Map egdes (ordered pairs of vertices) to facets
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> edge_map;

  for (std::size_t i = 0; i < facets.size(); i++)
  {
    const std::array<std::size_t, 3>& facet = facets[i];
    std::size_t prev = facet[facet.size()-1];
    for (auto vit = facet.begin(); vit != facet.end(); vit++)
    {
      const std::pair<std::size_t, std::size_t> e(prev, *vit);
      edge_map[e] = i;
      prev = *vit;
    }
    
    skip.insert(i);
  }

  // const std::size_t global_max = skip.size();

  std::set<std::size_t> visited;
  //std::vector<std::size_t> included;
  std::deque<std::size_t> queue;
  if (skip.count(start) > 0)
    queue.push_front(start);
  // else
  //   std::cout << "  Already added" << std::endl;

  std::size_t global_count = 0;
  while (!queue.empty())
  {
    // dolfin_assert(global_count <= global_max);

    std::size_t current = queue.front();
    queue.pop_front();

    const std::array<size_t, 3>& current_facet = facets[current];
    // std::cout << "-- Processing " << current << ", vertices: " << current_facet[0] << ", " << current_facet[1] << ", " << current_facet[2] << std::endl;

    visited.insert(current);

    // don't skip this facet
    if (skip.count(current) > 0)
    {
      // bool intersects = false;
      // for (auto it = included.begin(); it != included.end(); it++)
      // {
      //   if (triangles_intersect(vertices, facets, current, *it))
      //   {
      //     intersects = true;
      //     break;
      //   }
      // }

      // if (intersects)
      // {
      //   std::cout << " SKIPPING" << std::endl;
      //   {int tmp; std::cin >> tmp;}
      //   continue;
      // }

      skip.erase(current);
      //included.push_back(current);

      std::size_t prev = current_facet[2];
      for (auto fit = current_facet.begin(); fit != current_facet.end(); fit++)
      {
        const std::pair<std::size_t, std::size_t> opposite(*fit, prev);
        // std::cout << "  opposite: " << opposite.first << " " << opposite.second << std::endl;
        if (edge_map.count(opposite) > 0)
        {
          std::size_t opposite_facet = edge_map[opposite];
          if (visited.count(opposite_facet) == 0)
          {
            queue.push_back(opposite_facet);
            //std::cout << "  pushing: " << edge_map[opposite] << std::endl;
          }
        }
        prev = *fit;
      }

      // std::cout << "Size of skip: " << skip.size() << ", size of queue: " << queue.size() << ", current: " << current << std::endl;
      // {int tmp; std::cin >> tmp;}
      global_count++;
    }
  }
}
//-----------------------------------------------------------------------------
std::pair<std::unique_ptr<std::vector<std::array<double, 3> > >,
          std::unique_ptr<std::vector<std::array<std::size_t, 3> > > >
SurfaceConsistency::merge_close_vertices(const std::vector<std::array<std::size_t, 3> >& facets,
                                         const std::vector<std::array<double, 3> >& vertices)
{
  std::map<std::array<double, 3>, std::size_t, Point3FuzzyStrictlyLess<std::array<double, 3> > > point_map;
  std::vector<std::size_t> vertex_mapping;
  vertex_mapping.reserve(vertices.size());


  for (std::size_t i = 0; i < vertices.size(); i++)
  {
    const std::array<double, 3>& v = vertices[i];
    if (point_map.count(v) > 0)
    {
      vertex_mapping.push_back(point_map[v]);
      // std::cout << "Close vertices: (" << v[0] << ", " << v[1] << ", " << v[2] << ") (" << vertices[point_map[v]][0] << ", " << vertices[point_map[v]][1] << ", " << vertices[point_map[v]][2] << ")" << std::endl;
    }
    else
    {
      vertex_mapping.push_back(i);
      point_map[vertices[i]] = i;
    }
  }

  // std::cout << "Distinct vertices: " << point_map.size() << std::endl;

  std::unique_ptr<std::vector<std::array<double, 3> > > new_vertices(new std::vector<std::array<double, 3> >);
  std::unique_ptr<std::vector<std::array<std::size_t, 3> > > new_facets(new std::vector<std::array<std::size_t, 3> >);

  return std::make_pair(std::move(new_vertices), std::move(new_facets));
}
//-----------------------------------------------------------------------------
void SurfaceConsistency::orient_component(std::vector<std::array<std::size_t, 3> >& facets,
                                          std::size_t start)
{
  log(dolfin::TRACE, "Checking facet orientation");

  // Map from edge (pair of vertices) to two triangles
  std::map<std::pair<std::size_t, std::size_t>, std::pair<std::size_t, std::size_t> > edge_map;
  for (std::size_t i = 0; i < facets.size(); i++)
  {
    const std::array<std::size_t, 3>& f = facets[i];
    for (std::size_t j = 0; j < 3; j++)
    {
      auto edge = std::make_pair(std::min(f[j], f[(j+1)%3]), std::max(f[j], f[(j+1)%3]));
      if (edge_map.count(edge) > 0)
        edge_map[edge].second = i;
      else
        edge_map[edge] = std::make_pair(i, std::numeric_limits<std::size_t>::max());
    }
  }

  std::deque<std::size_t> queue;
  std::set<std::size_t> visited;
  queue.push_back(start);

  std::size_t counter = 0;
  std::size_t flipped = 0;
  while (!queue.empty())
  {
    std::size_t current = queue.front();
    const std::array<std::size_t, 3>& current_facet = facets[current];
    queue.pop_front();

    if (visited.count(current) == 0)
    {
      counter++;
      visited.insert(current);
      
      for (std::size_t j = 0; j < 3; j++)
      {
        auto edge = std::make_pair(std::min(current_facet[j], current_facet[(j+1)%3]), std::max(current_facet[j], current_facet[(j+1)%3]));
        dolfin_assert(edge_map.count(edge) > 0);

        auto facet_pair = edge_map[edge];
        const std::size_t opposite = facet_pair.first == current ? facet_pair.second : facet_pair.first;

        if (opposite != std::numeric_limits<std::size_t>::max() && visited.count(opposite) == 0)
        {
          dolfin_assert(opposite != current);
          
          std::array<std::size_t, 3>& opposite_facet = facets[opposite];
          for (std::size_t k = 0; k < 3; k++)
          {
            const std::size_t a = opposite_facet[k];
            const std::size_t b = opposite_facet[(k+1)%3];
            if ( (a == current_facet[0] && b == current_facet[1]) ||
                 (a == current_facet[1] && b == current_facet[2]) ||
                 (a == current_facet[2] && b == current_facet[0]))
            {
              flipped++;
              std::swap(opposite_facet[0], opposite_facet[1]);
              break;
            }
          }
          queue.push_back(opposite);
        }
      }
    }
  }
  log(dolfin::TRACE, "Flipped %u triangles", flipped);
}
}

