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

#include <dolfin/log/log.h>
#include <vector>
#include <map>
#include <deque>
#include <iostream>

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

void SurfaceConsistency::checkConnectivity(const std::vector<std::array<std::size_t, 3> >& facets,
                                           std::set<std::size_t>& duplicating, bool error)
{
  // Store all halfedges
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> halfedges;

  std::size_t facet_no = 0;
  for (auto it = facets.begin(); it != facets.end(); ++it)
  {
    // Check for (topologically) degenerate facets
    if ( (*it)[0] == (*it)[1] || (*it)[0] == (*it)[2] || (*it)[1] == (*it)[2] )
      dolfin::dolfin_error("SurfaceConsistency.cpp",
                           "confirm surface connectivity",
                           "Facet %d is degenerate", facet_no);


    for (int i = 0; i < 3; i++)
    {
      std::pair<std::size_t, std::size_t> e( (*it)[i], (*it)[(i+1)%3] );
      if (halfedges.count( e ) > 0 )
      {
        if (error)
        {
          dolfin::dolfin_error("SurfaceConsistency.cpp",
                               "confirm halfedge connectivity",
                               "Facet %d and %d share halfedge", halfedges[e], facet_no);
        }
        else
        {
          duplicating.insert(facet_no);
        }
      }
      else
      {
        halfedges[e] = facet_no;
      }
    }
    facet_no++;
  }

  // TODO: Check for border edges
}

void SurfaceConsistency::filterFacets(const std::vector<std::array<std::size_t, 3> >& facets,
                                      const std::vector<std::array<double, 3> >& vertices,
                                      std::size_t start, std::set<std::size_t>& skip)
{
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

  std::cout << "Size of skip: " << skip.size() << std::endl;
  const std::size_t global_max = skip.size();

  std::set<std::size_t> visited;
  //std::vector<std::size_t> included;
  std::deque<std::size_t> queue;
  if (skip.count(start) > 0)
    queue.push_front(start);
  else
    std::cout << "  Already added" << std::endl;

  std::size_t global_count = 0;
  while (!queue.empty())
  {
    dolfin_assert(global_count <= global_max);

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


}

