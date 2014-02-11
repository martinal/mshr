// Copyright (C) 2012 Benjamin Kehlet
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

#ifndef __MESH_CLEAN_H
#define __MESH_CLEAN_H

#include <CGAL/basic.h>
#include <CGAL/Kernel/global_functions.h>

namespace mshr
{

//-----------------------------------------------------------------------------
template<typename Polyhedron>
static inline double
get_edge_length(typename Polyhedron::Halfedge::Halfedge_handle halfedge)
{
  return CGAL::to_double((halfedge->vertex()->point() -
    halfedge->opposite()->vertex()->point()).squared_length());
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
static inline double get_triangle_area(typename Polyhedron::Facet_handle facet)
{
  const typename Polyhedron::Halfedge_handle edge = facet->halfedge();
  const typename Polyhedron::Point_3 a = edge->vertex()->point();
  const typename Polyhedron::Point_3 b = edge->next()->vertex()->point();
  const typename Polyhedron::Point_3 c
    = edge->next()->next()->vertex()->point();

  return CGAL::to_double(CGAL::cross_product(b-a, c-a).squared_length());
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
static inline double
get_min_edge_length(typename Polyhedron::Facet_handle facet)
{
  typename Polyhedron::Facet::Halfedge_around_facet_circulator half_edge
    = facet->facet_begin();
  double min_length = CGAL::to_double((half_edge->vertex()->point()
      - half_edge->opposite()->vertex()->point()).squared_length());

  half_edge++;
  min_length = std::min(min_length, get_edge_length<Polyhedron>(half_edge));

  half_edge++;
  min_length = std::min(min_length, get_edge_length<Polyhedron>(half_edge));

  return min_length;
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
bool facet_is_degenerate(typename Polyhedron::Facet_handle facet,
                         const double threshold)
{
  return get_min_edge_length<Polyhedron>(facet) < threshold
      || get_triangle_area<Polyhedron>(facet) < threshold;
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
static int number_of_degenerate_facets(Polyhedron& p, const double threshold)
{
  int count = 0;
  for (typename Polyhedron::Facet_iterator facet = p.facets_begin();
       facet != p.facets_end(); facet++)
  {
    dolfin_assert(facet->is_triangle());
    if ( facet_is_degenerate<Polyhedron>(facet, threshold) )
      count++;
  }
  return count;
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
static typename Polyhedron::Halfedge_handle
get_longest_edge(typename Polyhedron::Facet_handle facet)
{
  typename Polyhedron::Halfedge_handle edge = facet->halfedge();
  double length = get_edge_length<Polyhedron>(edge);

  {
    typename Polyhedron::Halfedge_handle e_tmp = edge->next();
    if (get_edge_length<Polyhedron>(e_tmp) > length)
    {
      length = get_edge_length<Polyhedron>(e_tmp);
      edge = e_tmp;
    }
  }

  {
    typename Polyhedron::Halfedge_handle e_tmp = edge->next()->next();
    if ( get_edge_length<Polyhedron>(e_tmp) > length )
    {
      length = get_edge_length<Polyhedron>(e_tmp);
      edge = e_tmp;
    }
  }

  return edge;
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
double shortest_edge(Polyhedron& p)
{
  double shortest = std::numeric_limits<double>::max();
  for (typename Polyhedron::Halfedge_iterator halfedge = p.halfedges_begin();
       halfedge != p.halfedges_end(); halfedge++)
  {
    const double length = get_edge_length<Polyhedron>(halfedge);
    shortest = std::min(shortest, length);
  }

  return shortest;
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
static void remove_edge(Polyhedron& p, typename
                        Polyhedron::Halfedge_handle& edge)
{

  // // FIXME: Is it possible to do this in a smarter way than a linear scan
  // for (csg::Polyhedron_3::Facet_iterator facet = p.facets_begin();
  //      facet != p.facets_end(); facet++)
  // {
  //   if ( facet_is_degenerate<csg::Polyhedron_3>(facet, threshold) )
  //   {
  //     //print_facet(facet);

  //     // Find a short edge
  //     csg::Polyhedron_3::Halfedge::Halfedge_handle shortest_edge = facet->facet_begin();
  //     csg::Polyhedron_3::Facet::Halfedge_around_facet_circulator current_edge = facet->facet_begin();
  //     double min_length = get_edge_length(current_edge);

  //     for (int i = 0; i < 2; i++)
  //     {
  // 	current_edge++;
  // 	if (get_edge_length(current_edge) < min_length)
  // 	{
  // 	  shortest_edge = current_edge;
  // 	  min_length = get_edge_length(current_edge);
  // 	}
  //     }

  // Join small triangles with neighbor facets
  edge = p.join_facet(edge->next());
  p.join_facet(edge->opposite()->prev());

  // The joined facets are now quads
  // Join the two close vertices
  p.join_vertex(edge);
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
static void remove_short_edges(Polyhedron& p, const double threshold)
{
  while (true)
  {
    bool removed = false;
    for (typename Polyhedron::Halfedge_iterator halfedge = p.halfedges_begin();
	 halfedge != p.halfedges_end(); halfedge++)
    {
      if (get_edge_length<Polyhedron>(halfedge) < threshold)
      {
	remove_edge<Polyhedron>(p, halfedge);
	removed = true;
	break;
      }
    }

    if (!removed)
      break;
  }
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
static typename Polyhedron::Point_3
facet_midpoint(typename Polyhedron::Facet_handle facet)
{
  typename Polyhedron::Point_3 p(CGAL::ORIGIN);

  typename Polyhedron::Facet::Halfedge_around_facet_circulator half_edge
    = facet->facet_begin();

  for (std::size_t i = 0; i < facet->facet_degree(); i++)
  {
    p = p + (half_edge->vertex()->point() - CGAL::ORIGIN);
    half_edge++;
  }

  p = CGAL::ORIGIN
    + (p - CGAL::ORIGIN)/static_cast<double>(facet->facet_degree());

  // std::dolfin::cout << "Center coordinates computed: " << p << std::dolfin::endl;

  // half_edge = facet->facet_begin();
  // for (std::size_t i = 0; i < facet->facet_degree(); i++)
  // {
  //   std::dolfin::cout << "Distance to point << " << half_edge->vertex()->point() << " = " << (half_edge->vertex()->point() - p).squared_length() << std::dolfin::endl;
  //   half_edge++;
  // }

  return p;
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
static void
remove_triangle(Polyhedron& p, typename Polyhedron::Facet_handle facet)
{
  dolfin_assert(facet->is_triangle());

  // dolfin::cout << "Removing triangle" << dolfin::endl;
  // print_facet<Polyhedron>(facet);

  // Find the longest edge
  typename Polyhedron::Halfedge_handle edge
    = get_longest_edge<Polyhedron>(facet);

  // dolfin::cout << "Longest edge" << dolfin::endl;
  // print_halfedge<Polyhedron>(edge);

  // dolfin::cout << "Opposite triangle" << dolfin::endl;
  // print_facet<Polyhedron>(edge->opposite()->facet());

  edge = p.join_facet(edge);
  // dolfin::cout << "Edge after join: " << dolfin::endl;
  // print_halfedge<Polyhedron>(edge);

  // dolfin::cout << "Facet after join" << dolfin::endl;
  // print_facet<Polyhedron>(edge->facet());

  typename Polyhedron::Point_3 new_center
    = facet_midpoint<Polyhedron>(edge->facet());

  edge = p.create_center_vertex(edge);

  edge->vertex()->point() = new_center;

  // std::dolfin::cout << "Center vertex: " << edge->vertex()->point() << std::dolfin::endl;

  // for (std::size_t i=0; i < 4; i++)
  // {
  //   print_facet<Polyhedron>(edge->facet());
  //   edge = edge->next()->opposite();
  // }
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
static void remove_small_triangles(Polyhedron& p, const double threshold)
{
  int n = number_of_degenerate_facets(p, threshold);

  while (n > 0)
  {
    for (typename Polyhedron::Facet_iterator facet = p.facets_begin();
	 facet != p.facets_end(); facet++)
    {
      dolfin_assert(facet->is_triangle());

      if (get_triangle_area<Polyhedron>(facet) < threshold)
      {
	// dolfin::cout << "Small triangle detected" << dolfin::endl;
	// print_facet<Polyhedron>(facet);
	remove_triangle<Polyhedron>(p, facet);
	n = number_of_degenerate_facets<Polyhedron>(p, threshold);
	break;
      }
    }
  }
}

//-----------------------------------------------------------------------------
template<typename Polyhedron>
static bool has_degenerate_facets(const Polyhedron& p,
                                  double threshold)
{
  for (typename Polyhedron::Facet_iterator facet = p.facets_begin();
       facet != p.facets_end(); facet++)
  {
    dolfin_assert(facet->is_triangle());
    if (facet_is_degenerate<Polyhedron>(facet, threshold))
      return true;
  }
  return false;
}
//-----------------------------------------------------------------------------
}

#endif
