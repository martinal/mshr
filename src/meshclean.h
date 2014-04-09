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

#ifndef __MESH_CLEAN_H
#define __MESH_CLEAN_H

#include <CGAL/basic.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Triangle_3.h>

//-----------------------------------------------------------------------------
template<typename Polyhedron>
inline double
get_edge_length(typename Polyhedron::Halfedge_const_handle halfedge)
{
  return CGAL::to_double((halfedge->vertex()->point() -
    halfedge->opposite()->vertex()->point()).squared_length());
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
inline double get_triangle_area(typename Polyhedron::Facet_handle facet)
{
  typedef typename Polyhedron::Traits::Triangle_3 Triangle;
  typename Polyhedron::Halfedge_const_handle h = facet->halfedge();
  Triangle t(h->vertex()->point(),
             h->next()->vertex()->point(),
             h->next()->next()->vertex()->point());
  return CGAL::to_double(t.squared_area());
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
inline typename Polyhedron::Halfedge_handle
  get_longest_edge(typename Polyhedron::Facet_handle facet)
{
  typename Polyhedron::Halfedge_handle current = facet->halfedge();
  typename Polyhedron::Halfedge_handle longest = current;
  double length = get_edge_length<Polyhedron>(current);

  current = current->next();
  if (get_edge_length<Polyhedron>(current) > length)
  {
    length = get_edge_length<Polyhedron>(current);
    longest = current;
  }

  current = current->next();
  if (get_edge_length<Polyhedron>(current) > length)
  {
    length = get_edge_length<Polyhedron>(current);
    longest = current;
  }

  return longest;
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
inline typename Polyhedron::Halfedge_const_handle
  get_longest_const_edge(typename Polyhedron::Facet_const_handle facet)
{
  typename Polyhedron::Halfedge_const_handle current = facet->halfedge();
  typename Polyhedron::Halfedge_const_handle longest = current;
  double length = get_edge_length<Polyhedron>(current);

  current = current->next();
  if (get_edge_length<Polyhedron>(current) > length)
  {
    length = get_edge_length<Polyhedron>(current);
    longest = current;
  }

  current = current->next();
  if (get_edge_length<Polyhedron>(current) > length)
  {
    length = get_edge_length<Polyhedron>(current);
    longest = current;
  }

  return longest;
}
//-----------------------------------------------------------------------------
// Compute the distance from the longest edge to the opposite vertex
template<typename Polyhedron>
inline double triangle_projection(typename Polyhedron::Facet_const_handle f)
{
  typename Polyhedron::Halfedge_const_handle longest = get_longest_const_edge<Polyhedron>(f);
  typename Polyhedron::Traits::Line_3 l(longest->vertex()->point(), longest->opposite()->vertex()->point());
  return CGAL::to_double( (l.projection(longest->next()->vertex()->point()) - longest->next()->vertex()->point()).squared_length());
}
//-----------------------------------------------------------------------------
// Print some info about a facet
// Meant for debugging
template<typename Polyhedron>
void print_facet(typename Polyhedron::Halfedge_const_handle h)
{
  typename Polyhedron::Halfedge_const_handle h2 = h->next(), h3 = h2->next();
  std::cout << "Vertices: " << std::endl;
  std::cout << "(" << h->vertex()->point() << ") "
            << "(" << h2->vertex()->point() << ") "
            << "(" << h3->vertex()->point() << ")" << std::endl;
  std::cout << "Edge lengths (squared):"
            << (h->vertex()->point() - h2->vertex()->point()).squared_length() << " "
            << (h2->vertex()->point() - h3->vertex()->point()).squared_length() << " "
            << (h3->vertex()->point() - h->vertex()->point()).squared_length() << std::endl;
  typename Polyhedron::Traits::Triangle_3 t(h->vertex()->point(),
                                            h2->vertex()->point(),
                                            h3->vertex()->point());
  std::cout << "Area: " << t.squared_area() << std::endl;
  
  typename Polyhedron::Halfedge_const_handle longest = get_longest_const_edge<Polyhedron>(h->facet());
  typename Polyhedron::Traits::Line_3 l(longest->vertex()->point(), longest->opposite()->vertex()->point());
  std::cout << "Colinearity: " << triangle_projection<Polyhedron>(h->facet()) << std::endl;
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
inline double
get_min_edge_length(typename Polyhedron::Facet_const_handle facet)
{
  typename Polyhedron::Facet::Halfedge_around_facet_const_circulator half_edge
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
inline double
get_max_edge_length(typename Polyhedron::Facet_const_handle facet)
{
  typename Polyhedron::Facet::Halfedge_around_facet_const_circulator h = facet->facet_begin();
  double max_length = get_edge_length<Polyhedron>(h);

  h++;
  max_length = std::max(max_length, get_edge_length<Polyhedron>(h));

  h++;
  max_length = std::max(max_length, get_edge_length<Polyhedron>(h));

  return max_length;
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
inline bool facet_is_degenerate(typename Polyhedron::Facet_const_handle facet,
                         const double tolerance)
{
  dolfin_assert(facet->is_triangle());

  return get_min_edge_length<Polyhedron>(facet) < tolerance 
    || triangle_projection<Polyhedron>(facet) < tolerance;
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
inline double shortest_edge(Polyhedron& p)
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
template<typename Polyhedron>
int number_of_degenerate_facets(const Polyhedron& p, const double tolerance)
{
  int count = 0;
  for (typename Polyhedron::Facet_const_iterator facet = p.facets_begin();
       facet != p.facets_end(); facet++)
  {
    dolfin_assert(facet->is_triangle());
    if ( facet_is_degenerate<Polyhedron>(facet, tolerance) )
      count++;
  }
  return count;
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
inline void collapse_edge(Polyhedron& p,
                          typename Polyhedron::Halfedge_handle& edge)
{
  // Join small triangles with neighbor facets
  edge = p.join_facet(edge->next());
  p.join_facet(edge->opposite()->prev());

  // The joined facets are now quads
  // Join the two close vertices
  p.join_vertex(edge);
  dolfin_assert(p.is_valid());
  dolfin_assert(p.is_pure_triangle());
}
//-----------------------------------------------------------------------------
// FIXME: Return the number of edges collapsed
template <typename Polyhedron>
void collapse_short_edges(Polyhedron& p, const double tolerance)
{
  bool removed;

  do
  {
    removed = false;

    for (typename Polyhedron::Halfedge_iterator halfedge = p.halfedges_begin();
	 halfedge != p.halfedges_end(); halfedge++)
    {
      if (get_edge_length<Polyhedron>(halfedge) < tolerance)
      {
	collapse_edge<Polyhedron>(p, halfedge);
        removed = true;
	break;
      }
    }
  } while (removed);
}
//-----------------------------------------------------------------------------
// FIXME: Return the number of facets fixed
template<typename Polyhedron>
void flip_edges(Polyhedron& p,
                double tolerance)
{
  typedef typename Polyhedron::Traits::Line_3 Line_3;
  typedef typename Polyhedron::Traits::Point_3 Point_3;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

  bool done;

  do
  {
    done = true;

    for (typename Polyhedron::Facet_iterator facet = p.facets_begin();
         facet != p.facets_end(); facet++)
    {
      dolfin_assert(facet->is_triangle());
      if (triangle_projection<Polyhedron>(facet) < tolerance)
      {     
        typename Polyhedron::Halfedge_handle longest 
          = get_longest_edge<Polyhedron>(facet);

        Halfedge_handle h = longest->next();
        
        Line_3 l(longest->vertex()->point(),
                 longest->opposite()->vertex()->point());
        Point_3 newpoint = l.projection(longest->next()->vertex()->point());
        Halfedge_handle flipped = p.flip_edge(longest);
        flipped->vertex()->point() = newpoint;
        
        // TODO: Check length of newly created edge and
        // collapse if necessary.
        collapse_short_edges(p, tolerance);
        done = false;
      }
    }
  } while (!done);
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
bool has_degenerate_facets(const Polyhedron& p,
                           double tolerance)
{
  for (typename Polyhedron::Facet_const_iterator facet = p.facets_begin();
       facet != p.facets_end(); facet++)
  {
    dolfin_assert(facet->is_triangle());
    if (facet_is_degenerate<Polyhedron>(facet, tolerance))
      return true;
  }
  return false;
}
//-----------------------------------------------------------------------------
// Remove degenerate facets of a triangular polyhedron by
// 1) Collapse edges with squared length less than tolerance
// 2) Remove (almost) colinear facets by flipping the longest edge of the 
//    colinear facet. Colinearity defined as the (squared) distance from a
//    vertex to the opposite edge being less than tolerance
template<typename Polyhedron>
std::size_t remove_degenerate(Polyhedron &p, double tolerance)
{
  dolfin_assert(p.is_pure_triangle());
  log(dolfin::TRACE, "Cleaning degenerate facets");

  log(dolfin::TRACE, "  Collapsing short edges");
  collapse_short_edges(p, tolerance);
  dolfin_assert(p.is_pure_triangle());

  // log(dolfin::TRACE, "Shortest edge: %f", shortest_edge());
  log(dolfin::TRACE, "  Removing colinear facets by edge flipping");
  flip_edges(p, tolerance);

  // Removal of facets should preserve the triangular structure
  // of the polyhedron
  dolfin_assert(p.is_pure_triangle());
  dolfin_assert(!has_degenerate_facets(p, tolerance));

  // FIXME: Return the number of facets removed.
  return 0;
}

#endif
