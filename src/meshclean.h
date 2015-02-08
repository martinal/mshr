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
//#include <CGAL/boost/graph/Euler_operations.h>

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
std::size_t get_vertex_id(const Polyhedron& p, typename Polyhedron::Vertex_const_handle v)
{
  return std::distance(p.vertices_begin(), v);
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
void print_edge(const Polyhedron& p, typename Polyhedron::Halfedge_const_handle h)
{
  std::cout << "(" << get_vertex_id(p, h->opposite()->vertex()) << ": " << h->opposite()->vertex()->point() << ", "
	    << get_vertex_id(p, h->vertex()) << ": " << h->vertex()->point() << ")" << std::endl;
}

// Print some info about a facet
// Meant for debugging
template<typename Polyhedron>
void print_facet(const Polyhedron& p, typename Polyhedron::Halfedge_const_handle h)
{
  std::cout << "Vertices: " << std::endl;
  typename Polyhedron::Halfedge_const_handle current = h;
  do
  {
    std::cout << "  " << get_vertex_id<Polyhedron>(p, current->vertex()) << ": (" << current->vertex()->point() << ") " << std::endl;
    current = current ->next();
  } while (current != h);
  
  std::cout << "Edge lengths (squared):";
  current = h;
  do
  {
    std::cout << "  " << (current->vertex()->point() - current->next()->vertex()->point()).squared_length() << std::endl;
    current = current->next();
  } while (current != h);
  
  typename Polyhedron::Traits::Triangle_3 t(h->vertex()->point(),
                                            h->next()->vertex()->point(),
                                            h->next()->next()->vertex()->point());
  std::cout << "Area (if triangle): " << t.squared_area() << std::endl;
  
  typename Polyhedron::Halfedge_const_handle longest = get_longest_const_edge<Polyhedron>(h->facet());
  typename Polyhedron::Traits::Line_3 l(longest->vertex()->point(), longest->opposite()->vertex()->point());
  std::cout << "Colinearity (if triangle): " << triangle_projection<Polyhedron>(h->facet()) << std::endl;
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
inline bool has_slivers(const Polyhedron& p)
{
  for (typename Polyhedron::Vertex_const_iterator vit = p.vertices_begin(); vit != p.vertices_end(); vit++)
  {
    if (vit->vertex_degree() < 3)
      return true;
  }

  for (typename Polyhedron::Halfedge_const_iterator it = p.halfedges_begin();
       it != p.halfedges_end(); it++)
  {
    if (it->next()->vertex() == it->opposite()->next()->vertex())
      return true;
  }  
  return false;
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
inline bool has_degree3_neighbors(const Polyhedron& p)
{
  for (typename Polyhedron::Halfedge_const_iterator it = p.halfedges_begin();
       it != p.halfedges_end(); it++)
  {
    if (it->vertex()->vertex_degree() < 4 &&
        it->opposite()->vertex()->vertex_degree() < 4)
      return true;
  }

  return false;
}
//-----------------------------------------------------------------------------
#define ASSERT_GOOD_STATE(p) do             \
{                                           \
  dolfin_assert(p.is_valid());              \
  dolfin_assert(p.is_pure_triangle());      \
  dolfin_assert(!has_slivers(p));           \
  dolfin_assert(p.size_of_vertices() < 5 || !has_degree3_neighbors(p)); \
} while(false);
//-----------------------------------------------------------------------------
template <typename Polyhedron>
inline void remove_degree3_center_vertex(Polyhedron& p,
                                         typename Polyhedron::Halfedge_handle h)
{
  // Remove center vertex, but assure the degree of the vertex is 3 and that at least
  // one of the sides of the incident triangles is short
  dolfin_assert(h->vertex()->vertex_degree() == 3);

  p.erase_center_vertex(h);
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
inline bool remove_degree3_with_short_edges(Polyhedron& p, double tolerance)
{
  bool removed = false;

  for (typename Polyhedron::Vertex_iterator it = p.vertices_begin();
       it != p.vertices_end(); it++)
  {
    if (it->vertex_degree() < 4)
    {
      dolfin_assert(it->is_trivalent());

      typename Polyhedron::Halfedge_handle h = it->halfedge();
      if (get_edge_length<Polyhedron>(h) < tolerance ||
          get_edge_length<Polyhedron>(h->next()) < tolerance ||
          get_edge_length<Polyhedron>(h->prev()) < tolerance)
      {

        // std::cout << "Removing degree 3 with short edges" << std::endl;

        p.erase_center_vertex(h);
        removed = true;
      }
    }
  }
  return removed;
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
inline void collapse_edge(Polyhedron& p,
                          typename Polyhedron::Halfedge_handle edge)
{
  std::cout << "--Collapse edge (" << get_vertex_id(p, edge->vertex()) << " (" << edge->vertex()->degree() << "), "
	    << get_vertex_id(p, edge->opposite()->vertex()) << "(" << edge->opposite()->vertex()->degree() << "))" << std::endl;

  ASSERT_GOOD_STATE(p);

  std::cout << "  OK" << std::endl;

  if (edge->vertex()->is_trivalent())
  {
    std::cout << "Edge trivalent" << std::endl;
    /*   std::cout << "Point: " << edge->vertex()->point() << ", degree: " << edge->vertex()->vertex_degree() << std::endl; */
    /*   typename Polyhedron::Halfedge_around_vertex_circulator start = edge->vertex()->vertex_begin(); */
    /*   typename Polyhedron::Halfedge_around_vertex_circulator current = start; */
    /*   do */
    /*   { */
    /*     std::cout << "Neighbor: " << current->opposite()->vertex()->point()  */
    /*               << ", length: " << (current->vertex()->point()-current->opposite()->vertex()->point()).squared_length() */
    /*               << ", degree: " << current->opposite()->vertex()->vertex_degree() << std::endl; */
    /*     current++; */
    /*   } while (current != start); */

    p.erase_center_vertex(edge);
    ASSERT_GOOD_STATE(p);
  }
  else if (edge->opposite()->vertex()->is_trivalent())
  {
    std::cout << "Opposite edge trivalent" << std::endl;
    p.erase_center_vertex(edge->opposite());
    ASSERT_GOOD_STATE(p);
  }
  else
  {
    std::cout << "  OK3" << std::endl;
  
    // Join small triangles with neighbor facets

    // Make sure we don't introduce slivers
    // (ie. vertices of degree 2)
    while (edge->next()->vertex()->is_trivalent())
    {
      std::cout << "Removing center vertex" << std::endl;
      remove_degree3_center_vertex(p, edge->next());
      ASSERT_GOOD_STATE(p);
    }

    std::cout << "  OK4" << std::endl;
  
    while (edge->opposite()->next()->vertex()->is_trivalent())
    {
      std::cout << "remove opposite vertex" << std::endl;
      remove_degree3_center_vertex(p, edge->opposite()->next());
      ASSERT_GOOD_STATE(p);
    }

    std::cout << "  OK5" << std::endl;
  
    if (edge->vertex()->is_trivalent())
    {
      std::cout << "Remove center vertex" << std::endl;
      p.erase_center_vertex(edge);
    }
    else if (edge->opposite()->vertex()->is_trivalent())
    {
      std::cout << "Remove opposite center vertex" << std::endl;
      p.erase_center_vertex(edge->opposite());
    }

    ASSERT_GOOD_STATE(p);

  /* std::cout << "Collapsing edge: ( " */
  /*           << "(" << get_vertex_id(p, edge->opposite()->vertex()) << ", " << edge->opposite()->vertex_degree() << ") --> " */
  /*           << "(" << get_vertex_id(p, edge->vertex()) << ", " << edge->vertex()->vertex_degree() << ") ), length: " */
  /*           << get_edge_length<Polyhedron>(edge) << std::endl; */

  /* std::cout << "First opposite: " */
  /*           << "(" << get_vertex_id(p, edge->next()->vertex()) << ", " << edge->next()->vertex_degree() << ")" << std::endl; */
  /* std::cout << "Second opposite: " */
  /*           << "(" << get_vertex_id(p, edge->opposite()->next()->vertex()) << ", " << edge->opposite()->next()->vertex()->vertex_degree() << ")" << std::endl; */

    std::cout << "  OK7" << std::endl;

    {
      std::ofstream outfile("Testout_pre.off");
      // outfile.precision(16);
      outfile << p;
    }
    
    std::cout << "Removing one:";
    print_edge(p, edge->next()); 
    
    edge = p.join_facet(edge->next());
    dolfin_assert(p.is_valid());
    dolfin_assert(!has_slivers(p));

    std::cout << "  OK8" << std::endl;
    std::cout << "Removing two: ";
    print_edge(p, edge->edge->opposite->prev());
    
    std::cout << "Joining facets: " << std::endl;
    print_facet<Polyhedron>(p, edge->opposite()->prev());
    print_facet<Polyhedron>(p, edge->opposite()->prev()->opposite());
    p.join_facet(edge->opposite()->prev());
    std::cout << "  OK8.5" << std::endl;
    dolfin_assert(p.is_valid());
    std::cout << "  OK8.6" << std::endl;
    {
      std::ofstream outfile("Testout.off");
      // outfile.precision(16);
      outfile << p;
    }
    std::cout << "Joining: " << "(" << edge->vertex()->point() << ", " << edge->opposite()->vertex()->point() << ")" << std::endl;

    dolfin_assert(edge->vertex() != edge->opposite()->vertex());
    
    // We can possibly have a sliver now

    // The joined facets are now quads
    // Join the two close vertices
    p.join_vertex(edge);
    std::cout << "  OK9" << std::endl;
    ASSERT_GOOD_STATE(p);
  }
}
//-----------------------------------------------------------------------------
// FIXME: Return the number of edges collapsed
template <typename Polyhedron>
bool collapse_short_edges(Polyhedron& p, const double tolerance)
{
  // Degree 3 vertices with short incident edges causes problems when collapsing
  // short edges. The very ad hoc solution that has shown to work is to remove
  // these center vertices (and by that the short edges)  before collapsing short
  // edges.

  bool edges_removed = remove_degree3_with_short_edges(p, tolerance);

  bool removed;

  do
  {
    std::cout << "Collapsing short edges. " << p.size_of_halfedges() << std::endl;
    removed = false;

    for (typename Polyhedron::Halfedge_iterator halfedge = p.halfedges_begin();
	 halfedge != p.halfedges_end(); halfedge++)
    {
      if (get_edge_length<Polyhedron>(halfedge) < tolerance)
      {
	std::cout << "Collapse edge" << std::endl;
	collapse_edge<Polyhedron>(p, halfedge);
	std::cout << "Done collapsing edge" << std::endl;
	//CGAL::Euler::collapse_edge<Polyhedron>(halfedge, p);
	
	std::cout << "remove degree 3" << std::endl;
        remove_degree3_with_short_edges(p, tolerance);
        removed = true;
        edges_removed = true;
	break;
      }
    }
  } while (removed);

  return edges_removed;
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
bool flip_edges(Polyhedron& p,
                double tolerance)
{
  typedef typename Polyhedron::Traits::Line_3 Line_3;
  typedef typename Polyhedron::Traits::Point_3 Point_3;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

  bool edge_flipped = false;
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

        Line_3 l(longest->vertex()->point(),
                 longest->opposite()->vertex()->point());
        Point_3 newpoint = l.projection(longest->next()->vertex()->point());
        Halfedge_handle flipped = p.flip_edge(longest);
        flipped->vertex()->point() = newpoint;
        
        // TODO: Check length of newly created edge and
        // collapse if necessary.
        collapse_short_edges(p, tolerance);
        done = false;
        edge_flipped = true;
      }
    }
  } while (!done);

  return edge_flipped;
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
bool remove_degenerate(Polyhedron &p, double tolerance)
{
  dolfin_assert(p.is_pure_triangle());
  ASSERT_GOOD_STATE(p);
  log(dolfin::TRACE, "Cleaning degenerate facets");

  log(dolfin::TRACE, "  Collapsing short edges");
  const bool collapsed = collapse_short_edges(p, tolerance);
  ASSERT_GOOD_STATE(p);

  // log(dolfin::TRACE, "Shortest edge: %f", shortest_edge());
  log(dolfin::TRACE, "  Removing colinear facets by edge flipping");
  const bool flipped = flip_edges(p, tolerance);
  ASSERT_GOOD_STATE(p);

  dolfin_assert(!has_degenerate_facets(p, tolerance));

  return collapsed || flipped;
}

#endif
