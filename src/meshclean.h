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

#ifndef __MESH_CLEAN_H
#define __MESH_CLEAN_H

#include <CGAL/basic.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Triangle_3.h>

namespace
{

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
typename Polyhedron::Halfedge_handle
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
template<typename Polyhedron>
inline void print_facet(typename Polyhedron::Facet_handle facet)
{
  typename Polyhedron::Halfedge_handle h1 = facet->halfedge(),
    h2 = h1->next(), h3 = h2->next();
  std::cout << "Vertices: " << std::endl;
  std::cout << "(" << h1->vertex()->point() << ") "
            << "(" << h2->vertex()->point() << ") "
            << "(" << h3->vertex()->point() << ")" << std::endl;
  std::cout << "Edge lengths (squared):"
            << (h1->vertex()->point()-h2->vertex()->point()).squared_length() << " "
            << (h2->vertex()->point()-h3->vertex()->point()).squared_length() << " "
            << (h3->vertex()->point()-h1->vertex()->point()).squared_length() << std::endl;
  typename Polyhedron::Traits::Triangle_3 t(h1->vertex()->point(),
                                            h2->vertex()->point(),
                                            h3->vertex()->point());
  std::cout << "Area: " << t.squared_area() << std::endl;
  
  typename Polyhedron::Halfedge_handle longest = get_longest_edge<Polyhedron>(facet);
  typename Polyhedron::Traits::Line_3 l(longest->vertex()->point(), longest->opposite()->vertex()->point());
  std::cout << "Colinear: " << (l.projection(longest->next()->vertex()->point())-longest->next()->vertex()->point()).squared_length() << std::endl;

  int tmp;
  std::cin >> tmp;


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
bool facet_is_degenerate(typename Polyhedron::Facet_const_handle facet,
                         const double threshold)
{
  // assume facet is a triangle
  typedef CGAL::Triangle_3<typename Polyhedron::Traits> Triangle;

  typename Polyhedron::Halfedge_const_handle h = facet->halfedge();
  Triangle t(h->vertex()->point(),
             h->next()->vertex()->point(),
             h->next()->next()->vertex()->point());

  return t.squared_area() < threshold || get_min_edge_length<Polyhedron>(facet) < threshold;
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
int number_of_degenerate_facets(const Polyhedron& p, const double threshold)
{
  int count = 0;
  for (typename Polyhedron::Facet_const_iterator facet = p.facets_begin();
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
void collapse_edge(Polyhedron& p,
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
template <typename Polyhedron>
void collapse_short_edges(Polyhedron& p, const double threshold)
{
  bool removed;

  do
  {
    removed = false;

    for (typename Polyhedron::Halfedge_iterator halfedge = p.halfedges_begin();
	 halfedge != p.halfedges_end(); halfedge++)
    {
      if (get_edge_length<Polyhedron>(halfedge) < threshold)
      {
        std::cout << "Short edge: " << halfedge->vertex()->point() << " " 
                  << halfedge->opposite()->vertex()->point() << ": " 
                  << (halfedge->vertex()->point() - halfedge->opposite()->vertex()->point()).squared_length() << std::endl;
	collapse_edge<Polyhedron>(p, halfedge);
        removed = true;
	break;
      }
    }
  } while (removed);

  std::cout << "Done removing short edges" << std::endl;
}

//-----------------------------------------------------------------------------
/* template<typename Polyhedron> */
/* void collapse_triangle(Polyhedron& P, typename Polyhedron::Facet_handle f) */
/* { */
/*   dolfin_assert(f->is_triangle()); */

/*   typedef typename Polyhedron::Halfedge_handle Halfedge_handle; */

/*   // Check if any of the vertices has degree 3, then we can just remove it as a */
/*   // center vertex */
/*   Halfedge_handle h_begin = f->halfedge(), h_current = h_begin; */
/*   do  */
/*   { */
/*     if (h_current->vertex()->degree() == 3) */
/*     { */
/*       P.erase_center_vertex(h_current); */
/*       return; */
/*     } */

/*     h_current = h_current->next(); */
/*   } while (h_current != h_begin); */


  //
// }
//-----------------------------------------------------------------------------
// Remove small triangles, ie. triangles with only short edges.
// Triangles are collapsed into a vertex.
/* template<typename Polyhedron> */
/* void collapse_small_triangles(Polyhedron& p, const double threshold) */
/* { */
/*   bool done; */
/*   do */
/*   { */
/*     done = true; */

/*     for (typename Polyhedron::Facet_iterator facet = p.facets_begin(); */
/* 	 facet != p.facets_end(); facet++) */
/*     { */
/*       dolfin_assert(facet->is_triangle()); */

/*       if (get_max_edge_length<Polyhedron>(facet) < threshold) */
/*       { */
/* 	// dolfin::cout << "Small triangle detected" << dolfin::endl; */
/* 	// print_facet<Polyhedron>(facet); */
/* 	collapse_triangle<Polyhedron>(p, facet); */
/*         done = false; */
/* 	break; */
/*       } */
/*     } */
/*   } while (!done); */
/* } */
//-----------------------------------------------------------------------------
template<typename Polyhedron>
void flip_edges(Polyhedron& p,
                double tolerance)
{
  std::cout << "Tolerance: " << tolerance << std::endl;

  bool done;

  do
  {
    done = true;

    for (typename Polyhedron::Facet_iterator facet = p.facets_begin();
         facet != p.facets_end(); facet++)
    {
      dolfin_assert(facet->is_triangle());
      if (get_triangle_area<Polyhedron>(facet) < tolerance)
      {     
        print_facet<Polyhedron>(facet);
        p.flip_edge(get_longest_edge<Polyhedron>(facet));
        dolfin_assert(p.is_valid());
        dolfin_assert(p.is_pure_triangle());
        done = false;
      }
    }
  } while (!done);

  std::cout << "Done flipping edges" << std::endl;
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
bool has_degenerate_facets(const Polyhedron& p,
                           double threshold)
{
  for (typename Polyhedron::Facet_const_iterator facet = p.facets_begin();
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
