// Copyright (C) 2014-2015 Benjamin Kehlet
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


#ifndef POLYHEDRON_UTILS_H__
#define POLYHEDRON_UTILS_H__

#include <dolfin/math/basic.h>

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Self_intersection_polyhedron_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/corefinement_operations.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Delaunay_mesher_no_edge_refinement_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <cmath>
#include <deque>
#include <fstream>

namespace mshr
{

class PolyhedronUtils
{
 public:

  //-----------------------------------------------------------------------------
  // Scans the vertices of the polyhedron the polyhedron and returns a
  // Polyhedron::Vertex_const_handle for each disconnected component.
  template <typename Polyhedron, typename OutputIterator>
  static void get_disconnected_components(const Polyhedron& p, OutputIterator it)
  {
    //typedef Polyhedron Polyhedron_t;
    typedef typename Polyhedron::Halfedge_around_vertex_const_circulator HV_const_circulator;
    typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;

    // store all vertices in a set
    std::set<Vertex_const_handle> v;
    for (typename Polyhedron::Vertex_const_iterator vit = p.vertices_begin();
         vit != p.vertices_end(); vit++)
      v.insert(vit);

    while (!v.empty())
    {
      // Add the component to the output
      typename std::set<Vertex_const_handle>::iterator start_it = v.begin();
      Vertex_const_handle start = *start_it;
      v.erase(start_it);

      *it = start;
      it++;

      // Remove rest of component
      std::deque<Vertex_const_handle> queue;
      queue.push_back(start);
      while (!queue.empty())
      {
        Vertex_const_handle current = queue.front();
        queue.pop_front();

        if (v.count(current) > 0)
        {
          v.erase(current);

          const HV_const_circulator h_start = current->vertex_begin();
          HV_const_circulator h_current = h_start;
          do
          {
            queue.push_back(h_current->opposite()->vertex());
            h_current++;
          } while (h_current != h_start);
        }
      }
    }
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static std::vector<typename Polyhedron::Halfedge_handle>
  get_holes(Polyhedron& P)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

    P.normalize_border();

    std::set<Halfedge_handle> border_edges;
    for (typename Polyhedron::Halfedge_iterator hit = P.border_halfedges_begin(); hit != P.halfedges_end(); hit++)
    {
      Halfedge_handle b = hit;
      dolfin_assert(b->is_border_edge());
      if (b->is_border())
        border_edges.insert(b);
      else if (b->opposite()->is_border())
        border_edges.insert(b->opposite());
      else
        dolfin_assert(false);
    }

    std::vector<Halfedge_handle> border_begins;
    while (!border_edges.empty())
    {
      Halfedge_handle current = *(border_edges.begin());
      border_begins.push_back(current);

      std::size_t counter = 0;
      const Halfedge_handle start = current;
      do
      {
        dolfin_assert(border_edges.count(current) == 1);
        border_edges.erase(current);
        counter++;
        current = current->next();
      } while (current != start);
    }

    return std::move(border_begins);
  }
  //-----------------------------------------------------------------------------
  // Count the number of edges in a facet or along a hole
  template<typename Halfedge_handle>
  static std::size_t edge_count(Halfedge_handle h)
  {
    Halfedge_handle current;
    std::size_t counter = 0;
    do
    {
      counter++;
      current = current->next();
    } while (current != h);

    return counter;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static typename Polyhedron::Traits::Triangle_3 get_facet_triangle(typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;

    dolfin_assert(h->facet()->is_triangle());

    return Triangle_3(h->vertex()->point(),
                      h->next()->vertex()->point(),
                      h->next()->next()->vertex()->point());
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
    static double cos_min_angle(typename Polyhedron::Traits::Vector_3 v,
                                typename Polyhedron::Halfedge_handle h1,
                                typename Polyhedron::Halfedge_handle h2)
  {
    v /= v.squared_length();

    double max = 0;
    typename Polyhedron::Halfedge_handle current = h1;
    do
    {
      typename Polyhedron::Traits::Vector_3 w(current->vertex()->point(),
                                              current->next()->vertex()->point());
      w /= std::sqrt(CGAL::to_double(w.squared_length()));
      max = std::max(std::abs(CGAL::to_double(w*v)));
    } while (current != h2);

    return max;
  }
  //-----------------------------------------------------------------------------
  // Given a facet and a vertex (assumed to be incident to the facet), find the
  // corresponding halfedge
  template<typename Polyhedron>
  static typename Polyhedron::Halfedge_handle
  find_edge(typename Polyhedron::Vertex_handle v,
            typename Polyhedron::Face_handle f)
  {
    dolfin_assert(v != typename Polyhedron::Vertex_handle());
    dolfin_assert(f != typename Polyhedron::Face_handle());

    typename Polyhedron::Halfedge_around_vertex_circulator start = v->vertex_begin();
    typename Polyhedron::Halfedge_around_vertex_circulator current = start;

    do
    {
      if (current->facet() == f)
      {
        dolfin_assert(current->vertex() == v && current->facet()  == f);
        return current;
      }

      current++;
    } while(current != start);

    dolfin_assert(false);
    return typename Polyhedron::Halfedge_handle();
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static typename Polyhedron::Facet_handle
  find_common_facet(typename Polyhedron::Vertex_handle v0,
                    typename Polyhedron::Vertex_handle v1,
                    typename Polyhedron::Vertex_handle v2)
  {
    typedef typename Polyhedron::Halfedge_around_vertex_circulator He_circulator;
    std::set<typename Polyhedron::Facet_handle> facets0;
    {
      He_circulator start = v0->vertex_begin();
      He_circulator current = start;
      do
      {
        if (!current->is_border())
          facets0.insert(current->facet());

        current++;
      } while (current != start);
    }

    std::set<typename Polyhedron::Facet_handle> facets1;
    {
      He_circulator start = v1->vertex_begin();
      He_circulator current = start;
      do
      {
        if (!current->is_border() && facets0.count(current->facet()) > 0)
          facets1.insert(current->facet());

        current++;
      } while (current != start);
    }

    std::set<typename Polyhedron::Facet_handle> facets2;
    {
      He_circulator start = v2->vertex_begin();
      He_circulator current = start;
      do
      {
        if (!current->is_border() && facets1.count(current->facet()) > 0)
          facets2.insert(current->facet());

        current++;
      } while (current != start);
    }

    dolfin_assert(facets2.size() < 2);
    dolfin_assert(facets2.size() > 0);

    return *(facets2.begin());
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void insert_edge(Polyhedron& P,
                          typename Polyhedron::Vertex_handle h,
                          typename Polyhedron::Vertex_handle g,
                          typename Polyhedron::Facet_handle f)
  {
    //std::pair<typename Polyhedron::Halfedge_handle, typename Polyhedron::Halfedge_handle>
    //edges = find_edges(h_edge, g_edge);
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

    Halfedge_handle h_edge, g_edge;
    const Halfedge_handle start = g->halfedge();
    Halfedge_handle current = start;
    do
    {
      if (current->vertex() == h)
        h_edge = current;
      else if(current->vertex() == g)
        g_edge == current;
    } while (current != start);

    dolfin_assert(h_edge != Halfedge_handle() && g_edge != Halfedge_handle());
    P.split_facet(h_edge, g_edge);
    dolfin_assert(P.is_valid());
  }
  //-----------------------------------------------------------------------------
  template<typename Triangle_3>
  static double get_triangle_cos_angle(Triangle_3 t1, 
                                       Triangle_3 t2)
  {
    typedef typename CGAL::Kernel_traits<Triangle_3>::Kernel::Vector_3 Vector_3;

    const Vector_3 v1 = t1.supporting_plane().orthogonal_vector();
    const Vector_3 v2 = t2.supporting_plane().orthogonal_vector();

    return CGAL::to_double((v1*v2)/std::sqrt(CGAL::to_double(v1.squared_length()*v2.squared_length())));  
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
    static double get_edge_cos_angle(typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;

    // std::cout << "get edge cos theta" << std::endl;
    
    Vector_3 h1_vec(h->vertex()->point(), h->prev()->vertex()->point());
    h1_vec = h1_vec/std::sqrt(CGAL::to_double(h1_vec.squared_length()));
    
    // std::cout << "h1_vec_normalized: " << h1_vec << std::endl;
    Vector_3 h2_vec(h->vertex()->point(), h->next()->vertex()->point());
    h2_vec = h2_vec/std::sqrt(CGAL::to_double(h2_vec.squared_length()));
    // std::cout << "h2_vec_normalized: " << h2_vec << std::endl;
    const double cos_theta = CGAL::to_double(h1_vec*h2_vec);
    // std::cout << "Cos theta: " << cos_theta << std::endl;
    return cos_theta;
  }
  //-----------------------------------------------------------------------------
  // Compute the plane fit quality of the vertices from h1 to h2 both included
  // Return fitting quality, max projection distance, max cos angle
  template<typename Polyhedron>
  static std::array<double, 3>
  get_plane_fit(const typename Polyhedron::Halfedge_handle h1,
                const typename Polyhedron::Halfedge_handle h2)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Plane_3 Plane_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;
    typedef typename Polyhedron::Traits::FT FT;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;
    typedef typename InexactKernel::Point_3 InexactPoint_3;
    //typedef typename InexactKernel::Segment_3 InexactSegment_3;
    //typedef typename InexactKernel::Vector_3 InexactVector_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;

    // std::cout << "Get plane fit" << std::endl;

    // std::cout << "Polygon ";
    std::vector<InexactPoint_3> points;
    //std::vector<InexactSegment_3> segments;
    Halfedge_handle current = h1;
    do
    {
      const Point_3& p = current->vertex()->point();
      // std::cout << ", " << p;
      points.push_back(InexactPoint_3(CGAL::to_double(p[0]),
                                      CGAL::to_double(p[1]),
                                      CGAL::to_double(p[2])));
      current = current->next();
    } while (current != h2);

    {
      const Point_3& p = h2->vertex()->point();
      // std::cout << ", " << p;
      points.push_back(InexactPoint_3(CGAL::to_double(p[0]),
                                      CGAL::to_double(p[1]),
                                      CGAL::to_double(p[2])));
    }

    dolfin_assert(points.size() > 2);
    // std::cout << "Size: " << points.size() << std::endl;
    //std::cout << std::endl;
    InexactPlane_3 fitting_plane_inexact;
    const double fit_quality = CGAL::linear_least_squares_fitting_3(points.begin(),
                                                                    points.end(),
                                                                    fitting_plane_inexact,
                                                                    CGAL::Dimension_tag<0>());
    Plane_3 fitting_plane(fitting_plane_inexact.a(),
                          fitting_plane_inexact.b(),
                          fitting_plane_inexact.c(),
                          fitting_plane_inexact.d());
    // std::cout << "Plane: " << fitting_plane << std::endl;
    // std::cout << "Length of normal: " << fitting_plane.orthogonal_vector().squared_length() << std::endl;
    const Vector_3 normal = fitting_plane.orthogonal_vector()/std::sqrt(CGAL::to_double(fitting_plane.orthogonal_vector().squared_length()));

    FT max_distance = (h1->vertex()->point()-fitting_plane.projection(h1->vertex()->point())).squared_length();
    FT max_angle = 0;

    Halfedge_handle prev = h1;
    current = h1->next();
    do
    {
      const Vector_3 v = current->vertex()->point()-prev->vertex()->point();
      const FT cos_angle = v/std::sqrt(CGAL::to_double(v.squared_length())) * normal;
      const Point_3 projection = fitting_plane.projection(current->vertex()->point());
      max_angle = std::max(max_angle, cos_angle);
      max_distance = std::max(max_distance, (current->vertex()->point()-projection).squared_length());

      prev = current;
      current = current->next();
    } while (prev != h2);

    // std::cout << "Fit quality: " << fit_quality << ", max distance: " << max_distance << ", cos_angle: " << max_angle << std::endl;
    return std::array<double, 3>{fit_quality, CGAL::to_double(max_distance), CGAL::to_double(max_angle)};
    // return -max_distance;
    //return CGAL::to_double(fit_quality - max_angle);
  }
  //-----------------------------------------------------------------------------
  // Compute the fit quality heuristic of the vertices from h1 to h2 both
  // included Some experiences:

  // The plane fit quality as returned from
  // * CGAL::linear_least_squares_fitting_3() is not suitable here. It measures
  // how distinct the best fitting plane is, rather than the actual quality of
  // the fit.
  // * The max angle between the segments and the normal work good is the hole
  // has no more than one "kink".
  // * The max distance from the points to their projection is also a pretty
  // good measure, but needs to be normalized in some clever way if not used
  // solely.
  template<typename Polyhedron>
  static double evaluate_hole_subdivision(const typename Polyhedron::Halfedge_handle h1,
                                          const typename Polyhedron::Halfedge_handle h2)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef CGAL::Simple_cartesian<double> InexactKernel;
    //typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;
    typedef typename InexactKernel::Point_3 InexactPoint_3;
    typedef typename InexactKernel::Segment_3 InexactSegment_3;
    typedef typename InexactKernel::Vector_3 InexactVector_3;

    // double plane1fit;
    double max_angle1 = 0;
    InexactPlane_3 fitting_plane1;
    double avg_distance_squared1 = 0;
    {
      std::vector<InexactSegment_3> segments;
      Halfedge_handle current = h1;
      do
      {
        const Point_3& p = current->vertex()->point();
        const Point_3& next = current->next()->vertex()->point();
        segments.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])),
                                            InexactPoint_3(CGAL::to_double(next[0]), CGAL::to_double(next[1]), CGAL::to_double(next[2]))));
        current = current->next();
      } while (current != h2);

      {
        const Point_3& p = h2->vertex()->point();
        const Point_3& next = h1->vertex()->point();
        segments.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])),
                                            InexactPoint_3(CGAL::to_double(next[0]), CGAL::to_double(next[1]), CGAL::to_double(next[2]))));
      }

      dolfin_assert(segments.size() > 2);

      /* plane1fit = */ CGAL::linear_least_squares_fitting_3(segments.begin(),
                                                             segments.end(),
                                                             fitting_plane1,
                                                             CGAL::Dimension_tag<1>());

      // std::cout << "  Plane 1: " << fitting_plane1 << ", " << plane1fit << std::endl;

      // Compute max angle between segment and fitting plane
      const InexactVector_3 normal = fitting_plane1.orthogonal_vector();
      // std::cout << "Length: " << normal.squared_length() << std::endl;
      dolfin_assert(dolfin::near(normal.squared_length(), 1, DOLFIN_EPS_LARGE));
      for (auto sit = segments.begin(); sit != segments.end(); sit++)
      {
        const InexactVector_3 v = InexactVector_3(*sit) / std::sqrt(sit->squared_length());
        // std::cout << "Length: " << (v.squared_length()-1) << std::endl;
        dolfin_assert(dolfin::near(v.squared_length(), 1, DOLFIN_EPS_LARGE));
        max_angle1 = std::max(max_angle1, std::abs(v*normal));
        // max_distance_squared1 = std::max(max_distance_squared1,
        avg_distance_squared1 += InexactVector_3(sit->source(), fitting_plane1.projection(sit->source())).squared_length();
      }
      avg_distance_squared1 /= segments.size();
    }

    /* double plane2fit; */
    InexactPlane_3 fitting_plane2;
    double max_angle2 = 0;
    //double max_distance_squared2 = 0;
    double avg_distance_squared2 = 0;
    {
      std::vector<InexactSegment_3> segments;
      Halfedge_handle current = h2;
      do
      {
        const Point_3& p = current->vertex()->point();
        const Point_3& next = current->next()->vertex()->point();
        segments.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])),
                                            InexactPoint_3(CGAL::to_double(next[0]), CGAL::to_double(next[1]), CGAL::to_double(next[2]))));
        current = current->next();
      } while (current != h1);

      const Point_3& p = h1->vertex()->point();
      const Point_3& next = h2->vertex()->point();

      segments.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])),
                                          InexactPoint_3(CGAL::to_double(next[0]), CGAL::to_double(next[1]), CGAL::to_double(next[2]))));

      dolfin_assert(segments.size() > 2);
      // std::cout << "  Size: " << segments.size() << std::endl;

      /* plane2fit = */ CGAL::linear_least_squares_fitting_3(segments.begin(),
                                                             segments.end(),
                                                             fitting_plane2,
                                                             CGAL::Dimension_tag<1>());

      // std::cout << "  Plane 1: " << fitting_plane1 << ", " << plane1fit << std::endl;

      // Compute max angle between plane and segments
      const InexactVector_3 normal = fitting_plane2.orthogonal_vector();
      dolfin_assert(dolfin::near(normal.squared_length(), 1, DOLFIN_EPS_LARGE));
      for (auto sit = segments.begin(); sit != segments.end(); sit++)
      {
        const InexactVector_3 v = InexactVector_3(*sit) / std::sqrt(sit->squared_length());
        dolfin_assert(dolfin::near(v.squared_length(), 1, DOLFIN_EPS_LARGE));
        max_angle2 = std::max(max_angle2, std::abs(v*normal));
        /* max_distance_squared2 = std::max(max_distance_squared2, */
        /*                                  InexactVector_3(sit->source(), fitting_plane2.projection(sit->source())).squared_length()); */
        avg_distance_squared2 += InexactVector_3(sit->source(), fitting_plane2.projection(sit->source())).squared_length();
      }
      avg_distance_squared2 /= segments.size();
    }

    // const double cos_angle = fitting_plane1.orthogonal_vector()*fitting_plane2.orthogonal_vector();
    // std::cout << "  Angle: " << cos_angle << "(" << acos(cos_angle)/(2*DOLFIN_PI)*360 << ")" << std::endl;
    // std::cout << "Max angles: " << max_angle1 << "(" << acos(max_angle1)/(2*DOLFIN_PI)*360 << "), " << max_angle2 << " (" << acos(max_angle2)/(2*DOLFIN_PI)*360 << ")" << std::endl;
    // return std::min(plane1fit, plane2fit); // - .04*cos_angle - .05*max_angle1 - .05*max_angle2;// + triangulation_extra;

    //return (-avg_distance_squared1 - avg_distance_squared2)/std::max(max_angle1, max_angle2); ///std::min(plane1fit, plane2fit);
    return -std::max(max_angle1, max_angle2);
  }
  //-----------------------------------------------------------------------------
  // This class is capable of filling a hole in a polyhedron by inserting
  // a 2D triangulation.
  // TODO: This needs to be fixed!
  template <class Polyhedron, class CDT>
  class Triangulation_inserter
    : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
  {
   public:
    typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
    typedef typename Polyhedron::Vertex Vertex;
    typedef typename Polyhedron::Halfedge_around_vertex_circulator He_circ;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename CDT::Geom_traits::Point_3 InexactPoint_3;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits Traits;
    typedef typename CGAL::Aff_transformation_3<typename CDT::Geom_traits> Aff_transformation_3;

    // Constructor
 Triangulation_inserter(CDT& cdt, Aff_transformation_3 tr, double z, Halfedge_handle h)
   : cdt(cdt), tr(tr), z(z), h(h)
    {}

    void operator()(HalfedgeDS& hds)
    {
      CGAL::HalfedgeDS_decorator<HalfedgeDS> decorator(hds);
      typedef typename HalfedgeDS::Halfedge Halfedge;

      std::map<std::pair<Vertex_handle, Vertex_handle>, Halfedge_handle> inserted_edges;

      // Collect the already existing edges on the boundary
      Halfedge_handle current = h;
      do
      {
        inserted_edges[std::make_pair(current->vertex(), current->next()->vertex())] = current->next();
        current = current->next();
      } while (current != h);

      std::size_t vcounter = 0;
      std::vector<Vertex_handle> inserted_vertices;

      for (typename CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
           vit != cdt.finite_vertices_end(); vit++)
      {
        if (vit->info() == Vertex_handle())
        {
          // std::cout << "Inserting new vertex" << std::endl;
          Vertex_handle v = decorator.vertices_push_back(Vertex());
          inserted_vertices.push_back(v);
          vit->info() = v;
          const InexactPoint_3 p(vit->point()[0], vit->point()[1], z);
          const InexactPoint_3 rotated = tr.transform(p);
          v->point() = Point_3(rotated[0], rotated[1], rotated[2]);
          vcounter++;
        }
      }

      // std::cout << vcounter << " vertices inserted" << std::endl;

      std::size_t fcounter = 0;
      for (typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
           fit != cdt.finite_faces_end(); fit++)
      {
        if (fit->is_in_domain())
        {
          // std::cout << "Inserting face" << std::endl;

          // Check if any of the edges are missing
          std::array<Vertex_handle, 3> vertices{fit->vertex(0)->info(),
                                                fit->vertex(1)->info(),
                                                fit->vertex(2)->info()} ;
          dolfin_assert(vertices[0] != Vertex_handle());
          dolfin_assert(vertices[1] != Vertex_handle());
          dolfin_assert(vertices[2] != Vertex_handle());

          // std::cout << "Creating edges" << std::endl;
          std::array<Halfedge_handle, 3> edges;
          for (std::size_t i = 0; i < 3; i++)
          {
            if (inserted_edges.count(std::make_pair(vertices[i], vertices[(i+1)%3])) == 0)
            {
              Halfedge_handle hnew = hds.edges_push_back(Halfedge(),
                                                         Halfedge());
              decorator.set_vertex(hnew, vertices[(i+1)%3]);
              decorator.set_vertex(hnew->opposite(), vertices[i]);
              inserted_edges[std::make_pair(vertices[i], vertices[(i+1)%3])] = hnew;
              inserted_edges[std::make_pair(vertices[(i+1)%3], vertices[i])] = hnew->opposite();
              edges[i] = hnew;
            }
            else
              edges[i] = inserted_edges[std::make_pair(vertices[i], vertices[(i+1)%3])];
          }

          // std::cout << "Now creating face" << std::endl;
          typedef typename Halfedge::Base HBase;
          edges[0]->HBase::set_next(edges[1]);
          decorator.set_prev(edges[1], edges[0]);
          edges[1]->HBase::set_next(edges[2]);
          decorator.set_prev(edges[2], edges[1]);
          edges[2]->HBase::set_next(edges[0]);
          decorator.set_prev(edges[0], edges[2]);

          decorator.fill_hole(edges[0]);

          fcounter++;
        }
      }

      // std::cout << "Added " << fcounter << " faces" << std::endl;
    }

    CDT& cdt;
    const CGAL::Aff_transformation_3<typename CDT::Geom_traits> tr;
    const double z;
    const Halfedge_handle h;
  };

  //-----------------------------------------------------------------------------
  // Compute the transformation that rotates a given vector (assumed to be of
  // unit length) into (0,0,1)
  template<typename Vector_3>
  static CGAL::Aff_transformation_3<typename CGAL::Kernel_traits<Vector_3>::Kernel>
    rotate_to_xy(Vector_3 a)
  {
    const typename CGAL::Kernel_traits<Vector_3>::Kernel::FT den = a[0]*a[0] + a[1]*a[1];
    return CGAL::Aff_transformation_3<typename CGAL::Kernel_traits<Vector_3>::Kernel>
      (1 - a[0]*a[0]*(1-a[2])/den, -a[0]*a[1]*(1-a[2])/den,   -a[0],
       - a[0]*a[1]*(1-a[2])/den,   1 -a[1]*a[1]*(1-a[2])/den, -a[1],
       a[0],                       a[1],                      1 + (-a[0]*a[0]-a[1]*a[1])*(1-a[2])/den);
  }
  //-----------------------------------------------------------------------------
  /// Attempts to triangulate a polygon in 3d by projecting vertices into the
  /// best fitting plane and triangulating in 2d.
  /// Return the new added edges.
  /// If the triangulation is not possible (the boundary self intersects in this 2d plane) then
  /// the return vector is empty
  template <typename Polyhedron>
  static bool triangulate_polygon_3d(Polyhedron& P,
                                     typename Polyhedron::Halfedge_handle h,
                                     bool check_for_intersections = true)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Segment_3 Segment_3;
    
    typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;
    typedef typename InexactKernel::Point_3 InexactPoint_3;
    typedef typename InexactKernel::Vector_3 InexactVector_3;
    typedef typename InexactKernel::Point_2 InexactPoint_2;
    typedef typename InexactKernel::Segment_2 InexactSegment_2;
    typedef typename InexactKernel::Segment_3 InexactSegment_3;

    typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_handle, InexactKernel> Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<InexactKernel> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
    typedef CGAL::No_intersection_tag Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<InexactKernel, TDS, Itag> CDT;
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Mesh_criteria_2;
    typedef CGAL::Delaunay_mesher_no_edge_refinement_2<CDT, Mesh_criteria_2>  CGAL_Mesher_2;

    dolfin_assert(P.is_valid());
    // std::cout << "Triangulating hole as 2d polygon" << std::endl;

    // Compute the best fitting plane of the points of the hole
    InexactPlane_3 fitting_plane;
    // double fit_quality;
    {
      std::vector<InexactSegment_3> boundary;
      Halfedge_handle current = h;
      do
      {
        const Point_3& p = current->vertex()->point();
        const Point_3& next = current->next()->vertex()->point();
        boundary.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])),
                                            InexactPoint_3(CGAL::to_double(next[0]), CGAL::to_double(next[1]), CGAL::to_double(next[2]))));

        current = current->next();
      } while (current != h);

      const double fit_quality = CGAL::linear_least_squares_fitting_3(boundary.begin(),
                                                                      boundary.end(),
                                                                      fitting_plane,
                                                                      CGAL::Dimension_tag<1>());
      double min_cos_normal_angle = 0;
      const InexactVector_3 normal = fitting_plane.orthogonal_vector();
      double max_squared_distance = 0;
      dolfin_assert(dolfin::near(normal.squared_length(), 1, DOLFIN_EPS_LARGE));

      InexactVector_3 prev = InexactVector_3(boundary[boundary.size()-1]) / std::sqrt(boundary[boundary.size()-1].squared_length());
      for (auto sit = boundary.begin(); sit != boundary.end(); sit++)
      {
        InexactVector_3 current = InexactVector_3(*sit) / std::sqrt(sit->squared_length());
        max_squared_distance = std::max(max_squared_distance, (sit->source()-fitting_plane.projection(sit->source())).squared_length());

        // std::cout << "Length: " << sit->squared_length() << ", " << current.squared_length() << ", " << (current*normal) << std::endl;
        dolfin_assert(dolfin::near(current.squared_length(), 1, DOLFIN_EPS_LARGE));
        min_cos_normal_angle = std::max(min_cos_normal_angle, std::abs(current*normal));
        prev = current;
      }

      // std::cout << "Max abs cos normal angle: " << min_cos_normal_angle << std::endl;
      // std::cout << "Plane quality: " << fit_quality << std::endl;
      // std::cout << "Max distance: " << std::sqrt(max_squared_distance) << std::endl;

      // FIXME: Improve this test
      if (fit_quality < .95 || min_cos_normal_angle > .3)
      {
        // std::cout << "Rejecting 2d triangulating" << std::endl;
        return false;
      }
    }

    // Compute rotation that will rotate the fitting plane to the xy plane
    const CGAL::Aff_transformation_3<InexactKernel> rotation = rotate_to_xy(fitting_plane.orthogonal_vector());
    double z = 0;//  = rotation.transform(InexactPoint_3(CGAL::to_double(h->vertex()->point()[0]),
    //                              CGAL::to_double(h->vertex()->point()[1]),
    //                                                 CGAL::to_double(h->vertex()->point()[2])))[2];

    // std::cout << "Rotate normal: " << rotation.transform(fitting_plane.orthogonal_vector()) << std::endl;
    dolfin_assert(dolfin::near(fitting_plane.orthogonal_vector().squared_length(), 1, DOLFIN_EPS_LARGE));

    CDT cdt;

    // std::cout << "Projected polygon" << std::endl;
    // std::cout << "Polygon";

    // Insert vertices into 2D triangulation
    std::vector<typename CDT::Vertex_handle> vertices;
    double max_squared_edge_length = 0;

    {
      Halfedge_handle current = h;
      Point_3 prev = current->prev()->vertex()->point();
      do
      {

        const Point_3& p = current->vertex()->point();
        max_squared_edge_length = std::max(max_squared_edge_length, CGAL::to_double(Segment_3(prev, p).squared_length()));
        const InexactPoint_3 rotated = rotation.transform(InexactPoint_3(CGAL::to_double(p[0]),
                                                                         CGAL::to_double(p[1]),
                                                                         CGAL::to_double(p[2])));
        z += rotated[2];
        // std::cout << " " << rotated << ", ";
      
        const InexactPoint_2 p_2d(rotated[0], rotated[1]);

        // std::cout << " " << p_2d << ", ";
      
        vertices.push_back(cdt.insert(p_2d));
        vertices.back()->info() = current->vertex();

        prev = p;
        current = current->next();
      } while (current != h);

      // std::cout << std::endl;
    }

    // std::cout << "Size of points: " << vertices.size() << std::endl;
    z /= vertices.size();

    // Check if any of the edges intersect (before actually adding the
    // constrained edges to the triangulation
    if (check_for_intersections)
    {
      for (std::size_t i = 0; i < vertices.size()-1; i++)
      {
        const Point_3& a = vertices[i]->info()->point(), b = vertices[+1]->info()->point();
        const InexactSegment_2 s(vertices[i]->point(), vertices[i+1]->point());
        const Segment_3 original(a, b);

        const InexactSegment_3 s2(fitting_plane.projection(InexactPoint_3(CGAL::to_double(a[0]),
                                                                          CGAL::to_double(a[1]),
                                                                          CGAL::to_double(a[2]))),
                                  fitting_plane.projection(InexactPoint_3(CGAL::to_double(b[0]),
                                                                          CGAL::to_double(b[1]),
                                                                          CGAL::to_double(b[2]))));

        for (std::size_t j = i+1; j < vertices.size(); j++)
        {
          InexactSegment_2 s2(vertices[j]->point(), vertices[(j+1)%vertices.size()]->point());
          
          const auto intersection = CGAL::intersection(s, s2);
          
          if (intersection)
          {
            if (const InexactPoint_2* intersection_point = boost::get<InexactPoint_2>(&*intersection))
            {
              if (j != i+1 && i != (j+1)%vertices.size())
              {
                // std::cout << "Non-neighbors (" << i << ", " << j << ")/" << vertices.size() << " intersect in single point" << std::endl;
                
                return false;
              }
            }
            else if (const InexactSegment_2* intersection_segment = boost::get<InexactSegment_2>(&*intersection))
            {
              // std::cout << "Intersects in segment" << std::endl;
              return false;
            }
            else
            {
              dolfin_assert(false);
              return false;
            }
          } // end if intersection
        } // end inner loop
      } // end outer loop
      
      // No edges intersect, so we can safely insert then as constraints to the
      // triangulation
    }

    // Insert the edges around the facet as constraints to the triangulation
    for (std::size_t i = 0; i < vertices.size(); i++)
    {
      cdt.insert_constraint(vertices[i], vertices[(i+1)%vertices.size()]);
    }

    // std::cout << "Done triangulating" << std::endl;
    // std::cout << "Num vertices: " << cdt.number_of_vertices() << std::endl;

    // Create mesher
    CGAL_Mesher_2 mesher(cdt);

    // Set shape and size criteria
    mesher.set_criteria(Mesh_criteria_2(.125, std::sqrt(max_squared_edge_length)));

    // Refine CGAL mesh/triangulation
    mesher.refine_mesh();

    // std::cout << "Done meshing. Num vertices: " << cdt.number_of_vertices() << std::endl;

    // Collecting faces inside the polygon
    std::set<typename CDT::Face_handle> faces_inside;
    std::size_t num_cells = 0;
    for (typename CDT::Finite_faces_iterator cgal_cell = cdt.finite_faces_begin();
         cgal_cell != cdt.finite_faces_end(); ++cgal_cell)
    {
      // Add cell if it is in the domain
      if (cgal_cell->is_in_domain())
      {
        faces_inside.insert(cgal_cell);
        num_cells++;
      }
    }

    // std::cout << "Collected " << num_cells << " faces inside" << std::endl;

    // Check if any of the triangles will intersect when transformed back to
    // the polyhedron.
    /* std::cout << "Checking if triangulation can be inserted without introducing self intersections" << std::endl; */
      
    /* std::vector<Triangle_3> triangle_set; */

    /* // Collect neighbor triangles */
    /* Halfedge_handle current = h; */
    /* do */
    /* { */
    /*   if (!current->opposite()->is_border() && current->opposite()->facet()->facet_degree() == 3) */
    /*   { */
    /*     triangle_set.push_back(get_facet_triangle<Polyhedron>(current->opposite())); */
    /*   } */
    /*   current = current->next(); */
    /* } while (current != h); */

    /* const auto rotate_back = rotation.inverse(); */

    /* for (auto fit = faces_inside.begin(); fit != faces_inside.end(); ++fit) */
    /* { */
    /*   Point_3 t0; */
    /*   if ((*fit)->vertex(0)->info() != Vertex_handle()) */
    /*     t0 = (*fit)->vertex(0)->info()->point(); */
    /*   else */
    /*   { */
    /*     InexactPoint_3 t_ = rotate_back.transform(InexactPoint_3( (*fit)->vertex(0)->point()[0], */
    /*                                                               (*fit)->vertex(0)->point()[1], */
    /*                                                               z)); */
    /*     t0 = Point_3(t_[0], t_[1], t_[2]); */
    /*   } */

    /*   Point_3 t1; */
    /*   if ((*fit)->vertex(1)->info() != Vertex_handle()) */
    /*     (*fit)->vertex(1)->info()->point(); */
    /*   else */
    /*   { */
    /*     InexactPoint_3 t_ = rotate_back.transform(InexactPoint_3( (*fit)->vertex(1)->point()[0], */
    /*                                                               (*fit)->vertex(1)->point()[1], */
    /*                                                               z)); */
    /*     t1 = Point_3(t_[0], t_[1], t_[2]); */
    /*   } */

    /*   Point_3 t2; */
    /*   if ((*fit)->vertex(2)->info() != Vertex_handle()) */
    /*     t2 = (*fit)->vertex(2)->info()->point(); */
    /*   else */
    /*   { */
    /*     InexactPoint_3 t_ = rotate_back.transform(InexactPoint_3( (*fit)->vertex(2)->point()[0], */
    /*                                                               (*fit)->vertex(2)->point()[1], */
    /*                                                               z)); */
    /*     t2 = Point_3(t_[0], t_[1], t_[2]); */
    /*   } */
      
    /*   Triangle_3 t(t0, t1, t2); */
    /*   if (triangle_set_intersect_triangle(t, triangle_set)) */
    /*   { */
    /*     std::cout << "2D triangulation cancelled! Will introduce self intersections" << std::endl; */
    /*     { int tmp; std::cin >> tmp; } */
    /*     return false; */
    /*   } */

    /*   triangle_set.push_back(t); */
    /* } */

    /* std::cout << "  Done" << std::endl; */
      
    // Triangles did not intersect any of the surrounding triangles

    dolfin_assert(cdt.is_valid());

    // Triangulation_inserter<Polyhedron, CDT> modifier(cdt, rotation.inverse(), z, h);
    // P.delegate(modifier);

    std::set<typename CDT::Face_handle> faces;
    for (typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
         fit != cdt.finite_faces_end(); fit++)
    {
      if (fit->is_in_domain())
        faces.insert(fit);
    }
    // std::cout << "Faces to be insered: " << faces.size() << std::endl;
    // std::size_t countdown = faces.size()+40;

    // TODO: This works, but is incredibly slow for large holes. Should be
    // possible to insert the new vertices and edges directly to the HDS of the
    // polyhedron using a modifier.
    
    while (faces.size() > 0)
    {
      // std::cout << "---- Face iteration (outer): " << h->facet()->facet_degree() << std::endl;
      auto fit = faces.begin();
      while (fit != faces.end())
      {
        dolfin_assert(P.is_valid());

        typename std::set<typename CDT::Face_handle>::iterator f = fit;
        fit++;

        dolfin_assert(!h->is_border());

        std::array<typename CDT::Vertex_handle, 3> vertices{ (*f)->vertex(0),
                                                             (*f)->vertex(1),
                                                             (*f)->vertex(2) };

        if (vertices[0]->info() != Vertex_handle() &&
            vertices[1]->info() != Vertex_handle() &&
            vertices[2]->info() != Vertex_handle())
        {
          // ear cut, all vertices already exists
          // std::cout << "All vertices present" << std::endl;
          const Facet_handle common_facet = find_common_facet<Polyhedron>(vertices[0]->info(),
                                                                          vertices[1]->info(),
                                                                          vertices[2]->info());
          if (common_facet->facet_degree() > 3)
          {
            dolfin_assert(common_facet == h->facet());
            std::array<Halfedge_handle, 3> halfedges { find_edge<Polyhedron>(vertices[0]->info(), h->facet()),
                                                       find_edge<Polyhedron>(vertices[1]->info(), h->facet()),
                                                       find_edge<Polyhedron>(vertices[2]->info(), h->facet()) };

            /* std::cout << "Vertices: " << std::distance(P.halfedges_begin(), halfedges[0]) << ", " */
            /*           << std::distance(P.halfedges_begin(), halfedges[1]) << ", " */
            /*           << std::distance(P.halfedges_begin(), halfedges[2]) << std::endl; */

            /* std::cout << "0: " << std::distance(P.halfedges_begin(), halfedges[0]->prev()) << " - " */
            /*           << std::distance(P.halfedges_begin(), halfedges[0]) << " - " */
            /*           << std::distance(P.halfedges_begin(), halfedges[0]->next()) << std::endl; */
            /* std::cout << "1: " << std::distance(P.halfedges_begin(), halfedges[1]->prev()) << " - " */
            /*           << std::distance(P.halfedges_begin(), halfedges[1]) << " - " */
            /*           << std::distance(P.halfedges_begin(), halfedges[1]->next()) << std::endl; */
            /* std::cout << "2: " << std::distance(P.halfedges_begin(), halfedges[2]->prev()) << " - " */
            /*           << std::distance(P.halfedges_begin(), halfedges[2]) << " - " */
            /*           << std::distance(P.halfedges_begin(), halfedges[2]->next()) << std::endl; */

            // Check how many edges are already present
            // insert edge (2 <--> 0), ie cut 1
            if ((halfedges[2]->next() == halfedges[1] &&
                 halfedges[2]->next()->next() == halfedges[0]) ||
                (halfedges[0]->next() == halfedges[1] &&
                 halfedges[0]->next()->next() == halfedges[2]))
            {
              // std::cout << "inserting (2-0)" << std::endl;
              h = halfedges[1]->next()->next();
              P.split_facet(halfedges[2], halfedges[0]);
            }
            // insert edge (2 <--> 1), ie cut 0
            else if ((halfedges[1]->next() == halfedges[0] &&
                        halfedges[1]->next()->next() == halfedges[2]) ||
                       (halfedges[2]->next() == halfedges[0] &&
                        halfedges[2]->next()->next() == halfedges[1]))
            {
              // std::cout << "inserting (2-1)" << std::endl;
              h = halfedges[0]->next()->next();
              P.split_facet(halfedges[2], halfedges[1]);
            }
            // insert edge (0 <--> 1), ie cut 2
            else if ((halfedges[0]->next() == halfedges[2] &&
                      halfedges[0]->next()->next() == halfedges[1]) ||
                     (halfedges[1]->next() == halfedges[2] &&
                      halfedges[1]->next()->next() == halfedges[0]))
            {
              // std::cout << "inserting (0-1)" << std::endl;
              h = halfedges[2]->next()->next();
              P.split_facet(halfedges[0], halfedges[1]);
            }
            else
            {
              // Note that even if all the vertices have been inserted, we may
              // miss edges
              continue;
            }
          }

          // print sometimes
          /* if (countdown-faces.size() > 100) */
          /* { */
          /*   std::cout << "Face iteration (inner): " << faces.size() << std::endl; */
          /*   countdown = faces.size(); */
          /* } */

          faces.erase(f);
          continue;
        }
        else
        {
          for (std::size_t i = 0; i < 3; i++)
          {
            if (vertices[i]->info() != Vertex_handle() &&
                vertices[(i+1)%3]->info() == Vertex_handle() &&
                vertices[(i+2)%3]->info() != Vertex_handle())
            {
              // std::cout << "vertex insertion (" << i << ")" << std::endl;
              Halfedge_handle h1 = find_edge<Polyhedron>(vertices[i]->info(), h->facet());
              Halfedge_handle h2 = find_edge<Polyhedron>(vertices[(i+2)%3]->info(), h->facet());

              if (h1->next() != h2)
              {
                Halfedge_handle tmp = h1;
                h1 = h2;
                h2 = tmp;
              }
              if (h1->next() != h2)
                break;

              dolfin_assert(h1->next() == h2);

              h = h2->next();

              // We are not allowed to introduce multiedges, so we can't do
              // split_facet(h1, h2)
              // Split facet w
              // std::cout << "adding diagonal" << std::endl;
              Halfedge_handle new_diagonal = P.split_facet(h2, h2->next()->next());
              dolfin_assert(new_diagonal->prev() == h2);
              // std::cout << "Adding new vertex" << std::endl;
              Halfedge_handle new_vertex = P.split_edge(new_diagonal);
              dolfin_assert(new_vertex->prev() == h2);
              // std::cout << "Rotating back" << std::endl;
              const auto rotate_back = rotation.inverse();
              InexactPoint_3 new_point = rotate_back.transform(InexactPoint_3(vertices[(i+1)%3]->point()[0],
                                                                              vertices[(i+1)%3]->point()[1],
                                                                              z));
              new_vertex->vertex()->point() = Point_3(new_point[0],
                                                      new_point[1],
                                                      new_point[2]);

              Vertex_handle v = new_vertex->vertex();
              typename CDT::Vertex_handle v_vdt = vertices[(i+1)%3];
              // std::cout << "Setting new vertex as info" << std::endl;

              vertices[(i+1)%3]->info() = new_vertex->vertex();
              // std::cout << "Splitting facet" << std::endl;
              P.split_facet(h1, new_vertex);

              dolfin_assert(h2->facet()->facet_degree() == 3);
              dolfin_assert(h2->prev()->vertex() == h1->vertex());
              // std::cout << "Joining facet" << std::endl;
              P.join_facet(new_diagonal);

              //dolfin_assert(edge_count(h2) == 3);
              // std::cout << "New facet size: " << h->facet()->facet_degree() << std::endl;
              // std::cout << "Erasing facet" << std::endl;
              faces.erase(f);
              // std::cout << "Breaking" << std::endl;
              break;
            }
          }
        }
      }
    }
    
    //P.normalize_border();
    // dolfin_assert(P.is_valid(false, 1));

    return true;
  }
  //-----------------------------------------------------------------------------
  /* template <typename Polyhedron> */
  /* void min_vertex_degree(const Polyhedron& p) */
  /* { */
  /*   std::size_t min_degree = std::numeric_limits<std::size_t>::max(); */
  /*   std::size_t min_degree_non_border = min_degree; */

  /*   for (typename Polyhedron::Vertex_const_iterator vit = p.vertices_begin(); */
  /*        vit != p.vertices_end(); vit++) */
  /*   { */
  /*     min_degree = std::min(min_degree, vit->vertex_degree); */
      
  /*   } */

  /*   std::cout << "Min vertex_degree: " << min_degree << std::endl; */
  /* } */
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void list_hole(typename Polyhedron::Halfedge_handle h)
  {
    std::size_t counter = 0;
    // std::cout << "Polygon";

    {
      typename Polyhedron::Halfedge_handle current = h;
      do
      {
        counter++;
        current = current->next();
      } while(current != h);
    }

    // if (counter < 250)
    // {
      typename Polyhedron::Halfedge_handle current = h;
      do
      {
        // std::cout << " " << current->vertex()->point() << ",";

        current = current->next();
      } while(current != h);
      // }
      // std::cout << std::endl;

      // std::cout << " size: " << counter << std::endl;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void print_triangle(typename Polyhedron::Halfedge_handle h)
  {
    std::cout << "Triangle "
              << h->prev()->vertex()->point() << ", " 
              << h->vertex()->point() << ", " 
              << h->next()->vertex()->point() << std::endl;
    std::cout << "Area: " << triangle_area<Polyhedron>(h) << std::endl;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static double triangle_area(typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;

    Triangle_3 t(h->prev()->vertex()->point(),
                 h->vertex()->point(),
                 h->next()->vertex()->point());

    return CGAL::to_double(t.squared_area());
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static typename Polyhedron::Vertex_handle get_common_vertex(typename Polyhedron::Facet_handle f1,
                                                              typename Polyhedron::Facet_handle f2)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;

    // Find common vertex
    Halfedge_handle h1 = f1->halfedge();
    Halfedge_handle current1 = h1;
    do
    {
      Halfedge_handle h2 = f2->halfedge();
      Halfedge_handle current2 = h2;
      do
      {
        if (current2->vertex() == current1->vertex())
          return current2->vertex();

        current2 = current2->next();
      } while (h2 != current2);

      current1 = current1->next();
    } while (h1 != current1);

    return Vertex_handle();
  }
  //----------------------------------------------------------------------------- 
  template <typename Polyhedron>
  static double facet_angle(typename Polyhedron::Halfedge_handle h)
  {
    dolfin_assert(h->is_border());
    dolfin_assert(h->next()->is_border());
    
    typedef typename Polyhedron::Traits::Line_3 Line_3;
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;

    const Line_3 l(h->prev()->vertex()->point(), 
                   h->vertex()->point());

    const Point_3& p11 = h->next()->vertex()->point();
    const Point_3  p12 = l.projection(p11);
    const Vector_3 v1(p12, p11);

    dolfin_assert(h->opposite()->facet()->is_triangle());

    const Point_3& p21 = h->opposite()->next()->vertex()->point();
    const Point_3  p22 = l.projection(p21);
    const Vector_3 v2(p22, p21);

    return CGAL::to_double((v1*v2)/std::sqrt(CGAL::to_double(v1.squared_length()*v2.squared_length())));
  }
  //-----------------------------------------------------------------------------
  template<typename Segment_3>
  static bool segment_intersects_triangle(const Segment_3& s,
                                          const typename CGAL::Kernel_traits<Segment_3>::Kernel::Triangle_3& t)
  {
    typedef typename CGAL::Kernel_traits<Segment_3>::Kernel::Point_3 Point_3;
    
    auto result = CGAL::intersection(s, t);
    if (!result)
      return false;

    if (const Point_3* p = boost::get<Point_3>(&*result))
    {
      if (*p == t[0] || *p == t[1] || *p == t[2])
        return false;
      else
        return true;
    }
    else if (const Segment_3* s_ = boost::get<Segment_3>(&*result))
    {
      if ( (s.source() == t[0] || s.source() == t[1] || s.source() == t[2]) &&
           (s.target() == t[0] || s.target() == t[1] || s.target() == t[2]) )
        return false;
      else
        return true;
    }

    dolfin_assert(false);
    return false;
  }
  //-----------------------------------------------------------------------------
  // Check if two triangles intersect.
  // Neighbor triangles (share a vertice or an edge) do not intersect
  // if t1 == t2 (geometrically), the triangles do not intersect 
  template<typename Triangle_3>
  static bool triangles_intersect(const Triangle_3& t1, const Triangle_3& t2)
  {
    typedef typename CGAL::Kernel_traits<Triangle_3>::Kernel::Point_3 Point_3;
    typedef typename CGAL::Kernel_traits<Triangle_3>::Kernel::Segment_3 Segment_3;

    if ( (t1[0] == t2[0] || t1[0] == t2[1] || t1[0] == t2[2]) &&
         (t1[1] == t2[0] || t1[1] == t2[1] || t1[1] == t2[2]) &&
         (t1[2] == t2[0] || t1[2] == t2[1] || t1[2] == t2[2]))
      return false;
    
    auto result = CGAL::intersection(t1, t2);

    if (!result)
      return false;
    
    if (const Point_3* p = boost::get<Point_3>(&*result))
    {
      if (t1[0] == t2[0] || t1[0] == t2[1] || t1[0] == t2[2] ||
          t1[1] == t2[0] || t1[1] == t2[1] || t1[1] == t2[2] ||
          t1[2] == t2[0] || t1[2] == t2[1] || t1[2] == t2[2])
        return false;
      else
        return true;
    }
    else if (const Segment_3* s = boost::get<Segment_3>(&*result))
    {
      std::size_t common_vertices = 0;
      if (t1[0] == t2[0] || t1[0] == t2[1] || t1[0] == t2[2])
        common_vertices++;

      if (t1[1] == t2[0] || t1[1] == t2[1] || t1[1] == t2[2])
        common_vertices++;

      if (t1[2] == t2[0] || t1[2] == t2[1] || t1[2] == t2[2])
        common_vertices++;

      if (common_vertices > 1)
        return false;
      else
        return true;
    }
    else if (const Triangle_3* t = boost::get<Triangle_3>(&*result))
      return true;
    else if (const std::vector<Point_3>* v = boost::get<std::vector<Point_3> >(&*result))
      return true;

    dolfin_assert(false);
    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Triangle_3>
  static bool triangle_set_intersects(const std::vector<Triangle_3>& t)
  {
    for (std::size_t i = 0; i < t.size(); i++)
    {
      for (std::size_t j = i+1; j < t.size(); j++)
      {
        if (triangles_intersect<Triangle_3>(t[i], t[j]))
          return true;
      }
    }
    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Halfedge_handle>
  static std::size_t vertex_count_halfedges(Halfedge_handle h)
  {
    std::size_t count = 0;
    Halfedge_handle current = h;
    do
    {
      ++count;
      current = current->next()->opposite();
    } while (current != h);

    dolfin_assert(count > 0);
    return count;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
    static std::size_t total_vertex_count_halfedges(const Polyhedron& P, std::map<typename Polyhedron::Vertex_const_handle, std::size_t>& m)
  {
    std::size_t total_count = 0;
    for (typename Polyhedron::Vertex_const_iterator it = P.vertices_begin(); it != P.vertices_end(); it++)
    {
      const std::size_t c = vertex_count_halfedges(it->halfedge());
      m[it] = c;
      // std::cout << " Halfedge count: " << c << std::endl;
      if (c == 0)
      {int tmp; std::cin >> tmp; }
      total_count += c;
    }
    return total_count;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool halfedge_is_in_polyhedron(const Polyhedron& P,
                                        typename Polyhedron::Halfedge_const_handle h)
  {
    for (auto hit = P.halfedges_begin(); hit != P.halfedges_end(); hit++)
    {
      if (hit == h)
        return true;
    }
    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool vertex_is_in_polyhedron(const Polyhedron& P,
                                      typename Polyhedron::Vertex_const_handle v)
  {
    for (auto vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
    {
      if (vit == v)
        return true;
    }

    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
    static bool check_vertex_consistency(const Polyhedron& P)
  {
    // std::cout << "Checking vertex consistency" << std::endl;
    std::size_t counter = 0;

    // Build a set of the list of vertices for faster lookup
    std::set<typename Polyhedron::Vertex_const_handle> vertex_set;
    for (auto vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
    {
      dolfin_assert(vertex_set.count(vit) == 0);
      vertex_set.insert(vit);
    }

    std::deque<typename Polyhedron::Halfedge_const_handle> queue;
    std::set<typename Polyhedron::Halfedge_const_handle> visited;

    queue.push_back(P.halfedges_begin());
    while (!queue.empty())
    {
      typename Polyhedron::Halfedge_const_handle current = queue.back();
      queue.pop_back();
      if (visited.count(current) == 0)
      {
        counter++;
        typename Polyhedron::Halfedge_const_handle start = current;
        visited.insert(current);

        // Walk around the facet (or hole). Check halfedges and queue opposites
        do
        {
          if (vertex_set.count(current->vertex()) == 0)
          {
            // std::cout << "Vertex not in vertex list: " << current->vertex()->point() << std::endl;
            return false;
          }

          // TODO: Add the opposite check: All vertices should be reachable via halfedges
          queue.push_back(current->opposite());
          current = current->next();
        } while(current != start);
      }
    }

    // std::cout << "  Checked " << counter << " halfedges" << std::endl;
    return true;
  }
  //-----------------------------------------------------------------------------
  template<typename Triangle_3>
  static bool triangle_set_intersect_triangle(Triangle_3 t,
                                              const std::vector<Triangle_3>& triangle_set)
  {
    for (auto tit = triangle_set.begin(); tit != triangle_set.end(); tit++)
    {
      if (triangles_intersect(*tit, t))
        return true;
    }
    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Segment_3>
    static bool segment_intersects_triangle_set(const Segment_3& s,
                                                const std::vector<typename CGAL::Kernel_traits<Segment_3>::Kernel::Triangle_3>& triangle_set)
  {
    for (auto tit = triangle_set.begin(); tit != triangle_set.end(); tit++)
    {
      if (segment_intersects_triangle(s, *tit))
        return true;
    }
    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool center_vertex_triangulation(Polyhedron& P,
                                          const typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

    // check if the triangulation with a center vertex intersects any of the
    // neighbor triangles
    // std::cout << "Attempting center vertex triangulation" << std::endl;
    dolfin_assert(P.is_valid(false, 0));
    dolfin_assert(halfedge_is_in_polyhedron(P, h));
    dolfin_assert(vertex_is_in_polyhedron(P, h->vertex()));
    
    const std::array<double, 3> plane_fit = get_plane_fit<Polyhedron>(h, h->prev());
    if (plane_fit[0] < .85)
    {
      // std::cout << "  Rejected. Not sufficiently planar: " << plane_fit[0] << std::endl;
      return false;
    }

    // Collect set of neighbor triangles and compute centroid
    std::vector<Triangle_3> triangles;
    Point_3 centroid = CGAL::ORIGIN;
    std::size_t counter = 0;
    {
      Halfedge_handle current = h;
      do
      {
        if (!current->opposite()->is_border() && current->opposite()->facet()->facet_degree() == 3)
        {
          triangles.push_back(get_facet_triangle<Polyhedron>(current->opposite()));
        }
        centroid = centroid + (current->vertex()->point()-CGAL::ORIGIN);
        counter++;
        
        current = current->next();
      } while (current != h);
    }
    // std::cout << "Number of triangles: " << triangles.size() << std::endl;
    centroid = CGAL::ORIGIN + (centroid-CGAL::ORIGIN)/counter;
    // std::cout << "Centroid: " << centroid << std::endl;
    
    Halfedge_handle current = h;
    do
    {
      Triangle_3 current_triangle(current->vertex()->point(), current->next()->vertex()->point(), centroid);
      for (auto tit = triangles.begin(); tit != triangles.end(); tit++)
      {
        if (triangles_intersect(current_triangle, *tit))
        {
          // std::cout << "No: Triangle " << current_triangle[0] << " " << current_triangle[1] << " " << current_triangle[2] << std::endl;
          // std::cout << "Triangle " << (*tit)[0] << " " << (*tit)[1] << " " << (*tit)[2] << std::endl;
          return false;
        }
      }
      
      triangles.push_back(current_triangle);
      
      current = current->next();
    } while (current != h);

    // std::cout << "Facet degree before center vertex: " << h->facet()->facet_degree() << std::endl;

    P.normalize_border();
    dolfin_assert(P.is_valid(false, 1));
    dolfin_assert(!h->is_border_edge());
    dolfin_assert(h->facet()->facet_degree() > 3);

    Halfedge_handle center = P.create_center_vertex(h);
    /* Halfedge_handle g = h->next()->next(); */
    /* dolfin_assert(check_vertex_consistency(P)); */
    /* dolfin_assert(halfedge_is_in_polyhedron(P, h)); */
    /* dolfin_assert(vertex_is_in_polyhedron(P, h->vertex())); */
    /* dolfin_assert(vertex_is_in_polyhedron(P, g->vertex())); */
    /* std::cout << "Splitting: Segment " << h->vertex()->point() << ", " << g->vertex()->point() << std::endl; */
    /* std::cout << "Vertex degrees: " << h->vertex()->vertex_degree() << ", " << g->vertex()->vertex_degree() << std::endl; */
    /* std::cout << "My count: " << vertex_count_halfedges(h) << ", " << vertex_count_halfedges(g) << std::endl; */

    /* std::cout << "Facet degree: " << h->facet()->facet_degree() << std::endl; */
    /* std::cout << "Opposite facet degree: " << h->opposite()->facet()->facet_degree() << std::endl; */
    /* std::cout << "Num halfedges: " << P.size_of_halfedges() << std::endl; */
    /* std::map<typename Polyhedron::Vertex_const_handle, std::size_t> vertex_degrees; */
    /* std::cout << "My total count: " << total_vertex_count_halfedges(P, vertex_degrees); */
    /* std::cout << "Mapped: " << vertex_degrees.at(h->vertex()) << ", " << vertex_degrees.at(g->vertex()) << std::endl; */
    /* std::cout << "Size of map: " << vertex_degrees.size() << std::endl; */
    /* std::cout << "-- Splitting facet --" << std::endl; */

    /* Halfedge_handle diagonal = P.split_facet(h, g); */

    /* std::cout << "Vertex degrees: " << h->vertex()->vertex_degree() << ", " << g->vertex()->vertex_degree() << std::endl; */
    /* std::cout << "My count: " << vertex_count_halfedges(h) << ", " << vertex_count_halfedges(g) << std::endl; */
    /* std::map<typename Polyhedron::Vertex_const_handle, std::size_t> vertex_degrees_after; */
    /* std::cout << "My total count: " << total_vertex_count_halfedges(P, vertex_degrees_after) << std::endl; */
    /* std::cout << "Size of map: " << vertex_degrees_after.size() << std::endl; */
    /* std::cout << "Num halfedges: " << P.size_of_halfedges() << std::endl; */
    /* std::cout << "Mapped: " << vertex_degrees_after.at(h->vertex()) << ", " << vertex_degrees_after.at(g->vertex()) << std::endl; */
    /* for (auto it = vertex_degrees_after.begin(); it != vertex_degrees_after.end(); ++it) */
    /* { */
    /*   if (vertex_degrees.at(it->first) != it->second) */
    /*   { */
    /*     std::cout << "DIFF!!!" << vertex_degrees[it->first] << " " << it->second << std::endl; */
    /*   } */

    /*   if (vertex_degrees.at(it->first) == 0 || it->second == 0) */
    /*   { */
    /*     std::cout << "zero degree " << it->first->point() << std::endl; */
    /*   } */
    /* } */
    /* dolfin_assert(P.is_valid(false)); */
    /* Halfedge_handle c = diagonal->vertex()->halfedge(); */
    /* Halfedge_handle start = c; */
    /* do */
    /* { */
    /*   if (c->opposite()->vertex() == diagonal->opposite()->vertex()) */
    /*     std::cout << "Yes!!!" << std::endl; */
    /*   c = c->opposite()->next(); */
    /* } while (c != start); */
    /* dolfin_assert(P.is_valid(false, 0)); */
    /* std::cout << "Splitting edge: Segment " << diagonal->opposite()->vertex()->point() << ", " << diagonal->vertex()->point() << std::endl; */
    /* Halfedge_handle center = P.split_edge(diagonal); */
    center->vertex()->point() = centroid;
    dolfin_assert(P.is_valid(false, 0));

    /* std::cout << "Splitting facet: " << diagonal->opposite()->vertex()->point() << ", " << diagonal->opposite()->prev()->prev()->vertex()->point() << std::endl; */
    /* P.split_facet(diagonal->opposite(), diagonal->opposite()->prev()->prev()); */

    /* do */
    /* { */
    /*   dolfin_assert(P.is_valid(false, 0)); */
    /*   std::cout << "adding diagonal: " << diagonal->next()->vertex()->point() << ", " << center->vertex()->point() << std::endl; */
    /*   diagonal = P.split_facet(diagonal->next(), center); */
    /*   diagonal = diagonal->opposite(); */

    /* } while (diagonal->next()->next() != center); */
    /* std::cout << "Center vertex degree: " << center->vertex()->vertex_degree() << std::endl; */
    /* //P.normalize_border(); */
    /* { */
    /*   std::ofstream ofile("center-vertex.off"); */
    /*   ofile << P; */
    /* } */

    /* dolfin_assert(P.is_valid()); */

    return true;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void triangulate_quad_facet(Polyhedron& P,
                             typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;

    // std::cout << "Tringulating quad" << std::endl;
    dolfin_assert(h->next()->next()->next()->next() == h);
    
    //P.fill_hole(h);

    if (get_triangle_cos_angle(Triangle_3(h->vertex()->point(),
                                          h->next()->vertex()->point(),
                                          h->next()->next()->vertex()->point()),
                               Triangle_3(h->vertex()->point(),
                                          h->next()->next()->vertex()->point(),
                                          h->prev()->vertex()->point())) <
        get_triangle_cos_angle(Triangle_3(h->next()->vertex()->point(),
                                          h->next()->next()->vertex()->point(),
                                          h->prev()->vertex()->point()),
                               Triangle_3(h->vertex()->point(),
                                          h->next()->vertex()->point(),
                                          h->prev()->vertex()->point())))
    {
      // The edge should go from h <--> h->next()->next()
      P.split_facet(h, h->next()->next());
    }
    else
    {
      // the edge should ho from h->next() <--> h->prev()
      P.split_facet(h->next(), h->prev());      
    }
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void triangulate_5_facet(Polyhedron& P,
                                  typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;
    
    double best_quality = 1;
    Halfedge_handle best_triangle;
    Halfedge_handle current = h;
    do
    {
      const double candidate_quality = get_triangle_cos_angle(Triangle_3(current->next()->vertex()->point(),
                                                                   current->next()->next()->vertex()->point(),
                                                                   current->next()->next()->next()->vertex()->point()),
                                                        Triangle_3(current->prev()->prev()->vertex()->point(),
                                                                   current->prev()->vertex()->point(),
                                                                   current->vertex()->point()));
      if (candidate_quality < best_quality)
      {
        best_triangle = current;
        best_quality = candidate_quality;
      }
      
      current = current->next();
    } while (current != h);

    //P.fill_hole(h);
    P.split_facet(best_triangle->next(),
                  best_triangle->next()->next()->next());
    P.split_facet(best_triangle->next()->next(), best_triangle);
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void close_hole(Polyhedron& P,
                         typename Polyhedron::Halfedge_handle h)
  {
    //typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

    // std::cout << "----------- Closing hole --------------------" << std::endl;
    dolfin_assert(P.is_valid(false, 0));
    dolfin_assert(P.is_pure_triangle());
    dolfin_assert(h->is_border());

    P.fill_hole(h);
    P.normalize_border();

    // std::cout << "Size of hole: " << h->facet()->facet_degree() << std::endl;
    list_hole<Polyhedron>(h);

    dolfin_assert(h->facet()->facet_degree() > 2);

    // Since the facet may be split, we push the facets to a fifo queue.
    std::deque<Halfedge_handle> queue;
    queue.push_back(h);

    while (!queue.empty())
    {
      // std::cout << "--- Popping facet from queue (" << queue.size() << ")" << std::endl;
      const Halfedge_handle current = queue.front();
      queue.pop_front();

      list_hole<Polyhedron>(current);
            
      dolfin_assert(P.is_valid(false, 0));
      dolfin_assert(halfedge_is_in_polyhedron(P, current));
      
      if (current->facet()->facet_degree() == 3)
      {
        //P.fill_hole(h);
        dolfin_assert(current->opposite()->facet()->facet_degree() != 3  ||
                      !triangles_intersect(get_facet_triangle<Polyhedron>(current),
                                           get_facet_triangle<Polyhedron>(current->opposite())));
        dolfin_assert(current->next()->opposite()->facet()->facet_degree() != 3  ||
                      !triangles_intersect(get_facet_triangle<Polyhedron>(current),
                                           get_facet_triangle<Polyhedron>(current->next()->opposite())));
        dolfin_assert(current->prev()->opposite()->facet()->facet_degree() != 3  ||
                      !triangles_intersect(get_facet_triangle<Polyhedron>(current),
                                           get_facet_triangle<Polyhedron>(current->prev()->opposite())));
      }
      else if (current->facet()->facet_degree() == 4)
      {
        triangulate_quad_facet(P, current);
      }
      /* else if (current->facet()->facet_degree() == 5) */
      /* { */
      /*   triangulate_5_facet(P, current); */
      /* } */
      else
      {
        // std::cout << "Attempting to triangulate in 2D" << std::endl;
        dolfin_assert(halfedge_is_in_polyhedron(P, current));
        if (!triangulate_polygon_3d(P, current, false))
        {
          dolfin_assert(halfedge_is_in_polyhedron(P, current));

          // if (!center_vertex_triangulation(P, current))
          // {
            //std::pair<Halfedge_handle, Halfedge_handle> facets = subdivide_facet(P, current);
            Halfedge_handle facet = subdivide_facet(P, current);
            // dolfin_assert(facets.first != Halfedge_handle());
            // dolfin_assert(facets.second != Halfedge_handle());

            queue.push_back(facet->opposite());
            queue.push_back(facet);
            //queue.push_back(facets.first);
            //queue.push_back(facets.second);
            // }
        }

        /* std::cout << "Points are not coplanar" << std::endl; */
        /* Polyhedron hole; */
        /* std::cout << "Compute convex hull" << std::endl; */
        /* CGAL::convex_hull_3(points.begin(), points.end(), hole); */
        /* {std::ofstream f("convex_hull.off"); f << hole; } */

        /* // Find intersection polyline */
        /* std::list<std::vector<Point_3> > intersection_polylines; */

        /* // Search for a common vertex */
        /* /\* bool found; *\/ */
        /* /\* Halfedge_handle current = h; *\/ */
        /* /\* do *\/ */
        /* /\* { *\/ */
          
        /* /\*   current = current->next(); *\/ */
        /* /\* } while (current != h); *\/ */
        
        /* CGALCSGOperator op; */

        /* std::cout << "Computing union" << std::endl; */
        /* op(P, hole, std::back_inserter(intersection_polylines), CGALCSGOperator::Join_tag); */
        /* std::cout << "Done computing union" << std::endl; */
      }
    }
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static typename Polyhedron::Halfedge_handle
  subdivide_facet(Polyhedron& P, typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Segment_3 Segment_3;
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;

    // Search for triangles that divide the hole, such that the dividing triangle
    // does not intersect triangles next to the hole

    // Store all triangles around the hole and compute max edge length
    std::vector<Triangle_3> border_triangles;
    double max_squared_edge_length = 0;
    {
      Halfedge_handle current = h;
      do
      {
        max_squared_edge_length = std::max(max_squared_edge_length,
                                           CGAL::to_double(Segment_3(current->prev()->vertex()->point(),
                                                                     current->vertex()->point()).squared_length()));
        if (current->opposite()->facet()->facet_degree() == 3)
        {
          border_triangles.push_back(get_facet_triangle<Polyhedron>(current->opposite()));
        }
        current = current->next();
      } while (current != h);
    }

    // std::cout << "Number of border triangles: " << border_triangles.size() << std::endl;
    // std::cout << "Plane fit: " << get_plane_fit<Polyhedron>(h, h)[0] << std::endl;
    dolfin_assert(border_triangles.size() > 4);

    // Search for the best dividing segment
    double best_quality = -1000;
    Halfedge_handle best_outer;
    Halfedge_handle best_inner;

    Halfedge_handle current_outer = h;

    // FIXME: This loop should run to h->prev()->prev(), but need
    // handle the specially
    const Halfedge_handle outer_end = h->prev()->prev();
    do
    {
      Halfedge_handle current_inner = current_outer->next()->next();
      const Halfedge_handle inner_end = h;
      do
      {
        if (current_inner->next() != current_outer &&
            current_inner->prev() != current_outer)
        {
        
          Segment_3 current_segment(current_outer->vertex()->point(),
                                    current_inner->vertex()->point());

          /* std::cout << "Checking segment: Segment " << current_segment.source() << ", " << current_segment.target() */
          /*           << " length: " << current_segment.squared_length() << std::endl; */
        
          // Check that this does not introduce an intersection
          if (!segment_intersects_triangle_set(current_segment, border_triangles))
          {
            //const double side1_quality = get_plane_fit<Polyhedron>(current_outer, current_inner);
            //const double side2_quality = get_plane_fit<Polyhedron>(current_inner, current_outer);
            //const double candidate_quality = std::min(side1_quality, side2_quality);
            const double candidate_quality = evaluate_hole_subdivision<Polyhedron>(current_inner, current_outer);
            //std::cout << "Plane qualities: " << side1_quality << ", " << side2_quality << std::endl;
            // std::cout << "Segment " << current_inner->vertex()->point() << ", " << current_outer->vertex()->point() << " : " << candidate_quality << std::endl;
            // std::cout << "Quality: " << candidate_quality << std::endl;

            if (candidate_quality > best_quality)
            {
              best_outer = current_outer;
              best_inner = current_inner;
              best_quality = candidate_quality;
            }
          }
        }
        
        current_inner = current_inner->next();
      } while (current_inner != inner_end);
      
      current_outer = current_outer->next();
    } while (current_outer != outer_end);

    dolfin_assert(best_outer != Halfedge_handle());
    dolfin_assert(best_inner != Halfedge_handle());

    // std::cout << "Found best subdivision: " << std::endl;

    // list_hole<Polyhedron>(best_outer);
    // std::cout << "Segment " << best_outer->vertex()->point()
    //           << ", " << best_inner->vertex()->point() << std::endl;
    // std::cout << "Quality: " << best_quality << std::endl;

    /* if (best_quality > .9) */
    /* { */
      dolfin_assert(P.is_valid(false, 0));
      const Halfedge_handle new_diagonal = P.split_facet(best_inner, best_outer);
      const Point_3& p = new_diagonal->opposite()->vertex()->point();

      const Vector_3 new_edge(new_diagonal->opposite()->vertex()->point(),
                              new_diagonal->vertex()->point());

      const int num_segments = static_cast<int>(sqrt(CGAL::to_double(new_edge.squared_length())/max_squared_edge_length)+.5);

      // std::cout << "Num segments: " << num_segments << std::endl;

      // Note: Don't use std::size_t as 0-1 becomes very large...
      for (int i = 1; i < num_segments; i++)
      {
        //std::cout << "Splitting segment" << std::endl;
        Halfedge_handle new_segment = P.split_edge(new_diagonal);
        new_segment->vertex()->point() = p + static_cast<double>(i)/num_segments * new_edge;
      }

      P.normalize_border();
      dolfin_assert(P.is_valid(false, 0));
      //return std::make_pair(new_diagonal, new_diagonal->opposite());
      return new_diagonal;
    /* } */
    /* else */
    /*   //return std::make_pair(Halfedge_handle(), Halfedge_handle()); */
    /*   return Halfedge_handle(); */

    /* Halfedge_handle v1, v2; */
    /* best_quality = -1000; */

    /* { */
    /*   // Candidate 1: best_inner, best_inner->next(), best_outer */
    /*   std::cout << "Candidate 1" << std::endl; */
    /*   // Check the four candidates where the chosen segment is an edge */
    /*   if (best_inner->next()->next() != best_outer) */
    /*   { */
    /*     const double candidate_quality= std::min(get_plane_fit<Polyhedron>(best_outer, best_inner), */
    /*                                              get_plane_fit<Polyhedron>(best_inner->next(), best_outer)); */

    /*     if (candidate_quality > best_quality) */
    /*     { */
    /*       v1 = best_outer; */
    /*       v2 = best_inner; */
    /*       best_quality = candidate_quality; */
    /*     } */
    /*   } */
    /* } */

    /* if (best_outer->next()->next() != best_inner) */
    /* { */
    /*   // Candidate 2: best_inner, best_outer, best_outer->next() */
    /*   double candidate_quality = std::min(get_plane_fit<Polyhedron>(best_outer->next(), best_inner), */
    /*                                       get_plane_fit<Polyhedron>(best_inner, best_outer)); */

    /*   if (candidate_quality > best_quality) */
    /*   { */
    /*     std::cout << "Candidate 2" << std::endl; */
    /*     v1 = best_inner; */
    /*     v2 = best_outer; */
    /*     best_quality = candidate_quality; */
    /*   } */
    /* } */

    /* if (best_outer->next()->next() != best_inner) */
    /* { */
    /*   // Candidate 3: best_inner->prev(), best_inner, best_outer */
    /*   const double candidate_quality = std::min(get_plane_fit<Polyhedron>(best_outer, best_inner->prev()), */
    /*                                             get_plane_fit<Polyhedron>(best_inner, best_outer)); */

    /*   if (candidate_quality > best_quality) */
    /*   { */
    /*     std::cout << "Candidate 3" << std::endl; */
    /*     v1 = best_outer; */
    /*     v2 = best_inner->prev(); */
    /*     best_quality = candidate_quality; */
    /*   } */
    /* } */

    /* if (best_inner->next()->next() != best_outer) */
    /* { */
    /*   // Candidate 4: best_inner, best_outer->prev(), best_outer */
    /*   const double candidate_quality = std::min(get_plane_fit<Polyhedron>(best_outer, best_inner), */
    /*                                             get_plane_fit<Polyhedron>(best_inner, best_outer->prev())); */

    /*   if (candidate_quality > best_quality) */
    /*   { */
    /*     std::cout << "Candidate 4" << std::endl; */
    /*     v1 = best_inner; */
    /*     v2 = best_outer->prev(); */
    /*   } */
    /* } */

    /* dolfin_assert(v1 != Halfedge_handle()); */
    /* dolfin_assert(v2 != Halfedge_handle()); */

    /* // Divide hole by chosen triangle (v1, v2, v2->next()) */
    /* Halfedge_handle f = P.fill_hole(h); */
    /* dolfin_assert(v1->facet() == v2->facet()); */

    /* Halfedge_handle hole1 = P.split_facet(v1, v2->next()); */
    /* Halfedge_handle hole2 = P.split_facet(v2, hole1->opposite()); */
    /* P.make_hole(hole1); */
    /* P.make_hole(hole2); */
    
    /* return true; */
  }
  //-----------------------------------------------------------------------------
  /* template<typename Polyhedron> */
  /* static double evaluate_heuristic(const Polyhedron& P, */
  /*                                  typename Polyhedron::Halfedge_handle h, */
  /*                                  double plane_fit) */
  /* { */
  /*   typedef typename Polyhedron::Traits::Triangle_3 Triangle_3; */
  /*   typedef typename Polyhedron::Traits::Vector_3 Vector_3; */
  /*   typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel; */
  /*   typedef typename InexactKernel::Plane_3 InexactPlane_3; */


  /*   // const double distance_to_plane_weight = 1.0; */
  /*   const double planarity_weight = 1.0; */
  /*   const double dihedral_weight  = 1.0; */
  /*   const double ear_angle_weight = 1.0; */

  /*   // Compute the planarity of the points excluding the current point */
  /*   InexactPlane_3 p; */
  /*   const double plane_fit_quality = get_plane_fit<Polyhedron>(h->next(), */
  /*                                                              h->prev(), */
  /*                                                              &p); */
  /*   // Compute the maximum of the dihedral angle to the neighbors */
  /*   const Triangle_3 candidate_triangle(h->prev()->vertex()->point(), */
  /*                                       h->vertex()->point(), */
  /*                                       h->next()->vertex()->point()); */
  /*   const double cos_dihedral = (std::min(get_triangle_cos_angle(candidate_triangle, */
  /*                                                                get_facet_triangle<Polyhedron>(h->opposite())), */
  /*                                         get_triangle_cos_angle(candidate_triangle, */
  /*                                                                get_facet_triangle<Polyhedron>(h->next()->opposite())))+1)/2; */

  /*   // Compute the angle of the cutted ear */
  /*   const Vector_3 v1(h->vertex()->point(), */
  /*                     h->prev()->vertex()->point()); */
  /*   const Vector_3 v2(h->vertex()->point(), */
  /*                     h->next()->vertex()->point()); */

  /*   const double cos_ear_angle = (CGAL::to_double((v1*v2)/std::sqrt(CGAL::to_double(v1.squared_length()*v2.squared_length())))+1)/2.0; */
  /*   const double ear_angle_quality = cos_ear_angle; */

  /*   std::cout << "Triangle " << candidate_triangle[0] */
  /*             << ", " << candidate_triangle[1]  */
  /*             << "," << candidate_triangle[2] << std::endl; */

  /*   std::cout << "Evaluate: planarity: " << (plane_fit_quality/plane_fit) */
  /*             << ", dihedral: " << cos_dihedral */
  /*             << ", ear angle: " << ear_angle_quality << std::endl; */

  /*   return planarity_weight*plane_fit_quality + dihedral_weight*cos_dihedral + ear_angle_quality*ear_angle_weight; */
    
  /*   /\* return planarity_weight*plane_fit_quality/plane_fit + *\/ */
  /*   /\*   dihedral_weight*cos_dihedral + *\/ */
  /*   /\*   (1-std::exp(-ear_angle_weight*cos_ear_angle))*cos_dihedral; *\/ */
  /* } */
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool vertex_is_border(typename Polyhedron::Vertex_const_handle v)
  {
    //typename Polyhedron::Vertex::Halfedge_around_vertex_circulator h_start = v->vertex_begin();
    auto h_start = v->vertex_begin();
    auto h_current = h_start;
    do 
    {
      if (h_current->is_border_edge())
        return true;
      h_current++;
    } while(h_current != h_start);

    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool facets_are_neighbors(typename Polyhedron::Facet_handle f1,
                                   typename Polyhedron::Facet_handle f2)
  {
    typename Polyhedron::Halfedge_handle h1 = f1->halfedge();
    typename Polyhedron::Halfedge_handle start1 = h1;
    do
    {
      typename Polyhedron::Halfedge_handle h2 = f2->halfedge();
      typename Polyhedron::Halfedge_handle start2 = h2;
      do 
      {
        if (h2->vertex() == h1->vertex())
          return true;

        h2 = h2->next();
      } while (h2 != start2);

      h1 = h1->next();
    } while (h1 != start1);
    
    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool segment_intersects_facets(typename Polyhedron::Vertex_handle v1,
                                        typename Polyhedron::Vertex_handle v2)
  {
    typedef typename Polyhedron::Halfedge_around_vertex_circulator Vertex_circulator;
    typedef typename Polyhedron::Traits::Segment_3 Segment_3;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;
        typedef typename Polyhedron::Traits::Point_3 Point_3;

    Segment_3 s(v1->point(), v2->point());
    
    Vertex_circulator start = v1->vertex_begin();
    Vertex_circulator current = start;
    do
    {
      if (!current->is_border() && current->facet()->is_triangle())
      {
        Triangle_3 t = get_facet_triangle<Polyhedron>(current);
        auto result = CGAL::intersection(t, s);
        dolfin_assert(result);

        if (const Point_3* p = boost::get<Point_3>(&*result))
        {
          dolfin_assert(*p == s.source() || *p == s.target());
        }
        else if (const Segment_3* s = boost::get<Segment_3>(&*result))
        {
          return true;
        }
        else
        {
          dolfin_assert(false);
        }
      }
      
      current++;
    } while (start != current);

    return false;
  }

  //-----------------------------------------------------------------------------
  /// Check if the triangle defined by (v, v2, v3) intersects any of the
  /// triangular facets adjacent to v
  template<typename Polyhedron>
  static bool facets_intersect(typename Polyhedron::Vertex_handle v,
                               typename Polyhedron::Vertex_handle v2,
                               typename Polyhedron::Vertex_handle v3)
  {
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;
    typedef typename Polyhedron::Traits::Segment_3 Segment_3;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Halfedge_around_vertex_circulator Vertex_circulator;

    Triangle_3 t(v->point(), v2->point(), v3->point());

    std::cout << "Facets_intersect()" << std::endl;
    std::cout << "--Vertex: " << v->point() << std::endl;
    std::cout << "--Trial triangle: " << t[0] << ", " << t[1] << ", " << t[2] << std::endl;

    dolfin_assert(v->point() == t[0] || v->point() == t[1] || v->point() == t[2]);

    Vertex_circulator start = v->vertex_begin();
    Vertex_circulator c = start;

    do
    {
      Halfedge_handle current = c;
      c++;
      // std::cout << "  Checking" << std::endl;
      if (!current->is_border())
      {
        dolfin_assert(current->facet()->is_triangle());

        Triangle_3 t_current(current->vertex()->point(),
                             current->next()->vertex()->point(),
                             current->next()->next()->vertex()->point());
        std::cout << "  Test triangle: " << t_current[0] << ", " << t_current[1] << ", " << t_current[2] << std::endl;

        auto result = CGAL::intersection(t, t_current);

        if (result)
        {
          if (const Point_3* p = boost::get<Point_3>(&*result))
          {
            dolfin_assert(*p == v->point());
            std::cout << "  Vertex point" << std::endl;
            std::cout << *p << std::endl;
            continue;
          }
          else if (const Segment_3* s = boost::get<Segment_3>(&*result))
          {
            std::cout << "  Segment intersection" << std::endl;

            // Check that they share a halfedge
            dolfin_assert(s->source() == v->point() || s->target() == v->point());
            std::array<Vertex_handle, 3> v_facet{v, v2, v3};
            std::array<Vertex_handle, 3> u_facet{current->vertex(), current->next()->vertex(), current->next()->next()->vertex()};

            bool neighbors_detected = false;
            for (std::size_t i = 0; i < 3; i++)
            {
              for (std::size_t j = 0; j < 3; j++)
              {
                if (v_facet[i] == u_facet[(j+1)%3] && v_facet[(i+1)%3] == u_facet[j])
                {
                  std::cout << "  Neighbors" << std::endl;
                  neighbors_detected = true;
                  break;
                }
              }
            }

            if (neighbors_detected)
              continue;
            /* if ( (s->source() == t_current[0] || s->source() == t_current[1] || s->source() == t_current[2]) && */
            /*      (s->target() == t_current[0] || s->target() == t_current[1] || s->target() == t_current[2]) ) */
            /* { */
            /*   std::cout << "  Edge segment" << std::endl; */
            /*   std::cout << s->source() << " " << s->target() << std::endl; */
            /*   continue; */
            /* } */
          }


          std::cout << "Facets intersect" << std::endl;

          if (const Point_3* p = boost::get<Point_3>(&*result))
          {
            std::cout << "    Intersection point: " << *p << std::endl;
            return true;
          }
          else if (const Segment_3* s = boost::get<Segment_3>(&*result))
          {
            std::cout << "    Intersection segment: " << s->source() << ", " << s->target() << std::endl;
            return true;
          }
          else if (const Triangle_3* tt = boost::get<Triangle_3>(&*result))
          {
            std::cout << "    Intersection triangle: " << (*tt)[0] << ", " << (*tt)[1] << ", " << (*tt)[2] << std::endl;
          }
          else if (const std::vector<Point_3>* v = boost::get<std::vector<Point_3> >(&*result))
          {
            std::cout << "    Intersection vector of points" << std::endl;
          }
          else
          {
            std::cout << "    Unknown intersection type" << std::endl;
            return true;
          }
        }
      }


    } while (c != start);

    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void list_self_intersections(Polyhedron& p)
  {
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Traits Polyhedron_traits;
    std::vector<std::pair<Facet_handle, Facet_handle> > intersections;
    CGAL::self_intersect<Polyhedron_traits>(p, std::back_inserter(intersections));

    for (auto iit = intersections.begin(); iit != intersections.end(); iit++)
    {
      std::cout << "Intersection (neighbors: " << (facets_are_neighbors<Polyhedron>(iit->first, iit->second) ? "Yes" : "No") << ")" << std::endl;
      print_triangle<Polyhedron>(iit->first->halfedge());
      print_triangle<Polyhedron>(iit->second->halfedge());
    }
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static double closest_vertices(const Polyhedron& p)
  {
    // TODO: Use a better optimal algorithm for closest pair problem
    std::cout << "Computing closest vertices" << std::endl;
    double min_distance = std::numeric_limits<double>::max();

    std::size_t counter = 0;
    std::cout << "Vertices: " << p.size_of_vertices() << std::endl;
    for (typename Polyhedron::Vertex_const_iterator v1 = p.vertices_begin();
         v1 != p.vertices_end(); v1++)
    {
      if (counter % 1000 == 0)
        std::cout << counter << std::endl;

      typename Polyhedron::Vertex_const_handle v2 = v1;
      v2++;
      std::size_t counter2 = 0;
      for (;v2 != p.vertices_end(); v2++)
      {
        min_distance = std::min(min_distance, 
                                CGAL::to_double(CGAL::squared_distance(v1->point(), v2->point())));

        if (min_distance == 0)
          return 0.0;

        counter2++;
      }

      counter++;
    }

    std::cout << "  Done computing closest" << std::endl;
    return min_distance;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static std::size_t min_vertex_degree(const Polyhedron& p)
  {
    std::size_t min_degree = std::numeric_limits<std::size_t>::max();
    for (typename Polyhedron::Vertex_const_iterator it = p.vertices_begin();
         it != p.vertices_end(); it++)
    {
      if (!vertex_is_border<Polyhedron>(it))
      {
        min_degree = std::min(min_degree, it->vertex_degree());
      }
    }
    return min_degree;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void remove_vertex(Polyhedron& P, typename Polyhedron::Vertex_handle v)
  {
    typedef typename Polyhedron::Halfedge_around_vertex_circulator Vertex_circulator;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    
    // std::cout << "Removing vertex" << std::endl;

    Vertex_circulator h = v->vertex_begin();
    Vertex_circulator start = h;

    std::vector<Halfedge_handle> to_be_removed;

    do
    {
      if (!h->is_border())
        to_be_removed.push_back(h);

      h++;
    } while (h != start);


    // std::cout << "Removing " << to_be_removed.size() << " halfedges" << std::endl;
    for (auto it = to_be_removed.begin(); it != to_be_removed.end(); it++)
    {
      P.erase_facet(*it);
    }

    // std::cout << "  done removing vertex" << std::endl;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static std::size_t remove_self_intersections(Polyhedron& P)
  {
    std::size_t removed = 0;

    typedef typename Polyhedron::Traits Polyhedron_traits;
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    std::vector<std::pair<Facet_handle, Facet_handle> > intersections;
    CGAL::self_intersect<Polyhedron_traits>(P, std::back_inserter(intersections));

    /* if (intersections.size() > 0) */
    /*   std::cout << "Removing self intersections" << std::endl; */

    while (intersections.size() > 0)
    {
      // std::cout << "  Removing pair" << std::endl;


      const typename Polyhedron::Facet_handle f1 = intersections.front().first;
      const typename Polyhedron::Facet_handle f2 = intersections.front().second;

      std::deque<Facet_handle> queue1;
      queue1.push_back(f1);

      std::deque<Facet_handle> queue2;
      queue2.push_back(f2);

      std::set<Facet_handle> to_be_removed1;
      to_be_removed1.insert(f1);

      std::set<Facet_handle> to_be_removed2;
      to_be_removed2.insert(f2);

      while (!queue1.empty())
      {
        {
          to_be_removed1.insert(queue1.front());
          Halfedge_handle start = queue1.front()->halfedge();
          queue1.pop_front();
          Halfedge_handle current = start;
          bool done = false;
          do
          {
            // std::cout << "Spreading out" << std::endl;
            if (!current->is_border_edge())
            {
              if (to_be_removed2.count(current->opposite()->facet()) > 0)
              {
                done = true;
                break;
              }
              else
              {
                queue1.push_back(current->opposite()->facet());
              }
            }
            current = current->next();
          } while (current != start);
          
          if (done)
            break;
        }

        {
          to_be_removed2.insert(queue2.front());
          Halfedge_handle start = queue2.front()->halfedge();
          queue2.pop_front();
          Halfedge_handle current = start;
          bool done = false;
          do
          {
            // std::cout << "Spreading out" << std::endl;
            if (!current->is_border_edge())
            {
              if (to_be_removed1.count(current->opposite()->facet()) > 0)
              {
                done = true;
                break;
              }
              else
              {
                queue2.push_back(current->opposite()->facet());
              }
            }
            current = current->next();
          } while (current != start);
          
          if (done)
            break;
        }
      }
      
      
      // std::cout << "To be removed 1: " << to_be_removed1.size() << std::endl;
      for (auto it = to_be_removed1.begin(); it != to_be_removed1.end(); it++)
      {
        P.erase_facet((*it)->halfedge());
      }

      // std::cout << "To be removed 2: " << to_be_removed2.size() << std::endl;
      for (auto it = to_be_removed2.begin(); it != to_be_removed2.end(); it++)
      {
        P.erase_facet((*it)->halfedge());
      }



      /* std::vector<std::pair<parent, Facet_handle> > tree; */

      /* tree.push_back(std::make_pair(0, f1)); */
      /* tree.push_back(std::make_pair(0, f1->halfedge()->opposite()->facet()); */
      
      /* bool done = false; */
      /* while (!done) */
      /* { */
      /*   std::size_t start = tree.size()/2; */
      /*   std::size_t level_size = tree.size()/2; */

      /*   std::cout << "Start: " << start << std::endl; */
      /*   std::cout << "Level size: " << level_size << std::endl; */

      /*   for (std::size_t i = 0; i < level_size; i++) */
      /*   { */
      /*     Facet_handle current = tree[start+i]; */
      /*     Facet_handle prev    = tree[(start-1)/2]; */

      /*     std::cout << "Her" << std::endl; */

      /*     Halfedge_handle h1 = current->halfedge(); */
      /*     if (h1->opposite()->facet() == prev) */
      /*       h1 = h1->next(); */
      /*     tree.push_back(h1->opposite()->facet()); */

      /*     h1 = h1->next(); */
      /*     if (h1->opposite()->facet() == prev) */
      /*       h1 = h1->next(); */

      /*     tree.push_back(h1->opposite()->facet()); */

      /*     if (tree[tree.size()-1] == f2) */
      /*     { */
      /*       done = true; */
      /*       break; */
      /*     } */
      /*   } */
      /* } */

      // std::cout << "Found route from f1 to f2" << std::endl;
      // break;

      intersections.clear();
      CGAL::self_intersect<Polyhedron_traits>(P, std::back_inserter(intersections));
    }

    // std::cout << "ok (" << facets.size() << " triangle pair(s))" << std::endl;

    return removed;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void filter_sharp_features(Polyhedron& P, int start_facet, double tolerance)
  {
    typedef typename Polyhedron::Facet_iterator Facet_iterator;
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;

    /* { */
    /*   std::ofstream outfile("before_filtering.off"); */
    /*   outfile << P; */
    /* } */
    
    /* std::cout << "Filter sharp features" << std::endl; */

    const double cos_tolerance = std::cos(tolerance);
    // std::cout << "tolerance: " << cos_tolerance << std::endl;
    Facet_iterator fit = P.facets_begin();
    for (int i = 0; i < start_facet; i++)
      fit++;

    std::deque<Facet_handle> queue;
    {
      Halfedge_handle h = fit->halfedge();
      // std::cout << "Starting facet: Triangle " << h->vertex()->point() << ", ";
      h = h->next();
      // std::cout << h->vertex()->point() << ", ";
      h = h->next();
      // std::cout << h->vertex()->point() << std::endl;
    }

    std::set<Facet_handle> visited;
    std::set<Facet_handle> to_be_removed;
    for (Facet_iterator fit = P.facets_begin(); fit != P.facets_end(); fit++)
      to_be_removed.insert(fit);

    // std::cout << "Number of facets: " << P.size_of_facets() << std::endl;

    queue.push_back(fit);
    while (!queue.empty())
    {
      // std::cout << "In queue" << std::endl;
      Facet_handle f = queue.front();
      queue.pop_front();

      if (visited.count(f) > 0)
      {
        // std::cout << "Already handled" << std::endl;
        continue;
      }
      
      visited.insert(f);
      to_be_removed.erase(f);

      const Halfedge_handle start = f->halfedge();
      Halfedge_handle current = start;
      do
      {
        // std::cout << "Exploring neighbor" << std::endl;
        if (!current->opposite()->is_border())
        { 
          
          Triangle_3 t1 = get_facet_triangle<Polyhedron>(current);
          Triangle_3 t2 = get_facet_triangle<Polyhedron>(current->opposite());
          if (get_triangle_cos_angle(t1, t2) > cos_tolerance)
          {
            queue.push_back(current->opposite()->facet());
          }
        }

        current = current->next();
      } while (current != start);
    }

    // std::cout << "Remove " << to_be_removed.size() << " facets" << std::endl;

    for (auto fit = to_be_removed.begin(); fit != to_be_removed.end(); fit++)
    {
      P.erase_facet( (*fit)->halfedge() );
    }

    /* { */
    /*   std::ofstream outfile("after_filtering.off"); */
    /*   outfile << P; */
    /* } */
  }
};

// Taken from demo/Polyhedron/Scene_nef_polyhedron_item.cpp in the
// CGAL source tree.
// Quick hacks to convert polyhedra from exact to inexact and
// vice-versa
template <class Polyhedron_input, class Polyhedron_output>
struct Copy_polyhedron_to
  : public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS>
{
  Copy_polyhedron_to(const Polyhedron_input& in_poly)
    : _in_poly(in_poly) {}

  void operator()(typename Polyhedron_output::HalfedgeDS& out_hds)
  {
    typedef typename Polyhedron_output::HalfedgeDS Output_HDS;
    //typedef typename Polyhedron_input::HalfedgeDS Input_HDS;

    CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

    typedef typename Polyhedron_input::Vertex_const_iterator Vertex_const_iterator;
    typedef typename Polyhedron_input::Facet_const_iterator  Facet_const_iterator;
    typedef typename Polyhedron_input::Halfedge_around_facet_const_circulator HFCC;

    builder.begin_surface(_in_poly.size_of_vertices(),
      _in_poly.size_of_facets(),
      _in_poly.size_of_halfedges());

    for(Vertex_const_iterator
      vi = _in_poly.vertices_begin(), end = _in_poly.vertices_end();
      vi != end ; ++vi)
    {
      typename Polyhedron_output::Point_3 p(::CGAL::to_double( vi->point().x()),
	::CGAL::to_double( vi->point().y()),
	::CGAL::to_double( vi->point().z()));
      builder.add_vertex(p);
    }

    typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
    Index index(_in_poly.vertices_begin(), _in_poly.vertices_end());

    for(Facet_const_iterator
      fi = _in_poly.facets_begin(), end = _in_poly.facets_end();
      fi != end; ++fi)
    {
      HFCC hc = fi->facet_begin();
      HFCC hc_end = hc;
      //     std::size_t n = circulator_size( hc);
      //     CGAL_assertion( n >= 3);
      builder.begin_facet ();
      do
      {
        builder.add_vertex_to_facet(index[hc->vertex()]);
        ++hc;
      } while( hc != hc_end);
      builder.end_facet();
    }
    builder.end_surface();
  } // end operator()(..)
private:
  const Polyhedron_input& _in_poly;
}; // end Copy_polyhedron_to<>

template <class Poly_A, class Poly_B>
void copy_to(const Poly_A& poly_a, Poly_B& poly_b)
{
  Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
  poly_b.delegate(modifier);
  // CGAL_assertion(poly_b.is_valid());
}
}
#endif
