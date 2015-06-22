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
//#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Self_intersection_polyhedron_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/corefinement_operations.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Delaunay_mesher_no_edge_refinement_2.h>

#include <cmath>
#include <deque>
#include <fstream>

#define TOLERANCE 1e-14

/* typedef CGAL::Exact_predicates_exact_constructions_kernel _K; */

/* namespace std */
/* { */
/*   // TODO: This is a hack to allow triangulation of exact polyhedrons.  A proper */
/*   // fix would be to provide a template specialization of compute_facet_normal */
/*   // for epeck. */
/*   inline _K::FT sqrt(_K::FT a) */
/*   { */
/*     return std::sqrt(CGAL::to_double(a)); */
/*   } */
/* } */

#include <CGAL/triangulate_polyhedron.h>

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
  static void insert_edge(Polyhedron& P,
                          typename Polyhedron::Vertex_handle h, 
                          typename Polyhedron::Vertex_handle g)
  {
    typedef typename Polyhedron::Halfedge_around_vertex_circulator Vertex_circulator;

    dolfin_assert(h != typename Polyhedron::Vertex_handle());
    dolfin_assert(g != typename Polyhedron::Vertex_handle());
    
    const Vertex_circulator h_start = h->vertex_begin();
    Vertex_circulator h_current = h_start;

    Vertex_circulator g_current;

    bool found = false;
    do
    {
      if (!h_current->is_border())
      {
        const Vertex_circulator g_start = g->vertex_begin();
        g_current = g_start;
        do
        {
          if (!g_current->is_border() && h_current->facet() == g_current->facet())
          {
            found = true;
            break;
          }
          g_current++;
        } while (g_start != g_current);

        if (found)
          break;
      }
      h_current++;
    } while (h_start != h_current);

    dolfin_assert(g_current != Vertex_circulator());
    dolfin_assert(found);
    dolfin_assert(h_current->facet() == g_current->facet());

    P.split_facet(h_current, g_current);
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

    std::cout << "get edge cos theta" << std::endl;
    
    Vector_3 h1_vec(h->vertex()->point(), h->prev()->vertex()->point());
    h1_vec = h1_vec/std::sqrt(CGAL::to_double(h1_vec.squared_length()));
    
    std::cout << "h1_vec_normalized: " << h1_vec << std::endl;
    Vector_3 h2_vec(h->vertex()->point(), h->next()->vertex()->point());
    h2_vec = h2_vec/std::sqrt(CGAL::to_double(h2_vec.squared_length()));
    std::cout << "h2_vec_normalized: " << h2_vec << std::endl;
    const double cos_theta = CGAL::to_double(h1_vec*h2_vec);
    std::cout << "Cos theta: " << cos_theta << std::endl;
    return cos_theta;
  }
  //-----------------------------------------------------------------------------
  // Compute the fit quality of the vertices from h1 to h2 both included
  template<typename Polyhedron>
  static double get_plane_fit(const typename Polyhedron::Halfedge_handle h1,
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
    //typedef typename InexactKernel::Vector_3 InexactVector_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;

    // std::cout << "Get plane fit" << std::endl;

    // std::cout << "Polygon ";
    std::vector<InexactPoint_3> points;
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
    std::cout << "Size: " << points.size() << std::endl;
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
    std::cout << "Plane: " << fitting_plane << std::endl;
    std::cout << "Length of normal: " << fitting_plane.orthogonal_vector().squared_length() << std::endl;
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

    std::cout << "Fit quality: " << fit_quality << ", max distance: " << max_distance << ", cos_angle: " << max_angle << std::endl;
    //return fit_quality;
    // return -max_distance;
    return CGAL::to_double(fit_quality - max_angle);
  }
  //-----------------------------------------------------------------------------
    // Compute the fit quality of the vertices from h1 to h2 both included
  template<typename Polyhedron>
  static double evaluate_hole_subdivision(const typename Polyhedron::Halfedge_handle h1,
                              const typename Polyhedron::Halfedge_handle h2)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    // typedef typename Polyhedron::Traits::Plane_3 Plane_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    // typedef typename Polyhedron::Traits::Vector_3 Vector_3;
    // typedef typename Polyhedron::Traits::FT FT;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;
    typedef typename InexactKernel::Point_3 InexactPoint_3;
    //typedef typename InexactKernel::Vector_3 InexactVector_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;

    // std::cout << "Get plane fit" << std::endl;

    // std::cout << "Polygon ";
    double plane1fit;
    InexactPlane_3 fitting_plane1;
    double max_distance_squared1 = 0;
    {
      std::vector<InexactPoint_3> points;
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

      const Point_3& p = h2->vertex()->point();
      // std::cout << ", " << p;
      points.push_back(InexactPoint_3(CGAL::to_double(p[0]),
                                      CGAL::to_double(p[1]),
                                      CGAL::to_double(p[2])));

      dolfin_assert(points.size() > 2);
      std::cout << "  Size: " << points.size() << std::endl;
      //std::cout << std::endl;
      plane1fit = CGAL::linear_least_squares_fitting_3(points.begin(),
                                                       points.end(),
                                                       fitting_plane1,
                                                       CGAL::Dimension_tag<0>());
      std::cout << "  Plane 1: " << fitting_plane1 << ", " << plane1fit << std::endl;
      std::cout << "  Length of normal: " << fitting_plane1.orthogonal_vector().squared_length() << std::endl;
      for (auto pit = points.begin(); pit != points.end(); pit++)
        max_distance_squared1 = std::max(max_distance_squared1, (*pit-fitting_plane1.projection(*pit)).squared_length());
      std::cout << "  Max squared distance: " << max_distance_squared1 << std::endl;
    }

    double plane2fit;
    InexactPlane_3 fitting_plane2;
    double max_distance_squared2 = 0;
    {
      std::vector<InexactPoint_3> points;
      Halfedge_handle current = h2;
      do
      {
        const Point_3& p = current->vertex()->point();
        // std::cout << ", " << p;
        points.push_back(InexactPoint_3(CGAL::to_double(p[0]),
                                        CGAL::to_double(p[1]),
                                        CGAL::to_double(p[2])));
        current = current->next();
      } while (current != h1);

      const Point_3& p = h1->vertex()->point();
      // std::cout << ", " << p;
      points.push_back(InexactPoint_3(CGAL::to_double(p[0]),
                                      CGAL::to_double(p[1]),
                                      CGAL::to_double(p[2])));

      dolfin_assert(points.size() > 2);
      std::cout << "  Size: " << points.size() << std::endl;
      //std::cout << std::endl;
      plane2fit = CGAL::linear_least_squares_fitting_3(points.begin(),
                                                       points.end(),
                                                       fitting_plane2,
                                                       CGAL::Dimension_tag<0>());
      std::cout << "  Plane 2: " << fitting_plane2 << ", fit: " << plane2fit << std::endl;
      std::cout << "  Length of normal: " << fitting_plane2.orthogonal_vector().squared_length() << std::endl;
      InexactPoint_3 max_distance_point;
      for (auto pit = points.begin(); pit != points.end(); pit++)
      {
        if ((*pit-fitting_plane2.projection(*pit)).squared_length() > max_distance_squared2)
        {
          max_distance_squared2 = (*pit-fitting_plane2.projection(*pit)).squared_length();
          max_distance_point = *pit;
        }
      }
      std::cout << "  Max squared segment: Segment " << max_distance_point << ", " << fitting_plane2.projection(max_distance_point) << std::endl;
      std::cout << "  Max squared distance: " << max_distance_squared2 << std::endl;
    }

    const double cos_angle = fitting_plane1.orthogonal_vector()*fitting_plane2.orthogonal_vector();
    std::cout << "  Angle: " << cos_angle << "(" << acos(cos_angle)/(2*DOLFIN_PI)*360 << ")" << std::endl;
    //return std::min(plane1fit, plane2fit) - cos_angle;
    return -std::max(max_distance_squared1, max_distance_squared2);
  }
  //-----------------------------------------------------------------------------
  template<typename Vector_3>
  static CGAL::Aff_transformation_3<typename CGAL::Kernel_traits<Vector_3>::Kernel>
    rotate_to_xy(Vector_3 a)
  {
    const typename CGAL::Kernel_traits<Vector_3>::Kernel::FT den = a[0]*a[0] + a[1]*a[1];
    return CGAL::Aff_transformation_3<typename CGAL::Kernel_traits<Vector_3>::Kernel>
      (1 - a[0]*a[0]*(1-a[2])/den, -a[0]*a[1]*(1-a[2])/den,   -a[0],
       - a[0]*a[1]*(1-a[2])/den,   1 -a[1]*a[1]*(1-a[2])/den, -a[1],
       a[0],                       a[1],                      1 +(-a[0]*a[0]+a[1]*a[1])*(1-a[2])/den);
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
    std::cout << "Triangulating hole as 2d polygon" << std::endl;
    list_hole<Polyhedron>(h);
    
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Segment_3 Segment_3;
    
    //typedef CGAL::Simple_cartesian<double> InexactKernel;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;
    typedef typename InexactKernel::Point_3 InexactPoint_3;
    // typedef typename SimpleCartesianKernel::Vector_3 InexactVector_3;
    typedef typename InexactKernel::Point_2 InexactPoint_2;
    typedef typename InexactKernel::Segment_2 InexactSegment_2;
    typedef typename InexactKernel::Segment_3 InexactSegment_3;

    typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_handle, InexactKernel> Vb;
    typedef CGAL::Constrained_triangulation_face_base_2<InexactKernel> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
    typedef CGAL::No_intersection_tag Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<InexactKernel, TDS, Itag> CDT;
    
    //dolfin_assert(vertices.size() > 3);
    //dolfin_assert(new_edges.empty());

    // Compute the best fitting plane of the points of the hole
    InexactPlane_3 fitting_plane;
    // double fit_quality;
    {
      std::vector<InexactPoint_3> boundary;
      Halfedge_handle current = h;
      do
      {
        const Point_3& p = current->vertex()->point();
        boundary.push_back(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])));
        current = current->next();
      } while (current != h);

      const double fit_quality = CGAL::linear_least_squares_fitting_3(boundary.begin(),
                                                                      boundary.end(),
                                                                      fitting_plane,
                                                                      CGAL::Dimension_tag<0>());
      std::cout << "Plane quality: " << fit_quality << std::endl;
    }

    const auto rotation = rotate_to_xy(fitting_plane.orthogonal_vector());

    std::cout << "Rotate normal: " << rotation.transform(fitting_plane.orthogonal_vector()) << std::endl;
    dolfin_assert(dolfin::near(fitting_plane.orthogonal_vector().squared_length(), 1));

    CDT cdt;

    std::cout << "Projected polygon" << std::endl;
    std::cout << "Polygon";

    // Insert vertices into triangulation
    std::vector<typename CDT::Vertex_handle> v;
    Halfedge_handle current = h;
    do
    {
      const Point_3& p = current->vertex()->point();
      const InexactPoint_3 rotated = rotation.transform(InexactPoint_3(CGAL::to_double(p[0]),
                                                                       CGAL::to_double(p[1]),
                                                                       CGAL::to_double(p[2])));
      std::cout << " " << rotated << ", ";
      
      /* const InexactPoint_2 p_2d = fitting_plane.to_2d(InexactPoint_3(CGAL::to_double(p[0]), */
      /*                                                                CGAL::to_double(p[1]),  */
      /*                                                                CGAL::to_double(p[2]))); */
      const InexactPoint_2 p_2d(rotated[0], rotated[1]);

      // std::cout << " " << p_2d << ", ";
      
      v.push_back(cdt.insert(p_2d));
      v.back()->info() = current->vertex();
      
      current = current->next();
    } while (current != h);

    std::cout << std::endl;

    std::cout << "Size of points: " << v.size() << std::endl;

    // Check if any of the edges intersect (before actually adding the
    // constrained edges to the triangulation
    if (check_for_intersections)
    {
      for (std::size_t i = 0; i < v.size()-1; i++)
      {
        const Point_3& a = v[i]->info()->point(), b = v[+1]->info()->point();
        const InexactSegment_2 s(v[i]->point(), v[i+1]->point());
        const Segment_3 original(a, b);

        const InexactSegment_3 s2(fitting_plane.projection(InexactPoint_3(CGAL::to_double(a[0]), CGAL::to_double(a[1]), CGAL::to_double(a[2]))),
                                  fitting_plane.projection(InexactPoint_3(CGAL::to_double(b[0]), CGAL::to_double(b[1]), CGAL::to_double(b[2]))));
        //std::cout << "Length (squared): " << s.squared_length() << ", projected: " << s2.squared_length() << ", original: " << original.squared_length() << ", relative: " << CGAL::to_double(s.squared_length()/original.squared_length()) << ", relative (proj): " << CGAL::to_double(s2.squared_length()/original.squared_length()) << std::endl;
      /* if (s.squared_length() < TOLERANCE) */
      /*   return false; */

        /*   // std::cout << "Outer " << i << ": " << s << ", length: " << s.squared_length() << std::endl; */

        for (std::size_t j = i+1; j < v.size(); j++)
        {
          InexactSegment_2 s2(v[j]->point(), v[(j+1)%v.size()]->point());
          // std::cout << "  Inner " << j << ": " << s2 << ", length: " << s2.squared_length() << std::endl;
          /* if (s2.squared_length() < TOLERANCE) */
          /*   return false; */

          /*     // If neighbor segments, check angle */
          /*     if (j == i+1 || (j+1)%v.size() == i) */
          /*     { */

          /*       const double cos_angle = (s.to_vector()*s2.to_vector())/std::sqrt(s.squared_length()*s2.squared_length()); */

          /*       /\* std::cout << "Angle: "  *\/ */
          /*       /\*           << (cos_angle+1) *\/ */
          /*       /\*           << std::endl; *\/ */
          /*       /\* std::cout << "\"Segment " << s.source() << ", " << s.target() << "\" \"Segment "  *\/ */
          /*       /\*           << s2.source() << ", " << s2.target() << "\"" << std::endl; *\/ */
          
          /*       if (cos_angle+1 < TOLERANCE) */
          /*         return false; */
          /*       else */
          /*         continue; */
          /*     } */

          /*     // std::cout << "  Intersecting " << j << ": Segment " << s2 << std::endl; */
          
          const auto intersection = CGAL::intersection(s, s2);
          
          if (intersection)
          {
            if (const InexactPoint_2* intersection_point = boost::get<InexactPoint_2>(&*intersection))
            {
              if (j != i+1 && i != (j+1)%v.size())
              {
                std::cout << "Non-neighbors (" << i << ", " << j << ")/" << v.size() << " intersect in single point" << std::endl;
                
                return false;
              }
            }
            else if (const InexactSegment_2* intersection_segment = boost::get<InexactSegment_2>(&*intersection))
            {
              std::cout << "Intersects in segment" << std::endl;
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
    
    for (std::size_t i = 0; i < v.size(); i++)
    {
      cdt.insert_constraint(v[i], v[(i+1)%v.size()]);
    }

    std::cout << "Done triangulating" << std::endl;
    std::cout << "Num vertices: " << cdt.number_of_vertices() << std::endl;
    std::set<std::pair<Vertex_handle, Vertex_handle> > edges_inside;
    {
      // Find an edge inside the polygon, starting from the infinite vertex
      typename CDT::Face_handle inside;
      {
        std::set<typename CDT::Face_handle> visited;
        std::deque<typename CDT::Face_handle> queue;
        queue.push_back(cdt.infinite_face());
        while (true)
        {
          typename CDT::Face_handle current = queue.front();
          queue.pop_front();
          if (visited.count(current) > 0)
            continue;
          visited.insert(current);

          bool found = false;
          for (std::size_t i = 0; i < 3; i++)
          {
            if (cdt.is_constrained(typename CDT::Edge(current, i)))
            {
              inside = current->neighbor(i);
              found = true;
              break;
            }
            else
            {
              queue.push_back(current->neighbor(i));
            }
          }
          
          if (found)
            break;
        }
      }

      // Collect faces inside the polygon
      std::set<typename CDT::Face_handle> visited;
      std::deque<typename CDT::Face_handle> queue;
      
      queue.push_back(inside);
      while (!queue.empty())
      {
        typename CDT::Face_handle current = queue.front();
        queue.pop_front();
        if (visited.count(current))
          continue;

        visited.insert(current);

        for (std::size_t i = 0; i < 3; i++)
        {
          if (!cdt.is_constrained(typename CDT::Edge(current, i)))
          {
            queue.push_back(current->neighbor(i));
          }
        }
      }
      
      // Write out to off file
      std::ofstream outfile("triangulation.off");
      outfile << "OFF\n" << cdt.number_of_vertices() << " " << visited.size() << " 0" << std::endl;
      std::size_t counter = 0;
      std::map<typename CDT::Vertex_handle, std::size_t> vertex_map;
      for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); vit++)
      {
        outfile << vit->point()[0] << " " << vit->point()[1] << " 0\n";
        vertex_map[vit] = counter;
        counter++;
      }
      
      for (auto fit = visited.begin(); fit != visited.end(); fit++)
      {
        typename CDT::Face_handle current = *fit;
        outfile << "3 " << vertex_map[current->vertex(0)]
                << " "  << vertex_map[current->vertex(1)]
                << " "  << vertex_map[current->vertex(2)] << "\n";

        for (std::size_t i = 0; i < 3; i++)
        {
          typename CDT::Edge e(current, i);
          if (!cdt.is_constrained(e))
          {
            const Vertex_handle a = e.first->vertex(cdt.cw(e.second))->info();
            const Vertex_handle b = e.first->vertex(cdt.ccw(e.second))->info();
            // We need to make sure the same edge (only opposite direction has been inserted already)
            if (edges_inside.count(std::make_pair(b, a)) == 0)
              edges_inside.insert(std::make_pair(a,b));
          }
        }
      }
    }
    
    dolfin_assert(cdt.is_valid());

    std::cout << "Edges to be added: " << edges_inside.size() << std::endl;

    // First make the entire hole a facet
    P.fill_hole(h);

    // Then add the edges as computed by the triangulation class
    std::size_t count = 0;
    for (auto eit = edges_inside.begin(); eit != edges_inside.end(); eit++)
    {
      std::cout << "Adding edge: Segment " << eit->first->point() << ", " << eit->second->point() << std::endl;
      count++;

      insert_edge(P, eit->first, eit->second);
    }

    std::cout << "Added " << count << " edges by triangulation" << std::endl;

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
  // Template specialization to allow triangulation of exact polyhedrons.
  // This uses an approximate normal (with sqrt computed in double precision)

  /* typedef CGAL::Polyhedron_3<_K>::Facet _K_Facet; */

  /* template <> */
  /* typename _K::Vector_3 compute_facet_normal<_K, _K_Facet>(const _K_Facet& f) */
  /* { */
  /*   typedef typename _K::Point_3 Point; */
  /*   typedef typename _K::Vector_3 Vector; */
  /*   typedef typename _K_Facet::Halfedge_around_facet_const_circulator HF_circulator; */
  /*   Vector normal = CGAL::NULL_VECTOR; */
  /*   /\* HF_circulator he = f.facet_begin(); *\/ */
  /*   /\* HF_circulator end = he; *\/ */
  /*   /\* CGAL_For_all(he,end) *\/ */
  /*   /\* { *\/ */
  /*   /\*   const Point& prev = he->prev()->vertex()->point(); *\/ */
  /*   /\*   const Point& curr = he->vertex()->point(); *\/ */
  /*   /\*   const Point& next = he->next()->vertex()->point(); *\/ */
  /*   /\*   Vector n = CGAL::cross_product(next-curr,prev-curr); *\/ */
  /*   /\*   normal = normal + n; *\/ */
  /*   /\* } *\/ */
  /*   /\* return normal / std::sqrt(normal * normal); *\/ */
  /*   return normal; */
  /* } */
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void list_hole(typename Polyhedron::Halfedge_handle h)
  {
    std::size_t counter = 0;
    std::cout << "Polygon";

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
        std::cout << " " << current->vertex()->point() << ",";

        current = current->next();
      } while(current != h);
      // }
    std::cout << std::endl;

    std::cout << " size: " << counter << std::endl;
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
  template<typename Triangle_3>
  static bool triangles_intersect(const Triangle_3& t1, const Triangle_3& t2)
  {
    typedef typename CGAL::Kernel_traits<Triangle_3>::Kernel::Point_3 Point_3;
    typedef typename CGAL::Kernel_traits<Triangle_3>::Kernel::Segment_3 Segment_3;
    
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
    static bool triangulate_center_vertex(Polyhedron& P,
                                          typename Polyhedron::Halfedge_handle h)
  {
    std::cout << "Triangulating with center vertex" << std::endl;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    
    // Compute centroid
    std::cout << "Checking if center vertex can be introduced safely" << std::endl;
    Point_3 centroid = h->vertex()->point();
    {
      Halfedge_handle current = h->next();
      std::size_t counter = 1;
      do
      {
        centroid = centroid + (current->vertex()->point()-CGAL::ORIGIN);
        counter++;
        current = current->next();
      } while (current != h);

      centroid = CGAL::ORIGIN + (centroid-CGAL::ORIGIN)/counter;
    }

    std::vector<Triangle_3> triangles;
    Halfedge_handle current = h;
    do
    {
      triangles.push_back(get_facet_triangle<Polyhedron>(current->opposite()));
      triangles.push_back(Triangle_3(centroid,
                                     current->vertex()->point(),
                                     current->next()->vertex()->point()));
      current = current->next();
    } while (current != h);

    if (triangle_set_intersects(triangles))
    {
      std::cout << "Couldn't triangulate by center vertex" << std::endl;
      return false;
    }

    std::cout << "adding center vertex" << std::endl;
    Halfedge_handle new_facet = P.fill_hole(h);
    Halfedge_handle center_vertex = P.create_center_vertex(new_facet);
    center_vertex->vertex()->point() = centroid;

    return true;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void fill_quad_hole(Polyhedron& P,
                             typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;

    dolfin_assert(h->next()->next()->next()->next() == h);
    
    P.fill_hole(h);

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
  static void fill_5hole(Polyhedron& P,
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

    P.fill_hole(h);
    P.split_facet(best_triangle->next(),
                  best_triangle->next()->next()->next());
    P.split_facet(best_triangle->next()->next(), best_triangle);
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void close_hole_agressive(Polyhedron& P,
                                   typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    // typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    
    std::vector<Point_3> points;
    typename Polyhedron::Halfedge_handle current = h;

    bool coplanar = true;
    do
    {
      if (!CGAL::coplanar(current->vertex()->point(),
                          current->next()->vertex()->point(),
                          current->next()->next()->vertex()->point(),
                          current->next()->next()->next()->vertex()->point()))
        coplanar = false;


      points.push_back(current->vertex()->point());
      current = current->next();
    } while (current != h);

    std::cout << "Size of hole: " << points.size() << std::endl;
    list_hole<Polyhedron>(h);

    dolfin_assert(points.size() > 2);
    if (points.size() == 3)
    {
      P.fill_hole(h);
      dolfin_assert(!triangles_intersect(get_facet_triangle<Polyhedron>(h),
                                         get_facet_triangle<Polyhedron>(h->opposite())));
      dolfin_assert(!triangles_intersect(get_facet_triangle<Polyhedron>(h),
                                         get_facet_triangle<Polyhedron>(h->next()->opposite())));
      dolfin_assert(!triangles_intersect(get_facet_triangle<Polyhedron>(h),
                                         get_facet_triangle<Polyhedron>(h->prev()->opposite())));
    }
    else if (points.size() == 4)
    {
      fill_quad_hole(P, h);
    }
    else if (points.size() == 5)
    {
      fill_5hole(P, h);
    }
    else
    {
      if (coplanar)
      {
        std::cout << "Points are coplanar" << std::endl;
        triangulate_polygon_3d(P, h, true);
      }
      else
      {
        // Try a center vertex
        if (!triangulate_center_vertex(P, h))
        {
          subdivide_hole(P, h);
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
  static bool subdivide_hole(Polyhedron& P, typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;
    // typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Segment_3 Segment_3;

    // Search for triangles that divide the hole, such that the dividing triangle
    // does not intersect triangles next to the hole

    // Store all triangles around the hole
    std::vector<Triangle_3> border_triangles;
    {
      Halfedge_handle current = h;
      do
      {
        border_triangles.push_back(get_facet_triangle<Polyhedron>(current->opposite()));
        current = current->next();
      } while (current != h);
    }

    std::cout << "Number of border triangles: " << border_triangles.size() << std::endl;
    std::cout << "Plane fit: " << get_plane_fit<Polyhedron>(h, h) << std::endl;
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
          std::cout << "Checking segment: Segment " << current_segment.source() << ", " << current_segment.target()
                    << " length: " << current_segment.squared_length() << std::endl;
        
          // Check that this does not introduce an intersection
          if (!segment_intersects_triangle_set(current_segment, border_triangles))
          {
            //const double side1_quality = get_plane_fit<Polyhedron>(current_outer, current_inner);
            //const double side2_quality = get_plane_fit<Polyhedron>(current_inner, current_outer);
            //const double candidate_quality = std::min(side1_quality, side2_quality);
            const double candidate_quality = evaluate_hole_subdivision<Polyhedron>(current_inner, current_outer);
            //std::cout << "Plane qualities: " << side1_quality << ", " << side2_quality << std::endl;
            std::cout << "Quality: " << candidate_quality << std::endl;

            if (candidate_quality > best_quality)
            {
              best_outer = current_outer;
              best_inner = current_inner;
              best_quality = candidate_quality;
              std::cout << "Currently best" << std::endl;
            }
            //{ int tmp; std::cin >> tmp; }
          }
        }
        
        current_inner = current_inner->next();
      } while (current_inner != inner_end);
      
      current_outer = current_outer->next();
    } while (current_outer != outer_end);

    dolfin_assert(best_outer != Halfedge_handle());
    dolfin_assert(best_inner != Halfedge_handle());

    std::cout << "Found best subdivision: " << std::endl;
    std::cout << "Segment " << best_outer->vertex()->point()
              << ", " << best_inner->vertex()->point() << std::endl;
    std::cout << "Quality: " << best_quality << std::endl;

    Halfedge_handle v1, v2;
    best_quality = -1000;

    {
      // Candidate 1: best_inner, best_inner->next(), best_outer
      std::cout << "Candidate 1" << std::endl;
      // Check the four candidates where the chosen segment is an edge
      if (best_inner->next()->next() != best_outer)
      {
        const double candidate_quality= std::min(get_plane_fit<Polyhedron>(best_outer, best_inner),
                                                 get_plane_fit<Polyhedron>(best_inner->next(), best_outer));

        if (candidate_quality > best_quality)
        {
          v1 = best_outer;
          v2 = best_inner;
          best_quality = candidate_quality;
        }
      }
    }

    if (best_outer->next()->next() != best_inner)
    {
      // Candidate 2: best_inner, best_outer, best_outer->next()
      double candidate_quality = std::min(get_plane_fit<Polyhedron>(best_outer->next(), best_inner),
                                          get_plane_fit<Polyhedron>(best_inner, best_outer));

      if (candidate_quality > best_quality)
      {
        std::cout << "Candidate 2" << std::endl;
        v1 = best_inner;
        v2 = best_outer;
        best_quality = candidate_quality;
      }
    }

    if (best_outer->next()->next() != best_inner)
    {
      // Candidate 3: best_inner->prev(), best_inner, best_outer
      const double candidate_quality = std::min(get_plane_fit<Polyhedron>(best_outer, best_inner->prev()),
                                                get_plane_fit<Polyhedron>(best_inner, best_outer));

      if (candidate_quality > best_quality)
      {
        std::cout << "Candidate 3" << std::endl;
        v1 = best_outer;
        v2 = best_inner->prev();
        best_quality = candidate_quality;
      }
    }

    if (best_inner->next()->next() != best_outer)
    {
      // Candidate 4: best_inner, best_outer->prev(), best_outer
      const double candidate_quality = std::min(get_plane_fit<Polyhedron>(best_outer, best_inner),
                                                get_plane_fit<Polyhedron>(best_inner, best_outer->prev()));

      if (candidate_quality > best_quality)
      {
        std::cout << "Candidate 4" << std::endl;
        v1 = best_inner;
        v2 = best_outer->prev();
      }
    }

    dolfin_assert(v1 != Halfedge_handle());
    dolfin_assert(v2 != Halfedge_handle());

    // Divide hole by chosen triangle (v1, v2, v2->next())
    Halfedge_handle f = P.fill_hole(h);
    dolfin_assert(v1->facet() == v2->facet());

    Halfedge_handle hole1 = P.split_facet(v1, v2->next());
    Halfedge_handle hole2 = P.split_facet(v2, hole1->opposite());
    P.make_hole(hole1);
    P.make_hole(hole2);
    
    return true;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void close_hole(Polyhedron& P, typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Traits Polyhedron_traits;
    // typedef typename Polyhedron_traits::Vector_3 Vector_3;
    typedef typename Polyhedron_traits::Triangle_3 Triangle_3;
    //typedef typename Polyhedron_traits::Plane_3 Plane_3;
    typedef typename Polyhedron_traits::Point_3 Point_3;

    //typedef CGAL::Simple_cartesian<double> InexactKernel
    typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;
    typedef typename InexactKernel::Point_3 InexactPoint_3;
    
    while (h->next()->next()->next() != h)
    {
      std::cout << "Closing hole" << std::endl;
      list_hole<Polyhedron>(h);
      
      // double max_cos_theta = std::numeric_limits<double>::lowest();
      // Halfedge_handle max_cos_theta_edge;

      // double max_cos_theta_intersections = std::numeric_limits<double>::lowest();
      // Halfedge_handle max_cos_theta_edge_intersections;

      //double max_facet_cos_angle = std::numeric_limits<double>::lowest();
      //Halfedge_handle max_facet_cos_angle_edge;

      double max_distance = std::numeric_limits<double>::lowest();
      Halfedge_handle max_distance_edge;
      
      //Halfedge_handle min_facet_angle_edge = h;
      Halfedge_handle start = h;

      std::vector<InexactPoint_3> boundary_points;
      do
      {
        const Point_3& p = h->vertex()->point();
        boundary_points.push_back(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])));
        h = h->next();
      } while (start != h);

      InexactPlane_3 fitting_plane;
      CGAL::linear_least_squares_fitting_3(boundary_points.begin(),
                                           boundary_points.end(),
                                           fitting_plane,
                                           CGAL::Dimension_tag<0>());
               
      do
      {
        dolfin_assert(h->is_border());
        // std::cout << "Vertex degree: " << h->vertex()->vertex_degree() << std::endl;
        
        if (h->vertex()->vertex_degree() > 2) 
        {
          /* Vector_3 v1(h->vertex()->point(), */
          /*             h->next()->vertex()->point()); */
          /* Vector_3 v2(h->vertex()->point(), */
          /*             h->prev()->vertex()->point()); */

          /* const double cos_theta = CGAL::to_double((v1*v2)/std::sqrt(CGAL::to_double(v1.squared_length()*v2.squared_length()))); */
          /* // std::cout << "Cos theta: " << cos_theta << ", " << max_cos_theta << std::endl; */
          /* if (cos_theta > max_cos_theta_intersections) */
          /* { */
          /*   max_cos_theta_intersections = cos_theta; */
          /*   max_cos_theta_edge_intersections = h; */
          /* } */

        /* Triangle_3 current(h->prev()->vertex()->point(), */
        /*                    h->vertex()->point(), */
        /*                    h->next()->vertex()->point()); */

        /* const double min_cos_angle = std::min(get_triangle_cos_angle<Polyhedron>(current, get_facet_triangle<Polyhedron>(h->opposite())), */
        /*                                       get_triangle_cos_angle<Polyhedron>(current, get_facet_triangle<Polyhedron>(h->next()->opposite()))); */


          const Point_3& p = h->vertex()->point();
          const InexactPoint_3 v(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2]));
          const double distance = (v-fitting_plane.projection(v)).squared_length();
          
          //if (min_cos_angle > max_facet_cos_angle &&
          if (distance > max_distance &&
            !facets_intersect<Polyhedron>(h->prev()->vertex(),
                                          h->vertex(),
                                          h->next()->vertex()) &&
            !facets_intersect<Polyhedron>(h->vertex(),
                                          h->next()->vertex(),
                                          h->prev()->vertex()) &&
            !facets_intersect<Polyhedron>(h->next()->vertex(),
                                          h->prev()->vertex(),
                                          h->vertex()))
          {
            // max_facet_cos_angle = min_cos_angle;
            // max_facet_cos_angle_edge = h;
            max_distance = distance;
            max_distance_edge = h;
          }
        }

        h = h->next();
      } while (h != start);
      
      //if (max_facet_cos_angle_edge == Halfedge_handle())
      if (max_distance_edge == Halfedge_handle())
      {
        std::cout << "Couldn't find vertex candidate without introducing self intersections" << std::endl;

        // FIXME: Pick the vertex more cleverly
        remove_vertex<Polyhedron>(P, h->vertex());
      }
      else
      {

        std::cout << "  Add facet to border, angle: " << max_distance << std::endl;

        /* Triangle_3 t(max_facet_cos_angle_edge->prev()->vertex()->point(), */
        /*              max_facet_cos_angle_edge->vertex()->point(), */
        /*              max_facet_cos_angle_edge->next()->vertex()->point()); */

        Triangle_3 t(max_distance_edge->prev()->vertex()->point(),
                     max_distance_edge->vertex()->point(),
                     max_distance_edge->next()->vertex()->point());

        
        std::cout << "Candidate: Triangle " << t[0] << ", " << t[1] << ", " << t[2] << std::endl;
        const bool intersects_simple = facets_intersect<Polyhedron>(max_distance_edge->prev()->vertex(),
                                                                    max_distance_edge->vertex(),
                                                                    max_distance_edge->next()->vertex()) ||
                                       facets_intersect<Polyhedron>(max_distance_edge->vertex(),
                                                                    max_distance_edge->next()->vertex(),
                                                                    max_distance_edge->prev()->vertex()) ||
                                       facets_intersect<Polyhedron>(max_distance_edge->next()->vertex(),
                                                                    max_distance_edge->prev()->vertex(),
                                                                    max_distance_edge->vertex());


        // Set h to a halfedge not affected by the ear clipping to ensure h is still
        // a border halfedge in the next iteration
        h = max_distance_edge->next()->next();

        print_triangle<Polyhedron>(max_distance_edge);
        dolfin_assert(max_distance_edge->prev()->is_border());
        dolfin_assert(max_distance_edge->next()->is_border());
        P.add_facet_to_border(max_distance_edge->prev(), max_distance_edge->next());

        std::vector<std::pair<Facet_handle, Facet_handle> > intersections;
        // CGAL::self_intersect<Polyhedron_traits>(P, std::back_inserter(intersections));

        if (intersections.size() > 0)
        {
          std::cout << "Self intersects (" << intersections.size() << "), remove common vertex" << std::endl;
          list_self_intersections(P);
          
          dolfin_assert(intersects_simple);
          dolfin_assert(false);

          for (auto iit = intersections.begin(); iit != intersections.end(); iit++)
          {
            Vertex_handle v = get_common_vertex<Polyhedron>(iit->first, iit->second);

            if (v == Vertex_handle())
              continue;

            // Ensure h is not one of the halfedges that will be removed
            if (h->vertex() == v)
              h = h->next();

            if (h->prev()->vertex() == v)
              h = h->next();

            std::cout << "Removing vertex: " << v->point() << std::endl;
            remove_vertex<Polyhedron>(P, v);
            break;
          }
        }
        else
        {
          dolfin_assert(!intersects_simple);
        }
      }
    }

    dolfin_assert(h->next()->next()->next() == h);

    P.fill_hole(h);
    std::cout << "Done filling" << std::endl;

    dolfin_assert(P.is_pure_triangle());
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void close_hole2(Polyhedron& P, typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

    std::cout << "Close hole 2" << std::endl;
    
    // Fill the hole, so it becomes a facet
    P.fill_hole(h);
    
    std::deque<Halfedge_handle> queue;
    queue.push_back(h);

    while (!queue.empty())
    {      
      Halfedge_handle facet_start = queue.front();
      queue.pop_front();

      std::cout << "Pop from queue" << std::endl;
      list_hole<Polyhedron>(facet_start);
      
      if (facet_start->facet()->is_triangle())
        continue;

      double best_fit = -1;
      std::pair<Halfedge_handle, Halfedge_handle> best_fit_edge;

      Halfedge_handle current = facet_start;
      std::size_t counter = 0;
      do
      {
        std::cout << "Outer counter: " << counter << std::endl;
        std::size_t inner_counter = 0;
        Halfedge_handle current2 = current->next()->next();
        do
        {
          std::cout << "Inner counter: " << inner_counter << std::endl;
          std::cout << "Checking segment:" << std::endl;
          std::cout << "Segment " << current->vertex()->point() << ", " << current2->vertex()->point() << std::endl;
          
          if (segment_intersects_facets<Polyhedron>(current->vertex(),
                                                    current2->vertex()))
          {
            std::cout << "Segment intersects" << std::endl;
          }
          else
          {
            const double fit1 = get_plane_fit<Polyhedron>(current, current2);
            const double fit2 = get_plane_fit<Polyhedron>(current2, current);
          
            if (std::min(fit1, fit2) > best_fit)
            {
              best_fit = std::min(fit1, fit2);
              best_fit_edge = std::make_pair(current, current2);
            }
          }

          // { int tmp; std::cin >> tmp; }
          inner_counter++;
          current2 = current2->next();
        } while (current2 != facet_start && current2->next() != current);

        counter++;
        current = current->next();
      } while (current->next() != facet_start);

      std::cout << "Polygon (again):" << std::endl;
      list_hole<Polyhedron>(best_fit_edge.first);
      std::cout << "Chosen split: " << std::endl;
      std::cout << "Segment " << best_fit_edge.first->vertex()->point() << ", " << best_fit_edge.second->vertex()->point() << std::endl;

      Halfedge_handle n = P.split_facet(best_fit_edge.first, best_fit_edge.second);
      queue.push_back(n);
      queue.push_back(n->opposite());

      { int tmp; std::cin >> tmp; }
    }
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static double evaluate_heuristic(const Polyhedron& P,
                                   typename Polyhedron::Halfedge_handle h,
                                   double plane_fit)
  {
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;


    // const double distance_to_plane_weight = 1.0;
    const double planarity_weight = 1.0;
    const double dihedral_weight  = 1.0;
    const double ear_angle_weight = 1.0;

    // Compute the planarity of the points excluding the current point
    InexactPlane_3 p;
    const double plane_fit_quality = get_plane_fit<Polyhedron>(h->next(),
                                                               h->prev(),
                                                               &p);
    // Compute the maximum of the dihedral angle to the neighbors
    const Triangle_3 candidate_triangle(h->prev()->vertex()->point(),
                                        h->vertex()->point(),
                                        h->next()->vertex()->point());
    const double cos_dihedral = (std::min(get_triangle_cos_angle(candidate_triangle,
                                                                 get_facet_triangle<Polyhedron>(h->opposite())),
                                          get_triangle_cos_angle(candidate_triangle,
                                                                 get_facet_triangle<Polyhedron>(h->next()->opposite())))+1)/2;

    // Compute the angle of the cutted ear
    const Vector_3 v1(h->vertex()->point(),
                      h->prev()->vertex()->point());
    const Vector_3 v2(h->vertex()->point(),
                      h->next()->vertex()->point());

    const double cos_ear_angle = (CGAL::to_double((v1*v2)/std::sqrt(CGAL::to_double(v1.squared_length()*v2.squared_length())))+1)/2.0;
    const double ear_angle_quality = cos_ear_angle;

    std::cout << "Triangle " << candidate_triangle[0]
              << ", " << candidate_triangle[1] 
              << "," << candidate_triangle[2] << std::endl;

    std::cout << "Evaluate: planarity: " << (plane_fit_quality/plane_fit)
              << ", dihedral: " << cos_dihedral
              << ", ear angle: " << ear_angle_quality << std::endl;

    return planarity_weight*plane_fit_quality + dihedral_weight*cos_dihedral + ear_angle_quality*ear_angle_weight;
    
    /* return planarity_weight*plane_fit_quality/plane_fit + */
    /*   dihedral_weight*cos_dihedral + */
    /*   (1-std::exp(-ear_angle_weight*cos_ear_angle))*cos_dihedral; */
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void close_hole3(Polyhedron& P, typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    /* typedef typename Polyhedron::Facet_handle Facet_handle; */
    /* typedef typename Polyhedron::Vertex_handle Vertex_handle; */
    /* typedef typename Polyhedron::Traits Polyhedron_traits; */
    // typedef typename Polyhedron_traits::Vector_3 Vector_3;
    // typedef typename Polyhedron_traits::Triangle_3 Triangle_3;
    //typedef typename Polyhedron_traits::Plane_3 Plane_3;
    // typedef typename Polyhedron_traits::Point_3 Point_3;

    /* typedef CGAL::Simple_cartesian<double> SimpleCartesianKernel; */
    /* typedef typename SimpleCartesianKernel::Plane_3 InexactPlane_3; */
    /* typedef typename SimpleCartesianKernel::Point_3 InexactPoint_3; */

    std::cout << "Close hole 3" << std::endl;
    dolfin_assert(h->is_border());
    list_hole<Polyhedron>(h);

    Halfedge_handle current = h;
    do
    {
      if (current->prev()->vertex()->vertex_degree() < 3 &&
          current->next()->vertex()->vertex_degree() < 3)
      {
        std::cout << "Cutting ear between degree 2 vertices" << std::endl;
        std::cout << "Segment " << h->prev()->vertex()->point() << ", " << h->next()->vertex()->point() << std::endl;
        if (segment_intersects_facets<Polyhedron>(current->prev()->vertex(),
                                                  current->next()->vertex()))
          std::cout << "Intersects!" << std::endl;

        if (segment_intersects_facets<Polyhedron>(current->next()->vertex(),
                                                  current->prev()->vertex()))
          std::cout << "Intersects!" << std::endl;

        
        Halfedge_handle new_facet = P.add_facet_to_border(current->prev(), 
                                                          current->next());
        current = new_facet->opposite();
        h = current->prev();
        dolfin_assert(current->is_border());
        dolfin_assert(h->is_border());
        const bool intersects = CGAL::self_intersect<typename Polyhedron::Traits, Polyhedron>(P);
        dolfin_assert(!intersects);
        { int tmp; std::cout << "Pausing..."; std::cin >> tmp; }
      }
      current = current->next();
    } while (h != current);

    while (h->next()->next()->next() != h)
    {
      std::cout << "Iteration" << std::endl;
      list_hole<Polyhedron>(h);
      //{ int tmp; std::cout << "Pausing..."; std::cin >> tmp; }

      std::cout << "Closing hole" << std::endl;
      /* if (triangulate_polygon_3d(P, h)) */
      /* { */
      /*   std::cout << "  Triangulated as 2d polygon" << std::endl; */
      /* } */
      /* else */
      {
        //std::cout << "Plane quality: " << get_plane_fit<Polyhedron>(h->next(), h->prev()) << std::endl;
        const double plane_fit = get_plane_fit<Polyhedron>(h, h->prev());
        double best_candidate_heuristic = -1;
        Halfedge_handle best_candidate;

        Halfedge_handle current = h;
        do
        {
          if (current->vertex()->vertex_degree() > 2 && 
              !segment_intersects_facets<Polyhedron>(current->prev()->vertex(), current->next()->vertex()) &&
              !segment_intersects_facets<Polyhedron>(current->next()->vertex(), current->prev()->vertex()))
          {
            const double heuristic_quality = evaluate_heuristic<Polyhedron>(P, current, plane_fit);
            if (heuristic_quality > best_candidate_heuristic)
            {
              best_candidate_heuristic = heuristic_quality;
              best_candidate = current;
            }
          }
        
          current = current->next();
        } while (current != h);

        if (best_candidate == Halfedge_handle())
        {
          std::cout << "Couldn't find non-intersecting candidate" << std::endl;
          exit(1);
        }

        std::cout << "Chosen triangle: " << std::endl;
        print_triangle<Polyhedron>(best_candidate);
        const double q = evaluate_heuristic<Polyhedron>(P, best_candidate, plane_fit);
        std::cout << "Quality: " << q << std::endl;
        
        h = best_candidate->next()->next();

        {
          std::ofstream outfile("before_clipping.off");
          outfile << P;
        }

        dolfin_assert(best_candidate->prev()->is_border());
        dolfin_assert(best_candidate->next()->is_border());
        P.add_facet_to_border(best_candidate->prev(), best_candidate->next());
      }
      
      // { int tmp; std::cin >> tmp; }      
      // const bool intersects = CGAL::self_intersect<typename Polyhedron::Traits, Polyhedron>(P);
      // dolfin_assert(!intersects);
    }

    dolfin_assert(h->next()->next()->next() == h);
    P.fill_hole(h);
  }
  //-----------------------------------------------------------------------------
  /* template<typename Polyhedron> */
  /* static void close_hole_advancing_front(Polyhedron& P, typename Polyhedron::Halfedge_handle h) */
  /* { */
  /*   typedef typename Polyhedron::Halfedge_handle Halfedge_handle; */
  /*   typedef typename Polyhedron::Traits::Point_3 Point_3; */
  /*   typedef typename Polyhedron::Traits::Vector_3 Vector_3; */

  /*   const double threshold = cos(5*DOLFIN_PI/12.); */
    
    
  /*   while (h->next()->next()->next() != h) */
  /*   { */
  /*     dolfin_assert(h->is_border()); */
      
  /*     list_hole<Polyhedron>(h); */
      
  /*     // Find vertex minimizing angle */
  /*     double max_cos_theta = -2; */
  /*     Halfedge_handle max_cos_theta_edge; */
      
  /*     Halfedge_handle current = h; */
  /*     do */
  /*     { */
  /*       if (current->vertex()->vertex_degree() > 2) */
  /*       { */
  /*         const double cos_theta = get_edge_cos_angle<Polyhedron>(current); */
  /*         if (cos_theta > max_cos_theta) */
  /*         { */
  /*           max_cos_theta = cos_theta; */
  /*           max_cos_theta_edge = current; */
  /*         } */
  /*       } */
  /*       current = current->next(); */
  /*     } while (current != h); */

  /*     dolfin_assert(max_cos_theta_edge != Halfedge_handle()); */

  /*     h = max_cos_theta_edge->next()->next(); */
  /*     std::cout << "Minimizing vertex: " << max_cos_theta_edge->vertex()->point() << ", theta: " << max_cos_theta << std::endl; */

  /*     if (max_cos_theta > threshold) */
  /*     { */
  /*       // Just cut ear */
  /*       std::cout << "Chosen Segment " << max_cos_theta_edge->prev()->vertex()->point() << ", " << max_cos_theta_edge->next()->vertex()->point() << std::endl; */
  /*       P.add_facet_to_border(max_cos_theta_edge->prev(), max_cos_theta_edge->next()); */
  /*       //{ int tmp; std::cin >> tmp; } */
  /*     } */
  /*     else */
  /*     { */
  /*       // Add one new vertex */
  /*       std::cout << "Chosen Segment " << max_cos_theta_edge->prev()->vertex()->point() << ", " << max_cos_theta_edge->next()->vertex()->point() << std::endl; */
  /*       std::cout << "Vertex: " << max_cos_theta_edge->vertex()->point() << std::endl; */
  /*       std::cout << "add vertex" << std::endl; */
  /*       get_edge_cos_angle<Polyhedron>(max_cos_theta_edge); */

  /*       const Halfedge_handle e = max_cos_theta_edge; */

  /*       const Point_3& prev = e->prev()->vertex()->point(); */
  /*       const Point_3& current = e->vertex()->point(); */
  /*       const Point_3& next = e->next()->vertex()->point(); */

  /*       dolfin_assert(e->is_border()); */
  /*       const Halfedge_handle new_vertex = P.add_vertex_and_facet_to_border(e->prev(), e->next()); */
  /*       dolfin_assert(!e->is_border()); */

  /*       // Compute new point */


  /*       const double length = (std::sqrt(CGAL::to_double((prev-current).squared_length())) + */
  /*                              std::sqrt(CGAL::to_double((next-current).squared_length())))/2; */
  /*       Vector_3 v = prev-current + (next-prev)/2; */
  /*       v = v*length/std::sqrt(CGAL::to_double(v.squared_length())); */
  /*       std::cout << "new vertex: " << (current+v) << std::endl; */
  /*       new_vertex->vertex()->point() = current+v; */

  /*       dolfin_assert(!P.is_pure_triangle()); */
  /*       P.split_facet(new_vertex, e); */
  /*       dolfin_assert(P.is_pure_triangle()); */
        
  /*       std::ofstream outfile("advancing_front.off"); */
  /*       outfile << P; */
  /*       //exit(1); */
  /*     } */
  /*     dolfin_assert(h->is_border()); */
  /*   }     */
  /* } */
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void close_holes(Polyhedron& P)
  {
    //typedef typename Polyhedron::Traits Kernel;
    std::cout << "Closing holes" << std::endl;
    std::cout << "Min vertex degree: " << min_vertex_degree(P) << std::endl;
    P.normalize_border();

    while (!P.is_closed())
    {
      {
        std::ofstream outfile("not-closed.off");
        outfile << P;
      }

      // close_hole_advancing_front(P, P.border_halfedges_begin()->opposite());
      close_hole3(P, P.border_halfedges_begin()->opposite());
      //const bool success = triangulate_polygon_3d(P, P.border_halfedges_begin()->opposite());

      dolfin_assert(P.is_pure_triangle());

      {
        std::ofstream outfile("closed.off");
        outfile << P;
      }

      /* if (!success) */
      /* { */
      /*   std::cout << "Couldn't triangulate by 2d projection" << std::endl; */
      /*   exit(1); */
                
      /* } */
      
      P.normalize_border();
      std::cout << "Closed hole" << std::endl;
      //const bool self_intersects = CGAL::self_intersect<Kernel, Polyhedron>(P);
      /* if (self_intersects) */
      /* { */
      /*   list_self_intersections(P); */
      /*   exit(1); */
      /* } */
      //dolfin_assert(!self_intersects);

    }

    dolfin_assert(P.is_closed());
  }
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
  static bool ear_is_cuttable(typename Polyhedron::Halfedge_handle h)
  {
    dolfin_assert(h->is_border());

    typename Polyhedron::Traits::Triangle_3 t(h->vertex()->point(),
                                              h->next()->vertex()->point(),
                                              h->prev()->vertex()->point());
    if (t.is_degenerate())
      return false;

    if (h->vertex()->vertex_degree() < 3)
    {
      std::cout << "Degree: " << h->vertex()->vertex_degree() << std::endl;
      dolfin_assert(h->next()->is_border());
      dolfin_assert(h->next()->vertex()->vertex_degree() > 2);
      dolfin_assert(h->prev()->is_border());
      if (h->prev()->vertex()->vertex_degree() < 3)
      {
        std::cout << "Prev degree: " << h->prev()->vertex()->vertex_degree() << std::endl;
        std::cout << "Facet degree: " << h->opposite()->facet_degree() << std::endl;
        std::cout << "Disconnected triangle: " << (h->next()->next()->next() == h ? "True" : "False") << std::endl;
      }
      dolfin_assert(h->prev()->vertex()->vertex_degree() > 2);
      return false;
    }
    

    return true;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void cut_ear_from_hole(Polyhedron& p, typename Polyhedron::Halfedge_handle h)
  {
    dolfin_assert(h->is_border());

    // count egdes along border
    std::size_t border_length = 0;
    typename Polyhedron::Halfedge_handle current = h;
    do 
    {
      border_length++;
      current = current->next();
    } while (current != h);

    dolfin_assert(border_length > 2);

    // Check if the remaining hole is triangular
    if (border_length == 3)
    {
      std::cout << "Filling hole" << std::endl;
      p.fill_hole(h);
    }
    else
    {

      typename Polyhedron::Halfedge_handle current = h;
      do
      {
        dolfin_assert(h->is_border());
        if (ear_is_cuttable<Polyhedron>(current))
        {
          typename Polyhedron::Halfedge_handle n = p.add_facet_to_border(current->prev(), current->next());

          break;
        }

        current = current->next();
        if (current->next() == h)
          dolfin::dolfin_error("Polyhedron_utils.h",
                               "triangulating hole",
                               "Couldn't triangulate by ear clipping");
      } while (current != h);
    }
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void cut_holes(Polyhedron& P)
  {
    dolfin_assert(P.is_pure_triangle());

    std::size_t global_count = 0;

    for (typename Polyhedron::Vertex_iterator vit = P.vertices_begin();
         vit != P.vertices_end(); vit++)
    {
      // Check if two consecutive edges around a vertex are both borders.
      // In that case, connect them
      std::size_t local_count = 0;
      typename Polyhedron::Halfedge_around_vertex_circulator start = vit->vertex_begin();
      bool prev_is_border = start->is_border_edge();
      start++;

      typename Polyhedron::Halfedge_around_vertex_circulator current = start;
      do
      {
        const bool current_is_border = current->is_border_edge();
        if (prev_is_border && current_is_border)
        {
          std::cout << "Two consecutive border edges" << std::endl;
          local_count++;
        }
        prev_is_border = current_is_border;
        current++;
      } while (current != start);
      if (local_count > 0)
        std::cout << "  Local: " << local_count << std::endl;
      global_count += local_count;
    }

    std::cout << "Global count: " << global_count << std::endl;
      int tmp;
    std::cin >> tmp;


    // Close holes by ear cutting
    std::cout << "Closing holes" << std::endl;
    P.normalize_border();

    while (!P.is_closed())
    {
      std::cout << "Cut ear from hole" << std::endl;
      cut_ear_from_hole(P, P.border_halfedges_begin()->opposite());
      P.normalize_border();
      dolfin_assert(P.is_pure_triangle());
    }

    dolfin_assert(P.is_closed());
    std::cout << "Done closing polyhedron" << std::endl;
    std::cout << "Min vertex degree: " << min_vertex_degree(P) << std::endl;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void triangulate_polyhedron(Polyhedron& p)
  {
    std::cout << "Triangulating!" << std::endl;

    typedef typename Polyhedron::Traits Kernel;

    while (!p.is_pure_triangle())
    {
      for (typename Polyhedron::Facet_iterator fit = p.facets_begin(); fit != p.facets_end(); fit++)
      {
        //typename Polyhedron::Facet_handle facet = fit;
        if (!fit->is_triangle())
          std::cout << "Triangulating facet of degree" << fit->facet_degree() << std::endl;


        std::size_t degree = fit->facet_degree();
        while (!fit->is_triangle())
        {
          std::cout << " degree: "<< fit->facet_degree() << std::endl;
          // Find an ear so that the newly created facet does not intersect any
          // neighbors

          typename Polyhedron::Halfedge_handle start = fit->halfedge(); 
          typename Polyhedron::Halfedge_handle current = start;
          do 
          {
            if (ear_is_cuttable<Polyhedron>(current))
            {
              
              typename Polyhedron::Halfedge_handle h = current->prev();
              typename Polyhedron::Halfedge_handle g = current->next();
              dolfin_assert(h != g);
              dolfin_assert(h->vertex() != g->vertex());
              typename Polyhedron::Halfedge_handle n = p.split_facet(h, g);
              {
                std::cout << "new facet" << std::endl;
                dolfin_assert(g->next()->next()->next() == g);
                typename Kernel::Triangle_3 t(g->vertex()->point(), 
                                              g->next()->vertex()->point(),
                                              g->prev()->vertex()->point());
                std::cout << "Area: " << t.squared_area() << std::endl;
                if (t.squared_area() == 0)
                {
                  int tmp;
                  std::cin >> tmp;
                }
              }
              break;
            }
              
            current = current->next();

            // We should never need to go the entire round
            dolfin_assert(current->next() != start);
          } while (current != start);

          assert(fit->facet_degree() == degree-1);
          degree = fit->facet_degree();
        }

        dolfin_assert(fit->is_triangle());

        {
          typename Polyhedron::Halfedge_handle g = fit->halfedge();
          std::cout << "original facet" << std::endl;
          dolfin_assert(g->next()->next()->next() == g);
          typename Kernel::Triangle_3 t(g->vertex()->point(), 
                                        g->next()->vertex()->point(),
                                              g->prev()->vertex()->point());
          std::cout << "Area: " << t.squared_area() << std::endl;
          if (t.squared_area() < 1e-10)
          {
            int tmp;
            std::cin >> tmp;
          }
        }
      }
    }
    std::cout << "Done triangulating" << std::endl;
  }
  //-----------------------------------------------------------------------------
  template<typename Point>
  bool point_x_less_than(const Point& a, const Point& b)
  {
    return a.x() < b.x();
  }

  template<typename Traits>
  typename Traits::FT closest_pair(const std::vector<typename Traits::FT>& P, std::size_t a, std::size_t b)
  {
    while (b-a > 3)
    {
      
    }
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
    
    std::cout << "Removing vertex" << std::endl;

    Vertex_circulator h = v->vertex_begin();
    Vertex_circulator start = h;

    std::vector<Halfedge_handle> to_be_removed;

    do
    {
      if (!h->is_border())
        to_be_removed.push_back(h);

      h++;
    } while (h != start);


    std::cout << "Removing " << to_be_removed.size() << " halfedges" << std::endl;
    for (auto it = to_be_removed.begin(); it != to_be_removed.end(); it++)
    {
      P.erase_facet(*it);
    }

    std::cout << "  done removing vertex" << std::endl;
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

    if (intersections.size() > 0)
      std::cout << "Removing self intersections" << std::endl;

    while (intersections.size() > 0)
    {
      std::cout << "  Removing pair" << std::endl;


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
      
      
      std::cout << "To be removed 1: " << to_be_removed1.size() << std::endl;
      for (auto it = to_be_removed1.begin(); it != to_be_removed1.end(); it++)
      {
        P.erase_facet((*it)->halfedge());
      }

      std::cout << "To be removed 2: " << to_be_removed2.size() << std::endl;
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

      std::cout << "Found route from f1 to f2" << std::endl;
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

    {
      std::ofstream outfile("before_filtering.off");
      outfile << P;
    }
    
    std::cout << "Filter sharp features" << std::endl;

    const double cos_tolerance = std::cos(tolerance);
    std::cout << "tolerance: " << cos_tolerance << std::endl;
    Facet_iterator fit = P.facets_begin();
    for (int i = 0; i < start_facet; i++)
      fit++;

    std::deque<Facet_handle> queue;
    {
      Halfedge_handle h = fit->halfedge();
      std::cout << "Starting facet: Triangle " << h->vertex()->point() << ", ";
      h = h->next();
      std::cout << h->vertex()->point() << ", ";
      h = h->next();
      std::cout << h->vertex()->point() << std::endl;
    }

    std::set<Facet_handle> visited;
    std::set<Facet_handle> to_be_removed;
    for (Facet_iterator fit = P.facets_begin(); fit != P.facets_end(); fit++)
      to_be_removed.insert(fit);

    std::cout << "Number of facets: " << P.size_of_facets() << std::endl;

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

    std::cout << "Remove " << to_be_removed.size() << " facets" << std::endl;

    for (auto fit = to_be_removed.begin(); fit != to_be_removed.end(); fit++)
    {
      P.erase_facet( (*fit)->halfedge() );
    }

    {
      std::ofstream outfile("after_filtering.off");
      outfile << P;
    }

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
  CGAL_assertion(poly_b.is_valid());
}
}
#endif
