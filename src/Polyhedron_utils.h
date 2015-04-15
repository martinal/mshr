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

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Self_intersection_polyhedron_3.h>

#include <cmath>
#include <deque>

typedef CGAL::Exact_predicates_exact_constructions_kernel _K;

namespace std
{
  // TODO: This is a hack to allow triangulation of exact polyhedrons.  A proper
  // fix would be to provide a template specialization of compute_facet_normal
  // for epeck.
  inline _K::FT sqrt(_K::FT a)
  {
    return std::sqrt(CGAL::to_double(a));
  }
}

#include <CGAL/triangulate_polyhedron.h>

namespace tanganyika
{

class PolyhedronUtils
{
 public:

    template<typename Polyhedron>
    static void recursive_remove(std::set<typename Polyhedron::Vertex_const_handle>& s,
                          typename Polyhedron::Vertex_const_handle h)
  {
    typedef typename Polyhedron::Halfedge_around_vertex_const_circulator HV_const_circulator;
    typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
    typedef Polyhedron Polyhedron_type;

    const HV_const_circulator start = h->vertex_begin();
    HV_const_circulator current = start;
    do
    {
      Vertex_const_handle current_vertex = current->opposite()->vertex();
      assert(current_vertex != h);
      if (s.count(current_vertex))
      {
        s.erase(current_vertex);
        recursive_remove<Polyhedron_type>(s, current_vertex);
      }
      current++;
    } while (current != start);
  }
  //-----------------------------------------------------------------------------
  // Scans the vertices of the polyhedron the polyhedron and returns a
  // Polyhedron::Vertex_const_handle for each disconnected component.
  template <typename Polyhedron, typename OutputIterator>
  static void get_disconnected_components(const Polyhedron& p, OutputIterator it)
  {
    typedef Polyhedron Polyhedron_t;
    typedef typename Polyhedron_t::Vertex_const_handle Vertex_const_handle;

    // store all vertices in a set
    std::set<Vertex_const_handle> v;
    for (typename Polyhedron_t::Vertex_const_iterator vit = p.vertices_begin();
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
      recursive_remove<Polyhedron_t>(v, start);
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
    static double get_triangle_cos_angle(typename Polyhedron::Traits::Triangle_3 t1, 
                                         typename Polyhedron::Traits::Triangle_3 t2)
  {
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;
    const Vector_3 v1 = t1.supporting_plane().orthogonal_vector();
    const Vector_3 v2 = t2.supporting_plane().orthogonal_vector();

    return CGAL::to_double((v1*v2)/std::sqrt(CGAL::to_double(v1.squared_length()*v2.squared_length())));  
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
    dolfin_assert(h->is_border());

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

    if (counter < 250)
    {
      typename Polyhedron::Halfedge_handle current = h;
      do
      {
        std::cout << " " << current->vertex()->point() << ",";

        current = current->next();
      } while(current != h);
    }
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
  template<typename Polyhedron>
  static void close_hole(Polyhedron& P, typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Traits Polyhedron_traits;
    // typedef typename Polyhedron_traits::Vector_3 Vector_3;
    typedef typename Polyhedron_traits::Triangle_3 Triangle_3;

    while (h->next()->next()->next() != h)
    {
      std::cout << "Closing hole" << std::endl;
      list_hole<Polyhedron>(h);
      
      // double max_cos_theta = std::numeric_limits<double>::lowest();
      // Halfedge_handle max_cos_theta_edge;

      // double max_cos_theta_intersections = std::numeric_limits<double>::lowest();
      // Halfedge_handle max_cos_theta_edge_intersections;

      double max_facet_cos_angle = std::numeric_limits<double>::lowest();
      Halfedge_handle max_facet_cos_angle_edge;

      //Halfedge_handle min_facet_angle_edge = h;
      Halfedge_handle start = h;

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

        Triangle_3 current(h->prev()->vertex()->point(),
                           h->vertex()->point(),
                           h->next()->vertex()->point());

        const double min_cos_angle = std::min(get_triangle_cos_angle<Polyhedron>(current, get_facet_triangle<Polyhedron>(h->opposite())),
                                              get_triangle_cos_angle<Polyhedron>(current, get_facet_triangle<Polyhedron>(h->next()->opposite())));

          
        if (min_cos_angle > max_facet_cos_angle &&
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
            max_facet_cos_angle = min_cos_angle;
            max_facet_cos_angle_edge = h;
          }
        }


        h = h->next();
      } while (h != start);

      if (max_facet_cos_angle_edge == Halfedge_handle())
      {
        std::cout << "Couldn't find vertex candidate without introducing self intersections" << std::endl;

        // FIXME: Pick the vertex more cleverly
        remove_vertex<Polyhedron>(P, h->vertex());
      }
      else
      {

        std::cout << "  Add facet to border, angle: " << max_facet_cos_angle << std::endl;

        Triangle_3 t(max_facet_cos_angle_edge->prev()->vertex()->point(),
                     max_facet_cos_angle_edge->vertex()->point(),
                     max_facet_cos_angle_edge->next()->vertex()->point());

        std::cout << "Candidate: Triangle " << t[0] << ", " << t[1] << ", " << t[2] << std::endl;
        const bool intersects_simple = facets_intersect<Polyhedron>(max_facet_cos_angle_edge->prev()->vertex(),
                                                                    max_facet_cos_angle_edge->vertex(),
                                                                    max_facet_cos_angle_edge->next()->vertex()) ||
                                       facets_intersect<Polyhedron>(max_facet_cos_angle_edge->vertex(),
                                                                    max_facet_cos_angle_edge->next()->vertex(),
                                                                    max_facet_cos_angle_edge->prev()->vertex()) ||
                                       facets_intersect<Polyhedron>(max_facet_cos_angle_edge->next()->vertex(),
                                                                    max_facet_cos_angle_edge->prev()->vertex(),
                                                                    max_facet_cos_angle_edge->vertex());


        // Set h to a halfedge not affected by the ear clipping to ensure h is still
        // a border halfedge in the next iteration
        h = max_facet_cos_angle_edge->next()->next();

        print_triangle<Polyhedron>(max_facet_cos_angle_edge);
        dolfin_assert(max_facet_cos_angle_edge->prev()->is_border());
        dolfin_assert(max_facet_cos_angle_edge->next()->is_border());
        P.add_facet_to_border(max_facet_cos_angle_edge->prev(), max_facet_cos_angle_edge->next());

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
  static void close_holes(Polyhedron& P)
  {
    typedef typename Polyhedron::Traits Kernel;
    std::cout << "Closing holes" << std::endl;
    std::cout << "Min vertex degree: " << min_vertex_degree(P) << std::endl;
    P.normalize_border();

    while (!P.is_closed())
    {
      close_hole(P, P.border_halfedges_begin()->opposite());
      dolfin_assert(P.is_pure_triangle());
      const bool self_intersects = CGAL::self_intersect<Kernel, Polyhedron>(P);
      dolfin_assert(!self_intersects);
      P.normalize_border();
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
        if (h2->vertex() != h1->vertex())
          return true;

        h2 = h2->next();
      } while (h2 != start2);

      h1 = h1->next();
    } while (h1 != start1);
    
    return false;
  }
  //-----------------------------------------------------------------------------
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
      std::cout << "Intersection " << (facets_are_neighbors<Polyhedron>(iit->first, iit->second) ? "True" : "False") << std::endl;
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
};
}
#endif
