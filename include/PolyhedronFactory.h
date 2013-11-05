// Copyright (C) 2013 Benjamin Kehlet
//
// This file is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2013-10-03
// Last changed: 2013-10-03

#ifndef PolyhedronFactoruy_h
#define PolyhedronFactoruy_h

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <iostream>

//-----------------------------------------------------------------------------
// Convenience routine to make debugging easier. Remove before releasing.
template<typename Polyhedron>
static void add_facet(CGAL::Polyhedron_incremental_builder_3<typename Polyhedron::HalfedgeDS>& builder,
		      std::vector<int>& vertices, bool print=false)
{
  static int facet_no = 0;

  if (print)
  {
    std::cout << "Begin facet " << facet_no << std::endl;
    if (!vertices.size())
    {
      std::cout << "No vertices in facet!" << std::endl;
      return;
    }

    // Print vertices
    for (std::vector<int>::iterator it=vertices.begin(); it != vertices.end(); it++)
      std::cout << "Vertex: " << (*it) << std::endl;

    if (builder.test_facet(vertices.begin(), vertices.end()))
      std::cout << "Facet ok, size: " << vertices.size() << std::endl;
    else
      std::cout << "Facet not ok" << std::endl;
  }

  builder.begin_facet();
  for (std::vector<int>::iterator it=vertices.begin(); it != vertices.end(); it++)
    builder.add_vertex_to_facet(*it);
  builder.end_facet();

  if (print)
    std::cout << "End facet" << std::endl;
  facet_no++;
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
static void add_verte(CGAL::Polyhedron_incremental_builder_3<typename Polyhedron::HalfedgeDS>& builder,
                      const typename Polyhedron::Traits::Point_3& point, 
                      bool print=false)
{
  static int vertex_no = 0;
  if (print)
    std::cout << "Adding vertex " << vertex_no << " at " << point << std::endl;

  builder.add_vertex(point);
  vertex_no++;
}
//-----------------------------------------------------------------------------
/* // Sphere */
/* //----------------------------------------------------------------------------- */
/* class Build_sphere : public CGAL::Modifier_base<csg::Exact_HalfedgeDS> */
/* { */
/*  public: */
/*   Build_sphere(const Sphere& sphere) : _sphere(sphere) {} */

/*   void operator()( csg::Exact_HalfedgeDS& hds ) */
/*   { */
/*     const std::size_t num_slices = _sphere._slices; */
/*     const std::size_t num_sectors = _sphere._slices*2 + 1; */

/*     const dolfin::Point top = _sphere.c + Point(_sphere.r, 0, 0); */
/*     const dolfin::Point bottom = _sphere.c - Point(_sphere.r, 0, 0); */
/*     const dolfin::Point axis = Point(1, 0, 0); */

/*     const int num_vertices = num_slices*num_sectors+2; */
/*     const int num_facets = num_sectors*2*num_slices; */

/*     CGAL::Polyhedron_incremental_builder_3<csg::Exact_HalfedgeDS> builder( hds, true ); */

/*     builder.begin_surface(num_vertices, num_facets); */

/*     const Point slice_rotation_axis(0, 1, 0); */

/*     for (std::size_t i = 0; i < num_slices; i++) */
/*     { */
/*       const Point sliced = axis.rotate(slice_rotation_axis, (i+1)*DOLFIN_PI/(num_slices+1)); */
/*       for (std::size_t j = 0; j < num_sectors; j++) */
/*       { */
/*         const Point direction = sliced.rotate(axis, j*2.0*DOLFIN_PI/num_sectors); */
/*         const Point v = _sphere.c + direction*_sphere.r; */
/*         add_vertex(builder, csg::Exact_Point_3 (v.x(), v.y(), v.z())); */
/*       } */
/*     } */

/*     // Add bottom has index num_vertices-1, top has index num_vertices-2 */
/*     add_vertex(builder, csg::Exact_Point_3(top.x(), top.y(), top.z())); */
/*     add_vertex(builder, csg::Exact_Point_3(bottom.x(), bottom.y(), bottom.z())); */

/*     // Add the side facets */
/*     for (std::size_t i = 0; i < num_slices-1; i++) */
/*     { */
/*       for (std::size_t j = 0; j < num_sectors; j++) */
/*       { */
/*         const std::size_t offset1 = i*num_sectors; */
/*         const std::size_t offset2 = (i+1)*num_sectors; */

/*         { */
/*           std::vector<int> f; */
/*           f.push_back(offset1 + j); */
/*           f.push_back(offset1 + (j+1)%num_sectors); */
/*           f.push_back(offset2 + j); */
/*           add_facet(builder, f); */
/*         } */

/*         { */
/*           std::vector<int> f; */
/*           f.push_back(offset2 + (j+1)%num_sectors); */
/*           f.push_back(offset2 + j); */
/*           f.push_back(offset1 + (j+1)%num_sectors); */
/*           add_facet(builder, f); */
/*         } */
/*       } */
/*     } */

/*     // Add the top and bottom facets */
/*     const std::size_t top_offset = num_sectors*(num_slices-1); */
/*     for (std::size_t i = 0; i < num_sectors; i++) */
/*     { */
/*       { */
/*         // Bottom facet */
/*         std::vector<int> f; */
/*         f.push_back( num_vertices-2 ); */
/*         f.push_back( (i+1)%num_sectors ); */
/*         f.push_back(i); */
/*         add_facet(builder, f); */
/*       } */

/*       { */
/*         // Top facet */
/*         std::vector<int> f; */
/*         f.push_back( num_vertices-1 ); */
/*         f.push_back( top_offset + (i%num_sectors) ); */
/*         f.push_back( top_offset + (i+1)%num_sectors ); */
/*         add_facet(builder, f); */
/*       } */
/*     } */
/*     builder.end_surface(); */
/*   } */

/*   private: */
/*   const Sphere& _sphere; */
/* }; */
/* //----------------------------------------------------------------------------- */
/* static void make_sphere(const Sphere* s, csg::Exact_Polyhedron_3& P) */
/* { */
/*   Build_sphere builder(*s); */
/*   P.delegate(builder); */
/*   dolfin_assert(P.is_valid()); */
/*   dolfin_assert(P.is_closed()); */
/* } */
/* //----------------------------------------------------------------------------- */
template <typename Polyhedron>
class Build_box : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
{
 public:
  typedef Polyhedron P;
  typedef typename P::Traits::Point_3 Point_3;

  Build_box(double x0, double x1, double x2, double y0, double y1, double y2) 
    : _x0(x0), _x1(x1), _x2(x2), _y0(y0), _y1(y1), _y2(y2) {}

  void operator()( typename P::HalfedgeDS& hds )
  {
    CGAL::Polyhedron_incremental_builder_3<typename P::HalfedgeDS> builder(hds, true);

    builder.begin_surface(8, 12);

    const double x0 = std::min(_x0, _y0);
    const double y0 = std::max(_x0, _y0);

    const double x1 = std::min(_x1, _y1);
    const double y1 = std::max(_x1, _y1);

    const double x2 = std::min(_x2, _y2);
    const double y2 = std::max(_x2, _y2);

    builder.add_vertex(Point_3(y0, x1, x2));
    builder.add_vertex(Point_3(x0, x1, y2));
    builder.add_vertex(Point_3(x0, x1, x2));
    builder.add_vertex(Point_3(x0, y1, x2));
    builder.add_vertex(Point_3(y0, x1, y2));
    builder.add_vertex(Point_3(x0, y1, y2));
    builder.add_vertex(Point_3(y0, y1, x2));
    builder.add_vertex(Point_3(y0, y1, y2));

    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(1); */
    /*   f.push_back(2); */
    /*   f.push_back(3); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(1);
    builder.add_vertex_to_facet(2);
    builder.add_vertex_to_facet(3);
    builder.end_facet();

    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(1); */
    /*   f.push_back(3); */
    /*   f.push_back(5); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(1);
    builder.add_vertex_to_facet(3);
    builder.add_vertex_to_facet(5);
    builder.end_facet();


    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(1); */
    /*   f.push_back(5); */
    /*   f.push_back(4); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(1);
    builder.add_vertex_to_facet(5);
    builder.add_vertex_to_facet(4);
    builder.end_facet();

    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(4); */
    /*   f.push_back(5); */
    /*   f.push_back(7); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(4);
    builder.add_vertex_to_facet(5);
    builder.add_vertex_to_facet(7);
    builder.end_facet();

    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(4); */
    /*   f.push_back(7); */
    /*   f.push_back(0); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(4);
    builder.add_vertex_to_facet(7);
    builder.add_vertex_to_facet(0);
    builder.end_facet();

    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(0); */
    /*   f.push_back(7); */
    /*   f.push_back(6); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(0);
    builder.add_vertex_to_facet(7);
    builder.add_vertex_to_facet(6);
    builder.end_facet();

    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(0); */
    /*   f.push_back(6); */
    /*   f.push_back(2); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(0);
    builder.add_vertex_to_facet(6);
    builder.add_vertex_to_facet(2);
    builder.end_facet();

    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(2); */
    /*   f.push_back(6); */
    /*   f.push_back(3); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(2);
    builder.add_vertex_to_facet(6);
    builder.add_vertex_to_facet(3);
    builder.end_facet();

    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(7); */
    /*   f.push_back(5); */
    /*   f.push_back(6); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(7);
    builder.add_vertex_to_facet(5);
    builder.add_vertex_to_facet(6);
    builder.end_facet();

    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(6); */
    /*   f.push_back(5); */
    /*   f.push_back(3); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(6);
    builder.add_vertex_to_facet(5);
    builder.add_vertex_to_facet(3);
    builder.end_facet();


    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(1); */
    /*   f.push_back(4); */
    /*   f.push_back(2); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(1);
    builder.add_vertex_to_facet(4);
    builder.add_vertex_to_facet(2);
    builder.end_facet();

    /* { */
    /*   std::vector<int> f; */
    /*   f.push_back(2); */
    /*   f.push_back(4); */
    /*   f.push_back(0); */
    /*   add_facet(builder, f); */
    /* } */

    builder.begin_facet();
    builder.add_vertex_to_facet(2);
    builder.add_vertex_to_facet(4);
    builder.add_vertex_to_facet(0);
    builder.end_facet();


    builder.end_surface();
  }

  const double _x0, _x1, _x2, _y0, _y1, _y2;
};
//-----------------------------------------------------------------------------
template <typename Polyhedron>
static void make_box(double x0, double x1, double x2, 
                     double y0, double y1, double y2, 
                     Polyhedron& P)
{
  Build_box<Polyhedron> builder(x0, x1, x2, y0, y1, y2);
  P.delegate(builder);
  assert(P.is_closed());
  assert(P.is_valid());
}
//-----------------------------------------------------------------------------
template <typename Polyhedron>
void make_tetrahedron(Polyhedron &P,
                      typename Polyhedron::Traits::Point_3 a,
                      typename Polyhedron::Traits::Point_3 b,
                      typename Polyhedron::Traits::Point_3 c,
                      typename Polyhedron::Traits::Point_3 d)
{
  P.make_tetrahedron(a, b, c, d);
}
//-----------------------------------------------------------------------------
/* class Build_cone : public CGAL::Modifier_base<csg::Exact_HalfedgeDS> */
/* { */
/*  public: */
/*   Build_cone(const Cone* cone) : _cone(cone) {} */

/*   void operator()(csg::Exact_HalfedgeDS& hds) */
/*   { */
/*     const dolfin::Point axis = (_cone->_top - _cone->_bottom)/(_cone->_top - _cone->_bottom).norm(); */
/*     dolfin::Point initial = generate_orthogonal(axis); */

/*     CGAL::Polyhedron_incremental_builder_3<csg::Exact_HalfedgeDS> builder(hds, true); */

/*     const int num_sides = _cone->_slices; */
/*     const bool top_degenerate = near(_cone->_top_radius, 0.0); */
/*     const bool bottom_degenerate = near(_cone->_bottom_radius, 0.0); */

/*     const int num_vertices = (top_degenerate || bottom_degenerate) ? num_sides+2 : num_sides*2+2; */

/*     builder.begin_surface(num_vertices, num_sides*4); */

/*     const double delta_theta = 2.0 * DOLFIN_PI / num_sides; */
/*     for (int i = 0; i < num_sides; ++i) */
/*     { */
/*       const double theta = i*delta_theta; */
/*       const Point rotated = initial.rotate(axis, theta); */
/*       if (!bottom_degenerate) */
/*       { */
/*         const Point p = _cone->_bottom + rotated*_cone->_bottom_radius; */
/*         const csg::Exact_Point_3 p_(p.x(), p.y(), p.z()); */
/*         add_vertex(builder, p_); */
/*       } */
/*       if (!top_degenerate) */
/*       { */
/*         const Point p = _cone->_top + rotated*_cone->_top_radius; */
/*         const csg::Exact_Point_3 p_(p.x(), p.y(), p.z()); */
/*         add_vertex(builder, p_); */
/*       } */
/*     } */

/*     // The top and bottom vertices */
/*     add_vertex(builder, csg::Exact_Point_3(_cone->_bottom.x(), _cone->_bottom.y(), */
/*                                            _cone->_bottom.z())); */
/*     add_vertex(builder, csg::Exact_Point_3(_cone->_top.x(), _cone->_top.y(), */
/*                                            _cone->_top.z())); */

/*     // bottom vertex has index num_vertices-2, top vertex has index num_vertices-1 */

/*     // Construct the facets on the side. */
/*     // Vertices must be sorted counter clockwise seen from inside. */
/*     for (int i = 0; i < num_sides; ++i) */
/*     { */
/*       if (top_degenerate) */
/*       { */
/*         std::vector<int> f; */
/*         f.push_back((i + 1)%num_sides); */
/*         f.push_back(i); */
/*         f.push_back(num_vertices - 1); */
/*         add_facet(builder, f); */
/*       } */
/*       else if (bottom_degenerate) */
/*       { */
/*         std::vector<int> f; */
/*         f.push_back( (i) ); */
/*         f.push_back( (i + 1) % num_sides); */
/*         f.push_back(num_vertices - 1); */
/*         add_facet(builder, f); */
/*       } */
/*       else */
/*       { */
/*         //Draw the sides as triangles. */
/*         const int vertex_offset = i*2; */

/*         // First triangle */
/*         std::vector<int> f; */
/*         f.push_back(vertex_offset); */
/*         f.push_back(vertex_offset + 1); */
/*         f.push_back((vertex_offset + 2) % (num_sides*2)); */
/*         add_facet(builder, f); */

/*         // Second triangle */
/*         std::vector<int> g; */
/*         g.push_back((vertex_offset + 3) % (num_sides*2)); */
/*         g.push_back((vertex_offset + 2) % (num_sides*2)); */
/*         g.push_back(vertex_offset + 1); */
/*         add_facet(builder, g); */
/*       } */
/*     } */

/*     // Construct the bottom facet. */
/*     if (!bottom_degenerate) */
/*     { */
/*       for (int i = num_sides-1; i >= 0; i -= 1) */
/*       { */
/*         std::vector<int> f; */
/*         if (!top_degenerate) */
/*         { */
/*           f.push_back(num_vertices-2); */
/*           f.push_back( i*2); */
/*           f.push_back( ( (i+1)*2) % (num_sides*2)); */
/*         } */
/*         else */
/*         { */
/*           f.push_back(num_vertices-2); */
/*           f.push_back(i); */
/*           f.push_back( (i+1)%num_sides ); */
/*         } */
/*         add_facet(builder, f); */
/*       } */
/*     } */

/*     // Construct the the top facet */
/*     if (!top_degenerate) */
/*     { */
/*       for (int i = 0; i < num_sides; i++) */
/*       { */
/*         if (!bottom_degenerate) */
/*         { */
/*           std::vector<int> f; */
/*           f.push_back(num_vertices-1); */
/*           f.push_back( ( (i+1)*2)%(num_sides*2) +1 ); */
/*           f.push_back( i*2 + 1 ); */
/*           add_facet(builder, f); */
/*         } */
/*         else */
/*         { */
/*           std::vector<int> f; */
/*           f.push_back(num_vertices-2); */
/*           f.push_back( (i+1)%num_sides); */
/*           f.push_back(i); */
/*           add_facet(builder, f); */
/*         } */
/*       } */
/*     } */

/*     builder.end_surface(); */
/*   } */
/* private: */
/*   const Cone* _cone; */
/* }; */
/* //----------------------------------------------------------------------------- */
/* static void make_cone(const Cone* c, csg::Exact_Polyhedron_3& P) */
/* { */
/*   Build_cone builder(c); */
/*   P.delegate(builder); */
/*   dolfin_assert(P.is_closed()); */
/*   dolfin_assert(P.is_valid()); */
/* } */
/* //----------------------------------------------------------------------------- */
/* static void make_surface3D(const Surface3D* s, csg::Exact_Polyhedron_3& P) */
/* { */
/*   dolfin_assert(s); */
/*   PolyhedronUtils::readSurfaceFile(s->_filename, P); */
/* } */
/* //----------------------------------------------------------------------------- */
#endif
