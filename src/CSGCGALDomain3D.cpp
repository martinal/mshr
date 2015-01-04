// Copyright (C) 2012-2014 Benjamin Kehlet
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

#include <mshr/CSGCGALDomain3D.h>
#include <mshr/CSGGeometry.h>
#include <mshr/CSGOperators.h>
#include <mshr/CSGPrimitives3D.h>
#include <mshr/STLFileReader.h>
#include <mshr/VTPFileReader.h>
#include <mshr/SurfaceConsistency.h>

#include "meshclean.h"
#include "triangulate_polyhedron.h"
#include "triangulation_refinement.h"
#include "Polyhedron_utils.h"

#include <dolfin/geometry/Point.h>
#include <dolfin/math/basic.h>
#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Origin.h>
#include <CGAL/Self_intersection_polyhedron_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>


#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>

#include <vector>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <set>
#include <cmath>
#include <memory>

using namespace mshr;

namespace
{

// Exact polyhedron
typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
typedef Exact_Kernel::Triangle_3                          Exact_Triangle_3;
typedef Exact_Kernel::Vector_3                            Exact_Vector_3;
typedef CGAL::Nef_polyhedron_3<Exact_Kernel>              Nef_polyhedron_3;
typedef CGAL::Polyhedron_3<Exact_Kernel>                  Exact_Polyhedron_3;
typedef Exact_Polyhedron_3::HalfedgeDS                    Exact_HalfedgeDS;
typedef Nef_polyhedron_3::Point_3                         Exact_Point_3;
typedef Exact_Kernel::Vector_3                            Vector_3;
typedef Exact_Kernel::Ray_3                               Ray_3;

// AABB tree primitives
typedef CGAL::AABB_face_graph_triangle_primitive<Exact_Polyhedron_3> Primitive;
typedef CGAL::AABB_traits<Exact_Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> AABB_Tree;


// Convenience routine to make debugging easier. Remove before releasing.
template<typename Builder>
inline void add_triangular_facet(Builder& builder,
                                 int v0, int v1, int v2)
{
  builder.begin_facet();
  builder.add_vertex_to_facet(v0);
  builder.add_vertex_to_facet(v1);
  builder.add_vertex_to_facet(v2);
  builder.end_facet();
}
//-----------------------------------------------------------------------------
template<typename Builder>
inline void add_vertex(Builder& builder,
                const Exact_Point_3& point)
{
  builder.add_vertex(point);
}

//-----------------------------------------------------------------------------
// Sphere
//-----------------------------------------------------------------------------
class Build_sphere : public CGAL::Modifier_base<Exact_HalfedgeDS>
{
 public:
  Build_sphere(const Sphere& sphere) : _sphere(sphere) {}

  void operator()( Exact_HalfedgeDS& hds )
  {
    std::vector<dolfin::Point> initial_vertices { dolfin::Point( 0.0, 0.0, 1.0 ),
                                                  dolfin::Point( 1.0, 0.0, 0.0 ),
                                                  dolfin::Point( 0.0,-1.0, 0.0 ),
                                                  dolfin::Point(-1.0, 0.0, 0.0 ),
                                                  dolfin::Point( 0.0, 1.0, 0.0 ),
                                                  dolfin::Point( 0.0, 0.0,-1.0 )};

    std::vector<std::array<std::size_t, 3> > initial_triangles { {0, 1, 2 },
                                                                 {0, 2, 3 },
                                                                 {0, 3, 4 },
                                                                 {0, 4, 1 },
                                                                 {5, 1, 4 },
                                                                 {5, 2, 1 },
                                                                 {5, 3, 2 },
                                                                 {5, 4, 3 } };

    std::vector<dolfin::Point> vertices;
    std::vector<std::array<std::size_t, 3> > triangles;
    if (_sphere._segments > 1 )
    {
      refine_triangulation(initial_vertices,
                           initial_triangles,
                           _sphere._segments,
                           vertices,
                           triangles);
    }
    else
    {
      vertices.reserve(initial_vertices.size());
      std::copy(initial_vertices.begin(), initial_vertices.end(), std::back_inserter(vertices));

      triangles.reserve(initial_triangles.size());
      std::copy(initial_triangles.begin(), initial_triangles.end(), std::back_inserter(triangles));
    }

    CGAL::Polyhedron_incremental_builder_3<Exact_HalfedgeDS> builder( hds, true );
    builder.begin_surface(vertices.size(), triangles.size());

    dolfin::Point center = _sphere.c;
    for (const dolfin::Point& p : vertices)
    {
      const double scaling = _sphere.r/std::sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z());
      add_vertex(builder, Exact_Point_3(center.x() + p.x()*scaling,
                                        center.y() + p.y()*scaling,
                                        center.z() + p.z()*scaling));
    }

    for (const std::array<std::size_t, 3>& t : triangles)
    {
      add_triangular_facet(builder, t[0], t[1], t[2]);
    }

    builder.end_surface();
  }

  private:
  const Sphere& _sphere;
};
//-----------------------------------------------------------------------------
void make_sphere(const Sphere* s, Exact_Polyhedron_3& P)
{
  Build_sphere builder(*s);
  P.delegate(builder);
  dolfin_assert(P.is_valid());
  dolfin_assert(P.is_closed());
}
//-----------------------------------------------------------------------------
class Build_box : public CGAL::Modifier_base<Exact_HalfedgeDS>
{
 public:
  Build_box(const Box* box) : _box(box) {}

  void operator()( Exact_HalfedgeDS& hds )
  {
    CGAL::Polyhedron_incremental_builder_3<Exact_HalfedgeDS> builder(hds, true);

    builder.begin_surface(8, 12);

    const double x0 = std::min(_box->a.x(), _box->b.x());
    const double y0 = std::max(_box->a.x(), _box->b.x());

    const double x1 = std::min(_box->a.y(), _box->b.y());
    const double y1 = std::max(_box->a.y(), _box->b.y());

    const double x2 = std::min(_box->a.z(), _box->b.z());
    const double y2 = std::max(_box->a.z(), _box->b.z());

    add_vertex(builder, Exact_Point_3(y0, x1, x2));
    add_vertex(builder, Exact_Point_3(x0, x1, y2));
    add_vertex(builder, Exact_Point_3(x0, x1, x2));
    add_vertex(builder, Exact_Point_3(x0, y1, x2));
    add_vertex(builder, Exact_Point_3(y0, x1, y2));
    add_vertex(builder, Exact_Point_3(x0, y1, y2));
    add_vertex(builder, Exact_Point_3(y0, y1, x2));
    add_vertex(builder, Exact_Point_3(y0, y1, y2));


    add_triangular_facet(builder, 1, 3, 2);
    add_triangular_facet(builder, 1, 5, 3);
    add_triangular_facet(builder, 1, 4, 5);
    add_triangular_facet(builder, 4, 7, 5);
    add_triangular_facet(builder, 4, 0, 7);
    add_triangular_facet(builder, 0, 6, 7);
    add_triangular_facet(builder, 0, 2, 6);
    add_triangular_facet(builder, 2, 3, 6);
    add_triangular_facet(builder, 7, 6, 5);
    add_triangular_facet(builder, 6, 3, 5);
    add_triangular_facet(builder, 1, 2, 4);
    add_triangular_facet(builder, 2, 0, 4);

    builder.end_surface();
  }

  const Box* _box;
};
//-----------------------------------------------------------------------------
void make_box(const Box* b, Exact_Polyhedron_3& P)
{
  Build_box builder(b);
  P.delegate(builder);
  dolfin_assert(P.is_closed());
  dolfin_assert(P.is_valid());
}
//-----------------------------------------------------------------------------
void make_tetrahedron(const Tetrahedron* b, Exact_Polyhedron_3& P)
{
  P.make_tetrahedron(Exact_Point_3(b->a.x(), b->a.y(), b->a.z()),
                     Exact_Point_3(b->b.x(), b->b.y(), b->b.z()),
                     Exact_Point_3(b->c.x(), b->c.y(), b->c.z()),
                     Exact_Point_3(b->d.x(), b->d.y(), b->d.z()));
}
//-----------------------------------------------------------------------------
// Return some vector orthogonal to a
dolfin::Point generate_orthogonal(const dolfin::Point& a)
{
  const dolfin::Point b(0, 1, 0);
  const dolfin::Point c(0, 0, 1);

  // Find a vector not parallel to a.
  const dolfin::Point d = (fabs(a.dot(b)) < fabs(a.dot(c))) ? b : c;
  return a.cross(d);
}
//-----------------------------------------------------------------------------
class Build_cylinder : public CGAL::Modifier_base<Exact_HalfedgeDS>
{
 public:
  Build_cylinder(const Cylinder* cylinder) : _cylinder(cylinder) {}

  void operator()(Exact_HalfedgeDS& hds)
  {
    const dolfin::Point axis = (_cylinder->_top - _cylinder->_bottom)/(_cylinder->_top - _cylinder->_bottom).norm();
    dolfin::Point initial = generate_orthogonal(axis);

    CGAL::Polyhedron_incremental_builder_3<Exact_HalfedgeDS> builder(hds, true);

    const int num_sides = _cylinder->_segments;
    const bool top_degenerate = dolfin::near(_cylinder->_top_radius, 0.0);
    const bool bottom_degenerate = dolfin::near(_cylinder->_bottom_radius, 0.0);

    const int num_vertices = (top_degenerate || bottom_degenerate) ? num_sides+2 : num_sides*2+2;

    builder.begin_surface(num_vertices, num_sides*4);

    const double delta_theta = 2.0 * DOLFIN_PI / num_sides;
    for (int i = 0; i < num_sides; ++i)
    {
      const double theta = i*delta_theta;
      const dolfin::Point rotated = initial.rotate(axis, theta);
      if (!bottom_degenerate)
      {
        const dolfin::Point p = _cylinder->_bottom + rotated*_cylinder->_bottom_radius;
        const Exact_Point_3 p_(p.x(), p.y(), p.z());
        add_vertex(builder, p_);
      }
      if (!top_degenerate)
      {
        const dolfin::Point p = _cylinder->_top + rotated*_cylinder->_top_radius;
        const Exact_Point_3 p_(p.x(), p.y(), p.z());
        add_vertex(builder, p_);
      }
    }

    // The top and bottom vertices
    add_vertex(builder, Exact_Point_3(_cylinder->_bottom.x(), _cylinder->_bottom.y(),
                                           _cylinder->_bottom.z()));
    add_vertex(builder, Exact_Point_3(_cylinder->_top.x(), _cylinder->_top.y(),
                                           _cylinder->_top.z()));

    // bottom vertex has index num_vertices-2, top vertex has index num_vertices-1

    // Construct the facets on the side.
    // Vertices must be sorted counter clockwise seen from inside.
    for (int i = 0; i < num_sides; ++i)
    {
      if (top_degenerate)
      {
        add_triangular_facet(builder, (i + 1)%num_sides, num_vertices - 1, i);
      }
      else if (bottom_degenerate)
      {
        add_triangular_facet(builder, i, num_vertices - 1, (i + 1) % num_sides);
      }
      else
      {
        //Draw the sides as triangles.
        const int vertex_offset = i*2;

        // First triangle
        add_triangular_facet(builder, vertex_offset, (vertex_offset + 2) % (num_sides*2), vertex_offset + 1);

        // Second triangle
        add_triangular_facet(builder, (vertex_offset + 3) % (num_sides*2), vertex_offset + 1, (vertex_offset + 2) % (num_sides*2));
      }
    }

    // Construct the bottom facet.
    if (!bottom_degenerate)
    {
      for (int i = num_sides-1; i >= 0; i -= 1)
      {
        if (!top_degenerate)
        {
          add_triangular_facet(builder, num_vertices-2,( (i+1)*2) % (num_sides*2), i*2);
        }
        else
        {
          add_triangular_facet(builder, num_vertices-2, (i+1)%num_sides, i);
        }
      }
    }

    // Construct the the top facet
    if (!top_degenerate)
    {
      for (int i = 0; i < num_sides; i++)
      {
        if (!bottom_degenerate)
        {
          add_triangular_facet(builder, num_vertices-1, i*2 + 1, ( (i+1)*2)%(num_sides*2) +1);
        }
        else
        {
          add_triangular_facet(builder, num_vertices-2, i, (i+1)%num_sides);
        }
      }
    }

    builder.end_surface();
  }
private:
  const Cylinder* _cylinder;
};
//-----------------------------------------------------------------------------
 void make_cylinder(const Cylinder* c, Exact_Polyhedron_3& P)
{
  Build_cylinder builder(c);
  P.delegate(builder);
  dolfin_assert(P.is_closed());
  dolfin_assert(P.is_valid());
}
//-----------------------------------------------------------------------------
template <class HDS>
class BuildFromFacetList : public CGAL::Modifier_base<HDS>
{
public:
  BuildFromFacetList(const std::vector<std::array<double, 3> >& vertices,
                     const std::vector<std::vector<std::size_t> >& facets)
    : vertices(vertices), facets(facets){}
  void operator()(HDS& hds)
  {
    CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);

    builder.begin_surface(vertices.size(), facets.size());
    
    for (std::vector<std::array<double, 3> >::const_iterator it = vertices.begin();
         it != vertices.end(); ++it)
      builder.add_vertex(Exact_Point_3( (*it)[0], (*it)[1], (*it)[2]));

    for (std::vector<std::vector<std::size_t> >::const_iterator it = facets.begin();
         it != facets.end(); ++it)
      builder.add_facet(it->begin(), it->end());

    builder.end_surface();

  }
  const std::vector<std::array<double, 3> > vertices;
  const std::vector<std::vector<std::size_t> > facets;
};
//-----------------------------------------------------------------------------
void make_surface3D(const Surface3D* s, Exact_Polyhedron_3& P)
{
  dolfin_assert(s);

  std::vector<std::array<double, 3> > vertices;
  std::vector<std::vector<std::size_t> > facets;

  boost::filesystem::path fpath(s->_filename);
  if (fpath.extension() == ".off")
  {
    std::ifstream infile(s->_filename);
    infile >> P;
    infile.close();

    if (P.size_of_vertices() == 0 || P.size_of_facets() == 0)
    {
      std::stringstream ss;
      ss << "File '" << s->_filename << "' was empty or contained syntactical errors";
      dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                           "read surface from off file",
                           ss.str());
    }
  }
  else
  {

    if (fpath.extension() == ".stl")
    {
      STLFileReader::read(s->_filename, vertices, facets);
    }
    else if (fpath.extension() == ".vtp")
    {
      // TODO: Only if vtk is installed
      VTPFileReader::read(s->_filename, vertices, facets);
    }
    else
    {
      dolfin::dolfin_error("PolyhedronUtils.cpp",
                           "open file to read 3D surface",
                           "Unknown file type");
    }

    SurfaceConsistency::checkConnectivity(facets);

    // Create the polyhedron
    BuildFromFacetList<Exact_HalfedgeDS> builder(vertices, facets);
    P.delegate(builder);
  }

  if (!P.is_valid())
  {
    dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                         "read surface from file",
                         "Polyhedron is not valid. If you are sure your file is valid, please file a bug report");
  }

  // Triangulate polyhedron
  triangulate_polyhedron(P);
  dolfin_assert (P.is_pure_triangle());

  if (!P.is_closed())
  {
    dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                         "read surface from file",
                         "Surface is not closed.");
  }

  if (s->degenerate_tolerance > 0)
  {
    if (remove_degenerate(P, s->degenerate_tolerance))
      log(dolfin::TRACE, "Removed degenerate facets from '%s'",
          s->_filename.c_str());
  }
}
//-----------------------------------------------------------------------------
void do_scaling(const CSGScaling& s, Nef_polyhedron_3& p)
{
  Exact_Kernel::Aff_transformation_3 transformation(CGAL::IDENTITY);
  if (s.translate)
  {
    Exact_Kernel::Aff_transformation_3 translation(CGAL::TRANSLATION,
                                                   Vector_3(-s.c.x(),
                                                            -s.c.y(),
                                                            -s.c.z()));
    transformation = translation * transformation;
  }

  Exact_Kernel::Aff_transformation_3 scaling(CGAL::SCALING, s.s);
  transformation = scaling * transformation;

  if (s.translate)
  {
    Exact_Kernel::Aff_transformation_3 translation(CGAL::TRANSLATION,
                                                   Vector_3(s.c.x(),
                                                            s.c.y(),
                                                            s.c.z()));
    transformation = translation * transformation;
  }

  p.transform(transformation);
}
//-----------------------------------------------------------------------------
void do_rotation(const CSGRotation& r, Nef_polyhedron_3& p)
{

  // Normalize rotation axis vector
  dolfin::Point axis = r.rot_axis/(r.rot_axis.norm());
  dolfin_assert(dolfin::near(axis.norm(), 1.0));

  Exact_Kernel::Aff_transformation_3 transformation(CGAL::IDENTITY);
  if (r.translate)
  {
    Exact_Kernel::Aff_transformation_3 translation(CGAL::TRANSLATION,
                                                   Vector_3(-r.c.x(),
                                                            -r.c.y(),
                                                            -r.c.z()));
    transformation = translation * transformation;
  }

  // The Euler-Rodrigues formula
  const double a = cos(r.theta/2);
  const double b = -axis.x()*sin(r.theta/2);
  const double c = -axis.y()*sin(r.theta/2);
  const double d = -axis.z()*sin(r.theta/2);

  Exact_Kernel::Aff_transformation_3 rotation(a*a+b*b-c*c-d*d,
                                              2*(b*c-a*d),
                                              2*(b*d+a*c),
                                              2*(b*c+a*d),
                                              a*a+c*c-b*b-d*d,
                                              2*(c*d-a*b),
                                              2*(b*d-a*c),
                                              2*(c*d+a*b),
                                              a*a+d*d-b*b-c*c);
  transformation = rotation * transformation;

  if (r.translate)
  {
    Exact_Kernel::Aff_transformation_3 translation(CGAL::TRANSLATION,
                                                   Vector_3(r.c.x(),
                                                            r.c.y(),
                                                            r.c.z()));
    transformation = translation * transformation;
  }

  p.transform(transformation);
}
//-----------------------------------------------------------------------------
std::shared_ptr<Nef_polyhedron_3>
convertSubTree(const CSGGeometry *geometry)
{
  switch (geometry->getType())
  {
    case CSGGeometry::Union :
    {
      const CSGUnion* u = dynamic_cast<const CSGUnion*>(geometry);
      dolfin_assert(u);
      std::shared_ptr<Nef_polyhedron_3> g0 = convertSubTree(u->_g0.get());
      std::shared_ptr<Nef_polyhedron_3> g1 = convertSubTree(u->_g1.get());
      (*g0) += (*g1);
      return g0;

      break;
    }
    case CSGGeometry::Intersection :
    {
      const CSGIntersection* u = dynamic_cast<const CSGIntersection*>(geometry);
      dolfin_assert(u);
      std::shared_ptr<Nef_polyhedron_3> g0 = convertSubTree(u->_g0.get());
      std::shared_ptr<Nef_polyhedron_3> g1 = convertSubTree(u->_g1.get());
      (*g0) *= (*g1);
      return g0;
      break;
    }
    case CSGGeometry::Difference :
    {
      const CSGDifference* u = dynamic_cast<const CSGDifference*>(geometry);
      dolfin_assert(u);
      std::shared_ptr<Nef_polyhedron_3> g0 = convertSubTree(u->_g0.get());
      std::shared_ptr<Nef_polyhedron_3> g1 = convertSubTree(u->_g1.get());
      (*g0) -= (*g1);
      return g0;
      break;
    }
    case CSGGeometry::Translation :
    {
      const CSGTranslation* t = dynamic_cast<const CSGTranslation*>(geometry);
      dolfin_assert(t);
      std::shared_ptr<Nef_polyhedron_3> g = convertSubTree(t->g.get());
      Exact_Kernel::Aff_transformation_3 translation(CGAL::TRANSLATION, Vector_3(t->t.x(), t->t.y(), t->t.z()));
      g->transform(translation);
      return g;
      break;
    }
    case CSGGeometry::Scaling :
    {
      const CSGScaling* t = dynamic_cast<const CSGScaling*>(geometry);
      dolfin_assert(t);
      std::shared_ptr<Nef_polyhedron_3> g = convertSubTree(t->g.get());
      do_scaling(*t, *g);
      return g;
      break;
    }
    case CSGGeometry::Rotation :
    {
      const CSGRotation* t = dynamic_cast<const CSGRotation*>(geometry);
      dolfin_assert(t);

      std::shared_ptr<Nef_polyhedron_3> g = convertSubTree(t->g.get());
      do_rotation(*t, *g);
      return g;
      break;
    }
    case CSGGeometry::Cylinder :
    {
      const Cylinder* c = dynamic_cast<const Cylinder*>(geometry);
      dolfin_assert(c);
      Exact_Polyhedron_3 P;
      make_cylinder(c, P);
      return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
      break;
    }
    case CSGGeometry::Sphere :
    {
      const Sphere* s = dynamic_cast<const Sphere*>(geometry);
      dolfin_assert(s);
      Exact_Polyhedron_3 P;
      make_sphere(s, P);
      return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
      break;
    }
    case CSGGeometry::Box :
    {
      const Box* b = dynamic_cast<const Box*>(geometry);
      dolfin_assert(b);
      Exact_Polyhedron_3 P;
      make_box(b, P);
      return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
      break;
    }

    case CSGGeometry::Tetrahedron :
    {
      const Tetrahedron* b = dynamic_cast<const Tetrahedron*>(geometry);
      dolfin_assert(b);
      Exact_Polyhedron_3 P;
      make_tetrahedron(b, P);
      return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
      break;
    }
    case CSGGeometry::Surface3D :
    {
      const Surface3D* b = dynamic_cast<const Surface3D*>(geometry);
      dolfin_assert(b);
      Exact_Polyhedron_3 P;
      make_surface3D(b, P);
      return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
      break;
    }
    default:
      dolfin::dolfin_error("GeometryToCGALConverter.cpp",
                           "converting geometry to cgal polyhedron",
                           "Unhandled primitive type");
  }

  // Make compiler happy.
  return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3);
}
//-----------------------------------------------------------------------------
void convert(const CSGGeometry& geometry,
             Exact_Polyhedron_3 &P)
{
  // If the tree has only one node, we don't have to convert to Nef
  // polyhedrons for csg manipulations
  if (!geometry.is_operator())
  {
    switch (geometry.getType())
    {

    case CSGGeometry::Cylinder :
    {
      const Cylinder* c = dynamic_cast<const Cylinder*>(&geometry);
      dolfin_assert(c);
      make_cylinder(c, P);
      break;
    }
    case CSGGeometry::Sphere :
    {
      const Sphere* s = dynamic_cast<const Sphere*>(&geometry);
      dolfin_assert(s);
      make_sphere(s, P);
      break;
    }
    case CSGGeometry::Box :
    {
      const Box* b = dynamic_cast<const Box*>(&geometry);
      dolfin_assert(b);
      make_box(b, P);
      break;
    }

    case CSGGeometry::Tetrahedron :
    {
      const Tetrahedron* b = dynamic_cast<const Tetrahedron*>(&geometry);
      dolfin_assert(b);
      make_tetrahedron(b, P);
      break;
    }
    case CSGGeometry::Surface3D :
    {
      const Surface3D* b = dynamic_cast<const Surface3D*>(&geometry);
      dolfin_assert(b);
      make_surface3D(b, P);
      break;
    }
    default:
      dolfin::dolfin_error("GeometryToCGALConverter.cpp",
                   "converting geometry to cgal polyhedron",
                   "Unhandled primitive type");
    }
  }
  else
  {
    log(dolfin::TRACE, "Convert to nef polyhedron");
    std::shared_ptr<Nef_polyhedron_3> cgal_geometry
      = convertSubTree(&geometry);
    dolfin_assert(cgal_geometry->is_valid());
    dolfin_assert(cgal_geometry->is_simple());
    cgal_geometry->convert_to_polyhedron(P);
  }

  if (P.size_of_facets() == 0)
    dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                         "Convert geometry to polyhedron",
                         "Geometry contains no facet");

  log(dolfin::TRACE, "Number of vertices: %d",  P.size_of_vertices());
  log(dolfin::TRACE, "Number of facets: %d", P.size_of_facets());
}
} //end unnamed namespace


namespace mshr
{

struct CSGCGALDomain3DImpl
{
  Exact_Polyhedron_3 p;
};
//-----------------------------------------------------------------------------
struct CSGCGALDomain3DQueryStructureImpl
{
  template <typename A>
  CSGCGALDomain3DQueryStructureImpl(A start, A end, const Exact_Polyhedron_3& p)
    : aabb_tree(start, end, p){}

  AABB_Tree aabb_tree;
};
//-----------------------------------------------------------------------------
CSGCGALDomain3DQueryStructure::CSGCGALDomain3DQueryStructure(std::unique_ptr<CSGCGALDomain3DQueryStructureImpl> impl)
{}
//-----------------------------------------------------------------------------
CSGCGALDomain3DQueryStructure::~CSGCGALDomain3DQueryStructure(){}
//-----------------------------------------------------------------------------
CSGCGALDomain3D::CSGCGALDomain3D()
: impl(new CSGCGALDomain3DImpl)
{
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
CSGCGALDomain3D::CSGCGALDomain3D(const mshr::CSGGeometry &csg)
: impl(new CSGCGALDomain3DImpl)
{
  parameters = default_parameters();

  if (csg.dim() != 3)
    dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                         "Creating polyhedral domain",
                         "Geometry has dimension %d, expected 3", csg.dim());

  convert(csg,
          impl->p);

}
//-----------------------------------------------------------------------------
CSGCGALDomain3D::~CSGCGALDomain3D(){}
//-----------------------------------------------------------------------------
std::size_t CSGCGALDomain3D::num_vertices() const 
{ return impl->p.size_of_vertices(); }
//-----------------------------------------------------------------------------
std::size_t CSGCGALDomain3D::num_facets() const 
{ return impl->p.size_of_facets(); }
//-----------------------------------------------------------------------------
std::size_t CSGCGALDomain3D::num_halfedges() const 
{ return impl->p.size_of_halfedges(); }
//-----------------------------------------------------------------------------
double CSGCGALDomain3D::volume() const
{
  double volume = .0;
  for (Exact_Polyhedron_3::Facet_iterator it = impl->p.facets_begin();
       it != impl->p.facets_end(); it++)
  {
    const Exact_Polyhedron_3::Halfedge_handle h = it->halfedge();
    const Vector_3 V0 = h->vertex()->point()-CGAL::ORIGIN;
    const Vector_3 V1 = h->next()->vertex()->point()-CGAL::ORIGIN;
    const Vector_3 V2 = h->next()->next()->vertex()->point()-CGAL::ORIGIN;

    volume += CGAL::to_double(V0*CGAL::cross_product(V1, V2));
  }

  return volume/6.0;
}
//-----------------------------------------------------------------------------
void CSGCGALDomain3D::get_vertices(std::vector<dolfin::Point> &v) const
{
  typedef typename Exact_Polyhedron_3::Vertex_const_iterator Vertex_const_iterator;

  v.reserve(impl->p.size_of_vertices());
  
  for(Vertex_const_iterator
        vi = impl->p.vertices_begin(), end = impl->p.vertices_end();
      vi != end ; ++vi)
  {
    v.push_back(dolfin::Point(::CGAL::to_double( vi->point().x()),
                              ::CGAL::to_double( vi->point().y()),
                              ::CGAL::to_double( vi->point().z())));
  }
}
//-----------------------------------------------------------------------------
void CSGCGALDomain3D::get_facets(std::vector< std::array<std::size_t, 3> > &f) const
{
  typedef typename Exact_Polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Exact_Polyhedron_3::Facet_const_iterator  Facet_const_iterator;
  typedef typename Exact_Polyhedron_3::Halfedge_around_facet_const_circulator HFCC;

  typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
  Index index(impl->p.vertices_begin(), impl->p.vertices_end());

  for(Facet_const_iterator
        fi = impl->p.facets_begin(), end = impl->p.facets_end();
      fi != end; ++fi)
  {
    f.push_back(std::array<std::size_t, 3>());
    std::array<std::size_t, 3> &v = f.back();

    HFCC hc = fi->facet_begin();
    v[0] = index[hc->vertex()];
    hc++;
    v[1] = index[hc->vertex()];
    hc++;
    v[2] = index[hc->vertex()];
  }
}
//-----------------------------------------------------------------------------
void CSGCGALDomain3D::get_points_in_holes(std::vector<dolfin::Point> h,
                                          std::shared_ptr<CSGCGALDomain3DQueryStructure> q) const
{
  std::vector<typename Exact_Polyhedron_3::Vertex_const_handle> parts;
  get_disconnected_components(impl->p,std::back_inserter(parts));

  for (std::vector<typename Exact_Polyhedron_3::Vertex_const_handle>::const_iterator it = parts.begin();
       it != parts.end(); it++)
  {
    typename Exact_Polyhedron_3::Halfedge_const_handle h = (*it)->halfedge();
    const Exact_Point_3 x1 = h->vertex()->point();
    const Exact_Point_3 x2 = h->next()->vertex()->point();
    const Exact_Point_3 x3 = h->next()->next()->vertex()->point();
    Exact_Point_3 p( (x1.x()+x2.x()+x3.x())/3,
                     (x1.y()+x2.y()+x3.y())/3,
                     (x1.z()+x2.z()+x3.z())/3);

    Vector_3 n = CGAL::cross_product(x2-x1, x3-x1);

    // shoot a ray from the facet and count the number of hits
    Ray_3 r(p, n);

  }
}
//-----------------------------------------------------------------------------
void CSGCGALDomain3D::remove_degenerate_facets(double tolerance) 
{
  remove_degenerate(impl->p, tolerance);
}
//-----------------------------------------------------------------------------
void CSGCGALDomain3D::ensure_meshing_preconditions()
{
  if (!impl->p.is_valid())
    dolfin::dolfin_error("CSGCGDomain3D.cpp",
                         "Checking meshing preconditions",
                         "Polyhedron is not valid");

  if (parameters["remove_degenerate"])
    remove_degenerate_facets(parameters["degenerate_tolerance"]);
}
//-----------------------------------------------------------------------------
std::shared_ptr<CSGCGALDomain3DQueryStructure> CSGCGALDomain3D::get_query_structure() const
{
  std::unique_ptr<CSGCGALDomain3DQueryStructureImpl> i(new CSGCGALDomain3DQueryStructureImpl(faces(impl->p).first,
                                                                                             faces(impl->p).second,
                                                                                             impl->p));
  return std::shared_ptr<CSGCGALDomain3DQueryStructure>(new CSGCGALDomain3DQueryStructure(std::move(i)));
}
//-----------------------------------------------------------------------------
bool CSGCGALDomain3D::is_insideout() const
{
  typedef Exact_Kernel::Ray_3 Exact_Ray_3;

  // Pick a point on a facet a. Shoot in the direction of the normal and count
  // the number of distinct hits (if the ray hit an edge, the same intersection
  // point will hit sevel facets). Excluding the hit on facet a the number
  // should be even for the polyhedron to be bounded.
  // TODO: Take QueryStructure argument.
  Exact_Polyhedron_3::Facet_iterator f = impl->p.facets_begin();

  Exact_Ray_3 ray;

  {
    Exact_Polyhedron_3::Halfedge_handle h = f->halfedge();
    Exact_Triangle_3 t(h->vertex()->point(),
                       h->next()->vertex()->point(),
                       h->next()->next()->vertex()->point());
    f++;

    //Compute the ray
    Exact_Vector_3 normal = CGAL::normal(t.vertex(0), t.vertex(1), t.vertex(2));
    Exact_Point_3  mid    = CGAL::centroid(t);

    ray = Exact_Ray_3(mid, normal);
  }

  // Collect the points to avoid double hits on facets.
  std::set<Exact_Point_3> points;
  for (; f != impl->p.facets_end(); ++f)
  {
    Exact_Polyhedron_3::Halfedge_handle h = f->halfedge();
    Exact_Triangle_3 t(h->vertex()->point(),
                       h->next()->vertex()->point(),
                       h->next()->next()->vertex()->point());


    auto result = intersection(t, ray);
    if(result)
    {
      if (const Exact_Point_3* p = boost::get<Exact_Point_3>(&*result))
      {
        points.insert(*p);
      }
    }
  }

  // std::cout << "Number of intersections: " << points.size() << std::endl;
  return points.size() % 2 != 0;
}
//-----------------------------------------------------------------------------
bool CSGCGALDomain3D::is_selfintersecting() const
{
  return CGAL::self_intersect<Exact_Kernel, Exact_Polyhedron_3>(impl->p);
}
//-----------------------------------------------------------------------------
std::size_t CSGCGALDomain3D::num_degenerate_facets(double threshold) const
{
  return number_of_degenerate_facets(impl->p, threshold);
}
//-----------------------------------------------------------------------------
void CSGCGALDomain3D::save_off(std::string filename) const
{
  {
    std::string message = "Writing to file: "+filename;
    log(dolfin::TRACE, message);
  }

  std::ofstream outfile(filename.c_str());
  // outfile.precision(16);

  outfile << impl->p;

  // std::map<Exact_Polyhedron_3::Vertex_const_handle, std::size_t> vertex_map;
  // std::size_t vertex_counter = 0;

  // // Write header
  // outfile << "OFF " << std::endl << impl->p.size_of_vertices() << " " << impl->p.size_of_facets() << " 0" << std::endl << std::endl;

  // // Write vertices
  // for (Exact_Polyhedron_3::Vertex_const_iterator vit = impl->p.vertices_begin(); vit != impl->p.vertices_end(); vit++)
  // {
  //   vertex_map[vit] = vertex_counter;
  //   outfile << std::fixed << std::setprecision(16) << CGAL::to_double(vit->point().x()) << " "
  //           << CGAL::to_double(vit->point().y()) << " "
  //           << CGAL::to_double(vit->point().z()) << std::endl;
  //   vertex_counter++;
  // }

  // for (Exact_Polyhedron_3::Facet_const_iterator fit = impl->p.facets_begin(); fit != impl->p.facets_end(); fit++)
  // {
  //   Exact_Polyhedron_3::Halfedge_const_iterator h = fit->halfedge();
  //   outfile << "3";
  //   outfile << vertex_map[h->vertex()] << " "
  //           << vertex_map[h->next()->vertex()] << " "
  //           << vertex_map[h->next()->next()->vertex()] << std::endl;
  // }

  outfile.close();
}
//-----------------------------------------------------------------------------
double CSGCGALDomain3D::shortest_edge() const
{
  Exact_Polyhedron_3::Edge_const_iterator it = impl->p.edges_begin();
  Exact_Kernel::FT shortest = (it->vertex()->point() - it->opposite()->vertex()->point()).squared_length();


  for (; it != impl->p.edges_end(); it++)
  {
    const Exact_Kernel::FT l = (it->vertex()->point() - it->opposite()->vertex()->point()).squared_length();
    if (l < shortest)
      shortest = l;
  }

  return CGAL::to_double(shortest);
}
//-----------------------------------------------------------------------------
std::string CSGCGALDomain3D::str(bool verbose) const
{
  std::stringstream ss;
  ss << "Triangular polyhedron with" << std::endl;
  ss << "  " << num_vertices() << " vertices," << std::endl;
  ss << "  " << num_facets() << " facets," << std::endl;
  ss << "  " << num_halfedges() << " halfedges." << std::endl;

  if (verbose)
  {
    ss << "Volume: " << volume() << std::endl;
    ss << "Shortest edge: " << shortest_edge() << std::endl;
    ss << "Degenerate facets: " << num_degenerate_facets(1e-12) << std::endl;
    ss << "Is inside out:        " << (is_insideout() ? "Yes" : "No") << std::endl;
    ss << "Is self-intersecting: " << (is_selfintersecting() ? "Yes" : "No");
  }

  return ss.str();
}
} // end namespace mshr
