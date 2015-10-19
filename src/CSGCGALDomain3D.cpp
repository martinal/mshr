// Copyright (C) 2012-2015 Benjamin Kehlet
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
#include <mshr/CSGCGALDomain2D.h>
#include <mshr/CSGCGALMeshGenerator2D.h>
#include <mshr/DolfinMeshUtils.h>

#include "meshclean.h"
#include "triangulation_refinement.h"
#include "Polyhedron_utils.h"

#include <dolfin/geometry/Point.h>
#include <dolfin/math/basic.h>
#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#ifndef MSHR_ENABLE_EXPERIMENTAL
#include <CGAL/Nef_polyhedron_3.h>
#else
#include <CGAL/corefinement_operations.h>
#endif
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Origin.h>
#include <CGAL/Self_intersection_polyhedron_3.h>


#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/convex_hull_3.h>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>

#include <vector>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <set>
#include <cmath>
#include <memory>

namespace
{

// Exact polyhedron
  typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
  typedef Exact_Kernel::Triangle_3                          Exact_Triangle_3;
  typedef Exact_Kernel::Triangle_2                          Exact_Triangle_2;
  typedef Exact_Kernel::Vector_3                            Exact_Vector_3;
  typedef CGAL::Polyhedron_3<Exact_Kernel>                  Exact_Polyhedron_3;
  typedef Exact_Polyhedron_3::HalfedgeDS                    Exact_HalfedgeDS;
  typedef Exact_Kernel::Point_3                             Exact_Point_3;
  typedef Exact_Kernel::Point_2                             Exact_Point_2;
  typedef Exact_Kernel::Vector_3                            Vector_3;
  typedef Exact_Kernel::Ray_3                               Ray_3;
  typedef Exact_Kernel::Aff_transformation_3                Aff_transformation_3;

#ifndef MSHR_ENABLE_EXPERIMENTAL
  typedef CGAL::Nef_polyhedron_3<Exact_Kernel>              Nef_polyhedron_3;
#endif

// AABB tree primitives
  typedef CGAL::AABB_face_graph_triangle_primitive<Exact_Polyhedron_3> Primitive;
  typedef CGAL::AABB_traits<Exact_Kernel, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> AABB_Tree;
}

namespace mshr
{
  struct CSGCGALDomain3DImpl
  {
    Exact_Polyhedron_3 p;
  };
}

namespace
{

  /*
  double get_polyline_squared_length(const std::vector<Exact_Point_3>& polyline)
  {
    double length = 0;
    std::vector<Exact_Point_3>::const_iterator it = polyline.begin();
    Exact_Point_3 prev = *it;
    it++;
    for (;it != polyline.end(); it++)
    {
      length += CGAL::to_double((*it-prev).squared_length());
    }
    length += CGAL::to_double((polyline.back()-polyline.front()).squared_length());

    return length;
  }
  */

  //-----------------------------------------------------------------------------
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
    Build_sphere(const mshr::Sphere& sphere) : _sphere(sphere) {}

    void operator()( Exact_HalfedgeDS& hds )
    {
      std::vector<dolfin::Point> initial_vertices { dolfin::Point( 0.0, 0.0, 1.0 ),
          dolfin::Point( 1.0, 0.0, 0.0 ),
          dolfin::Point( 0.0,-1.0, 0.0 ),
          dolfin::Point(-1.0, 0.0, 0.0 ),
          dolfin::Point( 0.0, 1.0, 0.0 ),
          dolfin::Point( 0.0, 0.0,-1.0 )};

      // Note: Some older compilers (eg. gcc on Ubuntu Precise) require std::array<std::size_t, 3> to be
      // given explicitly in the initializer list.
      std::vector<std::array<std::size_t, 3> > initial_triangles { std::array<std::size_t, 3>{{0, 2, 1}},
          std::array<std::size_t, 3>{{0, 3, 2 }},
            std::array<std::size_t, 3>{{0, 4, 3 }},
              std::array<std::size_t, 3>{{0, 1, 4,}},
                std::array<std::size_t, 3>{{5, 4, 1}},
                  std::array<std::size_t, 3>{{5, 1, 2 }},
                    std::array<std::size_t, 3>{{5, 2, 3 }},
                      std::array<std::size_t, 3>{{5, 3, 4 }} };

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
    const mshr::Sphere& _sphere;
  };
  //-----------------------------------------------------------------------------
  void make_sphere(const mshr::Sphere* s, Exact_Polyhedron_3& P)
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
    Build_box(const mshr::Box* box) : _box(box) {}

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

    const mshr::Box* _box;
  };
  //-----------------------------------------------------------------------------
  void make_box(const mshr::Box* b, Exact_Polyhedron_3& P)
  {
    Build_box builder(b);
    P.delegate(builder);
    dolfin_assert(P.is_closed());
    dolfin_assert(P.is_valid());
  }
  //-----------------------------------------------------------------------------
  void make_tetrahedron(const mshr::Tetrahedron* b, Exact_Polyhedron_3& P)
  {
    P.make_tetrahedron(Exact_Point_3(b->a.x(), b->a.y(), b->a.z()),
                       Exact_Point_3(b->b.x(), b->b.y(), b->b.z()),
                       Exact_Point_3(b->c.x(), b->c.y(), b->c.z()),
                       Exact_Point_3(b->d.x(), b->d.y(), b->d.z()));
  }
  //-----------------------------------------------------------------------------
  // Return some unit vector orthogonal to a
  dolfin::Point generate_orthogonal(const dolfin::Point& a)
  {
    const dolfin::Point b(0, 1, 0);
    const dolfin::Point c(0, 0, 1);

    // Find a vector not parallel to a.
    const dolfin::Point d = (fabs(a.dot(b)) < fabs(a.dot(c))) ? b : c;
    const dolfin::Point orthogonal = a.cross(d);
    return orthogonal/orthogonal.norm();
  }
  //-----------------------------------------------------------------------------
  class Build_cylinder : public CGAL::Modifier_base<Exact_HalfedgeDS>
  {
  public:
    Build_cylinder(const mshr::Cylinder* cylinder) : _cylinder(cylinder) {}

    void operator()(Exact_HalfedgeDS& hds)
    {
      const dolfin::Point axis = (_cylinder->_top - _cylinder->_bottom)/(_cylinder->_top - _cylinder->_bottom).norm();
      dolfin::Point initial = generate_orthogonal(axis);
      dolfin_assert(dolfin::near(initial.norm(), 1.0));

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
          add_triangular_facet(builder, (vertex_offset + 3) % (num_sides*2), vertex_offset + 1,
                               (vertex_offset + 2) % (num_sides*2));
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
    const mshr::Cylinder* _cylinder;
  };
  //-----------------------------------------------------------------------------
  void make_cylinder(const mshr::Cylinder* c, Exact_Polyhedron_3& P)
  {
    Build_cylinder builder(c);
    P.delegate(builder);
    dolfin_assert(P.is_closed());
    dolfin_assert(P.is_valid());
  }
  //-----------------------------------------------------------------------------
  class Build_ellipsoid : public CGAL::Modifier_base<Exact_HalfedgeDS>
  {
  public:
    Build_ellipsoid(const mshr::Ellipsoid& ellipsoid) : _ellipsoid(ellipsoid) {}

    void operator()(Exact_HalfedgeDS& hds)
    {
      std::vector<dolfin::Point> initial_vertices { dolfin::Point( 0.0, 0.0, 1.0 ),
          dolfin::Point( 1.0, 0.0, 0.0 ),
          dolfin::Point( 0.0,-1.0, 0.0 ),
          dolfin::Point(-1.0, 0.0, 0.0 ),
          dolfin::Point( 0.0, 1.0, 0.0 ),
          dolfin::Point( 0.0, 0.0,-1.0 )};

      // Note: Some older compilers (eg. gcc on Ubuntu Precise) require std::array<std::size_t, 3> to be
      // given explicitly in the initializer list.
      std::vector<std::array<std::size_t, 3> > initial_triangles { std::array<std::size_t, 3>{{0, 2, 1 }},
          std::array<std::size_t, 3>{{0, 3, 2 }},
            std::array<std::size_t, 3>{{0, 4, 3 }},
              std::array<std::size_t, 3>{{0, 1, 4 }},
                std::array<std::size_t, 3>{{5, 4, 1 }},
                  std::array<std::size_t, 3>{{5, 1, 2 }},
                    std::array<std::size_t, 3>{{5, 2, 3 }},
                      std::array<std::size_t, 3>{{5, 3, 4 }} };

      std::vector<dolfin::Point> vertices;
      std::vector<std::array<std::size_t, 3> > triangles;
      if (_ellipsoid._segments > 1 )
      {
        refine_triangulation(initial_vertices,
                             initial_triangles,
                             _ellipsoid._segments,
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

      dolfin::Point center = _ellipsoid.center;
      const double a = _ellipsoid.a, b = _ellipsoid.b, c = _ellipsoid.c;
      for (const dolfin::Point& p : vertices)
      {
        const double scaling = 1.0/std::sqrt( (p.x()*p.x())/(a*a) +
                                              (p.y()*p.y())/(b*b) +
                                              (p.z()*p.z())/(c*c) );
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

    const mshr::Ellipsoid& _ellipsoid;
  };
  //-----------------------------------------------------------------------------
  void make_ellipsoid(const mshr::Ellipsoid* e, Exact_Polyhedron_3& P)
  {
    Build_ellipsoid builder(*e);
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
                       const std::vector<std::array<std::size_t, 3> >& facets,
                       const std::set<std::size_t>& facets_to_be_skipped)
      : vertices(vertices), facets(facets), facets_to_be_skipped(facets_to_be_skipped){}
    void operator()(HDS& hds)
    {
      std::set<std::size_t> isolated_vertices;
      if (facets_to_be_skipped.size() > 0)
      {
        for (std::size_t i = 0; i < vertices.size(); i++)
          isolated_vertices.insert(i);

        for (std::size_t j = 0; j < facets.size(); j++)
        {
          if (facets_to_be_skipped.count(j) == 0)
          {
            for (auto it2 = facets[j].begin(); it2 != facets[j].end(); it2++)
              isolated_vertices.erase(*it2);
          }
        }
      }

      // std::cout << "Isolated vertices: " << isolated_vertices.size() << std::endl;

      CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);

      builder.begin_surface(vertices.size(), facets.size());

      for (std::size_t i = 0; i < vertices.size(); i++)
      {
        if (isolated_vertices.count(i) == 0)
          builder.add_vertex(Exact_Point_3( vertices[i][0], vertices[i][1], vertices[i][2]) );
      }

      const bool has_isolated_vertices = isolated_vertices.size() > 0;
      std::vector<std::size_t> iv(isolated_vertices.begin(), isolated_vertices.end());
      std::sort(iv.begin(), iv.end());
      for (std::size_t i = 0; i < facets.size(); i++)
      {
        if (facets_to_be_skipped.count(i) == 0)
        {
          if (has_isolated_vertices)
          {
            std::vector<std::size_t> f;
            f.reserve(facets[i].size());
            for (auto it = facets[i].begin(); it != facets[i].end(); it++)
            {
              f.push_back(*it - std::distance(iv.begin(), std::lower_bound(iv.begin(), iv.end(), *it)));
            }
            builder.add_facet(f.begin(), f.end());
          }
          else
            builder.add_facet(facets[i].begin(), facets[i].end());

          if (builder.error())
          {
            dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                                 "read surface from file",
                                 "error in polyhedron builder");
          }
        }
        else
        {
          // std::cout << "Skipping" << std::endl;
        }
      }

      builder.end_surface();

    }
    const std::vector<std::array<double, 3> >& vertices;
    const std::vector<std::array<std::size_t, 3> >& facets;
    const std::set<std::size_t>& facets_to_be_skipped;
  };
  //-----------------------------------------------------------------------------
  void make_surface3D(const mshr::Surface3D* s, Exact_Polyhedron_3& P)
  {
    dolfin_assert(s);

    std::vector<std::array<double, 3> > vertices;
    std::vector<std::array<std::size_t, 3> > facets;
    std::set<std::size_t> skip;

    if (s->_filename == "")
    {
      dolfin_assert(s->mesh);

      std::unique_ptr<dolfin::BoundaryMesh> b;

      if (s->use_cell_domain)
      {
        std::shared_ptr<dolfin::Mesh> m = mshr::DolfinMeshUtils::extract_subdomain(s->mesh, s->cell_domain);
        b.reset(new dolfin::BoundaryMesh(*m, "exterior", false));
      }
      else
      {
        // Extract global boundary of mesh, order with outward pointing normals
        b.reset(new dolfin::BoundaryMesh(*(s->mesh), "exterior", false));
      }

      for (dolfin::VertexIterator v(*b); !v.end(); ++v)
      {
        const dolfin::Point& p = v->point();
        vertices.push_back(std::array<double, 3>{{p[0], p[1], p[2]}});
      }

      for (dolfin::CellIterator c(*b); !c.end(); ++c)
      {
        const unsigned int* vertices = c->entities(0);
        facets.push_back(std::array<std::size_t, 3>{{vertices[0], vertices[1], vertices[2]}});
      }
    }
    else
    {
      boost::filesystem::path fpath(s->_filename);
      if (fpath.extension() == ".off")
      {
        std::ifstream infile(s->_filename);
        infile >> P;
        infile.close();

        if (infile.bad())
        {
          std::stringstream ss;
          ss << "Could not read polyhedral surface from '" << s->_filename << "'";
          dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                               "read surface from off file",
                               ss.str());
        }
      }
      else
      {
        if (fpath.extension() == ".stl")
        {
          mshr::STLFileReader::read(s->_filename, vertices, facets);
        }
        else if (fpath.extension() == ".vtp")
        {
          // TODO: Only if vtk is installed
          mshr::VTPFileReader::read(s->_filename, vertices, facets);
        }
        else
        {
          dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                               "open file to read 3D surface",
                               "Unknown file type");
        }

        log(dolfin::TRACE, "Done reading file");
      }

      if (s->flip_facets)
      {
        log(dolfin::TRACE, "Flipping facets");
        for (std::vector<std::array<std::size_t, 3> >::iterator it = facets.begin();
             it != facets.end(); it++)
        {
          std::array<std::size_t, 3>& t = *it;
          std::swap(t[1], t[2]);
        }
      }

      // std::pair<std::unique_ptr<std::vector<std::array<double, 3> > >,
      //           std::unique_ptr<std::vector<std::array<std::size_t, 3> > > > filtered =
      //   SurfaceConsistency::merge_close_vertices(facets, vertices);

      // log(dolfin::TRACE, "Checking connectivity");

      //mshr::closest_vertices(vertices);

      if (s->repair)
      {
        log(dolfin::TRACE, "Keep only connected component");
        mshr::SurfaceConsistency::orient_component(facets, 0);
        std::size_t start_facet = s->first_facet;
        std::cout << "Starting facet: " << start_facet << std::endl;

        std::set<std::size_t> duplicating;

        // FIXME: Commented out for now. Reintroduct before pushing
        // mshr::SurfaceConsistency::checkConnectivity(facets, duplicating, false);
        // log(dolfin::TRACE, "%u facets filtered out", duplicating.size());
        // skip.insert(duplicating.begin(), duplicating.end());


        // SurfaceConsistency::filterFacets(facets,
        // vertices,
        // start_facet,
        // skip);

        // std::vector<std::array<std::size_t, 3> > filtered_facets;
        // filtered_facets.reserve(facets.size()-skip.size());
        // for (std::size_t i = 0; i < facets.size(); i++)
        // {
        //   if (skip.count(i) == 0)
        //     filtered_facets.push_back(facets[i]);
        // }
      }

      // std::cout << "Duplicating: " << duplicating.size() << std::endl;
      // for (auto it = duplicating.begin(); it != duplicating.end(); it++)
      //   std::cout << *it << " ";
      // std::cout << std::endl;

    } // end read from file

    // Create the polyhedron
    BuildFromFacetList<Exact_HalfedgeDS> builder(vertices, facets, skip);
    P.delegate(builder);
    log(dolfin::TRACE, "Done creating polyhedron");

    if (!P.is_valid())
    {
      dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                           "read surface from file",
                           "Polyhedron is not valid. If you are sure your file is valid, please file a bug report");
    }

    P.normalize_border();

    // {
    //   std::size_t num_vertices = P.size_of_vertices();
    //   std::list<typename Exact_Polyhedron_3::Vertex_const_handle> components;
    //   mshr::PolyhedronUtils::get_disconnected_components(P, std::back_inserter(components));
    //   std::cout << "Number of components: " << components.size() << std::endl;
    //   const unsigned int deleted = P.keep_largest_connected_components(1);
    //   std::cout << "Deleted " << deleted << " disconnected components with " << num_vertices-P.size_of_vertices() << " vertices" << std::endl;
    //   std::cout << "Min vertex degree before closing: " << mshr::PolyhedronUtils::min_vertex_degree(P) << std::endl;
    //   std::cout << "Is pure triangular: " << (P.is_pure_triangle() ? "True" : "False") << std::endl;
    //   int tmp;
    //   std::cin >> tmp;
    // }

    // Triangulate polyhedron
    //mshr::PolyhedronUtils::triangulate_polyhedron(P);
    dolfin_assert (P.is_pure_triangle());

    // remove self-intersecting facets
    if (s->repair)
    {
      if (s->sharp_features_filter >= 0)
      {
        mshr::PolyhedronUtils::filter_sharp_features(P,
                                                     s->sharp_features_filter,
                                                     DOLFIN_PI/6.0);
      }

      // FIXME: Should this be enabled?
      // mshr::PolyhedronUtils::remove_self_intersections(P);
      // mshr::PolyhedronUtils::close_holes(P);
    }

    mshr::PolyhedronUtils::list_self_intersections(P);

    // if (s->degenerate_tolerance > 0)
    // {
    //   if (remove_degenerate(P, s->degenerate_tolerance))
    //     log(dolfin::TRACE, "Removed degenerate facets from '%s'",
    //         s->_filename.c_str());
    // }

    // if (!P.is_closed())
    // {
    //   if (s->repair)
    //   {
    //     mshr::PolyhedronUtils::cut_holes(P);

    //     dolfin_assert(P.is_closed());
    //     dolfin_assert(P.is_pure_triangle());
    //     remove_degenerate(P, s->degenerate_tolerance);
    //   }
    //   else
    //   {
    //     dolfin::dolfin_error("CSGCGALDomain3D.cpp",
    //                          "read surface from file",
    //                          "Surface is not closed.");
    //   }
    // }
  }
//-----------------------------------------------------------------------------
  template <class Polyhedron>
  struct Insert_polyhedron_to
    : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
  {
    Insert_polyhedron_to(const Polyhedron& in_poly)
      : _in_poly(in_poly) {}

    void operator()(typename Polyhedron::HalfedgeDS& hds)
    {
      std::cout << "Copying polyhedron" << std::endl;
      typedef typename Polyhedron::HalfedgeDS HDS;

      CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds);

      typedef typename Polyhedron::Vertex_const_iterator Vertex_const_iterator;
      typedef typename Polyhedron::Facet_const_iterator  Facet_const_iterator;
      typedef typename Polyhedron::Halfedge_around_facet_const_circulator HFCC;

      builder.begin_surface(_in_poly.size_of_vertices(),
                            _in_poly.size_of_facets(),
                            _in_poly.size_of_halfedges());

      for(Vertex_const_iterator
            vi = _in_poly.vertices_begin(), end = _in_poly.vertices_end();
          vi != end ; ++vi)
      {
        builder.add_vertex(vi->point());
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
      std::cout << "Done copying" << std::endl;
    } // end operator()(..)
  private:
    const Polyhedron& _in_poly;
  }; // end Copy_polyhedron_to<>
//-----------------------------------------------------------------------------
  template <class HDS>
  class BuildExtrude2D : public CGAL::Modifier_base<HDS>
  {
  public:
    BuildExtrude2D(const mshr::CSGGeometry& polygon, double z)
      : polygon(polygon), z(z){}
    void operator()(HDS& hds)
    {
      // Let the 2d mesh generator triangulate the polygon
      dolfin::Mesh mesh2d;
      {
        mshr::CSGCGALMeshGenerator2D generator;
        generator.parameters["mesh_resolution"] = 2.0;

        generator.generate(polygon, mesh2d);
      }

      CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);

      builder.begin_surface(0, 0);

      // Copy vertices to the new 3d polyhedron
      const std::vector<double>& vertices = mesh2d.coordinates();
      for (std::size_t i = 0; i < vertices.size()/2; i++)
      {
        builder.add_vertex(Exact_Point_3(vertices[2*i], vertices[2*i+1], 0));
        builder.add_vertex(Exact_Point_3(vertices[2*i], vertices[2*i+1], z));
      }

      // Add the triangles from the 2d mesh at z=0 and z=z
      for (dolfin::CellIterator c(mesh2d); !c.end(); ++c)
      {
        const unsigned int* v_indices = c->entities(0);
        const bool flip = Exact_Triangle_2(Exact_Point_2(vertices[2*v_indices[0]], vertices[2*v_indices[0]+1]),
                                           Exact_Point_2(vertices[2*v_indices[1]], vertices[2*v_indices[1]+1]),
                                           Exact_Point_2(vertices[2*v_indices[2]], vertices[2*v_indices[2]+1])).orientation() == CGAL::POSITIVE;
        builder.begin_facet();
        builder.add_vertex_to_facet(2*v_indices[0]);
        if (flip)
        {
          builder.add_vertex_to_facet(2*v_indices[2]);
          builder.add_vertex_to_facet(2*v_indices[1]);
        }
        else
        {
          builder.add_vertex_to_facet(2*v_indices[1]);
          builder.add_vertex_to_facet(2*v_indices[2]);
        }

        builder.end_facet();

        builder.begin_facet();
        builder.add_vertex_to_facet(2*v_indices[0]+1);
        if (!flip)
        {
          builder.add_vertex_to_facet(2*v_indices[2]+1);
          builder.add_vertex_to_facet(2*v_indices[1]+1);
        }
        else
        {
          builder.add_vertex_to_facet(2*v_indices[1]+1);
          builder.add_vertex_to_facet(2*v_indices[2]+1);
        }
        builder.end_facet();
      }

      // Connect the two polygons
      dolfin::BoundaryMesh bdr(mesh2d, "exterior", false);
      const dolfin::MeshFunction<std::size_t>& vertex_map = bdr.entity_map(0);
      for (dolfin::CellIterator cell(bdr); !cell.end(); ++cell)
      {
        const unsigned int* v_indices = cell->entities(0);
        builder.begin_facet();
        builder.add_vertex_to_facet(2*vertex_map[v_indices[0]]);
        builder.add_vertex_to_facet(2*vertex_map[v_indices[1]]);
        builder.add_vertex_to_facet(2*vertex_map[v_indices[0]]+1);
        builder.end_facet();

        builder.begin_facet();
        builder.add_vertex_to_facet(2*vertex_map[v_indices[1]]);
        builder.add_vertex_to_facet(2*vertex_map[v_indices[1]]+1);
        builder.add_vertex_to_facet(2*vertex_map[v_indices[0]]+1);
        builder.end_facet();

      }

      builder.end_surface();
    }

    const mshr::CSGGeometry& polygon;
    const double z;
  };
// ----
  void make_extrude2D(const mshr::Extrude2D* e, Exact_Polyhedron_3& P)
  {
    dolfin_assert(e);

    //CSGCGALDomain2D polygon(e->geometry_2d.get());
    BuildExtrude2D<Exact_HalfedgeDS> builder(*e->geometry_2d, e->z);
    P.delegate(builder);
  }
//-----------------------------------------------------------------------------
  Aff_transformation_3 get_scaling(const mshr::CSGScaling& s)
  {
    Aff_transformation_3 transformation(CGAL::IDENTITY);
    if (s.translate)
    {
      Aff_transformation_3 translation(CGAL::TRANSLATION,
                                       Vector_3(-s.c.x(),
                                                -s.c.y(),
                                                -s.c.z()));
      transformation = translation * transformation;
    }

    Aff_transformation_3 scaling(CGAL::SCALING, s.s);
    transformation = scaling * transformation;

    if (s.translate)
    {
      Aff_transformation_3 translation(CGAL::TRANSLATION,
                                       Vector_3(s.c.x(),
                                                s.c.y(),
                                                s.c.z()));
      transformation = translation * transformation;
    }

    return transformation;
  }
//-----------------------------------------------------------------------------
  Aff_transformation_3 get_rotation(const mshr::CSGRotation& r)
  {
    // Normalize rotation axis vector
    dolfin::Point axis = r.rot_axis/(r.rot_axis.norm());
    dolfin_assert(dolfin::near(axis.norm(), 1.0));

    Aff_transformation_3 transformation(CGAL::IDENTITY);
    if (r.translate)
    {
      Aff_transformation_3 translation(CGAL::TRANSLATION,
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

    Aff_transformation_3 rotation(a*a+b*b-c*c-d*d,
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
      Aff_transformation_3 translation(CGAL::TRANSLATION,
                                       Vector_3(r.c.x(),
                                                r.c.y(),
                                                r.c.z()));
      transformation = translation * transformation;
    }

    return transformation;
  }
//-----------------------------------------------------------------------------
#ifdef MSHR_ENABLE_EXPERIMENTAL
  typedef CGAL::Polyhedron_corefinement<Exact_Polyhedron_3> CGALCSGOperator;

  void convertSubTree(const mshr::CSGGeometry* geometry, Exact_Polyhedron_3& P)
  {
    switch (geometry->getType())
    {
    case mshr::CSGGeometry::Union :
    {
      const mshr::CSGUnion* u = dynamic_cast<const mshr::CSGUnion*>(geometry);
      dolfin_assert(u);
      convertSubTree(u->_g0.get(), P);
      Exact_Polyhedron_3 P2;
      convertSubTree(u->_g1.get(), P2);

      std::list<std::vector<Exact_Point_3> > intersection_polylines;
      CGALCSGOperator op;
      op(P, P2, std::back_inserter(intersection_polylines), CGALCSGOperator::Join_tag);

      // Check that intersection is not degenerate
      for (std::list<std::vector<Exact_Point_3> >::iterator it=intersection_polylines.begin();
           it != intersection_polylines.end(); it++)
      {
	if (get_polyline_squared_length(*it) < DOLFIN_EPS)
	{
	  dolfin::dolfin_error("CSGCGALDomain3D.cpp",
			       "union of csg geometries",
			       "degenerate intersection polyline (geometries meet in a single point?)");
	}
      }
      break;
    }
    case mshr::CSGGeometry::Intersection :
    {
      const mshr::CSGIntersection* u = dynamic_cast<const mshr::CSGIntersection*>(geometry);
      dolfin_assert(u);
      convertSubTree(u->_g0.get(), P);
      Exact_Polyhedron_3 P2;
      convertSubTree(u->_g1.get(), P2);

      std::list<std::vector<Exact_Point_3> > intersection_polylines;
      CGALCSGOperator op;
      op(P, P2, std::back_inserter(intersection_polylines), CGALCSGOperator::Intersection_tag);

      // Check that intersection is not degenerate
      for (std::list<std::vector<Exact_Point_3> >::iterator it=intersection_polylines.begin();
           it != intersection_polylines.end(); it++)
      {
	if (get_polyline_squared_length(*it) < DOLFIN_EPS)
	{
	  dolfin::dolfin_error("CSGCGALDomain3D.cpp",
			       "intersection of csg geometries",
			       "degenerate intersection polyline (geometries meet in a single point?)");
	}
      }

      break;
    }
    case mshr::CSGGeometry::Difference :
    {
      std::cout << "Difference operator" << std::endl;
      const mshr::CSGDifference* u = dynamic_cast<const mshr::CSGDifference*>(geometry);
      dolfin_assert(u);
      convertSubTree(u->_g0.get(), P);
      Exact_Polyhedron_3 P2;
      convertSubTree(u->_g1.get(), P2);

      std::list<std::vector<Exact_Point_3> > intersection_polylines;
      CGALCSGOperator op;
      op(P, P2, std::back_inserter(intersection_polylines), CGALCSGOperator::P_minus_Q_tag);

      std::cout << "Checking polyline" << std::endl;
      // Check that intersection is not degenerate
      for (std::list<std::vector<Exact_Point_3> >::iterator it=intersection_polylines.begin();
           it != intersection_polylines.end(); it++)
      {
	if (get_polyline_squared_length(*it) < DOLFIN_EPS)
	{
	  dolfin::dolfin_error("CSGCGALDomain3D.cpp",
			       "difference of csg geometries",
			       "degenerate intersection polyline (geometries meet in a single point?)");
	}
      }

      break;
    }
    case mshr::CSGGeometry::Translation :
    {
      const mshr::CSGTranslation* t = dynamic_cast<const mshr::CSGTranslation*>(geometry);
      dolfin_assert(t);
      convertSubTree(t->g.get(), P);
      Aff_transformation_3 translation(CGAL::TRANSLATION, Vector_3(t->t.x(), t->t.y(), t->t.z()));
      std::transform(P.points_begin(), P.points_end(), P.points_begin(), translation);
      break;
    }
    case mshr::CSGGeometry::Scaling :
    {
      const mshr::CSGScaling* t = dynamic_cast<const mshr::CSGScaling*>(geometry);
      dolfin_assert(t);
      convertSubTree(t->g.get(), P);

      Aff_transformation_3 scaling = get_scaling(*t);
      std::transform(P.points_begin(), P.points_end(), P.points_begin(), scaling);
      break;
    }
    case mshr::CSGGeometry::Rotation :
    {
      const mshr::CSGRotation* t = dynamic_cast<const mshr::CSGRotation*>(geometry);
      dolfin_assert(t);

      convertSubTree(t->g.get(), P);
      Aff_transformation_3 rotation = get_rotation(*t);
      std::transform(P.points_begin(), P.points_end(), P.points_begin(), rotation);
      break;
    }
    case mshr::CSGGeometry::Cylinder :
    {
      const mshr::Cylinder* c = dynamic_cast<const mshr::Cylinder*>(geometry);
      dolfin_assert(c);
      make_cylinder(c, P);
      break;
    }
    case mshr::CSGGeometry::Sphere :
    {
      const mshr::Sphere* s = dynamic_cast<const mshr::Sphere*>(geometry);
      dolfin_assert(s);
      make_sphere(s, P);
      break;
    }
    case mshr::CSGGeometry::Box :
    {
      const mshr::Box* b = dynamic_cast<const mshr::Box*>(geometry);
      dolfin_assert(b);
      make_box(b, P);
      break;
    }
    case mshr::CSGGeometry::Tetrahedron :
    {
      const mshr::Tetrahedron* b = dynamic_cast<const mshr::Tetrahedron*>(geometry);
      dolfin_assert(b);
      make_tetrahedron(b, P);
      break;
    }
    case mshr::CSGGeometry::Ellipsoid :
    {
      const mshr::Ellipsoid* b = dynamic_cast<const mshr::Ellipsoid*>(geometry);
      dolfin_assert(b);
      make_ellipsoid(b, P);
      break;
    }
    case mshr::CSGGeometry::Surface3D :
    {
      const mshr::Surface3D* b = dynamic_cast<const mshr::Surface3D*>(geometry);
      dolfin_assert(b);
      make_surface3D(b, P);
      break;
    }
    case mshr::CSGGeometry::TriPolyhedron :
    {
      const mshr::CSGCGALDomain3D* b = dynamic_cast<const mshr::CSGCGALDomain3D*>(geometry);
      dolfin_assert(b);
      Insert_polyhedron_to<Exact_Polyhedron_3> inserter(b->impl->p);
      P.delegate(inserter);
      dolfin_asert(P.is_valid());
      break;
    case CSGGeometry::Extrude2D :
    {
      const Extrude2D* e = dynamic_cast<const Extrude2D*>(geometry);
      dolfin_assert(e);
      make_extrude2D(e, P);
      break;
    }
    default:
      dolfin::dolfin_error("CSGCGALDomain.cpp",
                           "converting geometry to cgal polyhedron",
                           "Unhandled primitive type");
    }
    }
#else
    std::shared_ptr<Nef_polyhedron_3>
      convertSubTree(const mshr::CSGGeometry *geometry)
    {
      switch (geometry->getType())
      {
      case mshr::CSGGeometry::Union :
      {
        const mshr::CSGUnion* u = dynamic_cast<const mshr::CSGUnion*>(geometry);
        dolfin_assert(u);
        std::shared_ptr<Nef_polyhedron_3> g0 = convertSubTree(u->_g0.get());
        std::shared_ptr<Nef_polyhedron_3> g1 = convertSubTree(u->_g1.get());
        (*g0) += (*g1);
        return g0;

        break;
      }
      case mshr::CSGGeometry::Intersection :
      {
        const mshr::CSGIntersection* u = dynamic_cast<const mshr::CSGIntersection*>(geometry);
        dolfin_assert(u);
        std::shared_ptr<Nef_polyhedron_3> g0 = convertSubTree(u->_g0.get());
        std::shared_ptr<Nef_polyhedron_3> g1 = convertSubTree(u->_g1.get());
        (*g0) *= (*g1);
        return g0;
        break;
      }
      case mshr::CSGGeometry::Difference :
      {
        const mshr::CSGDifference* u = dynamic_cast<const mshr::CSGDifference*>(geometry);
        dolfin_assert(u);
        std::shared_ptr<Nef_polyhedron_3> g0 = convertSubTree(u->_g0.get());
        std::shared_ptr<Nef_polyhedron_3> g1 = convertSubTree(u->_g1.get());
        (*g0) -= (*g1);
        return g0;
        break;
      }
      case mshr::CSGGeometry::Translation :
      {
        const mshr::CSGTranslation* t = dynamic_cast<const mshr::CSGTranslation*>(geometry);
        dolfin_assert(t);
        std::shared_ptr<Nef_polyhedron_3> g = convertSubTree(t->g.get());
        Aff_transformation_3 translation(CGAL::TRANSLATION, Vector_3(t->t.x(), t->t.y(), t->t.z()));
        g->transform(translation);
        return g;
        break;
      }
      case mshr::CSGGeometry::Scaling :
      {
        const mshr::CSGScaling* t = dynamic_cast<const mshr::CSGScaling*>(geometry);
        dolfin_assert(t);
        std::shared_ptr<Nef_polyhedron_3> g = convertSubTree(t->g.get());
        Aff_transformation_3 scaling = get_scaling(*t);
        g->transform(scaling);
        return g;
        break;
      }
      case mshr::CSGGeometry::Rotation :
      {
        const mshr::CSGRotation* t = dynamic_cast<const mshr::CSGRotation*>(geometry);
        dolfin_assert(t);

        std::shared_ptr<Nef_polyhedron_3> g = convertSubTree(t->g.get());
        Aff_transformation_3 rotation = get_rotation(*t);
        g->transform(rotation);
        return g;
        break;
      }
      case mshr::CSGGeometry::Cylinder :
      {
        const mshr::Cylinder* c = dynamic_cast<const mshr::Cylinder*>(geometry);
        dolfin_assert(c);
        Exact_Polyhedron_3 P;
        make_cylinder(c, P);
        return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
        break;
      }
      case mshr::CSGGeometry::Sphere :
      {
        const mshr::Sphere* s = dynamic_cast<const mshr::Sphere*>(geometry);
        dolfin_assert(s);
        Exact_Polyhedron_3 P;
        make_sphere(s, P);
        return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
        break;
      }
      case mshr::CSGGeometry::Box :
      {
        const mshr::Box* b = dynamic_cast<const mshr::Box*>(geometry);
        dolfin_assert(b);
        Exact_Polyhedron_3 P;
        make_box(b, P);
        return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
        break;
      }
      case mshr::CSGGeometry::Tetrahedron :
      {
        const mshr::Tetrahedron* b = dynamic_cast<const mshr::Tetrahedron*>(geometry);
        dolfin_assert(b);
        Exact_Polyhedron_3 P;
        make_tetrahedron(b, P);
        return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
        break;
      }
      case mshr::CSGGeometry::Ellipsoid :
      {
        const mshr::Ellipsoid* b = dynamic_cast<const mshr::Ellipsoid*>(geometry);
        dolfin_assert(b);
        Exact_Polyhedron_3 P;
        make_ellipsoid(b, P);
        return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
        break;
      }
      case mshr::CSGGeometry::Surface3D :
      {
        const mshr::Surface3D* b = dynamic_cast<const mshr::Surface3D*>(geometry);
        dolfin_assert(b);
        Exact_Polyhedron_3 P;
        make_surface3D(b, P);
        return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
        break;
      }
      case mshr::CSGGeometry::TriPolyhedron :
      {
        const mshr::CSGCGALDomain3D* b = dynamic_cast<const mshr::CSGCGALDomain3D*>(geometry);
        dolfin_assert(b);
        return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(b->impl->p));
        break;
      }
      case mshr::CSGGeometry::Extrude2D :
      {
        const mshr::Extrude2D* e = dynamic_cast<const mshr::Extrude2D*>(geometry);
        dolfin_assert(e);
        Exact_Polyhedron_3 P;
        make_extrude2D(e, P);
        return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
        break;
      }
      default:
        dolfin::dolfin_error("CSGCGALDomain.cpp",
                             "converting geometry to cgal polyhedron",
                             "Unhandled primitive type");
      }
      // Make compiler happy.
      return std::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3);
    }
#endif
//-----------------------------------------------------------------------------
    void convert(const mshr::CSGGeometry& geometry,
                 Exact_Polyhedron_3 &P)
    {
      // If the tree has only one node, we don't have to convert to Nef
      // polyhedrons for csg manipulations
      if (!geometry.is_operator())
      {
        switch (geometry.getType())
        {

        case mshr::CSGGeometry::Cylinder :
        {
          const mshr::Cylinder* c = dynamic_cast<const mshr::Cylinder*>(&geometry);
          dolfin_assert(c);
          make_cylinder(c, P);
          break;
        }
        case mshr::CSGGeometry::Sphere :
        {
          const mshr::Sphere* s = dynamic_cast<const mshr::Sphere*>(&geometry);
          dolfin_assert(s);
          make_sphere(s, P);
          break;
        }
        case mshr::CSGGeometry::Box :
        {
          const mshr::Box* b = dynamic_cast<const mshr::Box*>(&geometry);
          dolfin_assert(b);
          make_box(b, P);
          break;
        }

        case mshr::CSGGeometry::Tetrahedron :
        {
          const mshr::Tetrahedron* b = dynamic_cast<const mshr::Tetrahedron*>(&geometry);
          dolfin_assert(b);
          make_tetrahedron(b, P);
          break;
        }
        case mshr::CSGGeometry::Ellipsoid :
        {
          const mshr::Ellipsoid* b = dynamic_cast<const mshr::Ellipsoid*>(&geometry);
          dolfin_assert(b);
          make_ellipsoid(b, P);
          break;
        }
        case mshr::CSGGeometry::TriPolyhedron :
        {
          const mshr::CSGCGALDomain3D* p = dynamic_cast<const mshr::CSGCGALDomain3D*>(&geometry);
          dolfin_assert(p);
          Insert_polyhedron_to<Exact_Polyhedron_3> inserter(p->impl->p);
          P.delegate(inserter);
          dolfin_assert(P.is_valid());
          break;
        }
        case mshr::CSGGeometry::Surface3D :
        {
          const mshr::Surface3D* b = dynamic_cast<const mshr::Surface3D*>(&geometry);
          dolfin_assert(b);
          make_surface3D(b, P);
          break;
        }
        case mshr::CSGGeometry::Extrude2D :
        {
          const mshr::Extrude2D* e = dynamic_cast<const mshr::Extrude2D*>(&geometry);
          dolfin_assert(e);
          make_extrude2D(e, P);
          break;
        }

        default:
          dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                               "converting geometry to cgal polyhedron",
                               "Unhandled primitive type");
        }
      }
      else
      {
#ifdef MSHR_ENABLE_EXPERIMENTAL
        convertSubTree(&geometry, P);
#else
        log(dolfin::TRACE, "Convert to nef polyhedron");
        std::shared_ptr<Nef_polyhedron_3> cgal_geometry
          = convertSubTree(&geometry);
        dolfin_assert(cgal_geometry->is_valid());
        dolfin_assert(cgal_geometry->is_simple());
        cgal_geometry->convert_to_polyhedron(P);
#endif
      }

      if (P.size_of_facets() == 0)
      {
        dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                             "Convert geometry to polyhedron",
                             "Geometry contains no facet");
      }

      log(dolfin::TRACE, "Number of vertices: %d",  P.size_of_vertices());
      log(dolfin::TRACE, "Number of facets: %d", P.size_of_facets());
    }
  } //end unnamed namespace


  namespace mshr
  {
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
      : impl(std::move(impl))
    {
      // Do nothing
    }
    //-----------------------------------------------------------------------------
    CSGCGALDomain3DQueryStructure::~CSGCGALDomain3DQueryStructure()
    {
      // Do nothing
    }
    //-----------------------------------------------------------------------------
    CSGCGALDomain3D::CSGCGALDomain3D() : impl(new CSGCGALDomain3DImpl)
    {
      parameters = default_parameters();
    }
    //-----------------------------------------------------------------------------
    CSGCGALDomain3D::CSGCGALDomain3D(const mshr::CSGGeometry &csg)
      : impl(new CSGCGALDomain3DImpl)
    {
      parameters = default_parameters();

      if (csg.dim() != 3)
      {
        dolfin::dolfin_error("CSGCGALDomain3D.cpp",
                             "Creating polyhedral domain",
                             "Geometry has dimension %d, expected 3", csg.dim());
      }

      convert(csg, impl->p);

    }
    //-----------------------------------------------------------------------------
    CSGCGALDomain3D::~CSGCGALDomain3D(){}
  //-----------------------------------------------------------------------------
    void CSGCGALDomain3D::insert(const CSGCGALDomain3D& p)
    {
      Insert_polyhedron_to<Exact_Polyhedron_3> inserter(p.impl->p);
      impl->p.delegate(inserter);
      CGAL_assertion(impl->p.is_valid());
    }
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
    std::shared_ptr<std::vector<double>> CSGCGALDomain3D::get_vertices() const
    {
      typedef typename Exact_Polyhedron_3::Vertex_const_iterator Vertex_const_iterator;

      //v.reserve(impl->p.size_of_vertices());
      std::shared_ptr<std::vector<double>> v(new std::vector<double>(impl->p.size_of_vertices()*3));

      std::size_t i = 0;
      for(Vertex_const_iterator
            vi = impl->p.vertices_begin(), end = impl->p.vertices_end();
          vi != end ; ++vi)
      {
        (*v)[i] = ::CGAL::to_double( vi->point().x() ); i++;
        (*v)[i] = ::CGAL::to_double( vi->point().y() ); i++;
        (*v)[i] = ::CGAL::to_double( vi->point().z() ); i++;
      }

      return v;
    }
//-----------------------------------------------------------------------------
    std::shared_ptr<std::vector<std::size_t>> CSGCGALDomain3D::get_facets() const
    {
      typedef typename Exact_Polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
      typedef typename Exact_Polyhedron_3::Facet_const_iterator  Facet_const_iterator;
      typedef typename Exact_Polyhedron_3::Halfedge_around_facet_const_circulator HFCC;

      typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
      Index index(impl->p.vertices_begin(), impl->p.vertices_end());

      std::shared_ptr<std::vector<std::size_t>> v(new std::vector<std::size_t>(impl->p.size_of_facets()*3));

      std::size_t i = 0;
      for(Facet_const_iterator
            fi = impl->p.facets_begin(), end = impl->p.facets_end();
          fi != end; ++fi)
      {

        HFCC hc = fi->facet_begin();
        (*v)[i] = index[hc->vertex()]; i++;
        hc++;
        (*v)[i] = index[hc->vertex()]; i++;
        hc++;
        (*v)[i] = index[hc->vertex()]; i++;
      }

      return v;
    }
//-----------------------------------------------------------------------------
    void CSGCGALDomain3D::get_points_in_holes(std::vector<dolfin::Point>& holes,
                                              std::shared_ptr<CSGCGALDomain3DQueryStructure> q) const
    {
      std::vector<typename Exact_Polyhedron_3::Vertex_const_handle> parts;
      mshr::PolyhedronUtils::get_disconnected_components(impl->p,std::back_inserter(parts));

      if (parts.size() > 1)
      {
        for (std::vector<typename Exact_Polyhedron_3::Vertex_const_handle>::const_iterator it = parts.begin();
             it != parts.end(); it++)
        {
          typename Exact_Polyhedron_3::Halfedge_const_handle h = (*it)->halfedge();
          const Exact_Point_3 x1 = h->vertex()->point();
          const Exact_Point_3 x2 = h->next()->vertex()->point();
          const Exact_Point_3 x3 = h->next()->next()->vertex()->point();
          const Exact_Point_3 origin( (x1.x()+x2.x()+x3.x())/3,
                                      (x1.y()+x2.y()+x3.y())/3,
                                      (x1.z()+x2.z()+x3.z())/3);

          Vector_3 n = CGAL::cross_product(x3-x1, x2-x1);

          // shoot a ray from the facet and count the number of hits
          Ray_3 ray(origin, n);

          std::list<AABB_Tree::Object_and_primitive_id> intersections;
          q->impl->aabb_tree.all_intersections(ray, std::back_inserter(intersections));
          // std::cout << "Number of intersections: " << intersections.size() << std::endl;

          // Collect unique intersection points not including the origin point
          std::set<Exact_Point_3> intersection_points;
          for (std::list<AABB_Tree::Object_and_primitive_id>::const_iterator
                 iit = intersections.begin();
               iit != intersections.end(); iit++)
          {
            if (const Exact_Point_3* p = CGAL::object_cast<Exact_Point_3>(&(iit->first)))
            {
              if (*p != origin)
                intersection_points.insert(*p);
            }
          }

          // std::cout << "Number of unique point intersections: " << intersection_points.size() << std::endl;
          if (intersection_points.size() % 2 == 0)
          {
            // This is a hole
            // Find a point strictly inside it

            std::set<Exact_Point_3>::const_iterator pit = intersection_points.begin();
            Exact_Point_3 closest_point = *pit;
            Exact_Kernel::FT distance_to_closest = (closest_point-origin).squared_length();
            pit++;

            for (;pit != intersection_points.end(); pit++)
            {
              Exact_Kernel::FT d = (*pit-origin).squared_length();
              if (d < distance_to_closest)
              {
                distance_to_closest = d;
                closest_point = *pit;
              }
            }

            const dolfin::Point h_point( (CGAL::to_double(origin.x())+CGAL::to_double(closest_point.x()))/2,
                                         (CGAL::to_double(origin.y())+CGAL::to_double(closest_point.y()))/2,
                                         (CGAL::to_double(origin.z())+CGAL::to_double(closest_point.z()))/2 );
            // std::cout << "Point in hole: " << h_point << std::endl;
            holes.push_back(h_point);
          }
        }
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
    bool CSGCGALDomain3D::is_selfintersecting(bool verbose) const
    {
      const bool selfintersects = CGAL::self_intersect<Exact_Kernel, Exact_Polyhedron_3>(impl->p);
      if (selfintersects && verbose)
        mshr::PolyhedronUtils::list_self_intersections(impl->p);
      return selfintersects;
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
      outfile << std::setprecision(16);

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
    std::size_t CSGCGALDomain3D::num_short_edges(double tolerance) const
    {
      std::size_t count = 0;
      for (Exact_Polyhedron_3::Edge_const_iterator it = impl->p.edges_begin(); it != impl->p.edges_end(); it++)
      {
        const Exact_Kernel::FT l = (it->vertex()->point() - it->opposite()->vertex()->point()).squared_length();
        if (l < tolerance)
          count++;
      }

      return count;
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
//-----------------------------------------------------------------------------
    std::shared_ptr<CSGCGALDomain3D> CSGCGALDomain3D::convex_hull(const CSGCGALDomain3D& c)
    {
      std::shared_ptr<CSGCGALDomain3D> res(new CSGCGALDomain3D);
      CGAL::convex_hull_3(c.impl->p.points_begin(), c.impl->p.points_end(), res->impl->p);

      return res;
    }
//-----------------------------------------------------------------------------
    void CSGCGALDomain3D::filter_facets(dolfin::Point start,
                                        double threshold,
                                        std::shared_ptr<CSGCGALDomain3DQueryStructure> q)
    {
      std::cout << "Filtering facets" << std::endl;

      typedef AABB_Tree::Point_and_primitive_id Point_and_primitive_id;
      typedef Exact_Polyhedron_3::Face_handle Face_handle;
      typedef Exact_Polyhedron_3::Halfedge_handle Halfedge_handle;
      Exact_Point_3 query_point(start[0], start[1], start[2]);
      Point_and_primitive_id pp = q->impl->aabb_tree.closest_point_and_primitive(query_point);
      std::cout << "Closest point: " << pp.first << std::endl;

      const double cos_threshold = cos(threshold);

      std::set<Face_handle> to_be_removed;
      {
        std::deque<Face_handle> queue;
        queue.push_back(pp.second);

        while (!queue.empty())
        {
          Face_handle f = queue.front();
          queue.pop_front();

          if (to_be_removed.count(f) > 0)
            continue;

          to_be_removed.insert(f);

          const Halfedge_handle start = f->halfedge();
          Halfedge_handle c = start;
          const Exact_Triangle_3 f_triangle = PolyhedronUtils::get_facet_triangle<Exact_Polyhedron_3>(f->halfedge());
          do
          {
            Halfedge_handle current = c;
            c = c->next();

            if (current->is_border_edge())
              continue;

            Face_handle f2 = current->opposite()->facet();
            if (PolyhedronUtils::get_triangle_cos_angle<Exact_Triangle_3>(f_triangle,
                                                                          PolyhedronUtils::get_facet_triangle<Exact_Polyhedron_3>(f2->halfedge())) > cos_threshold &&
                to_be_removed.count(f2) == 0)
            {
              queue.push_back(f2);
            }

          } while (c != start);
        }
      }

      std::cout << "Removing " << to_be_removed.size() << " facets" << std::endl;

      for (auto fit = to_be_removed.begin(); fit != to_be_removed.end(); fit++)
      {
        impl->p.erase_facet((*fit)->halfedge());
      }

      dolfin_assert(impl->p.is_valid());
    }
//-----------------------------------------------------------------------------
    void CSGCGALDomain3D::inside_out()
    {
      impl->p.inside_out();
    }
//-----------------------------------------------------------------------------
    std::size_t CSGCGALDomain3D::num_holes() const
    {
      impl->p.normalize_border();

      dolfin_assert(impl->p.is_valid(false, 0));

      const std::vector<Exact_Polyhedron_3::Halfedge_handle> holes = PolyhedronUtils::get_holes(impl->p);
      return holes.size();
    }
    //-----------------------------------------------------------------------------
    void CSGCGALDomain3D::close_holes(std::size_t max, std::size_t offset)
    {
      dolfin::warning("Hole closing is an experimental feature");
      dolfin_assert(impl->p.is_valid(false, 0));
      impl->p.normalize_border();

      const std::vector<Exact_Polyhedron_3::Halfedge_handle> holes = PolyhedronUtils::get_holes(impl->p);
      // std::cout << "Num holes: " << holes.size() << std::endl;

      std::size_t counter = 0;
      while (!impl->p.is_closed() && (max == 0 || counter < max) && offset+counter < holes.size())
      {
        // save_off("not_intersecting.off");
        PolyhedronUtils::close_hole(impl->p,
                                    holes[offset+counter]);
        impl->p.normalize_border();

        counter++;

        dolfin_assert(impl->p.is_valid(false, 0));
        dolfin_assert(impl->p.is_pure_triangle());
        // save_off("closed_hole_intersecting.off");
        dolfin_assert(!is_selfintersecting());
      }
    }
  } // end namespace mshr
