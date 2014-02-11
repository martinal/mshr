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

#include "meshclean.h"

#include <dolfin/geometry/Point.h>
#include <dolfin/math/basic.h>
#include <dolfin/log/LogStream.h>

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_polyhedron_3.h>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>

using namespace mshr;

namespace
{

// Exact polyhedron
typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
typedef Exact_Kernel::Triangle_3 Exact_Triangle_3;
typedef CGAL::Nef_polyhedron_3<Exact_Kernel> Nef_polyhedron_3;
typedef CGAL::Polyhedron_3<Exact_Kernel> Exact_Polyhedron_3;
typedef Exact_Polyhedron_3::HalfedgeDS Exact_HalfedgeDS;
typedef Nef_polyhedron_3::Point_3 Exact_Point_3;


// Convenience routine to make debugging easier. Remove before releasing.
void add_facet(CGAL::Polyhedron_incremental_builder_3<Exact_HalfedgeDS>& builder,
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
void add_vertex(CGAL::Polyhedron_incremental_builder_3<Exact_HalfedgeDS>& builder,
                const Exact_Point_3& point, bool print=false)
{
  static int vertex_no = 0;
  if (print)
    std::cout << "Adding vertex " << vertex_no << " at " << point << std::endl;

  builder.add_vertex(point);
  vertex_no++;
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
    const std::size_t num_slices = _sphere._slices;
    const std::size_t num_sectors = _sphere._slices*2 + 1;

    const dolfin::Point top = _sphere.c + dolfin::Point(_sphere.r, 0, 0);
    const dolfin::Point bottom = _sphere.c - dolfin::Point(_sphere.r, 0, 0);
    const dolfin::Point axis = dolfin::Point(1, 0, 0);

    const int num_vertices = num_slices*num_sectors+2;
    const int num_facets = num_sectors*2*num_slices;

    CGAL::Polyhedron_incremental_builder_3<Exact_HalfedgeDS> builder( hds, true );

    builder.begin_surface(num_vertices, num_facets);

    const dolfin::Point slice_rotation_axis(0, 1, 0);

    for (std::size_t i = 0; i < num_slices; i++)
    {
      const dolfin::Point sliced = axis.rotate(slice_rotation_axis, (i+1)*DOLFIN_PI/(num_slices+1));
      for (std::size_t j = 0; j < num_sectors; j++)
      {
        const dolfin::Point direction = sliced.rotate(axis, j*2.0*DOLFIN_PI/num_sectors);
        const dolfin::Point v = _sphere.c + direction*_sphere.r;
        add_vertex(builder, Exact_Point_3 (v.x(), v.y(), v.z()));
      }
    }

    // Add bottom has index num_vertices-1, top has index num_vertices-2
    add_vertex(builder, Exact_Point_3(top.x(), top.y(), top.z()));
    add_vertex(builder, Exact_Point_3(bottom.x(), bottom.y(), bottom.z()));

    // Add the side facets
    for (std::size_t i = 0; i < num_slices-1; i++)
    {
      for (std::size_t j = 0; j < num_sectors; j++)
      {
        const std::size_t offset1 = i*num_sectors;
        const std::size_t offset2 = (i+1)*num_sectors;

        {
          std::vector<int> f;
          f.push_back(offset1 + j);
          f.push_back(offset1 + (j+1)%num_sectors);
          f.push_back(offset2 + j);
          add_facet(builder, f);
        }

        {
          std::vector<int> f;
          f.push_back(offset2 + (j+1)%num_sectors);
          f.push_back(offset2 + j);
          f.push_back(offset1 + (j+1)%num_sectors);
          add_facet(builder, f);
        }
      }
    }

    // Add the top and bottom facets
    const std::size_t top_offset = num_sectors*(num_slices-1);
    for (std::size_t i = 0; i < num_sectors; i++)
    {
      {
        // Bottom facet
        std::vector<int> f;
        f.push_back( num_vertices-2 );
        f.push_back( (i+1)%num_sectors );
        f.push_back(i);
        add_facet(builder, f);
      }

      {
        // Top facet
        std::vector<int> f;
        f.push_back( num_vertices-1 );
        f.push_back( top_offset + (i%num_sectors) );
        f.push_back( top_offset + (i+1)%num_sectors );
        add_facet(builder, f);
      }
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

    const double x0 = std::min(_box->_x0, _box->_y0);
    const double y0 = std::max(_box->_x0, _box->_y0);

    const double x1 = std::min(_box->_x1, _box->_y1);
    const double y1 = std::max(_box->_x1, _box->_y1);

    const double x2 = std::min(_box->_x2, _box->_y2);
    const double y2 = std::max(_box->_x2, _box->_y2);

    add_vertex(builder, Exact_Point_3(y0, x1, x2));
    add_vertex(builder, Exact_Point_3(x0, x1, y2));
    add_vertex(builder, Exact_Point_3(x0, x1, x2));
    add_vertex(builder, Exact_Point_3(x0, y1, x2));
    add_vertex(builder, Exact_Point_3(y0, x1, y2));
    add_vertex(builder, Exact_Point_3(x0, y1, y2));
    add_vertex(builder, Exact_Point_3(y0, y1, x2));
    add_vertex(builder, Exact_Point_3(y0, y1, y2));

    {
      std::vector<int> f;
      f.push_back(1);
      f.push_back(2);
      f.push_back(3);
      add_facet(builder, f);
    }

    {
      std::vector<int> f;
      f.push_back(1);
      f.push_back(3);
      f.push_back(5);
      add_facet(builder, f);
    }

    {
      std::vector<int> f;
      f.push_back(1);
      f.push_back(5);
      f.push_back(4);
      add_facet(builder, f);
    }

    {
      std::vector<int> f;
      f.push_back(4);
      f.push_back(5);
      f.push_back(7);
      add_facet(builder, f);
    }

    {
      std::vector<int> f;
      f.push_back(4);
      f.push_back(7);
      f.push_back(0);
      add_facet(builder, f);
    }

    {
      std::vector<int> f;
      f.push_back(0);
      f.push_back(7);
      f.push_back(6);
      add_facet(builder, f);
    }

    {
      std::vector<int> f;
      f.push_back(0);
      f.push_back(6);
      f.push_back(2);
      add_facet(builder, f);
    }

    {
      std::vector<int> f;
      f.push_back(2);
      f.push_back(6);
      f.push_back(3);
      add_facet(builder, f);
    }

    {
      std::vector<int> f;
      f.push_back(7);
      f.push_back(5);
      f.push_back(6);
      add_facet(builder, f);
    }

    {
      std::vector<int> f;
      f.push_back(6);
      f.push_back(5);
      f.push_back(3);
      add_facet(builder, f);
    }

    {
      std::vector<int> f;
      f.push_back(1);
      f.push_back(4);
      f.push_back(2);
      add_facet(builder, f);
    }

    {
      std::vector<int> f;
      f.push_back(2);
      f.push_back(4);
      f.push_back(0);
      add_facet(builder, f);
    }

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
  P.make_tetrahedron(Exact_Point_3(b->_x0.x(), b->_x0.y(), b->_x0.z()),
                     Exact_Point_3(b->_x1.x(), b->_x1.y(), b->_x1.z()),
                     Exact_Point_3(b->_x2.x(), b->_x2.y(), b->_x2.z()),
                     Exact_Point_3(b->_x3.x(), b->_x3.y(), b->_x3.z()));
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
class Build_cone : public CGAL::Modifier_base<Exact_HalfedgeDS>
{
 public:
  Build_cone(const Cone* cone) : _cone(cone) {}

  void operator()(Exact_HalfedgeDS& hds)
  {
    const dolfin::Point axis = (_cone->_top - _cone->_bottom)/(_cone->_top - _cone->_bottom).norm();
    dolfin::Point initial = generate_orthogonal(axis);

    CGAL::Polyhedron_incremental_builder_3<Exact_HalfedgeDS> builder(hds, true);

    const int num_sides = _cone->_slices;
    const bool top_degenerate = dolfin::near(_cone->_top_radius, 0.0);
    const bool bottom_degenerate = dolfin::near(_cone->_bottom_radius, 0.0);

    const int num_vertices = (top_degenerate || bottom_degenerate) ? num_sides+2 : num_sides*2+2;

    builder.begin_surface(num_vertices, num_sides*4);

    const double delta_theta = 2.0 * DOLFIN_PI / num_sides;
    for (int i = 0; i < num_sides; ++i)
    {
      const double theta = i*delta_theta;
      const dolfin::Point rotated = initial.rotate(axis, theta);
      if (!bottom_degenerate)
      {
        const dolfin::Point p = _cone->_bottom + rotated*_cone->_bottom_radius;
        const Exact_Point_3 p_(p.x(), p.y(), p.z());
        add_vertex(builder, p_);
      }
      if (!top_degenerate)
      {
        const dolfin::Point p = _cone->_top + rotated*_cone->_top_radius;
        const Exact_Point_3 p_(p.x(), p.y(), p.z());
        add_vertex(builder, p_);
      }
    }

    // The top and bottom vertices
    add_vertex(builder, Exact_Point_3(_cone->_bottom.x(), _cone->_bottom.y(),
                                           _cone->_bottom.z()));
    add_vertex(builder, Exact_Point_3(_cone->_top.x(), _cone->_top.y(),
                                           _cone->_top.z()));

    // bottom vertex has index num_vertices-2, top vertex has index num_vertices-1

    // Construct the facets on the side.
    // Vertices must be sorted counter clockwise seen from inside.
    for (int i = 0; i < num_sides; ++i)
    {
      if (top_degenerate)
      {
        std::vector<int> f;
        f.push_back((i + 1)%num_sides);
        f.push_back(i);
        f.push_back(num_vertices - 1);
        add_facet(builder, f);
      }
      else if (bottom_degenerate)
      {
        std::vector<int> f;
        f.push_back( (i) );
        f.push_back( (i + 1) % num_sides);
        f.push_back(num_vertices - 1);
        add_facet(builder, f);
      }
      else
      {
        //Draw the sides as triangles.
        const int vertex_offset = i*2;

        // First triangle
        std::vector<int> f;
        f.push_back(vertex_offset);
        f.push_back(vertex_offset + 1);
        f.push_back((vertex_offset + 2) % (num_sides*2));
        add_facet(builder, f);

        // Second triangle
        std::vector<int> g;
        g.push_back((vertex_offset + 3) % (num_sides*2));
        g.push_back((vertex_offset + 2) % (num_sides*2));
        g.push_back(vertex_offset + 1);
        add_facet(builder, g);
      }
    }

    // Construct the bottom facet.
    if (!bottom_degenerate)
    {
      for (int i = num_sides-1; i >= 0; i -= 1)
      {
        std::vector<int> f;
        if (!top_degenerate)
        {
          f.push_back(num_vertices-2);
          f.push_back( i*2);
          f.push_back( ( (i+1)*2) % (num_sides*2));
        }
        else
        {
          f.push_back(num_vertices-2);
          f.push_back(i);
          f.push_back( (i+1)%num_sides );
        }
        add_facet(builder, f);
      }
    }

    // Construct the the top facet
    if (!top_degenerate)
    {
      for (int i = 0; i < num_sides; i++)
      {
        if (!bottom_degenerate)
        {
          std::vector<int> f;
          f.push_back(num_vertices-1);
          f.push_back( ( (i+1)*2)%(num_sides*2) +1 );
          f.push_back( i*2 + 1 );
          add_facet(builder, f);
        }
        else
        {
          std::vector<int> f;
          f.push_back(num_vertices-2);
          f.push_back( (i+1)%num_sides);
          f.push_back(i);
          add_facet(builder, f);
        }
      }
    }

    builder.end_surface();
  }
private:
  const Cone* _cone;
};
//-----------------------------------------------------------------------------
 void make_cone(const Cone* c, Exact_Polyhedron_3& P)
{
  Build_cone builder(c);
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
                     const std::vector<std::array<std::size_t, 3> >& facets)
    : vertices(vertices), facets(facets){}
  void operator()(HDS& hds)
  {
    CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);

    builder.begin_surface(vertices.size(), facets.size());
    
    for (std::vector<std::array<double, 3> >::const_iterator it = vertices.begin();
         it != vertices.end(); ++it)
      builder.add_vertex(Exact_Point_3( (*it)[0], (*it)[1], (*it)[2]));

    for (std::vector<std::array<std::size_t, 3> >::const_iterator it = facets.begin();
         it != facets.end(); ++it)
      builder.add_facet(it->begin(), it->end());

    builder.end_surface();

  }
  const std::vector<std::array<double, 3> > vertices;
  const std::vector<std::array<std::size_t, 3> > facets;
};


 void make_surface3D(const Surface3D* s, Exact_Polyhedron_3& P)
{
  dolfin_assert(s);

  std::vector<std::array<double, 3> > vertices;
  std::vector<std::array<std::size_t, 3> > facets;

  boost::filesystem::path fpath(s->_filename);
  if (fpath.extension() == ".stl")
  {
    STLFileReader::read(s->_filename, vertices, facets);
  }
  else if (fpath.extension() == ".vtp")
  {
    // TODO: Only if vtk is installed
    VTPFileReader::read(s->_filename, vertices, facets);
  }
  else if(fpath.extension() == ".off")
  {
    // TODO: Let cgal parse the file
  }
  else
  {
    dolfin::dolfin_error("PolyhedronUtils.cpp",
                         "open file to read 3D surface",
                         "Unknown file type");
  }

  // Create the polyhedron
  BuildFromFacetList<Exact_HalfedgeDS> builder(vertices, facets);
  P.delegate(builder);
}
//-----------------------------------------------------------------------------
boost::shared_ptr<Nef_polyhedron_3>
convertSubTree(const CSGGeometry *geometry)
{
  switch (geometry->getType())
  {
    case CSGGeometry::Union :
    {
      const CSGUnion* u = dynamic_cast<const CSGUnion*>(geometry);
      dolfin_assert(u);
      boost::shared_ptr<Nef_polyhedron_3> g0 = convertSubTree(u->_g0.get());
      boost::shared_ptr<Nef_polyhedron_3> g1 = convertSubTree(u->_g1.get());
      (*g0) += (*g1);
      return g0;

      break;
    }
    case CSGGeometry::Intersection :
    {
      const CSGIntersection* u = dynamic_cast<const CSGIntersection*>(geometry);
      dolfin_assert(u);
      boost::shared_ptr<Nef_polyhedron_3> g0 = convertSubTree(u->_g0.get());
      boost::shared_ptr<Nef_polyhedron_3> g1 = convertSubTree(u->_g1.get());
      (*g0) *= (*g1);
      return g0;
      break;
    }
    case CSGGeometry::Difference :
    {
      const CSGDifference* u = dynamic_cast<const CSGDifference*>(geometry);
      dolfin_assert(u);
      boost::shared_ptr<Nef_polyhedron_3> g0 = convertSubTree(u->_g0.get());
      boost::shared_ptr<Nef_polyhedron_3> g1 = convertSubTree(u->_g1.get());
      (*g0) -= (*g1);
      return g0;
      break;
    }
    case CSGGeometry::Cone :
    {
      const Cone* c = dynamic_cast<const Cone*>(geometry);
      dolfin_assert(c);
      Exact_Polyhedron_3 P;
      make_cone(c, P);
      return boost::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
      break;
    }
    case CSGGeometry::Sphere :
    {
      const Sphere* s = dynamic_cast<const Sphere*>(geometry);
      dolfin_assert(s);
      Exact_Polyhedron_3 P;
      make_sphere(s, P);
      return boost::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
      break;
    }
    case CSGGeometry::Box :
    {
      const Box* b = dynamic_cast<const Box*>(geometry);
      dolfin_assert(b);
      Exact_Polyhedron_3 P;
      make_box(b, P);
      return boost::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
      break;
    }

    case CSGGeometry::Tetrahedron :
    {
      const Tetrahedron* b = dynamic_cast<const Tetrahedron*>(geometry);
      dolfin_assert(b);
      Exact_Polyhedron_3 P;
      make_tetrahedron(b, P);
      return boost::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
      break;
    }
    case CSGGeometry::Surface3D :
    {
      const Surface3D* b = dynamic_cast<const Surface3D*>(geometry);
      dolfin_assert(b);
      Exact_Polyhedron_3 P;
      make_surface3D(b, P);
      return boost::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3(P));
      break;
    }
    default:
      dolfin::dolfin_error("GeometryToCGALConverter.cpp",
                           "converting geometry to cgal polyhedron",
                           "Unhandled primitive type");
  }

  // Make compiler happy.
  return boost::shared_ptr<Nef_polyhedron_3>(new Nef_polyhedron_3);
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

    case CSGGeometry::Cone :
    {
      const Cone* c = dynamic_cast<const Cone*>(&geometry);
      dolfin_assert(c);
      make_cone(c, P);
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
    dolfin::cout << "Convert to nef polyhedron" << dolfin::endl;
    boost::shared_ptr<Nef_polyhedron_3> cgal_geometry
      = convertSubTree(&geometry);
    dolfin_assert(cgal_geometry->is_valid());
    dolfin_assert(cgal_geometry->is_simple());
    cgal_geometry->convert_to_polyhedron(P);
  }

  dolfin::cout << "Number of vertices: " << P.size_of_vertices() << dolfin::endl;
  dolfin::cout << "Number of facets:   " << P.size_of_facets() << dolfin::endl;
}

} //end unnamed namespace


namespace mshr
{

struct CSGCGALDomain3DImpl
{
  Exact_Polyhedron_3 p;
};


CSGCGALDomain3D::CSGCGALDomain3D()
: impl(new CSGCGALDomain3DImpl)
{}
//-----------------------------------------------------------------------------
CSGCGALDomain3D::CSGCGALDomain3D(const mshr::CSGGeometry &csg)
: impl(new CSGCGALDomain3DImpl)
{
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
void CSGCGALDomain3D::remove_degenerated_facets(double threshold) 
{
  // int degenerate_facets = number_of_degenerate_facets(impl->p, threshold);

  // dolfin::cout << "Number of degenerate facets: " << degenerate_facets << dolfin::endl;

  // // FIXME: Use has_degenerate_facets() when in production code
  // if (degenerate_facets > 0)
  // {
  //   dolfin_assert(p.is_pure_triangle());

  //   shortest_edge(p);

  //   dolfin::cout << "Removing triangles with short edges" << dolfin::endl;
  //   remove_short_edges(p, threshold);

  //   dolfin::cout << "Number of degenerate facets: "
  //        << number_of_degenerate_facets(p, threshold) << dolfin::endl;

  //   dolfin::cout << "Removing small triangles" << dolfin::endl;
  //   remove_small_triangles(p, threshold);

  //   dolfin::cout << "Number of degenerate facets: "
  //        << number_of_degenerate_facets(p, threshold) << dolfin::endl;

  //   // Removal of triangles should preserve the triangular structure
  //   // of the polyhedron
  //   dolfin_assert(p.is_pure_triangle());
  // }

}
//-----------------------------------------------------------------------------

} // end namespace mshr

