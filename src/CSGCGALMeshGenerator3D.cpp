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
//
// Modified by Joachim B Haga 2012


#include <mshr/CSGCGALMeshGenerator3D.h>
#include <mshr/CSGGeometry.h>
#include <mshr/CSGCGALDomain3D.h>

#include <dolfin/log/LogStream.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <boost/scoped_ptr.hpp>


#define CGAL_NO_DEPRECATED_CODE
#define CGAL_MESH_3_VERBOSE
//#define PROTECTION_DEBUG

#define CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
#define CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

// Bounding sphere computation
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>

namespace mshr
{

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type MeshPolyhedron_3;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef K::Triangle_3 Triangle_3;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> PolyhedralMeshDomain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<PolyhedralMeshDomain>::type Tr;
// typedef CGAL::Mesh_complex_3_in_triangulation_3<
//  Tr,PolyhedralMeshDomain::Corner_index,PolyhedralMeshDomain::Curve_segment_index> C3t3;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;


// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
//-----------------------------------------------------------------------------
/*static*/ void build_dolfin_mesh(const C3t3& c3t3, dolfin::Mesh& mesh)
{
  typedef C3t3::Triangulation Triangulation;
  typedef Triangulation::Vertex_handle Vertex_handle;

  // CGAL triangulation
  const Triangulation& triangulation = c3t3.triangulation();

  // Clear mesh
  mesh.clear();

  // Count cells in complex
  std::size_t num_cells = 0;
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();
      cit != c3t3.cells_in_complex_end();
      ++cit)
  {
    num_cells++;
  }

  // Create and initialize mesh editor
  dolfin::MeshEditor mesh_editor;
  mesh_editor.open(mesh, 3, 3);
  mesh_editor.init_vertices(triangulation.number_of_vertices());
  mesh_editor.init_cells(num_cells);

  // Add vertices to mesh
  std::size_t vertex_index = 0;
  std::map<Vertex_handle, std::size_t> vertex_id_map;

  for (Triangulation::Finite_vertices_iterator
         cgal_vertex = triangulation.finite_vertices_begin();
       cgal_vertex != triangulation.finite_vertices_end(); ++cgal_vertex)
  {
    vertex_id_map[cgal_vertex] = vertex_index;

      // Get vertex coordinates and add vertex to the mesh
    dolfin::Point p(cgal_vertex->point()[0], cgal_vertex->point()[1], cgal_vertex->point()[2]);
    mesh_editor.add_vertex(vertex_index, p);

    ++vertex_index;
  }

  // Add cells to mesh
  std::size_t cell_index = 0;
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();
      cit != c3t3.cells_in_complex_end();
      ++cit)
  {
    mesh_editor.add_cell(cell_index,
                         vertex_id_map[cit->vertex(0)],
                         vertex_id_map[cit->vertex(1)],
                         vertex_id_map[cit->vertex(2)],
                         vertex_id_map[cit->vertex(3)]);

    ++cell_index;
  }

  // Close mesh editor
  mesh_editor.close();
}
//-----------------------------------------------------------------------------
struct Copy_polyhedron_to
  : public CGAL::Modifier_base<typename MeshPolyhedron_3::HalfedgeDS>
{
  Copy_polyhedron_to(const CSGCGALDomain3D& in_poly)
    : _in_poly(in_poly) {}

  void operator()(typename MeshPolyhedron_3::HalfedgeDS& out_hds)
  {
    typedef typename MeshPolyhedron_3::HalfedgeDS Output_HDS;
    CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

    builder.begin_surface(_in_poly.num_vertices(),
                          _in_poly.num_facets(),
                          _in_poly.num_halfedges());

    {
      std::vector<dolfin::Point> v;
      _in_poly.get_vertices(v);

      for(std::vector<dolfin::Point>::iterator it = v.begin();
          it != v.end(); it++)
      {
        typename MeshPolyhedron_3::Point_3 p( (*it)[0], (*it)[1], (*it)[2] );
        builder.add_vertex(p);
      }
    }

    {
      std::vector<std::array<std::size_t, 3> > f;
      _in_poly.get_facets(f);

      for (std::vector<std::array<std::size_t, 3> >::iterator it = f.begin();
           it != f.end(); it++)
      {
        builder.begin_facet ();
        builder.add_vertex_to_facet( (*it)[0] );
        builder.add_vertex_to_facet( (*it)[1] );
        builder.add_vertex_to_facet( (*it)[2] );
        builder.end_facet();
      }

    }
    builder.end_surface();
  }
private:
  const CSGCGALDomain3D& _in_poly;
}; 

static void convert_to_inexact(const CSGCGALDomain3D &exact_domain, 
                               MeshPolyhedron_3 &inexact_domain)
{
  Copy_polyhedron_to modifier(exact_domain);
  inexact_domain.delegate(modifier);
  CGAL_assertion(inexact_domain.is_valid());
}
//-----------------------------------------------------------------------------
/*static*/ double get_bounding_sphere_radius(const MeshPolyhedron_3& polyhedron)
{
  typedef CGAL::Min_sphere_of_spheres_d_traits_3<K, K::FT> MinSphereTraits;
  typedef CGAL::Min_sphere_of_spheres_d<MinSphereTraits> Min_sphere;
  typedef MinSphereTraits::Sphere Sphere;

  std::vector<Sphere> S;

  for (MeshPolyhedron_3::Vertex_const_iterator it=polyhedron.vertices_begin();
       it != polyhedron.vertices_end(); ++it)
  {
    S.push_back(Sphere(it->point(), 0.0));
  }

  Min_sphere ms(S.begin(), S.end());
  CGAL_assertion(ms.is_valid());

  return CGAL::to_double(ms.radius());
}
//-----------------------------------------------------------------------------
CSGCGALMeshGenerator3D::CSGCGALMeshGenerator3D(const CSGGeometry& geometry)
{
  boost::shared_ptr<const CSGGeometry> tmp = dolfin::reference_to_no_delete_pointer<const CSGGeometry>(geometry);
  _geometry = tmp;
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
CSGCGALMeshGenerator3D::CSGCGALMeshGenerator3D(boost::shared_ptr<const CSGGeometry> geometry)
  : _geometry(geometry)
{
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
CSGCGALMeshGenerator3D::~CSGCGALMeshGenerator3D() {}
//-----------------------------------------------------------------------------
void CSGCGALMeshGenerator3D::generate(dolfin::Mesh& mesh) const
{
  CSGCGALDomain3D exact_domain(*_geometry);
  if (parameters["remove_degenerated"])
    // TODO: Make the threshold a parameter
    exact_domain.remove_degenerated_facets(1e-6);

  // Create CGAL mesh domain
  MeshPolyhedron_3 p;
  convert_to_inexact(exact_domain, p);

  PolyhedralMeshDomain domain(p);

  if (parameters["detect_sharp_features"])
  {
    dolfin::cout << "Detecting sharp features" << dolfin::endl;
    //const int feature_threshold = parameters["feature_threshold"];
    domain.detect_features();
  }

  // Workaround, cgal segfaulted when assigning new mesh criterias
  // within the if-else blocks.
  boost::scoped_ptr<Mesh_criteria> criteria;

  int mesh_resolution = parameters["mesh_resolution"];
  if (mesh_resolution > 0)
  {
    // Try to compute reasonable parameters
    const double r = get_bounding_sphere_radius(p);
    //dolfin::cout << "Radius of bounding sphere: " << r << dolfin::endl;
    //dolfin::cout << "Mesh resolution" << mesh_resolution << dolfin::endl;
    const double cell_size = r/static_cast<double>(mesh_resolution)*2.0;
    // dolfin::cout << "Cell size: " << cell_size << dolfin::endl;

    criteria.reset(new Mesh_criteria(CGAL::parameters::edge_size = cell_size/5, // ???
                                          CGAL::parameters::facet_angle = 30.0,
                                          CGAL::parameters::facet_size = cell_size,
                                          CGAL::parameters::facet_distance = cell_size/10.0, // ???
                                          CGAL::parameters::cell_radius_edge_ratio = 3.0,
                                          CGAL::parameters::cell_size = cell_size));
  }
  else
  {
    // Mesh criteria
    criteria.reset(new Mesh_criteria(CGAL::parameters::edge_size = parameters["edge_size"],
                                          CGAL::parameters::facet_angle = parameters["facet_angle"],
                                          CGAL::parameters::facet_size = parameters["facet_size"],
                                          CGAL::parameters::facet_distance = parameters["facet_distance"],
                                          CGAL::parameters::cell_radius_edge_ratio = parameters["cell_radius_edge_ratio"],
                                          CGAL::parameters::cell_size = parameters["cell_size"]));
  }

  // Mesh generation
  dolfin::cout << "Generating mesh" << dolfin::endl;
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, *criteria,
                                                CGAL::parameters::no_perturb(),
                                                CGAL::parameters::no_exude());

  if (parameters["odt_optimize"])
  {
    dolfin::cout << "Optimizing mesh by odt optimization" << dolfin::endl;
    odt_optimize_mesh_3(c3t3, domain);
  }

  if (parameters["lloyd_optimize"])
  {
    dolfin::cout << "Optimizing mesh by lloyd optimization" << dolfin::endl;
    lloyd_optimize_mesh_3(c3t3, domain);
  }

  if (parameters["perturb_optimize"])
  {
    dolfin::cout << "Optimizing mesh by perturbation" << dolfin::endl;
    // TODO: Set time limit
    CGAL::perturb_mesh_3(c3t3, domain);
  }

  if (parameters["exude_optimize"])
  {
    dolfin::cout << "Optimizing mesh by sliver exudation" << dolfin::endl;
    exude_mesh_3(c3t3);
  }

  build_dolfin_mesh(c3t3, mesh);
}
//-----------------------------------------------------------------------------
void CSGCGALMeshGenerator3D::save_off(std::string filename) const
{
  CSGCGALDomain3D exact_domain(*_geometry);
  if (parameters["remove_degenerated"])
    // TODO: Make the threshold a parameter
    exact_domain.remove_degenerated_facets(1e-6);

  // Create CGAL mesh domain
  MeshPolyhedron_3 p;
  convert_to_inexact(exact_domain, p);

  dolfin::cout << "Writing to file " << filename << dolfin::endl;
  std::ofstream outfile(filename.c_str());

  outfile << p;
  outfile.close();
}
}
