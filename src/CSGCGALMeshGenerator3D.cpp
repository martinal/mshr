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


#include <mshr/CSGCGALMeshGenerator3D.h>
#include <mshr/CSGGeometry.h>
#include <mshr/CSGCGALDomain3D.h>

#include <dolfin/common/MPI.h>
#include <dolfin/log/LogStream.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshPartitioning.h>

#include <memory>

#define CGAL_NO_DEPRECATED_CODE
#define CGAL_MESH_3_VERBOSE
//#define PROTECTION_DEBUG

#define CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
#define CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include "Polyhedral_multicomponent_mesh_domain_with_features_3.h"
#include "make_multicomponent_mesh_3.h"

// Bounding sphere computation
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>

//#define NO_MULTICOMPONENT_DOMAIN

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type MeshPolyhedron_3;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef K::Triangle_3 Triangle_3;

#ifdef NO_MULTICOMPONENT_DOMAIN
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> PolyhedralMeshDomain;
#else
typedef Polyhedral_multicomponent_mesh_domain_with_features_3<K> PolyhedralMeshDomain;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<PolyhedralMeshDomain>::type Tr;
// typedef CGAL::Mesh_complex_3_in_triangulation_3<
//  Tr,PolyhedralMeshDomain::Corner_index,PolyhedralMeshDomain::Curve_segment_index> C3t3;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

namespace
{


//-----------------------------------------------------------------------------
void build_dolfin_mesh(const C3t3& c3t3, dolfin::Mesh& mesh)
{
  typedef C3t3::Triangulation Triangulation;
  typedef Triangulation::Vertex_handle Vertex_handle;

  // Collect the vertices that are part of the complex
  std::map<Vertex_handle, std::size_t> vertex_id_map;
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();
      cit != c3t3.cells_in_complex_end();
      ++cit)
  {
    for (std::size_t i = 0; i < 4; i++)
    {
      if (!vertex_id_map.count(cit->vertex(i)))
      {
        const std::size_t vertex_index = vertex_id_map.size();
        vertex_id_map[cit->vertex(i)] = vertex_index;
      }
    }
  }

  // Create and initialize mesh editor
  mesh.clear();
  dolfin::MeshEditor mesh_editor;
  mesh_editor.open(mesh, 3, 3);
  mesh_editor.init_vertices(vertex_id_map.size());
  mesh_editor.init_cells(c3t3.number_of_cells_in_complex());

  // Add vertices to mesh
  for (std::map<Vertex_handle, std::size_t>::const_iterator it = vertex_id_map.cbegin();
       it != vertex_id_map.cend(); it++)
  {
    // Get vertex coordinates and add vertex to the mesh
    dolfin::Point p(it->first->point()[0], it->first->point()[1], it->first->point()[2]);
    mesh_editor.add_vertex(it->second, p);
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
  Copy_polyhedron_to(const mshr::CSGCGALDomain3D& in_poly, bool flip)
    : _in_poly(in_poly), flip(flip) {}

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
        if (flip)
        {
          builder.add_vertex_to_facet( (*it)[2] );
          builder.add_vertex_to_facet( (*it)[1] );
        }
        else
        {
          builder.add_vertex_to_facet( (*it)[1] );
          builder.add_vertex_to_facet( (*it)[2] );
        }

        builder.end_facet();
      }
    }
    builder.end_surface();
  }
private:
  const mshr::CSGCGALDomain3D& _in_poly;
  const bool flip;
}; 

static void convert_to_inexact(const mshr::CSGCGALDomain3D &exact_domain,
                               MeshPolyhedron_3 &inexact_domain,
                               bool flip)
{
  Copy_polyhedron_to modifier(exact_domain, flip);
  inexact_domain.delegate(modifier);
  CGAL_assertion(inexact_domain.is_valid());
}
//-----------------------------------------------------------------------------
double get_bounding_sphere_radius(const MeshPolyhedron_3& polyhedron)
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
} //end anonymous namespace
//-----------------------------------------------------------------------------
namespace mshr
{
CSGCGALMeshGenerator3D::CSGCGALMeshGenerator3D()
{
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
CSGCGALMeshGenerator3D::~CSGCGALMeshGenerator3D() {}
//-----------------------------------------------------------------------------
void CSGCGALMeshGenerator3D::generate(std::shared_ptr<const CSGCGALDomain3D> csgdomain,
                                      dolfin::Mesh& mesh) const
{

  // Note that if not in parallel (ie. size() == 0)
  // then both receiver and broadcaster will return false
  if (!dolfin::MPI::is_receiver(mesh.mpi_comm()))
  {
    // Create CGAL mesh domain
    MeshPolyhedron_3 p;
    convert_to_inexact(*csgdomain, p, !csgdomain->is_insideout());

    // Reset the (memory consuming exact arithmetic) domain object
    // will be deleted if the pointer has been moved to the function
    csgdomain.reset();

    // Workaround, cgal segfaulted when assigning new mesh criterias
    // within the if-else blocks.
    std::unique_ptr<Mesh_criteria> criteria;
    double edge_size;

    const bool criteria_changed = parameters["edge_size"].change_count() > 0
      || parameters["facet_angle"].change_count() > 0
      || parameters["facet_size"].change_count() > 0
      || parameters["facet_distance"].change_count() > 0
      || parameters["cell_radius_edge_ratio"].change_count() > 0
      || parameters["cell_size"].change_count() > 0;

    if (parameters["mesh_resolution"].change_count() > 0 && criteria_changed)
      dolfin::warning("Attempt to set both mesh_resolution and other meshing criterias which are mutually exclusive");

    if (criteria_changed)
    {
      log(dolfin::TRACE, "Using user specified meshing criterias");

      // Mesh criteria
      criteria.reset(new Mesh_criteria(CGAL::parameters::edge_size = parameters["edge_size"],
                                       CGAL::parameters::facet_angle = parameters["facet_angle"],
                                       CGAL::parameters::facet_size = parameters["facet_size"], //  <-----------
                                       CGAL::parameters::facet_distance = parameters["facet_distance"],
                                       CGAL::parameters::cell_radius_edge_ratio = parameters["cell_radius_edge_ratio"],
                                       CGAL::parameters::cell_size = parameters["cell_size"])); // <--------------
      edge_size = parameters["edge_size"];
    }
    else
    {
      // Try to compute reasonable parameters
      const double mesh_resolution = parameters["mesh_resolution"];
      const double r = get_bounding_sphere_radius(p);
      const double cell_size = r/mesh_resolution*2.0;
      log(dolfin::TRACE, "Computing meshing criterias. Chose cell size %f", cell_size);

      criteria.reset(new Mesh_criteria(CGAL::parameters::edge_size = cell_size,
                                       CGAL::parameters::facet_angle = 30.0,
                                       CGAL::parameters::facet_size = cell_size,
                                       CGAL::parameters::facet_distance = cell_size/10.0, // ???
                                       CGAL::parameters::cell_radius_edge_ratio = 3.0,
                                       CGAL::parameters::cell_size = cell_size));
      edge_size = cell_size;
    }

    #ifdef NO_MULTICOMPONENT_DOMAIN
    PolyhedralMeshDomain domain(p);
    #else
    PolyhedralMeshDomain domain(p, edge_size);
    #endif

    if (parameters["detect_sharp_features"])
    {
      log(dolfin::TRACE, "Detecting sharp features");

      //const int feature_threshold = parameters["feature_threshold"];
      domain.detect_features();
    }

    // Mesh generation
    log(dolfin::TRACE, "Generating mesh");
    C3t3 c3t3;
    make_multicomponent_mesh_3_impl<C3t3>(c3t3,
                                          domain,
                                          *criteria,
                                          CGAL::parameters::no_exude(),
                                          CGAL::parameters::no_perturb(),
                                          CGAL::parameters::no_odt(),
                                          CGAL::parameters::no_lloyd(),
                                          true);



    if (parameters["odt_optimize"])
    {
      log(dolfin::TRACE, "Optimizing mesh by odt optimization");
      odt_optimize_mesh_3(c3t3, domain);
    }

    if (parameters["lloyd_optimize"])
    {
      log(dolfin::TRACE, "Optimizing mesh by lloyd optimization");
      lloyd_optimize_mesh_3(c3t3, domain);
    }

    if (parameters["perturb_optimize"])
    {
      log(dolfin::TRACE, "Optimizing mesh by perturbation");
      // TODO: Set time limit
      CGAL::perturb_mesh_3(c3t3, domain);
    }

    if (parameters["exude_optimize"])
    {
      log(dolfin::TRACE, "Optimizing mesh by sliver exudation");
      exude_mesh_3(c3t3);
    }

    build_dolfin_mesh(c3t3, mesh);
  }

  // Distribute the mesh (if in parallel)
  dolfin::MeshPartitioning::build_distributed_mesh(mesh);
}
//-----------------------------------------------------------------------------
std::shared_ptr<dolfin::Mesh>
CSGCGALMeshGenerator3D::generate(std::shared_ptr<const CSGCGALDomain3D> csgdomain) const
{
  std::shared_ptr<dolfin::Mesh> mesh(new dolfin::Mesh());
  generate(csgdomain, *mesh);
  return mesh;
}
}
