// Copyright (C) 2012 Johannes Ring, 2012-2014 Benjamin Kehlet
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


#include <vector>
#include <cmath>
#include <limits>

#include <CGAL/compiler_config.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/MPI.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshDomains.h>
#include <dolfin/mesh/MeshValueCollection.h>
#include <dolfin/mesh/MeshPartitioning.h>
#include <dolfin/log/log.h>

#include <mshr/CSGCGALMeshGenerator2D.h>
#include <mshr/CSGGeometry.h>
#include <mshr/CSGOperators.h>
#include <mshr/CSGPrimitives2D.h>
#include <mshr/CSGCGALDomain2D.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Inexact_Kernel;
//typedef CGAL::Exact_predicates_exact_constructions_kernel Inexact_Kernel;

typedef CGAL::Triangulation_vertex_base_2<Inexact_Kernel>  Vertex_base;
typedef CGAL::Delaunay_mesh_face_base_2<Inexact_Kernel> Face_base;

template <class Gt, class Fb >
class Enriched_face_base_2 : public Fb
{
 public:
  typedef Gt Geom_traits;
  typedef typename Fb::Vertex_handle Vertex_handle;
  typedef typename Fb::Face_handle Face_handle;

  template <typename TDS2>
  struct Rebind_TDS
  {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other Fb2;
    typedef Enriched_face_base_2<Gt,Fb2> Other;
  };

protected:
  int c;

public:
  Enriched_face_base_2(): Fb(), c(-1) {}

  Enriched_face_base_2(Vertex_handle v0,
                       Vertex_handle v1,
                       Vertex_handle v2)
    : Fb(v0,v1,v2), c(-1) {}

  Enriched_face_base_2(Vertex_handle v0,
                       Vertex_handle v1,
                       Vertex_handle v2,
                       Face_handle n0,
                       Face_handle n1,
                       Face_handle n2)
    : Fb(v0,v1,v2,n0,n1,n2), c(-1) {}

  inline void set_counter(int i) { c = i; }
  inline int counter() const { return c; }
};

typedef CGAL::Triangulation_vertex_base_2<Inexact_Kernel> Vb;
typedef CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Inexact_Kernel, Vb> Vbb;
typedef Enriched_face_base_2<Inexact_Kernel, Face_base> Fb;
typedef CGAL::Triangulation_data_structure_2<Vbb, Fb> TDS;
// typedef CGAL::Triangulation_data_structure_2<Vbb, Face_base> TDS;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Inexact_Kernel, TDS, Itag> CDT;
///typedef CGAL::Constrained_Delaunay_triangulation_2<Inexact_Kernel, TDS> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Mesh_criteria_2;
typedef CGAL::Delaunay_mesher_2<CDT, Mesh_criteria_2> CGAL_Mesher_2;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Face_handle Face_handle;
typedef CDT::All_faces_iterator All_faces_iterator;

typedef Inexact_Kernel::Point_2 Point_2;
typedef Inexact_Kernel::Segment_2 Segment_2;

namespace mshr
{

//-----------------------------------------------------------------------------
CSGCGALMeshGenerator2D::CSGCGALMeshGenerator2D()
{
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
CSGCGALMeshGenerator2D::~CSGCGALMeshGenerator2D() {}
//-----------------------------------------------------------------------------
void explore_subdomain(CDT &ct,
                       CDT::Face_handle start,
                       std::list<CDT::Face_handle>& other_domains)
{
  std::list<Face_handle> queue;
  queue.push_back(start);

  while (!queue.empty())
  {
    CDT::Face_handle face = queue.front();
    queue.pop_front();

    for(int i = 0; i < 3; i++)
    {
      Face_handle n = face->neighbor(i);
      if (ct.is_infinite(n))
        continue;

      const CDT::Edge e(face,i);

      if (n->counter() == -1)
      {
        if (ct.is_constrained(e))
        {
          // Reached a border
          other_domains.push_back(n);
        } 
        else
        {
          // Set neighbor interface to the same and push to queue
          n->set_counter(face->counter());
          queue.push_back(n);
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
// Set the member in_domain and counter for all faces in the cdt
void explore_subdomains(CDT& cdt,
                        const CSGCGALDomain2D& total_domain,
                        const std::vector<std::pair<std::size_t,
                        CSGCGALDomain2D> > &subdomain_geometries)
{
  // Set counter to -1 for all faces
  for (CDT::Finite_faces_iterator it = cdt.finite_faces_begin();
       it != cdt.finite_faces_end(); ++it)
  {
    it->set_counter(-1);
  }

  std::list<CDT::Face_handle> subdomains;
  subdomains.push_back(cdt.finite_faces_begin());

  //print_face(subdomains.front());
  //dolfin_assert(face_in_domain(subdomains.front(), total_domain));

  while (!subdomains.empty())
  {
    const CDT::Face_handle f = subdomains.front();
    subdomains.pop_front();

    if (f->counter() < 0)
    {
      // Set default marker (0)
      f->set_counter(0);

      const Point_2 p0 = f->vertex(0)->point();
      const Point_2 p1 = f->vertex(1)->point();
      const Point_2 p2 = f->vertex(2)->point();

      dolfin::Point p( (p0[0] + p1[0] + p2[0]) / 3.0,
                       (p0[1] + p1[1] + p2[1]) / 3.0 );

      for (int i = subdomain_geometries.size(); i > 0; --i)
      {
        if (subdomain_geometries[i-1].second.point_in_domain(p))
        {
          f->set_counter(subdomain_geometries[i-1].first);
          break;
        }
      }

      explore_subdomain(cdt, f, subdomains);
    }
  }
}
//-----------------------------------------------------------------------------
double shortest_constrained_edge(const CDT &cdt)
{
  double min_length = std::numeric_limits<double>::max();
  for (CDT::Finite_edges_iterator it = cdt.finite_edges_begin();
       it != cdt.finite_edges_end();
       it++)
  {
    if (!cdt.is_constrained(*it))
      continue;

    // An edge is an std::pair<Face_handle, int>
    // see CGAL/Triangulation_data_structure_2.h
    CDT::Face_handle f = it->first;
    const int i = it->second;

    CDT::Point p1 = f->vertex( (i+1)%3 )->point();
    CDT::Point p2 = f->vertex( (i+2)%3 )->point();

    min_length = std::min(CGAL::to_double((p1-p2).squared_length()),
                          min_length);
  }

  return min_length;
}
//-----------------------------------------------------------------------------
void CSGCGALMeshGenerator2D::generate(const CSGGeometry& geometry, dolfin::Mesh& mesh)
{
  // Note that if not in parallel (ie. size() == 0)
  // then both receiver and broadcaster will return false
  if (!dolfin::MPI::is_receiver(mesh.mpi_comm()))
  {

    std::vector<std::pair<std::size_t, CSGCGALDomain2D> >
      subdomain_geometries;

    log(dolfin::TRACE, "Converting geometry to CGAL polygon");
    CSGCGALDomain2D total_domain(&geometry);

    // Empty polygon, will be populated when traversing the subdomains
    CSGCGALDomain2D overlaying;

    // Add the subdomains to the PSLG. Traverse in reverse order to get the latest
    // added subdomain on top
    std::list<std::pair<std::size_t, std::shared_ptr<const CSGGeometry> > >::const_reverse_iterator it;

    if (!geometry.get_subdomains().empty())
      log(dolfin::TRACE, "Processing subdomains");

    for (it = geometry.get_subdomains().rbegin(); it != geometry.get_subdomains().rend();
         ++it)
    {
      const std::size_t current_index = it->first;
      std::cout << "Handling subdomain: " << current_index << std::endl;
      std::shared_ptr<const CSGGeometry> current_subdomain = it->second;

      CSGCGALDomain2D cgal_geometry(current_subdomain.get());
      std::cout << "Current subdomain: " << cgal_geometry.str(true) << std::endl;

      // Only the part inside the total domain
      cgal_geometry.intersect_inplace(total_domain);

      std::cout << "Subdomain intersecting total: " << cgal_geometry.str(true) << std::endl;

      // Only the part outside overlaying subdomains
      cgal_geometry.difference_inplace(overlaying);

      std::cout << "Subdomain minus overlays: " << cgal_geometry.str(true) << std::endl;

      subdomain_geometries.push_back(std::make_pair(current_index,
                                                    cgal_geometry));

      overlaying.join_inplace(cgal_geometry);
    }

    CSGCGALDomain2D remaining(total_domain);
    remaining.difference_inplace(overlaying);

    std::cout << "Remaining part of total domain: " << remaining.str(true) << std::endl;

    subdomain_geometries.push_back(std::make_pair(0, remaining));

    // const double pixel_size = parameters["pixel_size"];

    // Compute cell size
    double cell_size;
    const double mesh_resolution = parameters["mesh_resolution"];
    if (mesh_resolution > 0)
    {
      const double min_radius = total_domain.compute_boundingcircle_radius();
      cell_size = 2.0*min_radius/mesh_resolution;
    }
    else
    {
      cell_size = parameters["cell_size"];
    }

    log(dolfin::TRACE, "Request cell size: %f", cell_size);
    //const double truncate_tolerance = parameters["edge_truncate_tolerance"];

    std::pair<std::vector<dolfin::Point>,
              std::vector<std::pair<std::size_t, std::size_t>>> pslg =
      CSGCGALDomain2D::compute_pslg(subdomain_geometries);
    // PSLG pslg(subdomain_geometries,
    //           pixel_size,
    //           truncate_tolerance < 0 ? cell_size/200 : truncate_tolerance);

    // Create empty CGAL triangulation and copy data from the PSLG
    CDT cdt;

    {
      std::vector<CDT::Vertex_handle> vertices;
      {
        // Insert the vertices into the triangulation
        vertices.reserve(pslg.first.size());
        for (const dolfin::Point& vertex : pslg.first)
        {
          std::cout << "Inserting vertex: " << vertex.str() << std::endl;
          vertices.push_back(cdt.insert(Point_2(vertex.x(), vertex.y())));
        }
      }

      // Insert the edges as constraints
      for (const std::pair<std::size_t, size_t>& edge : pslg.second)
      {
        std::cout << "Inserting edge: (" << edge.first << ", " << edge.second << ")" << std::endl;
        cdt.insert_constraint(vertices[edge.first], vertices[edge.second]);
      }
    }
  

    // log(dolfin::TRACE, "Exploring subdomains in triangulation");
    // explore_subdomains(cdt, total_domain, subdomain_geometries);

    log(dolfin::TRACE, "Initializing mesh refinement");

    // Create mesher
    CGAL_Mesher_2 mesher(cdt);

    // Add seeds for all faces in the total domain
    std::list<Point_2> list_of_seeds;
    for(CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
        fit != cdt.finite_faces_end(); ++fit)
    {
      // Calculate center of triangle and add to list of seeds
      Point_2 p0 = fit->vertex(0)->point();
      Point_2 p1 = fit->vertex(1)->point();
      Point_2 p2 = fit->vertex(2)->point();
      const double x = (p0[0] + p1[0] + p2[0]) / 3;
      const double y = (p0[1] + p1[1] + p2[1]) / 3;

      if (total_domain.point_in_domain(dolfin::Point(x, y)))
      {
        list_of_seeds.push_back(Point_2(x, y));
      }
    }

    // log(dolfin::TRACE, "Added %d seed points", list_of_seeds.size());
    mesher.set_seeds(list_of_seeds.begin(), list_of_seeds.end(), true);

    // Set shape and size criteria
    const double shape_bound = parameters["triangle_shape_bound"];
    mesher.set_criteria(Mesh_criteria_2(shape_bound, cell_size));

    // Refine CGAL mesh/triangulation
    mesher.refine_mesh();
  
    // Make sure triangulation is valid
    dolfin_assert(cdt.is_valid());

    // Mark the subdomains
    log(dolfin::TRACE, "Exploring subdomains in mesh");
    explore_subdomains(cdt, total_domain, subdomain_geometries);

    // Clear mesh
    mesh.clear();

    const std::size_t gdim = cdt.finite_vertices_begin()->point().dimension();
    const std::size_t tdim = cdt.dimension();
    const std::size_t num_vertices = cdt.number_of_vertices();

    // Count valid cells
    std::size_t num_cells = 0;
    for (CDT::Finite_faces_iterator cgal_cell = cdt.finite_faces_begin();
         cgal_cell != cdt.finite_faces_end(); ++cgal_cell)
    {
      // Add cell if it is in the domain
      if (cgal_cell->is_in_domain())
      {
        num_cells++;
      }
    }

    log(dolfin::DBG, "Mesh with %d vertices and %d cells created", num_vertices, num_cells);

    // Create a MeshEditor and open
    dolfin::MeshEditor mesh_editor;
    mesh_editor.open(mesh, tdim, gdim);
    mesh_editor.init_vertices(num_vertices);
    mesh_editor.init_cells(num_cells);

    // Add vertices to mesh
    std::size_t vertex_index = 0;
    for (CDT::Finite_vertices_iterator cgal_vertex = cdt.finite_vertices_begin();
         cgal_vertex != cdt.finite_vertices_end(); ++cgal_vertex)
    {
      // Get vertex coordinates and add vertex to the mesh
      dolfin::Point p;
      p[0] = cgal_vertex->point()[0];
      p[1] = cgal_vertex->point()[1];

      // Add mesh vertex
      mesh_editor.add_vertex(vertex_index, p);

      // Attach index to vertex and increment
      cgal_vertex->info() = vertex_index++;
    }

    dolfin_assert(vertex_index == num_vertices);

    // Add cells to mesh and build domain marker mesh function
    dolfin::MeshDomains &domain_markers = mesh.domains();
    std::size_t cell_index = 0;
    const bool mark_cells = geometry.has_subdomains();
    for (CDT::Finite_faces_iterator cgal_cell = cdt.finite_faces_begin();
         cgal_cell != cdt.finite_faces_end(); ++cgal_cell)
    {
      // Add cell if it is in the domain
      if (cgal_cell->is_in_domain())
      {
        mesh_editor.add_cell(cell_index,
                             cgal_cell->vertex(0)->info(),
                             cgal_cell->vertex(1)->info(),
                             cgal_cell->vertex(2)->info());

        if (mark_cells)
        {
          domain_markers.set_marker(std::make_pair(cell_index,
                                                   cgal_cell->counter()), 2);
        }
        ++cell_index;
      }
    }
    dolfin_assert(cell_index == num_cells);

    // Close mesh editor
    mesh_editor.close();
  }

  // Distribute the mesh (if in parallel)
  dolfin::MeshPartitioning::build_distributed_mesh(mesh);
}
}
//-----------------------------------------------------------------------------
