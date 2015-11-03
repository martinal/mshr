// Copyright (C) 2015 Benjamin Kehlet
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

#include <mshr/SurfaceReconstruction.h>
#include <dolfin/log/log.h>

#define CGAL_EIGEN3_ENABLED true
#include <CGAL/trace.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>

// Based on CGAL/Surface_reconstruction_points_3/poisson_reconstruction_example.cpp
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::vector<Point_with_normal> PointList;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;
typedef typename C2t3::Triangulation Tr;
typedef typename Tr::Vertex_handle Vertex_handle;
typedef typename Tr::Edge Edge;
typedef typename Tr::Facet Facet;
typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;


namespace
{

template <class Vertex_handle>
int get_vertex_index(std::vector<std::array<double, 3>>& vertices,
                     Vertex_handle vh,
                     std::map<Vertex_handle,
                     int>& V,
                     int& inum)
{
  typedef typename std::map<Vertex_handle, int>::iterator map_iterator;
  std::pair<map_iterator,bool> insert_res = V.insert( std::make_pair(vh,inum) );
  if ( insert_res.second )
  {
    typename Tr::Point p = static_cast<typename Tr::Point>(vh->point());
    vertices.push_back({CGAL::to_double(p[0]),
          CGAL::to_double(p[1]),
          CGAL::to_double(p[2])});
    ++inum;
  }
  return insert_res.first->second;
}
//-----------------------------------------------------------------------------
void export_triangulation(const C2t3& c2t3,
                          std::vector<std::array<double, 3>>& vertices,
                          std::vector<std::array<std::size_t, 3>>& facets)
{
  const Tr& tr = c2t3.triangulation();
  const typename Tr::size_type number_of_facets = c2t3.number_of_facets();

  vertices.reserve(tr.number_of_vertices());
  facets.reserve(number_of_facets);

  // Finite vertices coordinates.
  Finite_facets_iterator fit = tr.finite_facets_begin();
  std::set<Facet> oriented_set;
  std::stack<Facet> stack;

  //CGAL_assertion_code(typename Tr::size_type nb_facets = 0; )

  while (oriented_set.size() != number_of_facets)
  {
    while ( fit->first->is_facet_on_surface(fit->second) == false ||
            oriented_set.find(*fit) != oriented_set.end() ||
            oriented_set.find(c2t3.opposite_facet(*fit)) != oriented_set.end() )
    {
      ++fit;
    }
    oriented_set.insert(*fit);
    stack.push(*fit);
    while(! stack.empty() )
    {
      Facet f = stack.top();
      stack.pop();
      for(int ih = 0 ; ih < 3 ; ++ih)
      {
        const int i1  = tr.vertex_triple_index(f.second, tr. cw(ih));
        const int i2  = tr.vertex_triple_index(f.second, tr.ccw(ih));
        if( c2t3.face_status(Edge(f.first, i1, i2)) == C2t3::REGULAR )
        {
          Facet fn = c2t3.neighbor(f, ih);
          if (oriented_set.find(fn) == oriented_set.end() &&
              oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
          {
            oriented_set.insert(fn);
            stack.push(fn);
          }
        } // end "if the edge is regular"
      } // end "for each neighbor of f"
    } // end "stack non empty"
  } // end "oriented_set not full"

  // Orients the whole mesh towards outside:
  // - find the facet with max z
  typename std::set<Facet>::const_iterator top_facet = oriented_set.begin();
  for(typename std::set<Facet>::const_iterator fit = oriented_set.begin();
      fit != oriented_set.end();
      ++fit)
  {
    double top_z =
      (top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 0))->point().z()
       + top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 1))->point().z()
       + top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 2))->point().z())/3.;
    double z =
      (fit->first->vertex(tr.vertex_triple_index(fit->second, 0))->point().z()
       + fit->first->vertex(tr.vertex_triple_index(fit->second, 1))->point().z()
       + fit->first->vertex(tr.vertex_triple_index(fit->second, 2))->point().z())/3.;
    if (top_z < z)
      top_facet = fit;
  }

  // - orient the facet with max z towards +Z axis
  Vertex_handle v0 = top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 0));
  Vertex_handle v1 = top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 1));
  Vertex_handle v2 = top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 2));
  Vector normal = cross_product(v1->point()-v0->point(), v2->point()-v1->point());
  const Vector Z(0, 0, 1);
  bool regular_orientation = (Z * normal >= 0);

  // used to set indices of vertices
  std::map<Vertex_handle, int> V;
  int inum = 0;

  for(typename std::set<Facet>::const_iterator fit =
        oriented_set.begin();
      fit != oriented_set.end();
      ++fit)
  {
    int indices[3];
    int index = 0;
    for (int i=0; i<3; i++)
    {
      indices[index++] = get_vertex_index(vertices,
                                          fit->first->vertex(tr.vertex_triple_index(fit->second, i)),
                                          V,
                                          inum);
    }

    facets.push_back({ indices[0],
          regular_orientation ? indices[1] : indices[2],
          regular_orientation ? indices[2] : indices[1] });
  }
  CGAL_assertion(nb_facets == number_of_facets);
}
}
//-----------------------------------------------------------------------------
void mshr::SurfaceReconstruction::reconstruct(const std::vector<double>& vertices,
                                              const std::vector<std::size_t>& facets,
                                              std::vector<std::array<double, 3>>& reconstructed_vertices,
                                              std::vector<std::array<std::size_t, 3>>& reconstructed_facets)
{
  // Poisson options
  FT sm_angle = 20.0; // Min triangle angle in degrees.
  FT sm_radius = 10; // Max triangle size w.r.t. point set average spacing.
  FT sm_distance = 0.10; // Surface Approximation error w.r.t. point set average spacing.

  // Reads the point set file in points[].
  // Note: read_xyz_points_and_normals() requires an iterator over points
  // + property maps to access each point's position and normal.
  // The position property map can be omitted here as we use iterators over Point_3 elements.
  PointList points;
  for (std::size_t i = 0; i < facets.size(); i += 3)
  {
    const Point a(vertices[facets[i]*3],   vertices[facets[i]*3+1],   vertices[facets[i]*3+2]);
    const Point b(vertices[facets[i+1]*3], vertices[facets[i+1]*3+1], vertices[facets[i+1]*3+2]);
    const Point c(vertices[facets[i+2]*3], vertices[facets[i+2]*3+1], vertices[facets[i+2]*3+2]);

    // compute normal
    const Vector normal = CGAL::cross_product(b-a, c-a);
    const Point centroid = CGAL::ORIGIN + ((a-CGAL::ORIGIN)+(b-CGAL::ORIGIN)+(c-CGAL::ORIGIN))/3.0;
    Point_with_normal pm( centroid,
                          normal/std::sqrt(normal.squared_length()));
    points.push_back(pm);
  }

  // Creates implicit function from the read points using the default solver.
  // Note: this method requires an iterator over points
  // + property maps to access each point's position and normal.
  // The position property map can be omitted here as we use iterators over Point_3 elements.
  log(dolfin::TRACE, "Construct implicit function");
  Poisson_reconstruction_function function(points.begin(), points.end(),
                                           CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()) );
  log(dolfin::TRACE, "Compute Poisson indicator function");
  // Computes the Poisson indicator function f()
  // at each vertex of the triangulation.
  if ( ! function.compute_implicit_function() )
  {
    dolfin::dolfin_error("SurfaceReconstruction.cpp",
                         "reconstruct surface",
                         "couldn't compute implicit function");
  }
  // Computes average spacing
  log(dolfin::TRACE, "Compute average spacing");
  FT average_spacing = CGAL::compute_average_spacing(points.begin(), points.end(),
                                                     6 /* knn = 1 ring */);
  // Gets one point inside the implicit surface
  // and computes implicit function bounding sphere radius.
  Point inner_point = function.get_inner_point();
  Sphere bsphere = function.bounding_sphere();
  FT radius = std::sqrt(bsphere.squared_radius());
  // Defines the implicit surface: requires defining a
  // conservative bounding sphere centered at inner point.
  FT sm_sphere_radius = 5.0 * radius;
  FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
  Surface_3 surface(function,
                    Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
                    sm_dichotomy_error/sm_sphere_radius);

  // Defines surface mesh generation criteria
  CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                      sm_radius*average_spacing,  // Max triangle size
                                                      sm_distance*average_spacing); // Approximation error
  // Generates surface mesh with manifold option
  log(dolfin::TRACE, "Invoke surface mesher");
  STr tr; // 3D Delaunay triangulation for surface mesh generation
  C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
  CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
                          surface,                              // implicit surface
                          criteria,                             // meshing criteria
                          CGAL::Manifold_with_boundary_tag());  // require manifold mesh
  if(tr.number_of_vertices() == 0)
    dolfin::dolfin_error("SurfaceReconstruction.cpp",
                         "reconstruct surface from point set",
                         "Couldn't reconstruct surface");

  export_triangulation(c2t3,
                       reconstructed_vertices,
                       reconstructed_facets);
}

