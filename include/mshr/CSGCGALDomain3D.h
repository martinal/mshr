// Copyright (C) 2013-2014 Benjamin Kehlet
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

#ifndef __MSHR_CSGCGAL_DOMAIN3D_H
#define __MSHR_CSGCGAL_DOMAIN3D_H

#include <mshr/CSGGeometry.h>
#include <dolfin/geometry/Point.h>

#include <memory>
#include <array>

namespace mshr
{

  // Forward declaration
  struct CSGCGALDomain3DImpl;
  struct CSGCGALDomain3DQueryStructureImpl;

class CSGCGALDomain3DQueryStructure
{
 public:
  CSGCGALDomain3DQueryStructure(std::unique_ptr<CSGCGALDomain3DQueryStructureImpl> impl);
  ~CSGCGALDomain3DQueryStructure();

  std::unique_ptr<CSGCGALDomain3DQueryStructureImpl> impl;
};

/// @brief A polyhedron meshing domain.
///
/// This class represents a polyhedral meshing domain which implements boolean
/// operations. It uses CGAL Polyhedron_3 as backend and utilizes CGAL
//  Nef_polyhedron for boolean operations.
class CSGCGALDomain3D : public dolfin::Variable
{
 public:
  /// @brief Create empty polyhedron
  CSGCGALDomain3D();

  /// @brief Construct polyhedron from CSG geometry
  CSGCGALDomain3D(const mshr::CSGGeometry &csg);

  /// @brief Destructor
  ~CSGCGALDomain3D();

  /// @brief Number of vertices in polyhedron
  std::size_t num_vertices() const;

  /// @brief Number of facets in polyhedron
  std::size_t num_facets() const;

  /// @brief Number of halfedges in polyhedron
  std::size_t num_halfedges() const;

  /// @brief Volume of polyhedron (experimental, use with care)
  double volume() const;

  /// @brief Is polyhedron  (experimental, use with care)
  bool is_insideout() const;

  /// @brief Count the number of degenerate facets wrt. the given tolerance
  std::size_t num_degenerate_facets(double threshold) const;

  /// @brief get length of shortest edge
  double shortest_edge() const;

  /// @brief Test if any facets intersects
  bool is_selfintersecting() const;

  /// @brief Save polyhedron to off file
  /// @param filename Filename to write to
  void save_off(std::string filename) const;

  // TODO: Define iterators to be more memory friendly

  /// @brief Output vertices in double precision
  void get_vertices(std::vector<dolfin::Point> &v) const;

  /// @brief Output facets as indices to the vertices array
  void get_facets(std::vector<std::array<std::size_t, 3> > &f) const;

  /// @brief get one point per hole, strictly inside the hole.
  /// @param holes the returned points
  /// @param q a query structure returned from get_query_structure()
  void get_points_in_holes(std::vector<dolfin::Point>& holes,
                           std::shared_ptr<CSGCGALDomain3DQueryStructure> q) const;

  /// @brief Attempt to remove degenerate facets.
  void remove_degenerate_facets(double tolerance);

  /// @brief Attempts to ensure that the preconditions
  /// for successfull meshing generation are fullfilled.
  void ensure_meshing_preconditions();

  /// @brief Return data structure that allows fast queries,
  /// essentially a wrapper for a CGAL AABB tree of the polyhedron triangles.
  /// When performing queries, the user is responsible for the query structure
  /// being in sync with the CSGCGALDomain3D object.
  std::shared_ptr<CSGCGALDomain3DQueryStructure> get_query_structure() const;

  /// Default parameter values
  static dolfin::Parameters default_parameters()
  {
    dolfin::Parameters p("csg_cgal_domain_3d");
    p.add("remove_degenerate", true);
    p.add("degenerate_tolerance", 1e-12);

    return p;
  }

  /// @brief Informal string representation
  /// @param verbose The verbosity level
  std::string str(bool verbose) const;

 private :
  std::unique_ptr<CSGCGALDomain3DImpl> impl;
};

}
#endif
