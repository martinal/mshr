// Copyright (C) 2012 Benjamin Kehlet
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

#ifndef __MSHR_CGAL_MESH_GENERATOR3D_H
#define __MSHR_CGAL_MESH_GENERATOR3D_H

#include <mshr/CSGCGALDomain3D.h>

#include <dolfin/common/Variable.h>
#include <memory>

namespace dolfin{ class Mesh; }

namespace mshr
{

  // Forward declaration
  class CSGGeometry;

  /// @brief Mesh generator for Constructive Solid Geometry (CSG)
  /// utilizing CGALs 3D Mesh generation package.
  ///
  /// This class gives access to the backend specific meshing
  /// criterias throught the parameter system.
  class CSGCGALMeshGenerator3D : public dolfin::Variable
  {
  public :

    /// @brief Create mesh generator
    CSGCGALMeshGenerator3D();

    /// @brief Destructor
    ~CSGCGALMeshGenerator3D();

    /// @brief Generate Dolfin mesh
    /// @param domain The polyhedral domain to be meshed
    /// @param mesh The mesh object to be filled. Will be cleared.
    void generate(std::shared_ptr<const CSGCGALDomain3D> domain, dolfin::Mesh& mesh) const;

    /// @brief Generate Dolfin mesh
    /// @param geometry The csg geometry to be meshed
    /// @param mesh The mesh object to be filled. Will be cleared.
    void generate(std::shared_ptr<const CSGGeometry> geometry, dolfin::Mesh& mesh) const;

    /// @brief Generate Dolfin mesh
    /// @param geometry The csg geometry to be meshed
    /// @param mesh The mesh object to be filled. Will be cleared.
    void generate(const CSGGeometry& geometry, dolfin::Mesh& mesh) const;

    /// @brief Default parameter values
    static dolfin::Parameters default_parameters()
    {
      dolfin::Parameters p("csg_cgal_meshgenerator");
      p.add("mesh_resolution", 64.0);
      p.add("perturb_optimize", false);
      p.add("exude_optimize", false);
      p.add("lloyd_optimize", false);
      p.add("odt_optimize", false);
      p.add("edge_size", 0.025);
      p.add("facet_angle", 25.0);
      p.add("facet_size", 0.05);
      p.add("facet_distance", 0.005);
      p.add("cell_radius_edge_ratio", 3.0);
      p.add("cell_size", 0.05);
      p.add("detect_sharp_features", true);

      return p;
    }
  };

}

#endif
