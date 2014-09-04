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

#ifndef __MSHR_TETGEN_MESH_GENERATOR3D_H
#define __MSHR_TETGEN_MESH_GENERATOR3D_H

#include <mshr/CSGCGALDomain3D.h>

#include <dolfin/common/Variable.h>
#include <memory>

namespace dolfin{ class Mesh; }

namespace mshr
{

  // Forward declaration
  class CSGGeometry;

  /// Mesh generator for Constructive Solid Geometry (CSG)
  /// utilizing CGALs boolean operation on Nef_polyhedrons.

  class TetgenMeshGenerator3D : public dolfin::Variable
  {
  public :
    TetgenMeshGenerator3D();
    ~TetgenMeshGenerator3D();

    void generate(std::shared_ptr<const CSGCGALDomain3D> domain, dolfin::Mesh& mesh) const;
    void generate(const CSGGeometry& geometry, dolfin::Mesh& mesh) const;

    /// Default parameter values
    static dolfin::Parameters default_parameters()
    {
      dolfin::Parameters p("tetgen_meshgenerator");
      p.add("mesh_resolution", 64.0);

      p.add("disable_quality_improvement", false);
      p.add("max_radius_edge_ratio", 2.0);
      p.add("min_dihedral_angle", 12.);

      // If set to a positive value, this will override "mesh_resolution"
      p.add("max_tet_volume", -1.0);

      return p;
    }
  };
}

#endif
