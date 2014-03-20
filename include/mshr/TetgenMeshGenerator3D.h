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
//
// Modified by Joachim B Haga 2012

#ifndef __MSHR_TETGEN_MESH_GENERATOR3D_H
#define __MSHR_TETGEN_MESH_GENERATOR3D_H

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
    TetgenMeshGenerator3D(const CSGGeometry& geometry);
    TetgenMeshGenerator3D(std::shared_ptr<const CSGGeometry> geometry);
    ~TetgenMeshGenerator3D();

    void generate(dolfin::Mesh& mesh) const;

    /// Default parameter values
    static dolfin::Parameters default_parameters()
    {
      dolfin::Parameters p("tetgen_meshgenerator");
      p.add("mesh_resolution", 64);

      return p;
    }

  private:
    std::shared_ptr<const CSGGeometry> _geometry;
  };

}

#endif
