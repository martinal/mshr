// Copyright (C) 2012 Anders Logg and 2012-2014 Benjamin Kehlet
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
// along with mshr. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Joachim B Haga, 2012


#ifndef __MSHR_MESH_GENERATOR_H
#define __MSHR_MESH_GENERATOR_H

#include "CSGGeometry.h"

namespace dolfin
{
  // Forward declarations
  class Mesh;
  class BoundaryMesh;
}

namespace mshr
{

  /// Generate mesh from CSG geometry
  void generate(dolfin::Mesh& mesh,
                const CSGGeometry& geometry,
                double resolution,
                std::string backend="cgal");

  /// Extract boundary of CSG geometry as dolfin::BoundaryMesh
  void get_boundary_mesh(dolfin::BoundaryMesh& mesh,
                         const CSGGeometry& geometry);
}

#endif
