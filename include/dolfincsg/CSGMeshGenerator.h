// Copyright (C) 2012 Anders Logg and 2012-2014 Benjamin Kehlet
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Joachim B Haga, 2012


#ifndef __DOLFINCSG_MESH_GENERATOR_H
#define __DOLFINCSG_MESH_GENERATOR_H

#include "CSGGeometry.h"
#include <boost/shared_ptr.hpp>

namespace dolfin
{
  // Forward declarations
  class Mesh;
  class BoundaryMesh;
}

namespace dolfincsg
{

  /// Mesh generator for Constructive Solid Geometry (CSG)

  class CSGMeshGenerator
  {
  public :

    /// Generate mesh from CSG geometry
    static void generate(dolfin::Mesh& mesh, 
                         const CSGGeometry& geometry,
                         std::size_t resolution);

    /// Generate boundary mesh from the surface of a CSG geometry
    static void generate(dolfin::BoundaryMesh& mesh, 
                         const CSGGeometry& geometry);
  };
}

#endif
