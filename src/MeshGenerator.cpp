// Copyright (C) 2012 Anders Logg, Benjamin Kehlet, Johannes Ring
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

#include <mshr/MeshGenerator.h>
#include <mshr/CSGGeometry.h>
#include <mshr/CSGCGALMeshGenerator2D.h>
#include <mshr/CSGCGALMeshGenerator3D.h>
#include <mshr/TetgenMeshGenerator3D.h>

#include <dolfin/log/log.h>
#include <dolfin/mesh/BoundaryMesh.h>


namespace mshr
{

//-----------------------------------------------------------------------------
void generate(dolfin::Mesh& mesh,
              const CSGGeometry& geometry,
              double resolution,
              std::string backend)
{
  if (geometry.dim() == 2)
  {
    CSGCGALMeshGenerator2D generator;
    generator.parameters["mesh_resolution"] = resolution;
    generator.generate(geometry, mesh);
  }
  else if (geometry.dim() == 3)
  {
    if (backend == "cgal")
    {
      CSGCGALMeshGenerator3D generator;
      generator.parameters["mesh_resolution"] = resolution;
      generator.generate(geometry, mesh);
    }
    else if (backend == "tetgen")
    {
      TetgenMeshGenerator3D generator(geometry);
      generator.parameters["mesh_resolution"] = resolution;
      generator.generate(mesh);

    }
    else
    {
      std::string e = "Unknown mesh generator backend: " + backend;
      dolfin::dolfin_error("MeshGenerator.cpp",
                           "Generator mesh of 3D geometry",
                           e);
    }
  }
  else
  {
    dolfin::dolfin_error("MeshGenerator.cpp",
                         "create mesh from CSG geometry",
                         "Unhandled geometry dimension %d", geometry.dim());
  }
}

}
