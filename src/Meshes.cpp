// Copyright (C) -2014 Benjamin Kehlet
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

#include <mshr/Meshes.h>
#include <mshr/TetgenMeshGenerator3D.h>
#include <mshr/CSGPrimitives3D.h>

#include <dolfin/mesh/MeshPartitioning.h>

namespace mshr
{
  UnitSphereMesh::UnitSphereMesh(std::size_t resolution)
    : dolfin::Mesh()
  {
    // Receive mesh according to parallel policy
    if (dolfin::MPI::is_receiver(this->mpi_comm()))
    {
      dolfin::MeshPartitioning::build_distributed_mesh(*this);
      return;
    }

    Sphere s(dolfin::Point(0,0,0), 1.0, resolution);

    TetgenMeshGenerator3D generator;
    const double max_volume = 1.0/std::pow(3.0, resolution);
    generator.parameters["max_tet_volume"] = max_volume;
    generator.parameters["max_radius_edge_ratio"] = 2.0;
    generator.parameters["min_dihedral_angle"] = 0.0;

    generator.generate(s, *this);

    // Broadcast mesh according to parallel policy
    if (dolfin::MPI::is_broadcaster(this->mpi_comm()))
    {
      dolfin::MeshPartitioning::build_distributed_mesh(*this);
      return;
    }
  }
}
