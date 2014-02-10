// Copyright (C) 2012 Anders Logg
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

#include <dolfin.h>
#include <mshr.h>

int main(int argc, char** argv)
{
  // Define 3D geometry
  mshr::Box box(0, 0, 0, 1, 1, 1);
  mshr::Sphere sphere(dolfin::Point(0, 0, 0), 0.3);
  mshr::Cone cone(dolfin::Point(0, 0, -1), dolfin::Point(0, 0, 1), .5, .5);

  const boost::shared_ptr<mshr::CSGGeometry> g3d = box + cone - sphere;

  // Test printing
  dolfin::info("\nCompact output of 3D geometry:");
  dolfin::info(*g3d);
  dolfin::info("\nVerbose output of 3D geometry:");
  dolfin::info(*g3d, true);

  // Plot geometry
  //dolfin::plot(g3d, "3D geometry (surface)");

  // Generate and plot mesh
  dolfin::Mesh mesh3d;

  mshr::CSGMeshGenerator::generate(mesh3d, *g3d, 24);
  dolfin::cout << "Done generating mesh" << dolfin::endl;
  dolfin::info(mesh3d);
  dolfin::plot(mesh3d, "3D mesh");

  dolfin::interactive();

  return 0;
}
