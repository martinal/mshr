// Copyright (C) 2014 Benjamin Kehlet
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

#ifndef __MSHR_DOLFIN_MESH_UTILS_H
#define __MSHR_DOLFIN_MESH_UTILS_H

#include <dolfin/mesh/Mesh.h>

class DolfinMeshUtils
{
 public:
  /// Compute the smallest and largest cell wrt. volume.
  /// @param m The mesh
  static std::pair<double, double> cell_volume_min_max(const dolfin::Mesh& m);
};

#endif
