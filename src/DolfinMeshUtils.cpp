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

#include <mshr/DolfinMeshUtils.h>

#include <dolfin/mesh/Cell.h>

#include <limits>

namespace mshr
{

std::pair<double, double> DolfinMeshUtils::cell_volume_min_max(const dolfin::Mesh& m)
{
  std::pair<double, double> res(std::numeric_limits<double>::max(), 0.0);

  for (dolfin::CellIterator cell(m); !cell.end(); ++cell)
  {
    const double v = cell->volume();
    res.first = std::min(res.first, v);
    res.second = std::max(res.second, v);
  }

  return res;
}
//-----------------------------------------------------------------------------
bool DolfinMeshUtils::has_isolated_vertices(const dolfin::Mesh& m)
{
  std::set<std::size_t> vertices;
  for (dolfin::CellIterator cit(m); !cit.end(); ++cit)
  {
    const unsigned int* v = cit->entities(0);
    for (std::size_t i = 0; i < cit->num_global_entities(0); i++)
    {
      vertices.insert(v[i]);
    }
  }

  bool isolated_vertices = false;
  for (std::size_t i = 0; i < m.num_vertices(); i++)
  {
    if (vertices.count(i) < 1)
    {
      log(dolfin::DBG, "Vertex %u has no incident cells", i);
      isolated_vertices = true;
    }
  }

  return isolated_vertices;
}
}
