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
#include <limits>
#include <dolfin/mesh/Cell.h>

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
