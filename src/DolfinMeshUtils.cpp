// Copyright (C) 2014-2015 Benjamin Kehlet
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
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/MeshEditor.h>

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
//-----------------------------------------------------------------------------
bool DolfinMeshUtils::check_mesh(const dolfin::Mesh& m)
{
  return !has_isolated_vertices(m);
}
//-----------------------------------------------------------------------------
std::shared_ptr<dolfin::Mesh>
   DolfinMeshUtils::extract_subdomain(std::shared_ptr<const dolfin::Mesh> mesh,
                                      std::size_t cell_domain)
{
  dolfin_assert(mesh->geometry().dim() == 3);
  dolfin_assert(mesh->topology().dim() == 3);

  // Collect all vertices incident to all marked cells
  std::map<std::size_t, std::size_t> collected_vertices;
  std::size_t num_cells = 0;
  for (const std::pair<std::size_t, std::size_t>& marker : mesh->domains().markers(3))
  {
    if (marker.second == cell_domain)
    {
      num_cells++;
      dolfin::Cell c(*mesh, marker.second);
      for (std::size_t i = 0; i < 4; i++)
      {
        const std::size_t s = collected_vertices.size();
        collected_vertices.insert(std::make_pair(c.entities(0)[i], s));
      }
    }
  }

  std::shared_ptr<dolfin::Mesh> outmesh(new dolfin::Mesh);
  dolfin::MeshEditor editor;
  editor.open(*outmesh, 3,3);

  editor.init_vertices(collected_vertices.size());
  for (std::pair<std::size_t, std::size_t> v : collected_vertices)
  {
    dolfin::Vertex existing_vertex(*mesh, v.first);
    editor.add_vertex(v.second, existing_vertex.point());
  }

  editor.init_cells(num_cells);
  std::size_t cell_counter = 0;
  for (const std::pair<std::size_t, std::size_t>& marker : mesh->domains().markers(3))
  {
    if (marker.second == cell_domain)
    {
      const dolfin::Cell c(*mesh, marker.second);
      const unsigned int* vertices = c.entities(0);
      editor.add_cell(cell_counter,
                      collected_vertices[vertices[0]],
                      collected_vertices[vertices[1]],
                      collected_vertices[vertices[2]],
                      collected_vertices[vertices[3]]);

    }
  }

  editor.close();
  return outmesh;
}
}
