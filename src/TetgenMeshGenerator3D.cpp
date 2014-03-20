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

#include <mshr/TetgenMeshGenerator3D.h>
#include <mshr/CSGCGALDomain3D.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/common/NoDeleter.h>
#include <tetgen.h>

namespace
{
//-----------------------------------------------------------------------------
static void build_dolfin_mesh(const tetgenio& tetgenmesh, dolfin::Mesh& dolfinmesh)
{
  // Clear mesh
  dolfinmesh.clear();

  // Create and initialize mesh editor
  dolfin::MeshEditor mesh_editor;
  mesh_editor.open(dolfinmesh, 3, 3);
  mesh_editor.init_vertices(tetgenmesh.numberofpoints);

  const int offset = tetgenmesh.firstnumber;

  for (int i = 0; i < tetgenmesh.numberofpoints; i++) 
  {

    dolfin::Point p(tetgenmesh.pointlist[i * 3], 
                    tetgenmesh.pointlist[i * 3 + 1], 
                    tetgenmesh.pointlist[i * 3 + 2]);
    mesh_editor.add_vertex(i, p);
  }


  mesh_editor.init_cells(tetgenmesh.numberoftetrahedra);
  dolfin_assert(tetgenmesh.numberofcorners == 4);

  for (int i = 0; i < tetgenmesh.numberoftetrahedra; i++) 
  {
    mesh_editor.add_cell(i,
                         tetgenmesh.tetrahedronlist[i*4 + 0]-offset,
                         tetgenmesh.tetrahedronlist[i*4 + 1]-offset,
                         tetgenmesh.tetrahedronlist[i*4 + 2]-offset,
                         tetgenmesh.tetrahedronlist[i*4 + 3]-offset);
  }

  // Close mesh editor
  mesh_editor.close();
}
} // end anonymous namespace
//-----------------------------------------------------------------------------
namespace mshr
{

TetgenMeshGenerator3D::TetgenMeshGenerator3D(const CSGGeometry& geometry)
{
  std::shared_ptr<const CSGGeometry> tmp = dolfin::reference_to_no_delete_pointer<const CSGGeometry>(geometry);
  _geometry = tmp;
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
TetgenMeshGenerator3D::TetgenMeshGenerator3D(std::shared_ptr<const CSGGeometry> geometry)
: _geometry(geometry)
{
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
TetgenMeshGenerator3D::~TetgenMeshGenerator3D()
{

}
//-----------------------------------------------------------------------------
void TetgenMeshGenerator3D::generate(dolfin::Mesh& mesh) const
{
  tetgenio in;

  std::vector<dolfin::Point> vertices;
  std::vector<std::array<std::size_t, 3> > facets;
  
  {
    CSGCGALDomain3D exact_domain(*_geometry);
    exact_domain.ensure_meshing_preconditions();
    exact_domain.get_vertices(vertices);
    exact_domain.get_facets(facets);
  }


  // Copy the vertices to the tetgen structure
  in.numberofpoints = vertices.size();
  in.pointlist = new REAL[in.numberofpoints * 3];

  int i = 0;
  for (std::vector<dolfin::Point>::const_iterator it = vertices.begin();
       it != vertices.end(); it++)
  {
    in.pointlist[i*3 + 0] = it->x();
    in.pointlist[i*3 + 1] = it->y();
    in.pointlist[i*3 + 2] = it->z();

    i++;
  }

  
  // Copy the facets
  in.numberoffacets = facets.size();
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  //in.facetmarkerlist = new int[in.numberoffacets];

  i = 0;
  for (std::vector<std::array<std::size_t, 3> >::const_iterator it = facets.begin();
       it != facets.end(); it++)
  {
    // Facet 1. The leftmost facet.
    tetgenio::facet *f = &in.facetlist[i];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    tetgenio::polygon *p = &f->polygonlist[0];
    p->numberofvertices = 3;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = (*it)[0];
    p->vertexlist[1] = (*it)[1];
    p->vertexlist[2] = (*it)[2];

    i++;
  }

  // Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
  //   do quality mesh generation (q) with a specified quality bound
  //   (1.414), and apply a maximum volume constraint (a0.1).

  std::string str("pq1.414a0.1");
  std::unique_ptr<char> writable(new char[str.size() + 1]);
  std::copy(str.begin(), str.end(), writable.get());
  writable.get()[str.size()] = '\0'; // terminating 0

  tetgenio out;
  tetrahedralize(writable.get(), &in, &out);

  build_dolfin_mesh(out, mesh);
}
//-----------------------------------------------------------------------------
} // end namespace mshr
