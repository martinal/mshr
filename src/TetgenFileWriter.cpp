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

#include <mshr/TetgenFileWriter.h>
#include <mshr/CSGCGALDomain3D.h>

#include <fstream>

namespace mshr
{

void TetgenFileWriter::write(const CSGCGALDomain3D &domain, std::string filename)
{
  std::vector<dolfin::Point> vertices;
  domain.get_vertices(vertices);

  std::vector<std::array<std::size_t, 3> > facets;
  domain.get_facets(facets);

  write(vertices, facets, filename);

}
//-----------------------------------------------------------------------------
void TetgenFileWriter::write(const std::vector<dolfin::Point> &vertices,
                             const std::vector<std::array<std::size_t, 3> > &facets,
                             std::string filename)
{
  std::ofstream ofile(filename.c_str());

  // Write the node list
  ofile << "# node count, 3 dim, no attribute, no boundary marker" << std::endl;
  ofile << vertices.size() << " 3 0 0 " << std::endl;
  for (auto it = vertices.begin(); it != vertices.end(); it++)
  {
    ofile << std::distance(vertices.begin(), it) << " "
          << it->x() << " " 
          << it->y() << " "
          << it->z() << std::endl;
  }

  // Write the facet list
  ofile << "# facet count, no boundary markers" << std::endl;
  ofile << facets.size() << " 0" << std::endl;
  for (auto it = facets.begin(); it != facets.end(); it++)
  {
    ofile << "1 # 1 polygon, no hole, no boundary marker" << std::endl;
    ofile << "3 " 
          << (*it)[0] << " " 
          << (*it)[1] << " " 
          << (*it)[2] << std::endl;
    
  }
}

}
