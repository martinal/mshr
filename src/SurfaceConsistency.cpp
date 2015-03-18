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

// OBS! Experimental code

#include <mshr/SurfaceConsistency.h>

#include <dolfin/log/log.h>
#include <vector>
#include <map>
#include <iostream>

namespace mshr
{

void SurfaceConsistency::checkConnectivity(const std::vector<std::vector<std::size_t> >& facets,
                                           std::set<std::size_t>& duplicating, bool error)
{
  // Store all halfedges
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> halfedges;

  std::size_t facet_no = 0;
  for (auto it = facets.begin(); it != facets.end(); ++it)
  {
    // Check for (topologically) degenerate facets
    if ( (*it)[0] == (*it)[1] || (*it)[0] == (*it)[2] || (*it)[1] == (*it)[2] )
      dolfin::dolfin_error("SurfaceConsistency.cpp",
                           "confirm surface connectivity",
                           "Facet %d is degenerate", facet_no);


    for (int i = 0; i < 3; i++)
    {
      std::pair<std::size_t, std::size_t> e( (*it)[i], (*it)[(i+1)%3] );
      if (halfedges.count( e ) > 0 )
      {
        if (error)
        {
          dolfin::dolfin_error("SurfaceConsistency.cpp",
                               "confirm halfedge connectivity",
                               "Facet %d and %d share halfedge", halfedges[e], facet_no);
        }
        else
        {
          duplicating.insert(facet_no);
        }
      }
      else
      {
        halfedges[e] = facet_no;
      }
    }
    facet_no++;
  }

  // TODO: Check for border edges
}
}

