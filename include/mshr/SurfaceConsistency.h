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

// OBS! Experimental

#ifndef _SURFACE_CONSISTENCY_H
#define _SURFACE_CONSISTENCY_H

#include <vector>
#include <set>
#include <array>

namespace mshr
{

class SurfaceConsistency
{
 public:

  /// Check that the connectivity of the facet is consistent, ie. all edges are
  /// shared by exactly two facets.
  /// If error is set to True, then an error will be thrown when a duplicated
  /// halfedges is encountered. Otherwise, the index of one of the facets will
  /// be stored in the set.
  static void checkConnectivity(const std::vector<std::array<std::size_t, 3> >& facets,
                                std::set<std::size_t>& duplicating, bool error);

  static void filterFacets(const std::vector<std::array<std::size_t, 3> >& facets,
                           const std::vector<std::array<double, 3> >& vertices,
                           std::size_t start, std::set<std::size_t>& skip);
};

}
#endif
