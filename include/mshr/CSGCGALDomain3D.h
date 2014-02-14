// Copyright (C) 2013-2014 Benjamin Kehlet
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

#ifndef __MSHR_CSGCGAL_DOMAIN3D_H
#define __MSHR_CSGCGAL_DOMAIN3D_H

#include <mshr/CSGGeometry.h>

#include <dolfin/geometry/Point.h>

#include <boost/scoped_ptr.hpp>

namespace mshr
{

  // Forward declaration
  struct CSGCGALDomain3DImpl;


// This class represents the polyhedral domain which implements boolean
// operations. It uses CGAL as backend
class CSGCGALDomain3D
{
 public:
  // Create empty polygon
  CSGCGALDomain3D();

  // Construct polyhedron from CSG geometry
  CSGCGALDomain3D(const mshr::CSGGeometry &csg);

  // Destructor
  ~CSGCGALDomain3D();

  std::size_t num_vertices() const;
  std::size_t num_facets() const;
  std::size_t num_halfedges() const;

  // Output in double precision
  // TODO: Define iterators to be more memory friendly
  void get_vertices(std::vector<dolfin::Point> &v) const;
  void get_facets(std::vector<std::array<std::size_t, 3> > &f) const;

  void remove_degenerated_facets(double threshold);

 private :
  boost::scoped_ptr<CSGCGALDomain3DImpl> impl;
};

}
#endif
