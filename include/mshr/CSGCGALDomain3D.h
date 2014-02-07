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

  void remove_degenerated();

 private :
  boost::scoped_ptr<CSGCGALDomain3DImpl> impl;
};

}
#endif
