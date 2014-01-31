// Copyright (C) 2012 Benjamin Kehlet
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

#ifndef __MSHR_GEOMETRY_TO_CGAL_CONVERTER_H
#define __MSHR_GEOMETRY_TO_CGAL_CONVERTER_H

#include "cgal_csg3d.h"

namespace mshr
{
  class CSGGeometry;

  // This class is capable of converting a 3D dolfin::CSGGeometry to a
  // CGAL::Polyhedron_3
  class GeometryToCGALConverter
  {
  public:
    static void convert(const CSGGeometry& geometry, csg::Polyhedron_3& p,
                        bool remove_degenerated=true);
  };
}

#endif
