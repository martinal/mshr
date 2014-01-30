// Copyright (C) 2012 Benjamin Kehlet
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __DOLFINCSG_GEOMETRY_TO_CGAL_CONVERTER_H
#define __DOLFINCSG_GEOMETRY_TO_CGAL_CONVERTER_H

#include "cgal_csg3d.h"

namespace dolfincsg
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
