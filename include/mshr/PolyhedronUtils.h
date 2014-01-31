// Copyright (C) 2012-2014 Benjamin Kehlet
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


#ifndef __MSHR_POLYHEDRON_UTILS_H
#define __MSHR_POLYHEDRON_UTILS_H

#include "cgal_csg3d.h"
#include "self_intersect.h"

// This file should be be included in any other header files to avoid
// CGAL leaking into the DOLFIN code.

namespace mshr
{

  // Some utilities for working with cgal polyhedrons
  class PolyhedronUtils
  {
  public:
    static void readSurfaceFile(std::string filename,
                                csg::Exact_Polyhedron_3& p);

    static void readSTLFile(std::string filename, csg::Exact_Polyhedron_3& p);
    static CGAL::Bbox_3 getBoundingBox(csg::Polyhedron_3& polyhedron);
    static double getBoundingSphereRadius(csg::Polyhedron_3& polyhedron);
    static bool has_degenerate_facets(csg::Exact_Polyhedron_3& p,
                                      double threshold);
    static void remove_degenerate_facets(csg::Exact_Polyhedron_3& p,
                                         const double threshold);

    template <typename Polyhedron>
    bool has_self_intersections(Polyhedron& p)
    {
      typedef typename Polyhedron::Triangle_3 Triangle;
      typedef typename std::back_insert_iterator<std::list<Triangle> > OutputIterator;

      std::list<Triangle> triangles; // intersecting triangles
      ::self_intersect<Polyhedron::Polyhedron_3, Polyhedron::Kernel,
          OutputIterator>(p, std::back_inserter(triangles));

      return triangles.size() > 0;
    }
  };
}

#endif
