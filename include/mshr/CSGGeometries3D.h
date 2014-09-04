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
//
// Modified by Johannes Ring, 2012

#ifndef __MSHR_GEOMETRIES_H
#define __MSHR_GEOMETRIES_H

#include "CSGGeometry.h"
#include <dolfin/geometry/Point.h>
#include <memory>

namespace mshr
{
  class CSGGeometries
  {
  public:

    // A standard LEGO brick starting at the point x with (n0, n1)
    // knobs and height n2. The height should be 1 for a thin brick or 3
    // for a regular brick.
    static std::shared_ptr<CSGGeometry> lego(std::size_t n0,
                                             std::size_t n1,
                                             std::size_t n2,
                                             dolfin::Point x = dolfin::Point(0,0,0));

    // A simple propeller with parameters r - radius of center body, R - length of blades,
    // w - width of blades and h - thicknes of blades
    static std::shared_ptr<CSGGeometry> propeller(double r=0.125,
                                                  double R=0.5,
                                                  double w=0.3,
                                                  double h=0.025,
                                                  bool include_tip=false);
  };
}

#endif
