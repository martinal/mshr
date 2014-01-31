// Copyright (C) 2012 Anders Logg
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
// Modified by Benjamin Kehlet, 2012-2013

#ifndef __MSHR_PRIMITIVE_H
#define __MSHR_PRIMITIVE_H

#include "CSGGeometry.h"

namespace mshr
{

  /// Base class for Constructive Solid Geometry (CSG) primitives.
  class CSGPrimitive : public CSGGeometry
  {
  public:
    virtual bool is_operator() const { return false; }
    virtual std::size_t dim() const = 0;
  };

}

#endif
