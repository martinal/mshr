// Copyright (C) 2012 Anders Logg
//
// This file is part of DolfinCSG.
//
// DolfinCSG is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// DolfinCSG is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with DolfinCSG.  If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Benjamin Kehlet, 2012

#include <dolfincsg/CSGPrimitives3D.h>

#include <dolfin/math/basic.h>
#include <dolfin/log/LogStream.h>

#include <sstream>

namespace dolfincsg
{

//-----------------------------------------------------------------------------
// Sphere
//-----------------------------------------------------------------------------
Sphere::Sphere(dolfin::Point center, double radius, std::size_t slices)
  : c(center), r(radius), _slices(slices)
{
  if (r < DOLFIN_EPS)
  {
    dolfin::dolfin_error("CSGPrimitives3D.cpp",
                         "Create sphere",
                         "Sphere with center (%f, %f, %f) has zero or negative radius", c.x(), c.y(), c.z());
  }

  if (slices < 1)
  {
    dolfin::dolfin_error("CSGPrimitives3D.cpp",
		 "Create sphere",
		 "Can't create sphere with zero slices");
  }
}
//-----------------------------------------------------------------------------
std::string Sphere::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Sphere at " << c << " "
      << "with radius " << r << ">";
  }
  else
    s << "Sphere(" << c << ", " << r << ")";

  return s.str();
}
//-----------------------------------------------------------------------------
// Box
//-----------------------------------------------------------------------------
Box::Box(double x0, double x1, double x2,
         double y0, double y1, double y2)
  : _x0(x0), _x1(x1), _x2(x2), _y0(y0), _y1(y1), _y2(y2)
{
  // FIXME: Check validity of coordinates here
  if (dolfin::near(x0, y0) || dolfin::near(x1, y2) || dolfin::near(x2, y2))
    dolfin::dolfin_error("CSGPrimitives3D.cpp",
                         "Create axis aligned box",
                         "Box with corner (%f, %f, %f) and (%f, %f, %f) degenerated", 
                         x0, x1, x2, y0, y1, y2);
}
//-----------------------------------------------------------------------------
std::string Box::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Box with first corner at (" << _x0 << ", " << _x1 << ", " << _x2 << ") "
      << "and second corner at (" << _y0 << ", " << _y1 << ", " << _y2 << ")>";
  }
  else
  {
    s << "Box("
      << _x0 << ", " << _x1 << ", " << _x2 << ", "
      << _y0 << ", " << _y1 << ", " << _y2 << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// Cone
//-----------------------------------------------------------------------------
Cone::Cone(dolfin::Point top, 
           dolfin::Point bottom, 
           double top_radius, 
           double bottom_radius,
           std::size_t slices)
  : _top(top), _bottom(bottom), _top_radius(top_radius),
    _bottom_radius(bottom_radius), _slices(slices)
{
  if (dolfin::near(top_radius, 0.0) && dolfin::near(bottom_radius, 0.0))
  {
    dolfin::dolfin_error("CSGPrimitives3D.cpp",
		   "Create cone",
		   "Cone with zero thickness");
  }

  if (top.distance(bottom) < DOLFIN_EPS)
  {
    dolfin::dolfin_error("CSGPrimitives3D.cpp",
		 "Create cone",
		 "Cone with zero length");
  }
}
//-----------------------------------------------------------------------------
std::string Cone::str(bool verbose) const
{
  std::stringstream s;
  if (verbose)
  {
    s << "<Cone with top at " << _top << ", top radius " << _top_radius
      << " and bottom at " << _bottom << ", bottom radius "
      << _bottom_radius << ", with " << _slices << " slices>";
  }
  else
  {
    s << "Cone( " << _top << ", " << _bottom << ", " << _top_radius
      << ", " << _bottom_radius << " )";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
Tetrahedron::Tetrahedron(dolfin::Point x0, 
                         dolfin::Point x1, 
                         dolfin::Point x2, 
                         dolfin::Point x3)
  : _x0(x0), _x1(x1), _x2(x2), _x3(x3)
{}
//-----------------------------------------------------------------------------
/// Informal string representation
std::string Tetrahedron::str(bool verbose) const
{
  std::stringstream s;
  if (verbose)
  {
    s << "<Tetrahedron with point at " << _x0 << ", " << _x1 << ", "
      << _x2 << ", " << _x3 << ">";
  }
  else
  {
    s << "Tetrahedron( " << _x0 << ", " << _x1 << ", " << _x2 << ", "
      << _x3 << ")";
  }
  return s.str();
}
//-----------------------------------------------------------------------------
Surface3D::Surface3D(std::string filename) : _filename(filename)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::string Surface3D::str(bool verbose) const
{
  return std::string("Surface3D from file ") + _filename;
}
//-----------------------------------------------------------------------------
}
