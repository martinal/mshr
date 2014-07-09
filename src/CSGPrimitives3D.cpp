// Copyright (C) 2012 Anders Logg, 2012, 2014 Benjamin Kehlet
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

#include <mshr/CSGPrimitives3D.h>

#include <dolfin/math/basic.h>
#include <dolfin/log/LogStream.h>
#include <sstream>

namespace mshr
{

//-----------------------------------------------------------------------------
// Sphere
//-----------------------------------------------------------------------------
Sphere::Sphere(dolfin::Point center, double radius, std::size_t segments)
  : c(center), r(radius), _segments(segments)
{
  if (r < DOLFIN_EPS)
  {
    dolfin::dolfin_error("CSGPrimitives3D.cpp",
                         "Create sphere",
                         "Sphere with center (%f, %f, %f) has zero or negative radius", c.x(), c.y(), c.z());
  }

  if (segments < 1)
  {
    dolfin::dolfin_error("CSGPrimitives3D.cpp",
		 "Create sphere",
		 "Can't create sphere with zero segments");
  }
}
//-----------------------------------------------------------------------------
std::string Sphere::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Sphere with center at " << c << " "
      << "withand radius " << r << ">";
  }
  else
    s << "Sphere(" << c << ", " << r << ")";

  return s.str();
}
//-----------------------------------------------------------------------------
// Box
//-----------------------------------------------------------------------------
Box::Box(dolfin::Point a, dolfin::Point b)
  : a(a), b(b)
{
  if (dolfin::near(a.x(), b.x()) || dolfin::near(a.y(), b.y()) || dolfin::near(a.z(), b.z()))
  {
    std::stringstream s;
    s << "Box with corner " << a.str() << " and " << b.str() << " degenerated";

    dolfin::dolfin_error("CSGPrimitives3D.cpp",
                         "Create axis aligned box",
                         s.str());
  }
}
//-----------------------------------------------------------------------------
std::string Box::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Box with first corner at (" << a.str(true) << ") "
      << "and second corner at (" << b.str(true) << ")>";
  }
  else
  {
    s << "Box(" << a.str(false) << ", " << b.str(false) << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// Cone
//-----------------------------------------------------------------------------
Cylinder::Cylinder(dolfin::Point top, 
                   dolfin::Point bottom, 
                   double top_radius, 
                   double bottom_radius,
                   std::size_t segments)
  : _top(top), _bottom(bottom), _top_radius(top_radius),
    _bottom_radius(bottom_radius), _segments(segments)
{
  if (dolfin::near(top_radius, 0.0) && dolfin::near(bottom_radius, 0.0))
  {
    dolfin::dolfin_error("CSGPrimitives3D.cpp",
		   "Create cylinder",
		   "Cylinder with zero thickness");
  }

  if (top.distance(bottom) < DOLFIN_EPS)
  {
    dolfin::dolfin_error("CSGPrimitives3D.cpp",
		 "Create cylinder",
		 "Cylinder with zero length");
  }
}
//-----------------------------------------------------------------------------
std::string Cylinder::str(bool verbose) const
{
  std::stringstream s;
  if (verbose)
  {
    s << "<Cylinder with top at " << _top << ", top radius " << _top_radius
      << " and bottom at " << _bottom << ", bottom radius "
      << _bottom_radius << ", with " << _segments << " segments>";
  }
  else
  {
    s << "Cylinder( " << _top << ", " << _bottom << ", " << _top_radius
      << ", " << _bottom_radius << " )";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
Tetrahedron::Tetrahedron(dolfin::Point a,
                         dolfin::Point b, 
                         dolfin::Point c, 
                         dolfin::Point d)
  : a(a), b(b), c(c), d(d)
{
  // TODO: Check validity of coordinates
}
//-----------------------------------------------------------------------------
/// Informal string representation
std::string Tetrahedron::str(bool verbose) const
{
  std::stringstream s;
  if (verbose)
  {
    s << "<Tetrahedron with points at " << a << ", " << b << ", "
      << c << ", " << d << ">";
  }
  else
  {
    s << "Tetrahedron( " << a << ", " << b << ", " << c << ", "
      << d << ")";
  }
  return s.str();
}
//-----------------------------------------------------------------------------
Surface3D::Surface3D(std::string filename, double degenerate_tolerance) 
 : _filename(filename), degenerate_tolerance(degenerate_tolerance)
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
