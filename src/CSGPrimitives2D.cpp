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
// Modified by Johannes Ring, 2012
// Modified by Benjamin Kehlet, 2012-2013


#include <dolfincsg/CSGPrimitives2D.h>
#include <dolfin/math/basic.h>
#include <dolfin/log/LogStream.h>

#include <sstream>
#include <limits>
#include <algorithm>


namespace dolfincsg
{

//-----------------------------------------------------------------------------
// Circle
//-----------------------------------------------------------------------------
Circle::Circle(double x0, double x1, double r, std::size_t fragments)
  : _x0(x0), _x1(x1), _r(r), _fragments(fragments)
{
  if (_r < DOLFIN_EPS)
  {
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create circle",
                 "Circle with center (%f, %f) has zero or negative radius",
                 _x0, _x1);
  }
  if (_fragments < 3)
  {
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create circle",
                 "Unable to create circle with less than 3 fragments");

  }
}
//-----------------------------------------------------------------------------
std::string Circle::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Circle at (" << _x0 << ", " << _x1 << ") with radius " << _r << ">";
  }
  else
  {
    s << "Circle("
      << _x0 << ", " << _x1 << ", " << _r << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// Ellipse
//-----------------------------------------------------------------------------
Ellipse::Ellipse(double x0, double x1, double a, double b,
                 std::size_t fragments)
  : _x0(x0), _x1(x1), _a(a), _b(b), _fragments(fragments)
{
  if (_a < DOLFIN_EPS || _b < DOLFIN_EPS)
  {
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create ellipse",
                 "Ellipse with center (%f, %f) has invalid semi-axis",
                 _x0, _x1);
  }
  if (_fragments < 3)
  {
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create ellipse",
                 "Unable to create ellipse with less than 3 fragments");

  }
}
//-----------------------------------------------------------------------------
std::string Ellipse::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Ellipse at (" << _x0 << ", " << _x1 << ") with horizontal semi-axis "
      << _a << " and vertical semi-axis " << _b << ">";
  }
  else
  {
    s << "Ellipse("
      << _x0 << ", " << _x1 << ", " << _a << ", " << _b << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// Rectangle
//-----------------------------------------------------------------------------
Rectangle::Rectangle(double x0, double x1, double y0, double y1)
  : _x0(x0), _x1(x1), _y0(y0), _y1(y1)
{
  if (dolfin::near(x0, y0) || dolfin::near(x1, y1))
  {
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create rectangle",
                 "Rectangle with corner (%f, %f) and (%f, %f) degenerated",
                 x0, x1, y0, y1);
  }
}
//-----------------------------------------------------------------------------
std::string Rectangle::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Rectangle with first corner at (" << _x0 << ", " << _x1 << ") "
      << "and second corner at (" << _y0 << ", " << _y1 << ")>";
  }
  else
  {
    s << "Rectangle("
      << _x0 << ", " << _x1 << ", " << _y0 << ", " << _y1 << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// Polygon
//-----------------------------------------------------------------------------
Polygon::Polygon(const std::vector<dolfin::Point>& vertices)
  : _vertices(vertices.begin(), vertices.end())
{
  if (_vertices.size() < 3)
  {
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create polygon",
                 "Polygon should have at least three vertices");
  }

  if (!ccw())
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create polygon",
                 "Polygonv vertices must be given in counter clockwise order");
}
//-----------------------------------------------------------------------------
std::string Polygon::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Polygon with vertices ";
    std::vector<dolfin::Point>::const_iterator p;
    for (p = _vertices.begin(); p != _vertices.end(); ++p)
    {
      s << "(" << p->x() << ", " << p->y() << ")";
      if ((p != _vertices.end()) && (p + 1 != _vertices.end()))
        s << ", ";
    }
    s << ">";
  }
  else
  {
    s << "Polygon(";
    std::vector<dolfin::Point>::const_iterator p;
    for (p = _vertices.begin(); p != _vertices.end(); ++p)
    {
      s << "(" << p->x() << ", " << p->y() << ")";
      if ((p != _vertices.end()) && (p + 1 != _vertices.end()))
        s << ", ";
    }
    s << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
bool Polygon::ccw() const
{
  double signed_area = 0.0;
  
  dolfin::Point prev = _vertices.back();
  for (std::vector<dolfin::Point>::const_iterator it = _vertices.begin(),
	 v_end = _vertices.end();
       it != v_end;
       ++it)
  {
    signed_area += (prev.x()*it->y())-(it->x()*prev.y());
    prev = *it;
  }

  return signed_area > 0;
}
}
