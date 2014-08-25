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
// Modified by Johannes Ring, 2012
// Modified by Benjamin Kehlet, 2012-2013


#include <mshr/CSGPrimitives2D.h>
#include <dolfin/math/basic.h>
#include <dolfin/log/LogStream.h>

#include <sstream>
#include <limits>
#include <algorithm>


namespace mshr
{
//-----------------------------------------------------------------------------
CSGPrimitive2D::CSGPrimitive2D()
{}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Circle
//-----------------------------------------------------------------------------
Circle::Circle(dolfin::Point(c), double r, std::size_t segments)
  : c(c), _r(r), _segments(segments)
{
  if (_r < DOLFIN_EPS)
  {
    std::stringstream s;
    s << "Circle with center " << c.str() << " has zero or negative radius";
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                         "create circle",
                         s.str());
  }

  if (_segments < 3)
  {
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                         "create circle",
                         "Unable to create circle with less than 3 segments");
  }
}
//-----------------------------------------------------------------------------
std::string Circle::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Circle at (" << c.str() << ") with radius " << _r << ">";
  }
  else
  {
    s << "Circle(" << c.str() << ", " << _r << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// Ellipse
//-----------------------------------------------------------------------------
Ellipse::Ellipse(dolfin::Point c, double a, double b,
                 std::size_t segments)
  : c(c), _a(a), _b(b), _segments(segments)
{
  if (_a < DOLFIN_EPS || _b < DOLFIN_EPS)
  {
    std::stringstream s;
    s << "Ellipse with center " << c.str() << " has invalid semi-axis";
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                         "create ellipse",
                         s.str());
  }

  if (_segments < 3)
  {
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create ellipse",
                 "Unable to create ellipse with less than 3 segmentss");
  }
}
//-----------------------------------------------------------------------------
std::string Ellipse::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Ellipse centered at (" << c.str() << ") with horizontal semi-axis "
      << _a << " and vertical semi-axis " << _b << ">";
  }
  else
  {
    s << "Ellipse(" << c.str() << ", " << _a << ", " << _b << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// Rectangle
//-----------------------------------------------------------------------------
Rectangle::Rectangle(dolfin::Point a, dolfin::Point b)
  : a(a), b(b)
{
  if (dolfin::near(a.x(), b.x()) || dolfin::near( a.y(), b.y()))
  {
    std::stringstream s;
    s << "Rectangle with corner " << a.str() << " and " << b.str() << " degenerated";
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create rectangle",
                         s.str());
  }
}
//-----------------------------------------------------------------------------
std::string Rectangle::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Rectangle with first corner at (" << a.str() << ") "
      << "and second corner at (" << b.str() << ")>";
  }
  else
  {
    s << "Rectangle( (" << a.str() << "), (" << b.str() << ") )";
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
