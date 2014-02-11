// Copyright (C) 2012 Anders Logg, 2012-2014 Benjamin Kehlet
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


#ifndef __MSHR_PRIMITIVES_2D_H
#define __MSHR_PRIMITIVES_2D_H

#include "CSGPrimitive.h"

#include <dolfin/geometry/Point.h>
#include <vector>



namespace mshr
{

  /// Base class for 2D primitives
  class CSGPrimitive2D : public CSGPrimitive
  {
  public:

    /// Return dimension of geometry
    std::size_t dim() const { return 2; }

  };

  /// This class describes a 2D circle which can be used to build
  /// geometries using Constructive Solid Geometry (CSG).
  class Circle : public CSGPrimitive2D
  {
  public:

    /// Create circle at x = (x0, x1) with radius r.
    ///
    /// *Arguments*
    ///     x0 (double)
    ///         x0-coordinate of center.
    ///     x1 (double)
    ///         x1-coordinate of center.
    ///     r (double)
    ///         radius.
    ///     fragments (std::size_t)
    ///         number of fragments.
    Circle(double x0, double x1, double r, std::size_t fragments=32);

    /// Informal string representation
    std::string str(bool verbose) const;
    Type getType() const { return CSGGeometry::Circle; }

    /// Return center of circle
    dolfin::Point center() const { return dolfin::Point(_x0, _x1); }

    /// Return radius of circle
    double radius() const { return _r; }

    /// Return number of fragments around the circle
    std::size_t fragments() const { return _fragments; }

  private:

    double _x0, _x1, _r;
    const std::size_t _fragments;

  };

  /// This class describes a 2D ellipse which can be used to build
  /// geometries using Constructive Solid Geometry (CSG).
  class Ellipse : public CSGPrimitive2D
  {
  public:

    /// Create ellipse at x = (x0, x1) with horizontal semi-axis a and
    /// vertical semi-axis b.
    ///
    /// *Arguments*
    ///     x0 (double)
    ///         x0-coordinate of center.
    ///     x1 (double)
    ///         x1-coordinate of center.
    ///     a (double)
    ///         horizontal semi-axis.
    ///     b (double)
    ///         vertical semi-axis.
    ///     fragments (std::size_t)
    ///         number of fragments.
    Ellipse(double x0, double x1, double a, double b, std::size_t fragments=32);

    /// Informal string representation
    std::string str(bool verbose) const;
    Type getType() const { return CSGGeometry::Ellipse; }

    /// Return center of ellipse
    dolfin::Point center() const { return dolfin::Point(_x0, _x1); }

    /// Return horizontal semi-axis
    double a() const { return _a; }

    /// Return vertical semi-axis
    double b() const { return _b; }

    /// Return number of fragments around the ellipse
    std::size_t fragments() const { return _fragments; }

  private:

    double _x0, _x1, _a, _b;
    const std::size_t _fragments;

  };

  /// This class describes a 2D rectangle which can be used to build
  /// geometries using Constructive Solid Geometry (CSG).
  class Rectangle : public CSGPrimitive2D
  {
  public:

    /// Create rectangle defined by two opposite corners
    /// x = (x0, x1) and y = (y0, y1).
    ///
    /// *Arguments*
    ///     x0 (double)
    ///         x0-coordinate of first corner.
    ///     x1 (double)
    ///         x1-coordinate of first corner.
    ///     y0 (double)
    ///         y0-coordinate of second corner.
    ///     y1 (double)
    ///         y1-coordinate of second corner.
    Rectangle(double x0, double x1, double y0, double y1);

    /// Informal string representation
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Rectangle; }

    /// Return first corner
    dolfin::Point first_corner() const { return dolfin::Point(_x0, _x1); }

    /// Return second corner
    dolfin::Point second_corner() const { return dolfin::Point(_y0, _y1); }

  private:

    double _x0, _x1, _y0, _y1;

  };

  /// This class describes a 2D polygon which can be used to build
  /// geometries using Constructive Solid Geometry (CSG).
  class Polygon : public CSGPrimitive2D
  {
  public:

    /// Create polygon defined by the given vertices.
    ///
    /// *Arguments*
    ///     vertices (std::vector<_dolfin::Point_>)
    ///         A vector of _dolfin::Point_ objects.
    ///         The points must be given in counter-clockwise order
    ///         (without repeating the first/last vertex) and the polygon
    //          must not self intersect.
    Polygon(const std::vector<dolfin::Point>& vertices);

    /// Informal string representation
    std::string str(bool verbose) const;
    Type getType() const { return CSGGeometry::Polygon; }

    // Check if vertices are counter clockwise oriented
    bool ccw() const;

    /// Return vertices in polygon
    const std::vector<dolfin::Point>& vertices() const { return _vertices; }

  private:

    const std::vector<dolfin::Point> _vertices;

  };
}

#endif