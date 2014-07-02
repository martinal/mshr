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

#ifndef __MSHR_PRIMITIVES_3D_H
#define __MSHR_PRIMITIVES_3D_H

#include <cstddef>
#include <dolfin/geometry/Point.h>
#include "CSGPrimitive.h"

namespace mshr
{

  /// Base class for 3D primitives
  class CSGPrimitive3D : public CSGPrimitive
  {
  public:

    /// Return dimension of geometry
    std::size_t dim() const { return 3; }

  };

  /// This class describes a 3D sphere which can be used to build
  /// geometries using Constructive Solid Geometry (CSG).
  class Sphere : public CSGPrimitive3D
  {
  public:

    /// Create sphere at c with radius r.
    ///
    /// *Arguments*
    ///     c  (dolfin::Point)
    ///         center.
    ///     r (double)
    ///         radius.
    Sphere(dolfin::Point center, double radius, std::size_t slices=16);

    /// Informal string representation
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Sphere; }

    const dolfin::Point c;
    const double r;
    const std::size_t _slices;

  };

  /// @brief This class describes a 3D box which can be used to build
  ///        geometries using Constructive Solid Geometry (CSG).
  class Box : public CSGPrimitive3D
  {
  public:

    /// @brief Create box defined by two opposite corners
    ///
    /// @param a The first corner
    /// @param b The second corner
    Box(dolfin::Point a, dolfin::Point b);

    /// @brief Informal string representation
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Box; }

    dolfin::Point a, b;
  };

  /// This class describes a 3D cone which can be used to build
  /// geometries using Constructive Solid Geometry (CSG).
  class Cone : public CSGPrimitive3D
  {
  public:

    /// Create cone defined by upper and lower center
    /// and radius respectively.
    ///
    /// *Arguments*
    ///     top (dolfin::Point)
    ///         Center at top of cone.
    ///     top_radius(double)
    ///         Radius bottom of cone.
    ///     bottom(dolfin::Point)
    ///         Center at top of cone.
    ///     bottom_radius (double)
    ///         radius at top of cone.
    ///     slices (std::size_t)
    ///         number of faces on the side when generating a
    ///         polyhedral approximation.
    Cone(dolfin::Point top, 
         dolfin::Point bottom, 
         double top_radius, 
         double bottom_radius,
         std::size_t slices=32);

    /// Informal string representation
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Cone; }

    const dolfin::Point _top, _bottom;
    const double _top_radius, _bottom_radius;
    const std::size_t _slices;
  };

  /// This class describes a 3D cylinder which can be used to build
  /// geometries using Constructive Solid Geometry (CSG). A cylinder
  /// is here just a special case of a cone.
  class Cylinder : public Cone
  {
  public:

    /// Create cylinder defined by upper and lower center
    /// and radius respectively.
    ///
    /// *Arguments*
    ///     top (dolfin::Point)
    ///         Center at top of cylinder.
    ///     bottom(dolfin::Point)
    ///         Center at top of cylinder.
    ///     r (double)
    ///         radius of cylinder.
    ///     slices (std::size_t)
    ///         number of faces on the side when generating a
    ///         polyhedral approximation.
    Cylinder(dolfin::Point top, dolfin::Point bottom, double r, std::size_t slices=32)
      : Cone(top, bottom, r, r, slices) {}
  };

  /// This class describes a Tetrahedron which can be used to build
  /// geometries using Constructive Solid Geometry (CSG).
  class Tetrahedron : public CSGPrimitive3D
  {
  public:
    /// Create tetrahedron defined by four corner points.
    ///
    /// *Arguments*
    ///     x0 (Point)
    ///         Point.
    ///     x1 (Point)
    ///         Point.
    ///     x2 (Point)
    ///         Point.
    ///     x3 (Point)
    ///         Point.
    Tetrahedron(dolfin::Point a,
                dolfin::Point b,
                dolfin::Point c,
                dolfin::Point d);

    /// Informal string representation
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Tetrahedron; }

    dolfin::Point a, b, c, d;
  };

  /// This class describes a 3D surface loaded from file.
  /// The supported file types
  class Surface3D : public CSGPrimitive3D
  {
  public:
    Surface3D(std::string filename);

    /// Informal string representation
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Surface3D; }

    std::string _filename;
  };
}
#endif
