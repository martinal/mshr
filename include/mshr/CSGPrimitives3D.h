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

  /// @brief Base class for 3D primitives
  class CSGPrimitive3D : public CSGPrimitive
  {
  public:

    /// @return get dimension of geometry
    std::size_t dim() const { return 3; }

  };

  /// @brief A 3D sphere
  ///
  /// { "small-icon" : "sphere-small.png" }
  class Sphere : public CSGPrimitive3D
  {
  public:

    /// @brief Create sphere at c with radius r.
    ///
    /// @param center center of sphere
    /// @param radius radius of sphere
    /// @param segments resolution when generating a polyhedral appoximation
    Sphere(dolfin::Point center, double radius, std::size_t segments=32);

    /// @brief Informal string representation
    /// @param verbose  Verbosity level
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Sphere; }

    const dolfin::Point c;
    const double r;
    const std::size_t _segments;
  };

  /// @brief A 3D axis aligned box
  ///
  /// { "small-icon" : "box-small.png" }
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

    Type getType() const { return CSGGeometry::Box; }

    dolfin::Point a, b;
  };

  /// @brief A 3D cylinder
  ///
  /// { "small-icon" : "cylinder-small.png" }
  class Cylinder : public CSGPrimitive3D
  {
  public:

    /// @brief Create cylinder defined by upper and lower center
    /// and radius respectively.
    ///
    /// @param top           Center at top of cylinder.
    /// @param bottom        Center at bottom of cylinder.
    /// @param top_radius    Radius top of cylinder.
    /// @param bottom_radius Radius at botoom of cylinder.  
    /// @param segments      number of faces on the side when generating a polyhedral approximation.
    Cylinder(dolfin::Point top, 
             dolfin::Point bottom, 
             double top_radius, 
             double bottom_radius,
             std::size_t segments=32);

    /// @brief  Informal string representation
    /// @param  verbose Verbosity level
    /// @return The description string
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Cylinder; }

    const dolfin::Point _top, _bottom;
    const double _top_radius, _bottom_radius;
    const std::size_t _segments;
  };

  /// @brief A 3D cone. 
  /// A cone is here just a special case of a cylinder.
  ///
  /// { "small-icon" : "cone-small.png" }
  class Cone : public Cylinder
  {
  public:

    /// @brief Create cone defined by upper and lower center and bottom radius respectively.
    /// @param top      Center at top of cone.
    /// @param bottom   Center at top of cone.
    /// @param r        bottom radius of cone.
    /// @param segments number of faces on the side when generating a polyhedral approximation.
    Cone(dolfin::Point top, dolfin::Point bottom, double r, std::size_t segments=32)
      : Cylinder(top, bottom, r, r, segments) {}
  };

  /// @brief A 3D tetrahedron
  ///
  /// { "small-icon" : "tetrahedron-small.png" }
  class Tetrahedron : public CSGPrimitive3D
  {
  public:
    /// @brief Create tetrahedron defined by four corner points.
    ///
    /// @param a Point
    /// @param b Point
    /// @param c Point
    /// @param d Point
    Tetrahedron(dolfin::Point a,
                dolfin::Point b,
                dolfin::Point c,
                dolfin::Point d);

    /// @brief Informal string representation
    /// @return The description string
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Tetrahedron; }

    dolfin::Point a, b, c, d;
  };

  /// @brief A triangular 3D surface read from file.
  ///
  /// { "small-icon" : "disk-small.png" }
  class Surface3D : public CSGPrimitive3D
  {
  public:
    Surface3D(std::string filename, double degenerate_tolerance=1e-12);

    /// @brief Informal string representation
    /// @return The description string
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Surface3D; }

    std::string _filename;
    double degenerate_tolerance;
  };
}
#endif
