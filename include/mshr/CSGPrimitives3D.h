// Copyright (C) 2012 Anders Logg, 2012-2015 Benjamin Kehlet
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
// Modified by Johannes Ring 2012

#ifndef __MSHR_PRIMITIVES_3D_H
#define __MSHR_PRIMITIVES_3D_H

#include "CSGPrimitive.h"
#include "CSGPrimitives2D.h"
#include <dolfin/mesh/Mesh.h>
#include <dolfin/geometry/Point.h>
#include <cstddef>

namespace mshr
{

  /// @brief Base class for 3D primitives
  class CSGPrimitive3D : public CSGPrimitive
  {
  protected:
    CSGPrimitive3D();

  public:
    /// @return get dimension of geometry
    std::size_t dim() const { return 3; }

  };

  // TODO: Make sphere a special case of ellipsoid.
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
    Sphere(dolfin::Point center, double radius, std::size_t segments=10);

    /// @brief Informal string representation
    /// @param verbose  Verbosity level
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Sphere; }

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

    const dolfin::Point a, b;
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
      : Cylinder(top, bottom, r, 0, segments) {}
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

    const dolfin::Point a, b, c, d;
  };

  /// @brief A triangular 3D surface read from file.
  ///
  /// { "small-icon" : "disk-small.png" }
  class Surface3D : public CSGPrimitive3D
  {
  public:
    Surface3D(std::string filename);

    // Create triangulated polyhedron from surface of mesh
    Surface3D(std::shared_ptr<const dolfin::Mesh> mesh);

    /// @brief Informal string representation
    /// @return The description string
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Surface3D; }

    const std::string _filename;
    std::shared_ptr<const dolfin::Mesh> mesh;

    /// @brief Tolerance when merging close vertices
    double vertex_tolerance;

    /// @brief Tolerance when removing degenerate facets
    double degenerate_tolerance;

    /// @brief Attempt to repair if surface is not topologically valid
    bool repair;
    
    // @brief Read only one connected_component. Only relevant if repair==true
    bool single_connected_component;

    int sharp_features_filter;
    
    /// @brief First facet, when reading only one connect component.
    std::size_t first_facet;
  };

  /// @brief An axis-aligned ellipsoid
  class Ellipsoid : public CSGPrimitive3D
  {
   public:
    /// @brief Create axis aligned ellipsoid
    ///
    /// @param center center of ellipsoid
    /// @param a semi-principal axis in x direction
    /// @param b semi-principal axis in y direction
    /// @param c semi-principal axis in z direction
    /// @param segments resolution when generating a polyhedral appoximation
    Ellipsoid(dolfin::Point center, double a, double b, double c, std::size_t segments=17);

    /// @brief Informal string representation
    /// @return The description string
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Ellipsoid; }

    const dolfin::Point center;
    const double a, b, c;
    const std::size_t _segments;
  };

  /// @brief A 2D polygon extruded alogn the z axis to 3D
  class Extrude2D : public CSGPrimitive3D
  {
   public :
    Extrude2D(std::shared_ptr<CSGGeometry>, double z);

    /// @brief Informal string representation
    /// @return The description string
    std::string str(bool verbose) const;

    Type getType() const
    { return CSGGeometry::Extrude2D; }

    std::shared_ptr<CSGGeometry> geometry_2d;
    const double z;
  };
}

#endif
