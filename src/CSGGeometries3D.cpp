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
// Modified by Johannes Ring 2014
// Modified by Anders Logg 2014

#include <mshr/CSGGeometries3D.h>
#include <mshr/CSGGeometry.h>
#include <mshr/CSGPrimitives3D.h>
#include <mshr/CSGOperators.h>

#include <cmath>

namespace mshr
{

// FIXME: Move somewhere else
// Convenience function for rotation
std::shared_ptr<CSGGeometry> rotate(std::shared_ptr<CSGGeometry> geometry,
                                     dolfin::Point rot_axis,
                                     dolfin::Point rot_center,
                                     double theta)
{
  return std::shared_ptr<CSGGeometry>(new CSGRotation(geometry,
                                                      rot_axis,
                                                      rot_center,
                                                      theta));
}

//-----------------------------------------------------------------------------
std::shared_ptr<CSGGeometry> CSGGeometries::lego(std::size_t n0,
                                                 std::size_t n1,
                                                 std::size_t n2,
                                                 dolfin::Point x)
{
  // Standard dimensions for LEGO bricks / m
  const double P = 8.0 * 0.001;
  const double h = 3.2 * 0.001;
  const double D = 5.0 * 0.001;
  const double b = 1.7 * 0.001;
  const double d = 0.2 * 0.001;

  // Create brick
  std::shared_ptr<CSGGeometry>
    lego(new Box(x + dolfin::Point(0.5*d, 0.5*d, 0),
                 x + dolfin::Point(n0*P - 0.5*d, n1*P - 0.5*d, n2*h)));

  // Add knobs
  for (std::size_t i = 0; i < n0; i++)
  {
    for (std::size_t j = 0; j < n1; j++)
    {
      const dolfin::Point knop_bottom = x + dolfin::Point( (i + 0.5)*P,
                                                           (j + 0.5)*P,
                                                           0);

      std::shared_ptr<CSGGeometry>
        knob(new Cylinder(knop_bottom,
                          knop_bottom + dolfin::Point(0, 0, n2*h + b),
                          0.5*D,
                          0.5*D));

      lego = lego + knob;
    }
  }

  return lego;
}
//-----------------------------------------------------------------------------
std::shared_ptr<CSGGeometry> CSGGeometries::propeller(double r,
                                                      double R,
                                                      double w,
                                                      double h,
                                                      bool rotate_blades,
                                                      bool include_tip)
{
  // Parameters
  const double v = 0.3;     // rotation of blades
  const double l = 0.8*w;  // length of cone

  // // Create blades
  std::shared_ptr<CSGGeometry>
    blade_0(new Box(dolfin::Point(0.8*r, -0.5*h, -0.5*w),
                    dolfin::Point(    R,  0.5*h,  0.5*w)));
  std::shared_ptr<CSGGeometry>
    blade_1(new Box(dolfin::Point(    -R, -0.5*h, -0.5*w),
                    dolfin::Point(-0.8*r,  0.5*h,  0.5*w)));
  std::shared_ptr<CSGGeometry>
    blade_2(new Box(dolfin::Point(-0.5*h, 0.8*r, -0.5*w),
                    dolfin::Point( 0.5*h,     R,  0.5*w)));
  std::shared_ptr<CSGGeometry>
    blade_3(new Box(dolfin::Point(-0.5*h,     -R, -0.5*w),
                    dolfin::Point( 0.5*h, -0.8*r,  0.5*w)));

  // Rotate blades
  if (rotate_blades)
  {
    blade_0 = rotate(blade_0, dolfin::Point(1, 0, 0), dolfin::Point(0, 0, 0), v);
    blade_1 = rotate(blade_1, dolfin::Point(1, 0, 0), dolfin::Point(0, 0, 0), -v);
    blade_2 = rotate(blade_2, dolfin::Point(0, 1, 0), dolfin::Point(0, 0, 0), v);
    blade_3 = rotate(blade_3, dolfin::Point(0, 1, 0), dolfin::Point(0, 0, 0), -v);
  }

  // Create blade tips
  std::shared_ptr<CSGGeometry>
    blade_tip_0(new Cylinder(dolfin::Point( R, -0.5*h, 0),
                             dolfin::Point( R,  0.5*h, 0), 0.5*w, 0.5*w));
  std::shared_ptr<CSGGeometry>
    blade_tip_1(new Cylinder(dolfin::Point(-R, -0.5*h, 0),
                             dolfin::Point(-R,  0.5*h, 0), 0.5*w, 0.5*w));
  std::shared_ptr<CSGGeometry>
    blade_tip_2(new Cylinder(dolfin::Point(-0.5*h,  R, 0),
                             dolfin::Point( 0.5*h,  R, 0), 0.5*w, 0.5*w));
  std::shared_ptr<CSGGeometry>
    blade_tip_3(new Cylinder(dolfin::Point(-0.5*h, -R, 0),
                             dolfin::Point( 0.5*h, -R, 0), 0.5*w, 0.5*w));

  // Rotate blades
  if (rotate_blades)
  {
    blade_tip_0 = rotate(blade_tip_0, dolfin::Point(1, 0, 0), dolfin::Point(0, 0, 0), v);
    blade_tip_1 = rotate(blade_tip_1, dolfin::Point(1, 0, 0), dolfin::Point(0, 0, 0), -v);
    blade_tip_2 = rotate(blade_tip_2, dolfin::Point(0, 1, 0), dolfin::Point(0, 0, 0), v);
    blade_tip_3 = rotate(blade_tip_3, dolfin::Point(0, 1, 0), dolfin::Point(0, 0, 0), -v);
  }

  // Add blade tips
  blade_0 = blade_0 + blade_tip_0;
  blade_1 = blade_1 + blade_tip_1;
  blade_2 = blade_2 + blade_tip_2;
  blade_3 = blade_3 + blade_tip_3;

  // // Add blades
  std::shared_ptr<CSGGeometry>
    blades = blade_0 + blade_1 + blade_2 + blade_3;

  // Create outer cylinder
  std::shared_ptr<CSGGeometry>
    cylinder_outer(new Cylinder(dolfin::Point(0, 0, -0.5*w),
                                dolfin::Point(0, 0, 0.5*w),
                                r, r));

  // Create inner cylinder
  std::shared_ptr<CSGGeometry>
    cylinder_inner(new Cylinder(dolfin::Point(0, 0, -0.5*w),
                                dolfin::Point(0, 0, 0.5*w),
                                0.5*r, 0.5*r));

  // Create center cone
  std::shared_ptr<CSGGeometry>
    cone(new Cylinder(dolfin::Point(0, 0, -0.5*w),
                      dolfin::Point(0, 0, -0.5*w - l), r, h));

  // Create sphere for tip of cone
  const double k = (r - h) / l;
  const double rc = h*sqrt(1 + k*k);
  const double xc = k*h;
  std::shared_ptr<CSGGeometry>
    tip(new Sphere(dolfin::Point(0, 0, -0.5*w - l + xc), rc));

  // Build propeller from parts
  if (include_tip)
    return cylinder_outer - cylinder_inner + blades + cone + tip;
  else
    return cylinder_outer - cylinder_inner + blades + cone;
}
//-----------------------------------------------------------------------------

}
