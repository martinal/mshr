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
// Modified by Benjamin Kehlet, 2013

#include <mshr/CSGOperators.h>

#include <dolfin/common/utils.h>
#include <dolfin/log/log.h>
#include <dolfin/common/constants.h>

#include <sstream>

namespace mshr
{

//-----------------------------------------------------------------------------
// CSGUnion
//-----------------------------------------------------------------------------
CSGOperator::CSGOperator()
{}

//-----------------------------------------------------------------------------
// CSGUnion
//-----------------------------------------------------------------------------
CSGUnion::CSGUnion(std::shared_ptr<CSGGeometry> g0,
                   std::shared_ptr<CSGGeometry> g1)
  : _g0(g0), _g1(g1)
{
  assert(g0);
  assert(g1);

  // Check dimensions
  if (g0->dim() != g1->dim())
  {
    dolfin::dolfin_error("CSGOperators.cpp",
                         "create union of CSG geometries",
                         "Dimensions of geomestries don't match (%d vs %d)",
                         g0->dim(), g1->dim());
  }

  dim_ = g0->dim();
}
//-----------------------------------------------------------------------------
std::string CSGUnion::str(bool verbose) const
{
  assert(_g0);
  assert(_g1);

  std::stringstream s;

  if (verbose)
  {
    s << "<Union>\n"
      << "{\n"
      << dolfin::indent(_g0->str(true))
      << "\n"
      << dolfin::indent(_g1->str(true))
      << "\n}";
  }
  else
  {
    s << "(" << _g0->str(false) << " + " << _g1->str(false) << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// CSGDifference
//-----------------------------------------------------------------------------
CSGDifference::CSGDifference(std::shared_ptr<CSGGeometry> g0,
			     std::shared_ptr<CSGGeometry> g1)
  : _g0(g0), _g1(g1)
{
  assert(g0);
  assert(g1);

  // Check dimensions
  if (g0->dim() != g1->dim())
  {
    dolfin::dolfin_error("CSGOperators.cpp",
                         "create difference of CSG geometries",
                         "Dimensions of geomestries don't match (%d vs %d)",
                         g0->dim(), g1->dim());
  }

  dim_ = g0->dim();
}
//-----------------------------------------------------------------------------
std::string CSGDifference::str(bool verbose) const
{
  assert(_g0);
  assert(_g1);

  std::stringstream s;

  if (verbose)
  {
    s << "<Difference>\n"
      << "{\n"
      << dolfin::indent(_g0->str(true))
      << "\n"
      << dolfin::indent(_g1->str(true))
      << "\n}";
  }
  else
  {
    s << "(" << _g0->str(false) << " - " << _g1->str(false) << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// CSGIntersection
//-----------------------------------------------------------------------------
CSGIntersection::CSGIntersection(std::shared_ptr<CSGGeometry> g0,
                                 std::shared_ptr<CSGGeometry> g1)
  : _g0(g0), _g1(g1)
{
  assert(g0);
  assert(g1);

  // Check dimensions
  if (g0->dim() != g1->dim())
  {
    dolfin::dolfin_error("CSGOperators.cpp",
                         "create intersection of CSG geometries",
                         "Dimensions of geomestries don't match (%d vs %d)",
                         g0->dim(), g1->dim());
  }

  dim_ = g0->dim();
}
//-----------------------------------------------------------------------------
std::string CSGIntersection::str(bool verbose) const
{
  assert(_g0);
  assert(_g1);

  std::stringstream s;

  if (verbose)
  {
    s << "<Intersection>\n"
      << "{\n"
      << dolfin::indent(_g0->str(true))
      << "\n"
      << dolfin::indent(_g1->str(true))
      << "\n}";
  }
  else
  {
    s << "(" << _g0->str(false) << " * " << _g1->str(false) << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// CSGTranslation
//-----------------------------------------------------------------------------
CSGTranslation::CSGTranslation(std::shared_ptr<CSGGeometry> g,
                               dolfin::Point t)
  : g(g), t(t) 
{
  assert(g);

  dim_ = g->dim();
}
//-----------------------------------------------------------------------------
std::string CSGTranslation::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Translation>\n"
      << "{\n"
      << dolfin::indent(g->str(true) + "\nby\n" + t.str(true))
      << "\n}";
  }
  else
  {
    s << "(" << g->str(false) << " + " << t.str(false) << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
// CSGScaling
//-----------------------------------------------------------------------------
CSGScaling::CSGScaling(std::shared_ptr<CSGGeometry> g,
                       dolfin::Point c,
                       double s)
  : g(g), c(c), s(s), translate(true)
{
  assert(g);

  dim_ = g->dim();
}
//-----------------------------------------------------------------------------
CSGScaling::CSGScaling(std::shared_ptr<CSGGeometry> g,
                       double s)
  : g(g), c(0,0,0), s(s), translate(false)
{
  assert(g);

  dim_ = g->dim();
}
//-----------------------------------------------------------------------------
std::string CSGScaling::str(bool verbose) const
{
  std::stringstream ss;

  if (verbose)
  {
    ss << "<Scaling>\n"
      << "{\n"
      << dolfin::indent(g->str(true) + "\nby\n" + std::to_string(s));

      if (translate)
        ss << "\naround " << c.str(true);

      ss << "\n}";
  }
  else
  {
    ss << "(" << g->str(false) << " * " << std::to_string(s);
    if (translate)
      ss << "(" << c.str(true) << ")";
    ss << ")";
  }

  return ss.str();
}
//-----------------------------------------------------------------------------
// CSGRotation
//-----------------------------------------------------------------------------
CSGRotation::CSGRotation(std::shared_ptr<CSGGeometry> g,
                         double theta)
  : g(g), rot_axis(.0,.0), c(.0,.0), theta(theta), translate(false)
{
  if (g->dim() > 2)
    dolfin::dolfin_error("CSGOperators.cpp",
                         "Constructing CSG rotation",
                         "Rotation axis must be given in 3D");
}
//-----------------------------------------------------------------------------
CSGRotation::CSGRotation(std::shared_ptr<CSGGeometry> g,
                         dolfin::Point v,
                         double theta)
  : g(g),
    rot_axis(v),
    c(v),
    theta(theta),
    translate(dim_ == 2 ? true : false)
{
  assert(g);

  dim_ = g->dim();

  // if (dim_ == 2)
  //   translate = true;
  // else
  //   translate = false;
}
//-----------------------------------------------------------------------------
CSGRotation::CSGRotation(std::shared_ptr<CSGGeometry> g,
                         dolfin::Point rot_axis,
                         dolfin::Point rot_center,
                         double theta)
  : g(g),
    rot_axis(rot_axis),
    c(rot_center),
    theta(theta),
    translate(true)
{
  assert(g);

  dim_ = g->dim();

  if (dim_ < 3)
    dolfin::dolfin_error("CSGOperators.cpp",
                 "Constructing CSG rotation",
                 "Can't give rotation axis for dimension < 3");
}
//-----------------------------------------------------------------------------
std::string CSGRotation::str(bool verbose) const
{
  std::stringstream ss;

  if (verbose)
  {
    ss << "<Rotation>\n"
      << "{\n"
      << dolfin::indent(g->str(true)
                        + (translate ? "\naround "+rot_axis.str(true) : "")
                        + "\nby " + std::to_string(theta/DOLFIN_PI) + " PI");

      ss << "\n}";
  }
  else
  {
    ss << "rotate(" << g->str(false)
       << ", " << std::to_string(theta/DOLFIN_PI) << " PI";

    if (translate)
      ss << ", " << rot_axis.str(true);

    ss << ")";
  }

  return ss.str();
}
//-----------------------------------------------------------------------------
}
