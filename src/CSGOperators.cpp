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

#include <sstream>

namespace mshr
{

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
  : g(g), c(c), s(s), translate(true), scale_around_center(false)
{
  assert(g);

  dim_ = g->dim();
}
//-----------------------------------------------------------------------------
CSGScaling::CSGScaling(std::shared_ptr<CSGGeometry> g,
                       double s,
                       bool translate)
  : g(g), c(0,0,0), s(s), translate(translate), scale_around_center(true)
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

      if (!scale_around_center)
        ss << "\naround " << c.str(true);

      ss << "\n}";
  }
  else
  {
    ss << "(" << g->str(false) << " * " << std::to_string(s);
    if (!scale_around_center)
      ss << "(" << c.str(true) << ")";
    ss << ")";
  }

  return ss.str();
}
//-----------------------------------------------------------------------------

}
