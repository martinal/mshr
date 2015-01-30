// Copyright (C) 2012 Anders Logg and 2013-2014 Benjamin Kehlet
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


#ifndef __MSHR_GEOMETRY_H
#define __MSHR_GEOMETRY_H

#include <memory>
#include <cstddef>
#include <vector>
#include <list>

#include <dolfin/common/Variable.h>

namespace mshr
{

  /// Geometry described by Constructive Solid Geometry (CSG)
  class CSGGeometry : public dolfin::Variable
  {
   protected:
    /// Constructor
    CSGGeometry();

    /// Destructor
    virtual ~CSGGeometry();

   public:
    /// Return dimension of geometry
    virtual std::size_t dim() const = 0;

    /// Informal string representation
    virtual std::string str(bool verbose) const = 0;

    /// Define subdomain. This feature is 2D only.
    /// The subdomain is itself a CSGGeometry and the corresponding
    /// cells in the resulting will be marked with i
    /// If subdomains overlap, the latest added will take precedence.
    void set_subdomain(std::size_t i, std::shared_ptr<CSGGeometry> s);

    /// Define subdomain. This feature is 2D only.
    /// The subdomain is itself a CSGGeometry and the corresponding
    /// cells in the resulting will be marked with i
    /// If subdomains overlap, the latest added will take precedence.
    void set_subdomain(std::size_t i, CSGGeometry& s);

    /// @brief Has subdomains been set
    bool has_subdomains() const;

    /// @brief Return const list of subdomain geometries
    const std::list<std::pair<std::size_t, std::shared_ptr<const CSGGeometry> > >& get_subdomains() const { return subdomains; }

    enum Type { Box, Sphere, Cylinder, Tetrahedron, Ellipsoid, Surface3D, Circle, Ellipse, Rectangle, Polygon, Union, Intersection, Difference, Translation, Scaling, Rotation };

    virtual Type getType() const = 0;
    virtual bool is_operator() const = 0;

   private :
    std::list<std::pair<std::size_t, std::shared_ptr<const CSGGeometry> > > subdomains;
  };
}

#endif
