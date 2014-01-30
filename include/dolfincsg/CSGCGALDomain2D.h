// Copyright (C) 2013-2014 Benjamin Kehlet
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

#ifndef __DOLFINCSG_CSGCGAL_DOMAIN2D_H
#define __DOLFINCSG_CSGCGAL_DOMAIN2D_H

#include <dolfincsg/CSGGeometry.h>

#include <dolfin/geometry/Point.h>
#include <boost/scoped_ptr.hpp>


namespace dolfincsg
{
  // Forward declaration
  struct CSGCGALDomain2DImpl;

class CSGCGALDomain2D
{
 public:
  // Create empty polygon
  CSGCGALDomain2D();

  // Construct polygon from Dolfin CSG geometry
  CSGCGALDomain2D(const dolfincsg::CSGGeometry *csg);

  // Destructor
  ~CSGCGALDomain2D();

  // Copy constructor
  CSGCGALDomain2D(const CSGCGALDomain2D &other);
  CSGCGALDomain2D &operator=(const CSGCGALDomain2D &other);

  // Boolean operators
  void join_inplace(const CSGCGALDomain2D& other);
  void intersect_inplace(const CSGCGALDomain2D& other);
  void difference_inplace(const CSGCGALDomain2D& other);

  bool point_in_domain(dolfin::Point p) const;
  double compute_boundingcircle_radius() const ;
  
  // TODO: Replace this with a more C++-ish
  // implementation, ie, take an outputiterator as arugment
  // or define iterator
  void get_vertices(std::list<std::vector<dolfin::Point> >& v, 
                    double truncate_threshold) const;

  void get_holes(std::list<std::vector<dolfin::Point> >& h, 
                 double truncate_threshold) const;

  boost::scoped_ptr<CSGCGALDomain2DImpl> impl;

};
}

#endif
