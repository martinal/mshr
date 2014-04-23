// Copyright (C) 2013-2014 Benjamin Kehlet
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

#ifndef __MSHR_CSGCGAL_DOMAIN2D_H
#define __MSHR_CSGCGAL_DOMAIN2D_H

#include <mshr/CSGGeometry.h>

#include <dolfin/geometry/Point.h>
#include <memory>

namespace mshr
{
// Forward declaration
struct CSGCGALDomain2DImpl;

class CSGCGALDomain2D : public dolfin::Variable
{
 public:
  // Create empty polygon
  CSGCGALDomain2D();

  // Construct polygon from Dolfin CSG geometry
  CSGCGALDomain2D(const mshr::CSGGeometry *csg);

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
                    double truncate_tolerance) const;

  void get_holes(std::list<std::vector<dolfin::Point> >& h, 
                 double truncate_tolerance) const;

  std::string str(bool verbose) const;

 private:
  std::unique_ptr<CSGCGALDomain2DImpl> impl;

};
}

#endif
