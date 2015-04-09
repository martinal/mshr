// Copyright (C) 2014 Benjamin Kehlet
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

#ifndef __MSHR_POINT3_FUZZY_STRICTLY_LESS_H_
#define __MSHR_POINT3_FUZZY_STRICTLY_LESS_H_

template<typename Point>
class Point3FuzzyStrictlyLess
{
 public:
  typedef Point Point_;

  Point3FuzzyStrictlyLess(double squared_distance=1e-5)
  : tol(squared_distance){}

  // strictly less
  inline bool operator()(const  Point_ a, const Point_ b) const
  {
    // check distance
    if ( (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]) < tol )
      return false;

    if (a[0] != b[0])
      return a[0] < b[0];

    if (a[1] != b[1])
      return a[1] < b[1];

    //assert a[2] != b[2];
    
    return a[2] < b[2];
  }

  private:
    const double tol;
};

#endif
