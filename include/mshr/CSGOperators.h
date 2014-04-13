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

#ifndef __MSHR_OPERATORS_H
#define __MSHR_OPERATORS_H

#include "CSGGeometry.h"
#include <dolfin/geometry/Point.h>
#include <dolfin/common/NoDeleter.h>

#include <memory>

namespace mshr
{

  //--- Operator classes (nodes) ---
  class CSGOperator : public CSGGeometry
  {
   public:
    virtual bool is_operator() const { return true; }
    std::size_t dim() const { return dim_;}

   protected :
    std::size_t dim_;
  };

  /// Union of CSG geometries
  class CSGUnion : public CSGOperator
  {
  public:

    /// Create union of two geometries
    CSGUnion(std::shared_ptr<CSGGeometry> g0,
             std::shared_ptr<CSGGeometry> g1);

    /// Informal string representation
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Union; }

    std::shared_ptr<CSGGeometry> _g0;
    std::shared_ptr<CSGGeometry> _g1;
  };

  /// Difference of CSG geometries
  class CSGDifference : public CSGOperator
  {
  public:

    /// Create union of two geometries
    CSGDifference(std::shared_ptr<CSGGeometry> g0,
             std::shared_ptr<CSGGeometry> g1);

    /// Informal string representation
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Difference; }

    std::shared_ptr<CSGGeometry> _g0;
    std::shared_ptr<CSGGeometry> _g1;
  };


  /// Intersection of CSG geometries
  class CSGIntersection : public CSGOperator
  {
  public:

    /// Create intersection of two geometries
    CSGIntersection(std::shared_ptr<CSGGeometry> g0,
                    std::shared_ptr<CSGGeometry> g1);

    /// Informal string representation
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Intersection; }

    std::shared_ptr<CSGGeometry> _g0;
    std::shared_ptr<CSGGeometry> _g1;

  };

  // Translate CSG geometry
  class CSGTranslation : public CSGOperator
  {
    public:

    /// Translate geoemtry by the vector represented by point
    CSGTranslation(std::shared_ptr<CSGGeometry> g,
                   dolfin::Point t);

    /// Informal string representation
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Translation; }

    std::shared_ptr<CSGGeometry> g;
    dolfin::Point t;
  };


  //--- Union operators ---

  /// Create union of two geometries
  inline std::shared_ptr<CSGUnion> operator+(std::shared_ptr<CSGGeometry> g0,
                                        std::shared_ptr<CSGGeometry> g1)
  {
    return std::shared_ptr<CSGUnion>(new CSGUnion(g0, g1));
  }

  /// Create union of two geometries
  inline std::shared_ptr<CSGUnion> operator+(CSGGeometry& g0,
                                        std::shared_ptr<CSGGeometry> g1)
  {
    return reference_to_no_delete_pointer(g0) + g1;
  }

  /// Create union of two geometries
  inline std::shared_ptr<CSGUnion> operator+(std::shared_ptr<CSGGeometry> g0,
                                        CSGGeometry& g1)
  {
    return g0 + reference_to_no_delete_pointer(g1);
  }

  /// Create union of two geometries
  inline std::shared_ptr<CSGUnion> operator+(CSGGeometry& g0,
                                        CSGGeometry& g1)
  {
    return reference_to_no_delete_pointer(g0) + reference_to_no_delete_pointer(g1);
  }

  //--- Difference operators ---

  /// Create difference of two geometries
  inline  std::shared_ptr<CSGDifference> operator-(std::shared_ptr<CSGGeometry> g0,
					     std::shared_ptr<CSGGeometry> g1)
  {
    return std::shared_ptr<CSGDifference>(new CSGDifference(g0, g1));
  }

  /// Create difference of two geometries
  inline std::shared_ptr<CSGDifference> operator-(CSGGeometry& g0,
					     std::shared_ptr<CSGGeometry> g1)
  {
    return reference_to_no_delete_pointer(g0) - g1;
  }

  /// Create union of two geometries
  inline std::shared_ptr<CSGDifference> operator-(std::shared_ptr<CSGGeometry> g0,
					     CSGGeometry& g1)
  {
    return g0 - reference_to_no_delete_pointer(g1);
  }

  /// Create difference of two geometries
  inline std::shared_ptr<CSGDifference> operator-(CSGGeometry& g0,
					     CSGGeometry& g1)
  {
    return reference_to_no_delete_pointer(g0) - reference_to_no_delete_pointer(g1);
  }


  //--- Intersection operators ---

  /// Create intersection  of two geometries
  inline std::shared_ptr<CSGIntersection> operator*(std::shared_ptr<CSGGeometry> g0,
                                               std::shared_ptr<CSGGeometry> g1)
  {
    return std::shared_ptr<CSGIntersection>(new CSGIntersection(g0, g1));
  }

  /// Create intersection of two geometries
  inline std::shared_ptr<CSGIntersection> operator*(CSGGeometry& g0,
                                               std::shared_ptr<CSGGeometry> g1)
  {
    return reference_to_no_delete_pointer(g0) * g1;
  }

  /// Create intersection of two geometries
  inline std::shared_ptr<CSGIntersection> operator*(std::shared_ptr<CSGGeometry> g0,
                                               CSGGeometry& g1)
  {
    return g0 * reference_to_no_delete_pointer(g1);
  }

  /// Create intersection of two geometries
  inline std::shared_ptr<CSGIntersection> operator*(CSGGeometry& g0,
                                               CSGGeometry& g1)
  {
    return reference_to_no_delete_pointer(g0) * reference_to_no_delete_pointer(g1);
  }

  //--- Translation operators ---
  inline std::shared_ptr<CSGTranslation> operator+(std::shared_ptr<CSGGeometry> g,
                                                   dolfin::Point t)
  {
    return std::shared_ptr<CSGTranslation>(new CSGTranslation(g, t));
  }
}

#endif
