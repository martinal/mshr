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

  /// @brief Base class for csg operators
  /// CSGOperator object are internal (non-leaf) nodes in the CSG tree.
  class CSGOperator : public CSGGeometry
  {
   protected:
    CSGOperator();
   public:
    virtual bool is_operator() const { return true; }
    std::size_t dim() const { return dim_;}

   protected :
    std::size_t dim_;
  };

  /// @brief Union of CSG geometries
  class CSGUnion : public CSGOperator
  {
  public:

    /// @brief Create union of two geometries
    /// @param g0 a CSG geometry
    /// @param g1 a CSG geometry
    CSGUnion(std::shared_ptr<CSGGeometry> g0,
             std::shared_ptr<CSGGeometry> g1);

    /// @brief get informal string representation
    /// @param verbose vervosity level
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Union; }

    const std::shared_ptr<CSGGeometry> _g0;
    const std::shared_ptr<CSGGeometry> _g1;
  };

  /// @brief Difference of CSG geometries
  class CSGDifference : public CSGOperator
  {
  public:

    /// @brief Create union of two geometries
    /// @param g0 a CSG geometry
    /// @param g1 a CSG CSG geometry
    CSGDifference(std::shared_ptr<CSGGeometry> g0,
                  std::shared_ptr<CSGGeometry> g1);

    /// @brief get informal string representation
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Difference; }

    const std::shared_ptr<CSGGeometry> _g0;
    const std::shared_ptr<CSGGeometry> _g1;
  };


  /// @brief Intersection of CSG geometries
  class CSGIntersection : public CSGOperator
  {
  public:

    /// @brief Create intersection of two geometries
    /// @param g0 a CSG geometry
    /// @param g1 a CSG geometry
    CSGIntersection(std::shared_ptr<CSGGeometry> g0,
                    std::shared_ptr<CSGGeometry> g1);

    /// @brief get informal string representation
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Intersection; }

    const std::shared_ptr<CSGGeometry> _g0;
    const std::shared_ptr<CSGGeometry> _g1;

  };

  /// @brief Translate CSG geometry by vector
  ///
  /// { 'small-icon' : 'translation-small.png' }
  class CSGTranslation : public CSGOperator
  {
    public:

    /// @brief create translation 
    /// @param g a CSG geometry
    /// @param t the translation vector
    CSGTranslation(std::shared_ptr<CSGGeometry> g,
                   dolfin::Point t);

    /// @brief get informal string representation
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Translation; }

    const std::shared_ptr<CSGGeometry> g;
    const dolfin::Point t;
  };

  /// @brief Scale CSG geometry
  ///
  /// { 'small-icon' : 'scaling-small.png' }
  class CSGScaling : public CSGOperator
  {
   public:

    /// @brief Create scaling
    /// @param g a CSG geometry
    /// @param scale_factor 
    CSGScaling(std::shared_ptr<CSGGeometry> g,
               double scale_factor);

    /// @brief Scale (translated) geometry. Geometry will be translated, scaled and translated back
    /// @param g a CSG geometry
    /// @param c translation 
    /// @param scale_factor the scale factor
    CSGScaling(std::shared_ptr<CSGGeometry> g,
               dolfin::Point c,
               double scale_factor);



    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Scaling; }

    const std::shared_ptr<CSGGeometry> g;
    const dolfin::Point c;
    const double s;
    const bool translate;
  };

  /// @brief Rotate CSG geometry
  ///
  /// { 'small-icon' : 'rotation-small.png' }
  class CSGRotation : public CSGOperator
  {
   public:
    /// @brief Create 2D rotation (2D only)
    /// @param g A CSG geometry
    /// @param theta Rotate by theta
    CSGRotation(std::shared_ptr<CSGGeometry> g,
                             double theta);

    /// @brief Create rotation
    /// @param g A CSG geometry
    /// @param v In 2D: the rotation center. In 3D: the rotation axis.
    /// @param theta Radians to rotate.
    CSGRotation(std::shared_ptr<CSGGeometry> g,
                dolfin::Point v,
                double theta);

    /// @brief create 3D rotation 
    /// @param g A CSG geometry
    /// @param rot_axis The rotation axis
    /// @param rot_center The rotation center
    /// @param theta Radians to rotate
    CSGRotation(std::shared_ptr<CSGGeometry> g,
                dolfin::Point rot_axis,
                dolfin::Point rot_center,
                double theta);
    Type getType() const { return CSGGeometry::Rotation; }
    std::string str(bool verbose) const;

    const std::shared_ptr<CSGGeometry> g;
    const dolfin::Point rot_axis;
    const dolfin::Point c;
    const double theta;
    const bool translate;
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

  //--- Scaling operators ---
  inline std::shared_ptr<CSGScaling> operator*(std::shared_ptr<CSGGeometry> g,
                                                   double s)
  {
    return std::shared_ptr<CSGScaling>(new CSGScaling(g, s, false));
  }

}

#endif
