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

// This file must be included to get the compiler flags
// They should idellay have been added via the command line,
// but since CGAL configure time is at mshr compile time, we don't
// have access to CGALConfig.cmake at configure time...
#include <CGAL/compiler_config.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

#include <mshr/CSGCGALDomain2D.h>
#include <mshr/CSGPrimitives2D.h>
#include <mshr/CSGOperators.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/LogStream.h>

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_set_2.h>

#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>

#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>

// Polygon typedefs
//typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
typedef CGAL::Quotient<CGAL::MP_Float>                    Quotient;
typedef CGAL::Cartesian<Quotient>                         Exact_Kernel;

typedef Exact_Kernel::Point_2                             Point_2;
typedef Exact_Kernel::Vector_2                            Vector_2;
typedef Exact_Kernel::Segment_2                           Segment_2;
typedef CGAL::Polygon_2<Exact_Kernel>                     Polygon_2;
typedef Polygon_2::Vertex_const_iterator                  Vertex_const_iterator;
typedef CGAL::Polygon_with_holes_2<Exact_Kernel>          Polygon_with_holes_2;
typedef Polygon_with_holes_2::Hole_const_iterator         Hole_const_iterator;
typedef CGAL::Polygon_set_2<Exact_Kernel>                 Polygon_set_2;

typedef CGAL::Snap_rounding_traits_2<Exact_Kernel>        snap_rounding_traits;
typedef std::list<Point_2>                                Polyline_2;
typedef std::list<Polyline_2>                             Polyline_list_2;

// Min enclosing circle typedefs
typedef CGAL::Min_circle_2_traits_2<Exact_Kernel>  Min_Circle_Traits;
typedef CGAL::Min_circle_2<Min_Circle_Traits>      Min_circle;
typedef CGAL::Circle_2<Exact_Kernel> CGAL_Circle;

namespace mshr
{

struct CSGCGALDomain2DImpl
{
  Polygon_set_2 polygon_set;

  CSGCGALDomain2DImpl(){}
  CSGCGALDomain2DImpl(const Polygon_set_2& p)
    : polygon_set(p) {}
};
//-----------------------------------------------------------------------------
Polygon_2 make_circle(const Circle* c)
{
  std::vector<Point_2> pts;
  pts.reserve(c->segments());

  for (std::size_t i = 0; i < c->segments(); i++)
  {
    const double phi = (2*DOLFIN_PI*i) / c->segments();
    const double x = c->center().x() + c->radius()*cos(phi);
    const double y = c->center().y() + c->radius()*sin(phi);
    pts.push_back(Point_2(x, y));
  }

  return Polygon_2(pts.begin(), pts.end());
}
//-----------------------------------------------------------------------------
Polygon_2 make_ellipse(const Ellipse* e)
{
  std::vector<Point_2> pts;

  for (std::size_t i = 0; i < e->segments(); i++)
  {
    const double phi = (2*DOLFIN_PI*i) / e->segments();
    const double x = e->center().x() + e->a()*cos(phi);
    const double y = e->center().y() + e->b()*sin(phi);
    pts.push_back(Point_2(x, y));
  }

  return Polygon_2(pts.begin(), pts.end());
}
//-----------------------------------------------------------------------------
Polygon_2 make_rectangle(const Rectangle* r)
{
  const double x0 = std::min(r->first_corner().x(), r->second_corner().x());
  const double y0 = std::min(r->first_corner().y(), r->second_corner().y());

  const double x1 = std::max(r->first_corner().x(), r->second_corner().x());
  const double y1 = std::max(r->first_corner().y(), r->second_corner().y());

  std::vector<Point_2> pts;
  pts.push_back(Point_2(x0, y0));
  pts.push_back(Point_2(x1, y0));
  pts.push_back(Point_2(x1, y1));
  pts.push_back(Point_2(x0, y1));

  Polygon_2 p(pts.begin(), pts.end());
  
  return p;
}
//-----------------------------------------------------------------------------
Polygon_2 make_polygon(const Polygon* p)
{
  std::vector<Point_2> pts;
  std::vector<dolfin::Point>::const_iterator v;
  for (v = p->vertices().begin(); v != p->vertices().end(); ++v)
    pts.push_back(Point_2(v->x(), v->y()));

  return Polygon_2(pts.begin(), pts.end());
}
//-----------------------------------------------------------------------------
std::unique_ptr<CSGCGALDomain2DImpl> do_transformation(const Polygon_set_2& p, Exact_Kernel::Aff_transformation_2 t)
{
  std::unique_ptr<CSGCGALDomain2DImpl> result(new CSGCGALDomain2DImpl);

  std::list<Polygon_with_holes_2> polygon_list;
  p.polygons_with_holes(std::back_inserter(polygon_list));

  std::list<Polygon_with_holes_2>::const_iterator pit;
  for (pit = polygon_list.begin(); pit != polygon_list.end(); ++pit)
  {
    const Polygon_with_holes_2& pwh = *pit;

    // Transform outer boundary
    Polygon_with_holes_2 transformed(CGAL::transform(t, pwh.outer_boundary()));

    // Transform holes
    for (Hole_const_iterator hit = pwh.holes_begin(); hit != pwh.holes_end(); hit++)
    {
      transformed.add_hole(CGAL::transform(t, *hit));
    }

    result->polygon_set.insert(transformed);
  }

  return result;
}
//-----------------------------------------------------------------------------
CSGCGALDomain2D::CSGCGALDomain2D()
  : impl(new CSGCGALDomain2DImpl)
{
  
}
//-----------------------------------------------------------------------------
CSGCGALDomain2D::~CSGCGALDomain2D()
{
}
//-----------------------------------------------------------------------------
CSGCGALDomain2D::CSGCGALDomain2D(const CSGGeometry *geometry)
: impl(new CSGCGALDomain2DImpl)
{

  if (geometry->dim() != 2)
    dolfin::dolfin_error("CSGCGALDomain2D.cpp",
                         "Creating polygonal domain",
                         "Geometry has dimension %d, expected 2", geometry->dim());


  switch (geometry->getType()) 
  {
    case CSGGeometry::Union:
    {
      const CSGUnion *u = dynamic_cast<const CSGUnion*>(geometry);
      dolfin_assert(u);

      CSGCGALDomain2D a(u->_g0.get());
      CSGCGALDomain2D b(u->_g1.get());

      impl.swap(a.impl);
      impl->polygon_set.join(b.impl->polygon_set);    
      break;
    }
    case CSGGeometry::Intersection:
    {
      const CSGIntersection* u = dynamic_cast<const CSGIntersection*>(geometry);
      dolfin_assert(u);

      CSGCGALDomain2D a(u->_g0.get());
      CSGCGALDomain2D b(u->_g1.get());
      
      impl.swap(a.impl);
      impl->polygon_set.intersection(b.impl->polygon_set);
      break;
    }
    case CSGGeometry::Difference:
    {
      const CSGDifference* u = dynamic_cast<const CSGDifference*>(geometry);
      dolfin_assert(u);
      CSGCGALDomain2D a(u->_g0.get());
      CSGCGALDomain2D b(u->_g1.get());
      
      impl.swap(a.impl);
      impl->polygon_set.difference(b.impl->polygon_set);
      break;
    }
    case CSGGeometry::Translation :
    {
      const CSGTranslation* t = dynamic_cast<const CSGTranslation*>(geometry);
      dolfin_assert(t);
      CSGCGALDomain2D a(t->g.get());
      Exact_Kernel::Aff_transformation_2 translation(CGAL::TRANSLATION, Vector_2(t->t.x(), t->t.y()));
      std::unique_ptr<CSGCGALDomain2DImpl> transformed = do_transformation(a.impl->polygon_set, translation);
      impl.swap(transformed);
      break;
    }
    case CSGGeometry::Scaling :
    {
      const CSGScaling* t = dynamic_cast<const CSGScaling*>(geometry);
      dolfin_assert(t);
      CSGCGALDomain2D a(t->g.get());
      Exact_Kernel::Aff_transformation_2 tr(CGAL::IDENTITY);

      // Translate if requested
      if (t->translate)
        tr = Exact_Kernel::Aff_transformation_2 (CGAL::TRANSLATION,
                                                 Vector_2(-t->c.x(), -t->c.y())) * tr;

      // Do the scaling
      tr = Exact_Kernel::Aff_transformation_2(CGAL::SCALING, t->s) * tr;

      if (t->translate)
        tr = Exact_Kernel::Aff_transformation_2(CGAL::TRANSLATION,
                                                Vector_2(t->c.x(), t->c.y())) * tr;

      std::unique_ptr<CSGCGALDomain2DImpl> transformed = do_transformation(a.impl->polygon_set,
                                                                           tr);
      impl.swap(transformed);
      break;
    }
    case CSGGeometry::Rotation :
    {
      const CSGRotation* t = dynamic_cast<const CSGRotation*>(geometry);
      dolfin_assert(t);
      CSGCGALDomain2D a(t->g.get());
      Exact_Kernel::Aff_transformation_2 tr(CGAL::IDENTITY);

      // Translate if requested
      if (t->translate)
        tr = Exact_Kernel::Aff_transformation_2 (CGAL::TRANSLATION,
                                                 Vector_2(-t->c.x(), -t->c.y())) * tr;

      // Do the rotation
      tr = Exact_Kernel::Aff_transformation_2(CGAL::ROTATION, sin(t->theta), cos(t->theta)) *tr;

      if (t->translate)
        tr = Exact_Kernel::Aff_transformation_2(CGAL::TRANSLATION,
                                                Vector_2(t->c.x(), t->c.y())) * tr;

      std::unique_ptr<CSGCGALDomain2DImpl> transformed = do_transformation(a.impl->polygon_set,
                                                                           tr);
      impl.swap(transformed);
      break;
    }
    case CSGGeometry::Circle:
    {
      const Circle* c = dynamic_cast<const Circle*>(geometry);
      dolfin_assert(c);
      impl->polygon_set.insert(make_circle(c));
      break;
    }
    case CSGGeometry::Ellipse:
    {
      const Ellipse* c = dynamic_cast<const Ellipse*>(geometry);
      dolfin_assert(c);
      impl->polygon_set.insert(make_ellipse(c));
      break;
    }
    case CSGGeometry::Rectangle:
    {
      const Rectangle* r = dynamic_cast<const Rectangle*>(geometry);
      dolfin_assert(r);
      impl->polygon_set.insert(make_rectangle(r));
      break;
    }
    case CSGGeometry::Polygon:
    {
      const Polygon* p = dynamic_cast<const Polygon*>(geometry);
      dolfin_assert(p);
      impl->polygon_set.insert(make_polygon(p));
      break;
    }
    default:
      dolfin::dolfin_error("CSGCGALMeshGenerator2D.cpp",
                           "converting geometry to cgal polyhedron",
                           "Unhandled primitive type");
  }
}
//-----------------------------------------------------------------------------
CSGCGALDomain2D::CSGCGALDomain2D(const CSGCGALDomain2D &other)
 : impl(new CSGCGALDomain2DImpl(other.impl->polygon_set))
{
}
//-----------------------------------------------------------------------------
CSGCGALDomain2D &CSGCGALDomain2D::operator=(const CSGCGALDomain2D &other)
{
  std::unique_ptr<CSGCGALDomain2DImpl> tmp(new CSGCGALDomain2DImpl(other.impl->polygon_set));
  
  impl.swap(tmp);

  return *this;
}
//-----------------------------------------------------------------------------
double CSGCGALDomain2D::compute_boundingcircle_radius() const
{
  std::list<Polygon_with_holes_2> polygon_list;
  impl->polygon_set.polygons_with_holes(std::back_inserter(polygon_list));

  std::vector<Point_2> points;

  for (std::list<Polygon_with_holes_2>::const_iterator pit = polygon_list.begin();
       pit != polygon_list.end(); ++pit)
    for (Polygon_2::Vertex_const_iterator vit = pit->outer_boundary().vertices_begin(); 
         vit != pit->outer_boundary().vertices_end(); ++vit)
      points.push_back(*vit);

  Min_circle min_circle (points.begin(),
                         points.end(),
                         true); //randomize point order

  return sqrt(CGAL::to_double(min_circle.circle().squared_radius()));
}
//-----------------------------------------------------------------------------
void CSGCGALDomain2D::join_inplace(const CSGCGALDomain2D& other)
{
  impl->polygon_set.join(other.impl->polygon_set);
}
//-----------------------------------------------------------------------------
void CSGCGALDomain2D::difference_inplace(const CSGCGALDomain2D& other)
{
  impl->polygon_set.difference(other.impl->polygon_set);
}
//-----------------------------------------------------------------------------
void CSGCGALDomain2D::intersect_inplace(const CSGCGALDomain2D &other)
{
  impl->polygon_set.intersection(other.impl->polygon_set);
}
//-----------------------------------------------------------------------------
bool CSGCGALDomain2D::point_in_domain(dolfin::Point p) const
{
  const Point_2 p_(p.x(), p.y());
  return impl->polygon_set.oriented_side(p_) == CGAL::ON_POSITIVE_SIDE;
}
//-----------------------------------------------------------------------------
std::string CSGCGALDomain2D::str(bool verbose) const
{
  std::stringstream ss;
  ss << "<Polygonal domain with" << std::endl;

  {
    std::list<Polygon_with_holes_2> polygon_list;
    impl->polygon_set.polygons_with_holes(std::back_inserter(polygon_list));
    ss << "  " << polygon_list.size() << " polygons" << std::endl;
  }

  // std::list<std::vector<dolfin::Point> > vertices;
  // get_vertices(vertices, 0);

  // if (verbose)
  // {
  //   for (std::list<std::vector<dolfin::Point> >::const_iterator lit = vertices.begin();
  //        lit != vertices.end(); lit++)
  //   {
  //     ss << "  Polygon" << std::endl;
  //     for (std::vector<dolfin::Point>::const_iterator vit = lit->begin();
  //          vit != lit->end(); vit++)
  //     {
  //       ss << "    " << vit->str(false) << std::endl;
  //     }
  //   }
  // }
  return ss.str();
}
//-----------------------------------------------------------------------------
struct PSLGImpl
{
 public:
  std::vector<Point_2> vertices;
  std::vector<std::pair<std::size_t, std::size_t> > edges;
};
//-----------------------------------------------------------------------------
static inline void add_simple_polygon(std::vector<Segment_2>& segments, 
                                      const Polygon_2& p, 
                                      const Quotient & truncate_tolerance_squared)
{
  Polygon_2::Vertex_const_iterator first = p.vertices_begin(); 
  Polygon_2::Vertex_const_iterator prev = first;
  Polygon_2::Vertex_const_iterator current = prev;
  current++;

  while (current != p.vertices_end())
  {
    Segment_2 s(*prev, *current);

    // Don't add segment if shorter than tolerance
    if (s.squared_length() < truncate_tolerance_squared)
    {
      current++;
      continue;
    }

    segments.push_back(s);
    prev = current;
    current++;

  }
  segments.push_back(Segment_2(*prev, *first));
}
//-----------------------------------------------------------------------------
static inline std::size_t get_vertex_index(std::map<Point_2, std::size_t>& m,
                                           std::vector<dolfin::Point>& v,
                                           Point_2 p)
{
  std::size_t i;
  std::map<Point_2, std::size_t>::const_iterator it = m.find(p);
  if (it == m.end())
  {
    v.push_back(dolfin::Point(CGAL::to_double(p.x()), CGAL::to_double(p.y())));
    m[p] = v.size()-1;
    i = v.size()-1;
  }
  else
  {
    i = it->second;
  }

  return i;
}
//-----------------------------------------------------------------------------
PSLG::PSLG(const std::vector<std::pair<std::size_t, CSGCGALDomain2D> >& domains, 
           double pixel_size, 
           double edge_truncate_tolerance)
{
  const Quotient truncate_tolerance(edge_truncate_tolerance);
  const Quotient truncate_tolerance_squared(truncate_tolerance*truncate_tolerance);


  if (domains.size() > 1 && pixel_size > 0)
  {
    // Collect all segments from all domains to send to snap rounding
    std::vector<Segment_2> segments;

    for (auto it = domains.begin(); it != domains.end(); it++)
    {
      const Polygon_set_2& p = it->second.impl->polygon_set;

      std::list<Polygon_with_holes_2> polygon_list;
      p.polygons_with_holes(std::back_inserter(polygon_list));

      for (std::list<Polygon_with_holes_2>::const_iterator pit = polygon_list.begin();
           pit != polygon_list.end(); ++pit)
      {
        add_simple_polygon(segments, pit->outer_boundary(), truncate_tolerance_squared);

        // Add holes
        Hole_const_iterator hit;
        for (hit = pit->holes_begin(); hit != pit->holes_end(); ++hit)
        {
          add_simple_polygon(segments, *hit, truncate_tolerance_squared);
        }
      }
    }

    log(dolfin::TRACE, "Snap rounding PSLG");
    Polyline_list_2 snapped_polylines;

    CGAL::snap_rounding_2<snap_rounding_traits,
                          std::vector<Segment_2>::const_iterator,
                          Polyline_list_2>
      (segments.begin(), segments.end(),  // input
       snapped_polylines,                 // output
       pixel_size,                        // pixel size
       false,                             // do iterated snap rounding
       false,                             // output as integers
       1);                                // number of kd-trees

    std::map<Point_2, std::size_t> point_to_index;

    for (Polyline_list_2::const_iterator iter1 = snapped_polylines.begin();
         iter1 != snapped_polylines.end(); ++iter1)
    {
      Polyline_2::const_iterator iter2 = iter1->begin();
      const std::size_t first = get_vertex_index(point_to_index,
                                                 vertices,
                                                 *iter2);
      std::size_t prev = first;
      iter2++;

      while(iter2 != iter1->end())
      {
        std::size_t current = get_vertex_index(point_to_index,
                                               vertices,
                                               *iter2);

        edges.push_back(std::make_pair(prev, current));
        prev = current;
        ++iter2;
      }
    }
  }
  else
  {
    log(dolfin::TRACE, "Skipping snap rounding");

    // There is only one domain, so the segment list is just a set of disjoint
    // simple polygons possibly with holes

    const Polygon_set_2& p = domains.front().second.impl->polygon_set;

    std::list<Polygon_with_holes_2> polygon_list;
    p.polygons_with_holes(std::back_inserter(polygon_list));

    for (std::list<Polygon_with_holes_2>::const_iterator pit = polygon_list.begin();
         pit != polygon_list.end(); ++pit)
    {
      {
        // Add the outer boundary
        const Polygon_2& p = pit->outer_boundary();
        Polygon_2::Vertex_const_iterator current = p.vertices_begin();
        Polygon_2::Vertex_const_iterator prev = current;

        const std::size_t vertex_index_start = vertices.size();
        vertices.push_back(dolfin::Point(CGAL::to_double(current->x()),
                                         CGAL::to_double(current->y())));

        std::size_t vertex_index = vertex_index_start;

        current++;

        while (current != p.vertices_end())
        {
          Segment_2 s(*prev, *current);

          // Don't add segment if shorter than tolerance
          if (s.squared_length() < truncate_tolerance_squared)
          {
            current++;
            continue;
          }

          vertices.push_back(dolfin::Point(CGAL::to_double(current->x()),
                                           CGAL::to_double(current->y())));

          edges.push_back(std::make_pair(vertex_index, vertex_index+1));

          prev = current;
          current++;
          vertex_index++;
        }

        edges.push_back(std::make_pair(vertex_index, vertex_index_start));
      }

      // Add holes
      Hole_const_iterator hit;
      for (hit = pit->holes_begin(); hit != pit->holes_end(); ++hit)
      {
        const Polygon_2& p = *hit;
        Polygon_2::Vertex_const_iterator current = p.vertices_begin();
        Polygon_2::Vertex_const_iterator prev = current;

        const std::size_t vertex_index_start = vertices.size();
        vertices.push_back(dolfin::Point(CGAL::to_double(current->x()),
                                         CGAL::to_double(current->y())));
        std::size_t vertex_index = vertex_index_start;

        current++;

        while (current != p.vertices_end())
        {
          Segment_2 s(*prev, *current);

          // Don't add segment if shorter than tolerance
          if (s.squared_length() < truncate_tolerance_squared)
          {
            current++;
            continue;
          }

          vertices.push_back(dolfin::Point(CGAL::to_double(current->x()),
                                           CGAL::to_double(current->y())));

          edges.push_back(std::make_pair(vertex_index, vertex_index+1));

          prev = current;
          current++;
          vertex_index++;
        }

        edges.push_back(std::make_pair(vertex_index, vertex_index_start));
      }
    }
  }
}
//-----------------------------------------------------------------------------
PSLG::~PSLG(){}
//-----------------------------------------------------------------------------
}
