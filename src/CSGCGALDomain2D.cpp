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

// Polygon typedefs
typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
typedef Exact_Kernel::Point_2                             Point_2;
typedef Exact_Kernel::Vector_2                            Vector_2;
typedef Exact_Kernel::Segment_2                           Segment_2;
typedef CGAL::Polygon_2<Exact_Kernel>                     Polygon_2;
typedef Polygon_2::Vertex_const_iterator                  Vertex_const_iterator;
typedef CGAL::Polygon_with_holes_2<Exact_Kernel>          Polygon_with_holes_2;
typedef Polygon_with_holes_2::Hole_const_iterator         Hole_const_iterator;
typedef CGAL::Polygon_set_2<Exact_Kernel>                 Polygon_set_2;

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
  pts.reserve(c->fragments());

  for (std::size_t i = 0; i < c->fragments(); i++)
  {
    const double phi = (2*DOLFIN_PI*i) / c->fragments();
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

  for (std::size_t i = 0; i < e->fragments(); i++)
  {
    const double phi = (2*DOLFIN_PI*i) / e->fragments();
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
// TODO: Consider if the 2D snap rounding package in CGAL can be used for this
class PSLGImpl
{
 public:
  std::vector<Point_2> vertices;
  std::vector<std::pair<std::size_t, std::size_t> > edges;
  //-----------------------------------------------------------------------------
  // Add point if not present. Return index of point
  std::size_t add_point(Point_2 p)
  {
    for (std::size_t i = 0; i < vertices.size(); i++)
    {
      if (vertices[i] == p)
        return i;
    }

    vertices.push_back(p);
    return vertices.size()-1;
  }
  //-----------------------------------------------------------------------------
  // Adds a new line segment from two existing vertices, possibly introduing new
  // vertices if the new line segment intersects existing ones
  // Recursive function that
  void add_linesegment(std::size_t p1, std::size_t p2)
  {
    Segment_2 s(vertices[p1], vertices[p2]);
    std::cout << "Adding segment: " << s << std::endl;
    for (auto it = edges.begin(); it != edges.end(); it++)
    {
      Segment_2 current(vertices[it->first], vertices[it->second]);
      std::cout << "  Intersecting with: " << current << std::endl;
      if (CGAL::do_intersect(s, current))
      {
        std::cout << "Intersects!";
        auto intersection = CGAL::intersection(s, current);
        Point_2 vertex;
        Segment_2 segment;
        if (CGAL::assign(vertex, intersection))
        {
          std::cout << "New vertex: " << vertex << std::endl;
        }
        else if(CGAL::assign(segment, intersection))
        {
          std::cout << "Intersection segment: " << segment << std::endl;
        }
        else
        {
          dolfin::dolfin_error("CSGCGALDomain2.cpp",
                               "add subdomain to PSLG",
                               "intersection if not point or segment");
        }
      }
    }

    // Did not intersect any line segments, so just insert it
    // (otherwise the function would have returned earlier
    edges.push_back(std::make_pair(p1, p2));
  }
  //-----------------------------------------------------------------------------
  void add_simple_polygon(const Polygon_2 & p)
  {
    Polygon_2::Vertex_const_iterator prev = p.vertices_begin(); 
    std::cout << "First vertex:" << *prev << std::endl;
    const std::size_t first_id = add_point(*prev);

    Polygon_2::Vertex_const_iterator current = prev;
    current++;

    // Assume that polygon has more than one vertex 
    // (are polygon with exactly one vertex allowed?)
    std::size_t prev_id = first_id;
    do
    {
      std::size_t current_id = add_point(*current);
      std::cout << "Pushing: " << *current << std::endl;
      add_linesegment(prev_id, current_id);
      
      prev_id = current_id;
      current++;
    } while (current != p.vertices_end());

    add_linesegment(prev_id, first_id);
  }
};
//-----------------------------------------------------------------------------
PSLG::PSLG(const CSGCGALDomain2D& domain)
  : impl(new PSLGImpl)
{

  // TODO: This can be optimized since we don't need to check for
  // existing vertices and intersections.
  add_subdomin(domain);

  // TODO: This can be optimized to avoid some copying
  std::list<Polygon_with_holes_2> polygon_list;
  domain.impl->polygon_set.polygons_with_holes(std::back_inserter(polygon_list));

  // Add the vertices from the outer boundary
  std::list<Polygon_with_holes_2>::const_iterator pit;
  for (pit = polygon_list.begin(); pit != polygon_list.end(); ++pit)
  {
    const Polygon_2 &outer = pit->outer_boundary();
    impl->add_simple_polygon(outer);

    // Add holes
    Hole_const_iterator hit;
    for (hit = pit->holes_begin(); hit != pit->holes_end(); ++hit)
    {
      impl->add_simple_polygon(*hit);
    }
  }
}
//-----------------------------------------------------------------------------
PSLG::~PSLG(){}
//-----------------------------------------------------------------------------
void PSLG::add_subdomain(const CSGCGALDomain2D& subdomain)
{
  // TODO: This can be optimized to avoid some copying
  std::list<Polygon_with_holes_2> polygon_list;
  subdomain.impl->polygon_set.polygons_with_holes(std::back_inserter(polygon_list));

  // Add the vertices from the outer boundary
  std::list<Polygon_with_holes_2>::const_iterator pit;
  for (pit = polygon_list.begin(); pit != polygon_list.end(); ++pit)
  {
    log(dolfin::TRACE, "Add polygon");
    const Polygon_2 &outer = pit->outer_boundary();

    impl->add_simple_polygon(outer);

    // Add holes
    for (Hole_const_iterator hit = pit->holes_begin(); hit != pit->holes_end(); ++hit)
    {
      impl->add_simple_polygon(*hit);
    }
  }
}
//-----------------------------------------------------------------------------
void PSLG::collapse_short_edges(double tolerance)
{

}
//-----------------------------------------------------------------------------
void PSLG::get_vertices(std::vector<dolfin::Point>& v) const
{
  v.clear();
  for (auto it = impl->vertices.begin(); it != impl->vertices.end(); it++)
  {
    v.push_back(dolfin::Point(CGAL::to_double(it->x()), CGAL::to_double(it->y())));
    dolfin::cout << v.back() << dolfin::endl;
  }
}
//-----------------------------------------------------------------------------
void PSLG::get_edges(std::vector<std::pair<std::size_t, std::size_t> >& e) const
{
  std::copy(impl->edges.begin(), impl->edges.end(), std::back_inserter(e));
  std::cout << "Copied edges" << impl->edges.size() << " " << e.size() << std::endl;
}
}
