// compute self-intersection of a CGAL triangle polyhedron mesh
// original code from Lutz Kettner
#ifndef intersection_segments_h
#define intersection_segments_h

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/intersections.h>

#include <boost/variant/get.hpp>
#include <boost/tuple/tuple.hpp>

#include <vector>
#include <list>
#include <deque>
#include <iterator>
#include <type_traits>

/// Compute all intersecting facets of two polyhedron. 
/// For efficiency the biggest polyhedron (in terms of facets) should be given
/// as first argument
template<typename Polyhedron, typename OutputIterator>
void compute_intersections(const Polyhedron &biggest, const Polyhedron &smallest, OutputIterator out)
{
  typedef typename Polyhedron::Traits K;
  typedef typename K::Point_3 Point;
  typedef typename K::Plane_3 Plane;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Triangle_3 Triangle;
  typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;
  typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
  typedef CGAL::AABB_face_graph_triangle_primitive<const Polyhedron> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;
  typedef typename Tree::Object_and_primitive_id Object_and_primitive_id;
  typedef typename Tree::Primitive_id Primitive_id;

  // constructs AABB tree of the smallest polyhedron
  Tree tree( smallest.facets_begin(), smallest.facets_end(), smallest);

  for (Facet_const_iterator it = biggest.facets_begin();
       it != biggest.facets_end(); ++it)
  {
    Triangle t(it->halfedge()->vertex()->point(),
               it->halfedge()->next()->vertex()->point(),
               it->halfedge()->next()->next()->vertex()->point());

    std::cout << "Testing triangle: " << t << std::endl;

    std::list<Object_and_primitive_id> intersections;

    tree.all_intersections(t, std::back_inserter(intersections));
    std::cout << "Number of intersections: " << intersections.size() << std::endl;

    Facet_const_handle f_biggest = it;

    for (typename std::list<Object_and_primitive_id>::iterator it = intersections.begin();
         it != intersections.end(); ++it)
    {
      std::cout << "Id: " << std::distance(it->second, smallest.facets_begin()) << std::endl;
      Facet_const_handle f_smallest = it->second;
      Triangle t(f_smallest->halfedge()->vertex()->point(),
                 f_smallest->halfedge()->next()->vertex()->point(),
                 f_smallest->halfedge()->next()->next()->vertex()->point());

      Segment s;
      if (CGAL::assign(s, it->first))
      {
        std::cout << "Is segment" << std::endl;
        out = boost::make_tuple(f_biggest, f_smallest, s);
      } 
      else
      {
        std::cout << "Not segment" << std::endl;
      }
    }
    std::cout << std::endl;
  }
}
//-----------------------------------------------------------------------------
// Container is assumed to be a container of boost::tuple<facet_handle,
// facet_handle, Segment>
// The container must not invalidate iterators when elements are removed
template<typename Polyhedron> //, typename Container>
void split_facets(Polyhedron &a, Polyhedron &b,
                  std::list<boost::tuple<typename Polyhedron::Facet_const_handle, typename Polyhedron::Facet_const_handle, typename Polyhedron::Traits::Segment_3> > &intersections)
{
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
  typedef typename Kernel::Triangle_3 Triangle;
  typedef typename Kernel::Segment_3 Segment;
  typedef typename Kernel::Point_3 Point_3;
  typedef std::list<boost::tuple<typename Polyhedron::Facet_const_handle, typename Polyhedron::Facet_const_handle, typename Polyhedron::Traits::Segment_3> > IntersectionList;

  typedef boost::tuple<Facet_const_handle, Facet_const_handle, Segment> IntersectionType;

  std::cout << "Sorting segments" << std::endl;
  
  std::map<Point_3, typename IntersectionList::iterator> map_head;
  std::map<Point_3, typename IntersectionList::iterator> map_trail;

  // Store all intersections in map for quick lookup based on their points
  for (typename IntersectionList::iterator it = intersections.begin();
       it != intersections.end(); ++it)
  {
    const Segment &s = it->get<2>();
    map_head[s.source()]  = it;
    map_trail[s.target()] = it;
  }

  std::list<std::deque<IntersectionType> > intersection_polylines;

  while(intersections.size() > 0)
  {
    // A new deque containing the polyline
    intersection_polylines.push_back(std::deque<IntersectionType>());
    std::deque<IntersectionType> &polyline = intersection_polylines.back();

    // Insert one segment into polyline
    IntersectionType e = intersections.front();
    intersections.pop_front();
    map_head.erase(e.get<2>().source());
    map_trail.erase(e.get<2>().target());
      
    polyline.push_back(e);

    while (map_head.count(polyline.back().get<2>().target()))
    {
      typename IntersectionList::iterator e = map_head[polyline.back().get<2>().target()];
      polyline.push_front(*e);
      map_head.erase(e->get<2>().target());
      map_trail.erase(e->get<2>().source());        
    }
    
    while(map_trail.count(polyline.front().get<2>().source()))
    {
      typename IntersectionList::iterator e = map_trail[polyline.front().get<2>().source()];
      polyline.push_back(*e);
      map_trail.erase(e->get<2>().source());
      map_head.erase(e->get<2>().target());
    }
  }
}

#endif 
