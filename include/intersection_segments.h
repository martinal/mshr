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

/// Sort to a continuous list of segments
template <typename Polyhedron, typename Container>
void sort_segments(const Polyhedron &a, const Polyhedron &b, 
                   const Container &segments)
{
  
}
//-----------------------------------------------------------------------------
// Container is assumed to be a container of boost::tuple<facet_handle,
// facet_handle, Segment>
template<typename Polyhedron, typename Container>
void split_facets(Polyhedron &a, Polyhedron &b,
                  Container &intersections)
{
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Kernel::Triangle_3 Triangle;
  typedef typename Kernel::Segment_3 Segment;
  typedef typename Kernel::Point_3 Point_3;

  typedef boost::tuple<Facet_handle, Facet_handle, Segment> IntersectionType;

  std::cout << "Sorting segments" << std::endl;
  
  std::map<Point_3, typename Container::iterator> map_head;
  std::map<Point_3, typename Container::iterator> map_trail;

  // Store all intersections in map for quick lookup based on their points
  for (typename Container::iterator it = intersections.begin();
       it != intersections.end(); ++it)
  {
    const Segment &s = it->get<2>();
    map_head[s.source()]  = it;
    map_trail[s.target()] = it;
  }

  std::list<std::deque<IntersectionType> > intersection_polylines;

  while(intersections.size() > 0)
  {
    intersection_polylines.push_back(std::deque<IntersectionType>);
    

  }
}

#endif 
