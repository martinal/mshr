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

#include <vector>
#include <iterator>
#include <type_traits>



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

    std::list<Primitive_id> l;

    tree.all_intersected_primitives(t, std::back_inserter(l));
    std::cout << "Number of intersections: " << l.size() << std::endl;

    Facet_const_handle f_biggest = it;

    for (typename std::list<Primitive_id>::iterator it = l.begin();
         it != l.end(); ++it)
    {
      std::cout << "Id: " << std::distance(*it, smallest.facets_begin()) << std::endl;
      Facet_const_handle f_smallest = *it;
      Triangle t(f_smallest->halfedge()->vertex()->point(),
                 f_smallest->halfedge()->next()->vertex()->point(),
                 f_smallest->halfedge()->next()->next()->vertex()->point());


      std::cout << "Triangle: " << t << std::endl;

      out = std::make_pair(f_biggest, f_smallest);

    }
    std::cout << std::endl;
  }
}

template<typename Polyhedron, typename Container>
void split_facets(Polyhedron &a, Polyhedron &b,
                 Container &intersections)
{
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Kernel::Triangle_3 Triangle;
  typedef typename Kernel::Segment_3 Segment;

  // Assume intersections is a container of std::pair<facet_handle, facet_handle>

  for (typename Container::iterator it = intersections.begin(); it != intersections.end(); ++it)
  {
    Triangle t1(it->first->halfedge()->vertex()->point(),
               it->first->halfedge()->next()->vertex()->point(),
               it->first->halfedge()->next()->next()->vertex()->point());

    Triangle t2(it->second->halfedge()->vertex()->point(),
               it->second->halfedge()->next()->vertex()->point(),
               it->second->halfedge()->next()->next()->vertex()->point());


    auto result = intersection(t1, t2);

    if (Segment *s = boost::get<Segment>(&*result))
    {
      std::cout << "Segment: " << *s << std::endl;

    }
    else
    {
      std::cout << "Intersection is not a segment" << std::endl;
    }
  }
}

#endif 
