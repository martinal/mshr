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

#include <vector>

template<typename Polyhedron, typename OutputIterator>
void compute_intersections(Polyhedron &a, const Polyhedron &b, OutputIterator out)
{
  typedef typename Polyhedron::Traits K;
  typedef typename K::Point_3 Point;
  typedef typename K::Plane_3 Plane;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Triangle_3 Triangle;
  typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;
  typedef CGAL::AABB_triangle_primitive<K, typename std::vector<Triangle>::iterator> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;
  typedef typename Tree::Object_and_primitive_id Object_and_primitive_id;
  typedef typename Tree::Primitive_id Primitive_id;

  std::vector<Triangle> triangles;

  const Polyhedron &biggest = a.size_of_facets() > b.size_of_facets() ? a : b;
  const Polyhedron &smallest = a.size_of_facets() > b.size_of_facets() ? b : a;

  //Building the AABB-tree on the smallest polyhedron
  for(Facet_const_iterator it = smallest.facets_begin(); it != smallest.facets_end(); ++it)
    triangles.push_back(Triangle(it->halfedge()->vertex()->point(), 
                                 it->halfedge()->next()->vertex()->point(),
                                 it->halfedge()->next()->next()->vertex()->point()));


  // constructs AABB tree of the smallest polyhedron
  Tree tree(triangles.begin(), triangles.end());

  std::list<Object_and_primitive_id> l;

  for (Facet_const_iterator it = biggest.facets_begin();
       it != biggest.facets_end(); ++it)
  {
    Triangle t(it->halfedge()->vertex()->point(),
               it->halfedge()->next()->vertex()->point(),
               it->halfedge()->next()->next()->vertex()->point());

    tree.all_intersections(t, std::back_inserter(l));
    std::cout << "Number of intersections: " << l.size() << std::endl;
  }

}

#endif 
