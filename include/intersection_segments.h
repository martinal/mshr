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

//-----------------------------------------------------------------------------
/// Some helper functions
template<typename Polyhedron>
inline typename Polyhedron::Halfedge_handle edge_with_point_on(typename Polyhedron::Facet_handle f, typename Polyhedron::Traits::Point_3 p)
{
  // This assumes that the point is actually on an edge.
  typedef typename Polyhedron::Traits::Segment_3 Segment;
  typename Polyhedron::Halfedge_handle h = f->halfedge();

  assert(Segment(h->vertex()->point(), h->opposite()->vertex()->point()).has_on(p) ||
         Segment(h->next()->vertex()->point(), h->next()->opposite()->vertex()->point()).has_on(p) ||
         Segment(h->next()->next()->vertex()->point(), h->next()->next()->opposite()->vertex()->point()).has_on(p));

  typename Polyhedron::Halfedge_handle res = Segment(h->vertex()->point(), h->opposite()->vertex()->point()).has_on(p) ? h :
    (Segment(h->next()->vertex()->point(), h->next()->opposite()->vertex()->point()).has_on(p) ? h->next() :
     h->next()->next());

  return res;
}

template<typename Polyhedron>
inline bool point_on_edge(typename Polyhedron::Facet_handle f, typename Polyhedron::Traits::Point_3 p)
{
  typedef typename Polyhedron::Traits::Segment_3 Segment;
  typename Polyhedron::Halfedge_handle h = f->halfedge();
  return  Segment(h->vertex()->point(), h->next()->vertex()->point()).has_on(p) ||
    Segment(h->next()->vertex()->point(), h->next()->next()->vertex()->point()).has_on(p) ||
    Segment(h->next()->next()->vertex()->point(), h->next()->next()->next()->vertex()->point()).has_on(p);
}
template <typename Polyhedron>
inline void print_facet(typename Polyhedron::Facet_handle f)
{
  std::cout << "Facet(" << &(*f) << ") : " ;
  typename Polyhedron::Halfedge_handle start = f->halfedge();
  typename Polyhedron::Halfedge_handle current = start;
  do 
  {
    std::cout << "(" << current->vertex()->point() << ") "; 
    current = current->next();
  } while (current != start);
  std::cout << std::endl;
}
//-----------------------------------------------------------------------------
/// Compute all intersecting facets of two polyhedron. 
/// For efficiency the biggest polyhedron (in terms of facets) should be given
/// as first argument
template<typename Polyhedron, typename OutputIterator>
void compute_intersections(Polyhedron &biggest, Polyhedron &smallest, OutputIterator out)
{
  typedef typename Polyhedron::Traits K;
  // typedef typename K::Point_3 Point;
  // typedef typename K::Plane_3 Plane;
  // typedef typename K::Vector_3 Vector;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Triangle_3 Triangle;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;
  typedef typename Tree::Object_and_primitive_id Object_and_primitive_id;
  // typedef typename Tree::Primitive_id Primitive_id;

  assert(biggest.is_pure_triangle());
  assert(smallest.is_pure_triangle());

  // constructs AABB tree of the smallest polyhedron
  Tree tree( smallest.facets_begin(), smallest.facets_end(), smallest);

  for (Facet_iterator it = biggest.facets_begin(); it != biggest.facets_end(); ++it)
  {
    Triangle t(it->halfedge()->vertex()->point(),
               it->halfedge()->next()->vertex()->point(),
               it->halfedge()->next()->next()->vertex()->point());

    std::cout << "Testing triangle: " << t << std::endl;

    std::list<Object_and_primitive_id> intersections;

    tree.all_intersections(t, std::back_inserter(intersections));
    std::cout << "Number of intersections: " << intersections.size() << std::endl;

    Facet_handle facet_biggest = it;

    for (typename std::list<Object_and_primitive_id>::iterator it = intersections.begin();
         it != intersections.end(); ++it)
    {
      std::cout << "Id: " << std::distance(it->second, smallest.facets_begin()) << std::endl;
      Facet_handle facet_smallest = it->second;
      Triangle t(facet_smallest->halfedge()->vertex()->point(),
                 facet_smallest->halfedge()->next()->vertex()->point(),
                 facet_smallest->halfedge()->next()->next()->vertex()->point());

      Segment s;
      if (CGAL::assign(s, it->first))
      {
        std::cout << "Is segment" << std::endl;
        out = boost::make_tuple(facet_biggest, facet_smallest, s);
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
template<typename Polyhedron>
void split_facets(Polyhedron &a, Polyhedron &b,
                  std::list<boost::tuple<typename Polyhedron::Facet_handle, typename Polyhedron::Facet_handle, typename Polyhedron::Traits::Segment_3> > &intersections)
{
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  // typedef typename Polyhedron::Vertex_handle Vertex_handle;
  // typedef typename Kernel::Triangle_3 Triangle;
  typedef typename Kernel::Segment_3 Segment;
  // typedef typename Kernel::Point_3 Point_3;

  typedef boost::tuple<Facet_handle, Facet_handle, Segment> IntersectionType;

  typedef std::list<IntersectionType> IntersectionList;


  std::cout << "Sorting segments" << std::endl;

  for (typename IntersectionList::iterator it = intersections.begin();
       it != intersections.end(); ++it)
    std::cout << it->get<2>() << std::endl;

  std::cout << "Total: " << intersections.size() << " segments" << std::endl;
  
  std::list<std::list<IntersectionType> > polylines;

  while(intersections.size())
  {
    bool found = false;
    IntersectionType e = intersections.front();
    intersections.pop_front();

    std::cout << "Looking for : " << e.get<2>() << std::endl;

    for (typename std::list<IntersectionList>::iterator p = polylines.begin();
         p != polylines.end(); ++p)
    {

      found = true;
      if (e.get<2>().target() == p->front().get<2>().source() ||
          e.get<2>().target() == p->front().get<2>().target() ||
          e.get<2>().source() == p->front().get<2>().source() ||
          e.get<2>().source() == p->front().get<2>().target())
        p->push_front(e);
      else if (e.get<2>().target() == p->back().get<2>().source() ||
               e.get<2>().target() == p->back().get<2>().target() ||
               e.get<2>().source() == p->back().get<2>().source() ||
               e.get<2>().source() == p->back().get<2>().target())
        p->push_back(e);
      else
        found = false;

      if (found)
      {
        std::cout << "  Found: " << p->back().get<2>() << std::endl;

        // Now check if we can connect some of the polylines
        for (typename std::list<IntersectionList>::iterator p2 = polylines.begin();
             p2 != polylines.end(); ++p2)
        {
          if (p2 == p)
            continue;

          if (p2->back().get<2>().source() == p->front().get<2>().source() ||
              p2->back().get<2>().target() == p->front().get<2>().source() ||
              p2->back().get<2>().source() == p->front().get<2>().target() ||
              p2->back().get<2>().target() == p->front().get<2>().target())
          {
            p->splice(p->begin(), *p2);
            polylines.erase(p2);
            break;
          } 
          else if (p2->front().get<2>().source() == p->front().get<2>().source() ||
                   p2->front().get<2>().target() == p->front().get<2>().source() ||
                   p2->front().get<2>().source() == p->front().get<2>().target() ||
                   p2->front().get<2>().target() == p->front().get<2>().target())
          {
            p2->reverse();
            p->splice(p->begin(), *p2);
            polylines.erase(p2);
            break;
          } 
          else if (p2->back().get<2>().source() == p->back().get<2>().source() ||
                   p2->back().get<2>().target() == p->back().get<2>().source() ||
                   p2->back().get<2>().source() == p->back().get<2>().target() ||
                   p2->back().get<2>().target() == p->back().get<2>().target())
          {
            p2->reverse();
            p->splice(p->end(), *p2);
            polylines.erase(p2);
            break;
          } 
          else if (p2->front().get<2>().source() == p->back().get<2>().source() ||
                   p2->front().get<2>().target() == p->back().get<2>().source() ||
                   p2->front().get<2>().source() == p->back().get<2>().target() ||
                   p2->front().get<2>().target() == p->back().get<2>().target())
          {
            p2->reverse();
            p->splice(p->end(), *p2);
            polylines.erase(p2);
            break;
          }
        }
        break;
      }
    }
    if (!found)
    {
      polylines.push_back(IntersectionList());
      polylines.back().push_back(e);
    }
  }

  // Swap the segments if needed
  for (typename std::list<std::list<IntersectionType> >::iterator it = polylines.begin();
       it != polylines.end(); ++it)
  {
    std::cout << "---Flipping segments" << std::endl;

    typename std::list<IntersectionType>::iterator it2 = it->begin();
    Segment prev_segment = it2->get<2>();

    {
      // Check if we need to flip the first segment
      typename std::list<IntersectionType>::iterator second = it2;
      second++;
      if (prev_segment.target() != second->get<2>().source() &&
          prev_segment.target() != second->get<2>().target())
      {
        it2->get<2>() = prev_segment.opposite();
      }
    }

    it2++;
    for (; it2 != it->end(); it2++)
    {
      Segment &current_segment = it2->get<2>();

      std::cout << "Prev: " << prev_segment << std::endl;
      std::cout << "Current: " << current_segment << std::endl;

      if (current_segment.source() != prev_segment.target())
      {
        current_segment = current_segment.opposite();
        std::cout << "  Flipping" << std::endl;
      }
      std::cout <<  std::endl;

      prev_segment = current_segment;
    }
  }  

  std::cout << "Done sorting" << std::endl;
  std::cout << "Found " << polylines.size() << " polylines" << std::endl;

  // Split facets so all intersections exists as edges
  for (typename std::list<std::list<IntersectionType> >::iterator it = polylines.begin();
       it != polylines.end(); ++it)
  {
    std::cout << "polyline" << std::endl;

    // Find element that intersects an edge if possible
    bool polyline_intersects_edges = false;
    typename std::list<IntersectionType>::iterator start = it->begin();

    {
      typename std::list<IntersectionType>::iterator prev = it->end();
      prev--;

      for (; start != it->end(); start++)
      {
        if (start->get<0>() != prev->get<0>())
        {
          polyline_intersects_edges = true;
          break;
        }
        prev = start;
      }
    }
    
    if (polyline_intersects_edges)
    {
      // Bring the segment intersecting an edge to the front of the list
      if (start != it->begin())
        it->splice(it->begin(), *it, start, it->end());

      assert(point_on_edge<Polyhedron>(it->begin()->get<0>(), it->begin()->get<2>().source()));
      assert(point_on_edge<Polyhedron>(it->back().get<0>(), it->back().get<2>().target()));

      std::cout << "Polyline: " << std::endl;
      for (typename std::list<IntersectionType>::iterator it2 = it->begin(); it2 != it->end(); it2++)
      {
        std::cout << "(" << &(*it2->get<0>()) << ")  " << it2->get<2>() << std::endl;
      }
      std::cout << std::endl;

      // Add the new vertex
      std::cout << "Add new start vertex" << std::endl;
      print_facet<Polyhedron>(it->begin()->get<0>());
      Halfedge_handle new_start_edge = a.split_edge(edge_with_point_on<Polyhedron>(it->begin()->get<0>(), it->begin()->get<2>().source()));
      new_start_edge->vertex()->point() = it->begin()->get<2>().source();

      // Save the edge where the polyline ends
      Halfedge_handle polyline_end_vertex = new_start_edge->opposite()->prev();

      for (typename std::list<IntersectionType>::iterator it2 = it->begin(); it2 != it->end(); it2++)
      {
        std::cout << "Adding to facet" << std::endl;
        // Do all job on a particular facet
        typename std::list<IntersectionType>::iterator start = it2;
        typename std::list<IntersectionType>::iterator end = start;
        
        assert(start->get<0>() == new_start_edge->facet());

        // Advance in list to find last segment on this facet
        {
          typename std::list<IntersectionType>::iterator ahead = start;
          ahead++;
          while (ahead != it->end() && ahead->get<0>() == start->get<0>())
          {
            std::cout << "  ---Advancing" << std::endl;
            end++;
            ahead++;
          }
        }
        assert(end->get<0>() == new_start_edge->facet());

        std::cout << "Getting/creating last vertex" << std::endl;

        // Add the last vertex
        // TODO: Handle the case where startpoint and endpoint coinside
        Halfedge_handle new_end_edge;
        if (polyline_end_vertex->facet() == end->get<0>())
        {
          new_end_edge = polyline_end_vertex;
          assert(end->get<2>().target() == polyline_end_vertex->vertex()->point());
        }
        else
        {
          std::cout << "Creating last vertex on facet" << std::endl;
          print_facet<Polyhedron>(end->get<0>());
          std::cout << "Point: " << end->get<2>().target() << std::endl;
          Halfedge_handle end_edge = edge_with_point_on<Polyhedron>(end->get<0>(), end->get<2>().target());
          new_end_edge = a.split_edge(end_edge);
          new_end_edge->vertex()->point() = end->get<2>().target();
        }

        Halfedge_handle inserted_edge = a.split_facet(new_start_edge, new_end_edge);
        
        // And finally: Insert the interior points from the polyline on the new edge
        for (typename std::list<IntersectionType>::iterator current = start; current != end; current++)
        {
          Halfedge_handle i = a.split_edge(inserted_edge);
          i->vertex()->point() = current->get<2>().target();
        }

        // save the start edge for the next facet
        new_start_edge = new_end_edge;

      }
    } 
    else
    {
      // The entire intersection polyline is within one facet
      std::cout << "Not implemented yet" << std::endl;

    }
  }
}
//-----------------------------------------------------------------------------
#endif 
