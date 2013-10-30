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
  // Add an assert?
  typedef typename Polyhedron::Traits::Segment_3 Segment;
  typename Polyhedron::Halfedge_handle h = f->halfedge();
  return  Segment(h->vertex()->point(), h->next()->vertex()->point()).has_on(p) ? h :
    (Segment(h->next()->vertex()->point(), h->next()->next()->vertex()->point()).has_on(p) ? h->next() :
     h->next()->next());
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
//-----------------------------------------------------------------------------
/// Compute all intersecting facets of two polyhedron. 
/// For efficiency the biggest polyhedron (in terms of facets) should be given
/// as first argument
template<typename Polyhedron, typename OutputIterator>
void compute_intersections(Polyhedron &biggest, Polyhedron &smallest, OutputIterator out)
{
  typedef typename Polyhedron::Traits K;
  typedef typename K::Point_3 Point;
  typedef typename K::Plane_3 Plane;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Triangle_3 Triangle;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;
  typedef typename Tree::Object_and_primitive_id Object_and_primitive_id;
  typedef typename Tree::Primitive_id Primitive_id;

  // constructs AABB tree of the smallest polyhedron
  Tree tree( smallest.facets_begin(), smallest.facets_end(), smallest);

  for (Facet_iterator it = biggest.facets_begin();
       it != biggest.facets_end(); ++it)
  {
    Triangle t(it->halfedge()->vertex()->point(),
               it->halfedge()->next()->vertex()->point(),
               it->halfedge()->next()->next()->vertex()->point());

    std::cout << "Testing triangle: " << t << std::endl;

    std::list<Object_and_primitive_id> intersections;

    tree.all_intersections(t, std::back_inserter(intersections));
    std::cout << "Number of intersections: " << intersections.size() << std::endl;

    Facet_handle f_biggest = it;

    for (typename std::list<Object_and_primitive_id>::iterator it = intersections.begin();
         it != intersections.end(); ++it)
    {
      std::cout << "Id: " << std::distance(it->second, smallest.facets_begin()) << std::endl;
      Facet_handle f_smallest = it->second;
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
template<typename Polyhedron>
void split_facets(Polyhedron &a, Polyhedron &b,
                  std::list<boost::tuple<typename Polyhedron::Facet_handle, typename Polyhedron::Facet_handle, typename Polyhedron::Traits::Segment_3> > &intersections)
{
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Kernel::Triangle_3 Triangle;
  typedef typename Kernel::Segment_3 Segment;
  typedef typename Kernel::Point_3 Point_3;

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

  for (typename std::list<std::list<IntersectionType> >::iterator it = polylines.begin();
       it != polylines.end(); ++it)
  {
    for (typename std::list<IntersectionType>::iterator it2 = it->begin(); it2 != it->end(); it2++)
    {
      std::cout << "  " << it2->get<2>() << std::endl;
    }
    std::cout << std::endl;
  }

  // Split facets so all intersections are exists as edges
  for (typename std::list<std::list<IntersectionType> >::iterator it = polylines.begin();
       it != polylines.end(); ++it)
  {
    std::cout << "polyline" << std::endl;

    // Rotate list to ensure that the first element intersects an edge
    // if possible
    // TODO: This is untested
    typename std::list<IntersectionType>::iterator it2 = it->begin();
    bool polyline_intersects_edges = false;
    for (;it2 != it->end(); it2++)
    {
      std::cout << "  Looking" << std::endl;
      Facet_handle f = it2->get<1>();
      assert(f->is_triangle());
      Halfedge_handle h = f->halfedge();
      Segment &s = it2->get<2>();

      // TODO: Avoid expensive geometric test here by comparing facet_handles
      if (Segment(h->vertex()->point(), h->next()->vertex()->point()).has_on(s.source()) ||
          Segment(h->next()->vertex()->point(), h->next()->next()->vertex()->point()).has_on(s.source()) ||
          Segment(h->next()->next()->vertex()->point(), h->vertex()->point()).has_on(s.source()))
      {
        std::cout << "    Breaking" << std::endl;
        polyline_intersects_edges = true;
        break;
      }
    }
    
    if (polyline_intersects_edges)
    {
      // Bring the segment intersecting an edge to the front of the list
      if (it2 != it->begin())
        it->splice(it->begin(), *it, it2, it->end());



      typename std::list<IntersectionType>::iterator it2 = it->begin();
      while (it2 != it->end())
      {
        
        typename std::list<IntersectionType>::iterator start = it2;
        // Add the new vertex
        Halfedge_handle start_edge = edge_with_point_on<Polyhedron>(start->get<0>(), start->get<2>().source());
        Halfedge_handle new_start_edge = a.split_edge(start_edge);
        new_start_edge->vertex()->point() = start->get<2>().source();
        
        // Advance in list to find last segment on this facet
        // TODO: This is cumbersome. Put in a separate function?
        typename std::list<IntersectionType>::iterator end = start;
        {
          typename std::list<IntersectionType>::iterator ahead = start;
          ahead++;
          while (ahead != it->end() || ahead->get<0>() == start->get<0>())
          {
            end++;
            ahead++;
          }
        }

        // Add the last vertex
        // TODO: Handle the case where startpoint and endpoint coinside
        Halfedge_handle end_edge = edge_with_point_on<Polyhedron>(end->get<0>(), end->get<2>().target());
        Halfedge_handle new_end_edge = a.split_edge(end_edge);
        new_end_edge->vertex()->point() = end->get<2>().target();

        Halfedge_handle new_edge = a.split_facet(new_start_edge, new_end_edge);

        // And finally: Insert the interior points from the polyline on the new edge
        for (typename std::list<IntersectionType>::iterator current = start; current != end; current++)
        {
          Halfedge_handle inserted = a.split_edge(new_edge);
          inserted->vertex()->point() = current->get<2>().target();
        }
      }      
    } 
    else
    {
      // The entire intersection polyline is within one facet
      

    }
  }
}
//-----------------------------------------------------------------------------
#endif 
