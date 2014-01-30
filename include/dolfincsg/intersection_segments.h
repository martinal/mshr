// compute self-intersection of a CGAL triangle polyhedron mesh
// original code from Lutz Kettner
#ifndef dolfincsg_intersection_segments_h
#define dolfincsg_intersection_segments_h

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
#include <set>
#include <fstream>
#include <functional>

#define EXPENSIVE_CHECKS true

//-----------------------------------------------------------------------------
/// Some helper functions
template<typename Polyhedron>
inline bool point_on_edge(typename Polyhedron::Facet_handle f, typename Polyhedron::Traits::Point_3 p)
{
  typedef typename Polyhedron::Traits::Segment_3 Segment;
  typename Polyhedron::Halfedge_handle start = f->halfedge();
  typename Polyhedron::Halfedge_handle h = start;

  do
  {
    if (Segment(h->vertex()->point(), h->opposite()->vertex()->point()).has_on(p))
      return true;
    
    h = h->next();
  } while (h != start);
  return false;
}

template<typename Polyhedron>
inline typename Polyhedron::Halfedge_handle edge_with_point_on(typename Polyhedron::Facet_handle f, typename Polyhedron::Traits::Point_3 p)
{
  // Assume that the point is actually on an edge.
  typedef typename Polyhedron::Traits::Segment_3 Segment;

  typedef typename Polyhedron::Traits::Segment_3 Segment;
  typename Polyhedron::Halfedge_handle start = f->halfedge();
  typename Polyhedron::Halfedge_handle h = start;

  do
  {
    if (Segment(h->vertex()->point(), h->opposite()->vertex()->point()).has_on(p))
    {
      return h;
    }
    
    h = h->next();
  } while (h != start);

  assert(false);

  // Should never happen
  return h;
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

template <typename Polyhedron>
typename Polyhedron::Traits::Triangle_3 get_facet_triangle(typename Polyhedron::Facet_handle f)
{
  typename Polyhedron::Halfedge_handle h = f->halfedge();
  return typename Polyhedron::Traits::Triangle_3(f->vertex()->point(),
                                                 f->next()->vertex()->point(),
                                                 f->next()->next()->vertex()->point());
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
  // typedef typename K::Plane_3 Plane;
  // typedef typename K::Vector_3 Vector;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Triangle_3 Triangle;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
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
    Halfedge_handle f = it->halfedge();
    Triangle t(f->vertex()->point(),
               f->next()->vertex()->point(),
               f->next()->next()->vertex()->point());

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
      print_facet<Polyhedron>(facet_smallest);
      Segment s;
      Triangle t;
      Point p;
      if (CGAL::assign(s, it->first))
      {
        std::cout << "Is segment" << std::endl;
        out = boost::make_tuple(facet_biggest, facet_smallest, s);
      }
      else if (CGAL::assign(t, it->first))
      {
        std::cout << "Is Triangle" << std::endl;
        std::cout << t << std::endl;
      } 
      else if (CGAL::assign(p, it->first))
      {
        std::cout << "Is point" << std::endl;
        std::cout << p << std::endl;
      }
      else
        std::cout << "Not known primitive" << std::endl;
    }
    std::cout << std::endl;
  }
}
//-----------------------------------------------------------------------------
template<typename Polyhedron>
void sort_polylines(Polyhedron &a, Polyhedron &b,
                    std::list<boost::tuple<typename Polyhedron::Facet_handle, typename Polyhedron::Facet_handle, typename Polyhedron::Traits::Segment_3> > &intersections,
                    std::list<std::list<boost::tuple<typename Polyhedron::Facet_handle, typename Polyhedron::Facet_handle, typename Polyhedron::Traits::Segment_3> > > &polylines)
{
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  //typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Kernel::Segment_3 Segment;
  typedef boost::tuple<Facet_handle, Facet_handle, Segment> IntersectionType;
  typedef std::list<IntersectionType> IntersectionList;


  std::cout << "Sorting segments" << std::endl;

  for (typename IntersectionList::iterator it = intersections.begin();
       it != intersections.end(); ++it)
    std::cout << it->get<2>() << std::endl;

  std::cout << "Total: " << intersections.size() << " segments" << std::endl;
  
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
      {
        std::cout << "  Found (front): " << p->front().get<2>() << std::endl;
        p->push_front(e);
      }
      else if (e.get<2>().target() == p->back().get<2>().source() ||
               e.get<2>().target() == p->back().get<2>().target() ||
               e.get<2>().source() == p->back().get<2>().source() ||
               e.get<2>().source() == p->back().get<2>().target())
      {
        std::cout << "  Found (back): " << p->back().get<2>() << std::endl;
        p->push_back(e);
      }
      else
        found = false;

      if (found)
      {
        std::cout << "Try to connect polylines" << std::endl;
        //bool connected = false;

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
}
//-----------------------------------------------------------------------------
// Split facets so all intersection segments are represented as edges. Return
// TODO: Consider using std::vector instead. Might be easier to work with
// indices that list iterators
template<typename Polyhedron, int N>
void split_facets(Polyhedron &a,
                  const std::list<std::list<boost::tuple<typename Polyhedron::Facet_handle, typename Polyhedron::Facet_handle, typename Polyhedron::Traits::Segment_3> > >  &polylines,
                  std::list<std::vector<typename Polyhedron::Halfedge_handle> > &intersection_list)
{
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Kernel::Segment_3 Segment;
  typedef typename Kernel::Point_3 Point_3;
  typedef boost::tuple<Facet_handle, Facet_handle, Segment> IntersectionType;
  typedef std::list<IntersectionType> IntersectionList;

  std::vector<typename std::list<IntersectionList>::const_iterator> postponed;

  // the list [1] is the polyline, the iterator [2] is the starting point, corresponding to the halfedge [0]
  std::vector<boost::tuple<Halfedge_handle, std::reference_wrapper<const IntersectionList>, typename IntersectionList::const_iterator> > polylines_with_startvertex;

  // Split edges to create a starting point for each polyline
  for (typename std::list<IntersectionList>::const_iterator it = polylines.begin();
       it != polylines.end(); ++it)
  {
    // Find element that intersects an edge if possible
    bool polyline_intersects_edges = false;
    typename IntersectionList::const_iterator start = it->begin();

    {
      typename IntersectionList::const_iterator prev = it->end();
      prev--;

      for (; start != it->end(); start++)
      {
        if (start->get<N>() != prev->get<N>())
        {
          polyline_intersects_edges = true;
          break;
        }
        prev = start;
      }
    }
    
    if (!polyline_intersects_edges)
    {
      // Need to create new edges, postpone
      postponed.push_back(it);
    }
    else
    {
      std::cout << "Intersects edges" << std::endl;

      // Removed spline here!

      std::cout << "polyline" << std::endl;
      for (typename std::list<IntersectionType>::const_iterator it2 = it->begin(); it2 != it->end(); it2++)
      {
        std::cout << "(" << &(*it2->get<N>()) << ")  " << it2->get<2>() << std::endl;
      }
      std::cout << std::endl;
      
      assert(point_on_edge<Polyhedron>(start->get<N>(), start->get<2>().source()));

      // Add the new vertex
      // TODO: Check that this point isn't already a vertex

      std::cout << "Add new start vertex (" << it->begin()->get<2>().source() << ")" << std::endl;
      print_facet<Polyhedron>(it->begin()->get<N>());

      const Halfedge_handle edge_to_be_splitted = edge_with_point_on<Polyhedron>(it->begin()->get<N>(),
                                                                                 it->begin()->get<2>().source());
      std::cout << "Edge: " << edge_to_be_splitted->vertex()->point() << ", " << edge_to_be_splitted->opposite()->vertex()->point() << std::endl;
      Halfedge_handle start_vertex = a.split_edge(edge_to_be_splitted);
      start_vertex->vertex()->point() = it->begin()->get<2>().source();
      std::cout << "Startvertex: " << start_vertex->vertex()->point() << ", prev: " << start_vertex->opposite()->vertex()->point() << 
        ", next: " << start_vertex->next()->vertex()->point() << std::endl;
      assert(start_vertex->facet() == it->begin()->get<N>());

      //boost::tuple<Halfedge_handle, std::reference_wrapper<IntersectionList>, typename IntersectionList::const_iterator>
      std::reference_wrapper<const IntersectionList> tmp = std::cref(*it);
      polylines_with_startvertex.push_back(boost::make_tuple(start_vertex, tmp /*std::ref(*it)*/, start));
    }
  }

  // TODO: Handle the postponed polylines (those not intersecting edges)
  assert(postponed.size() == 0);

  std::cout << " --- Splitting facets" << std::endl;

  std::set<Halfedge_handle> intersection_facets;

  // Now split the facets so the intersection polyline is represented as edges
  // in the polyhedron
  for (typename std::vector<boost::tuple<Halfedge_handle, 
                                         std::reference_wrapper<const IntersectionList>, 
                                         typename IntersectionList::const_iterator> >::iterator polyline_it = polylines_with_startvertex.begin();
       polyline_it != polylines_with_startvertex.end(); ++polyline_it)
  {
    std::cout << "Starting polyline" << std::endl;

    // Unpack tuple for convenience
    const Halfedge_handle polyline_start_vertex = polyline_it->get<0>();
    const IntersectionList &the_list = polyline_it->get<1>();
    typename IntersectionList::const_iterator circulator_begin = polyline_it->get<2>();

    // Add the new vector to return list
    intersection_list.push_back(std::vector<Halfedge_handle>(the_list.size()));
    std::vector<Halfedge_handle> &current_intersection_list = intersection_list.back();

    Halfedge_handle current_facet_start = polyline_start_vertex;

    // Temporarily store the interior points on current facet
    std::vector<Point_3> interior_points;

    // TODO: Handle the case where a facet has been splitted, so segment
    // doesn't end on current facet

    typename IntersectionList::const_iterator circulator_current = circulator_begin;
    const int polyline_offset = std::distance(the_list.begin(), circulator_begin);
    int facet_offset = 0;
    //int facet_segment_counter = 0;

    std::cout << "Circulating polyline" << std::endl;
    std::cout << "offset: " << polyline_offset << std::endl;
    // Circulate around the polyline from the start point
    do
    {
      const Segment &current_segment = circulator_current->get<2>();
      std::cout << "Current segment: " << current_segment << std::endl;

      // check if this is the last segment on this facet
      // Avoid geometric test here
      if (point_on_edge<Polyhedron>(current_facet_start->facet(), current_segment.target()))
      {
        std::cout << "Edge point" << std::endl;
        Halfedge_handle current_facet_end;

        // TODO: Reenable this test
        if (//current_facet_start->opposite()->facet() == polyline_start_vertex->facet() &&
            polyline_start_vertex->vertex()->point() == current_segment.target())
        {
          current_facet_end = polyline_start_vertex->opposite()->prev();
        }
        else
        {
          std::cout << "Getting/creating last vertex on facet" << std::endl;
          const Halfedge_handle end_vertex = edge_with_point_on<Polyhedron>(current_facet_start->facet(),
                                                                            current_segment.target());
          current_facet_end = a.split_edge(end_vertex);
          current_facet_end->vertex()->point() = current_segment.target();
          std::cout << "Done" << std::endl;
        }

        // CGAL does not allow multiedges (not even temporarily)
        // Split the edge (and join it afterwards) as a workaround
        // TODO
        Halfedge_handle tmp_intermediate;
        bool has_tmp_intermediate = false;
        if (current_facet_start->next() == current_facet_end)
        {
          tmp_intermediate = a.split_edge(current_facet_start->next());
          has_tmp_intermediate = true;
        }
        else if (current_facet_start->prev() == current_facet_end)
        {
          tmp_intermediate = a.split_edge(current_facet_start);
          has_tmp_intermediate = true;
        }

        // Insert the end vertex
        const Halfedge_handle inserted_edge = a.split_facet(current_facet_start, current_facet_end);
        intersection_facets.insert(inserted_edge);
        std::cout << "Inserting at position " << (polyline_offset + facet_offset + interior_points.size()) << std::endl;
        current_intersection_list[polyline_offset + facet_offset + interior_points.size()] = inserted_edge;
        

        // Insert all interior points
        for (uint c = 0; c < interior_points.size(); ++c)
        {
          Halfedge_handle i = a.split_edge(inserted_edge);
          i->vertex()->point() = interior_points[c];
          intersection_facets.insert(i);
          current_intersection_list[polyline_offset+facet_offset+c] = i;
          std::cout << "Inserting interior at position " << (polyline_offset+facet_offset+c) << std::endl;
        }

        // Clear intermediate edge (if it has been created)
        if (has_tmp_intermediate)
        {
          a.join_vertex(tmp_intermediate->opposite());
        }

        // Prepare for next segment
        facet_offset += interior_points.size()+1;
        interior_points.clear();
        current_facet_start = current_facet_end->opposite()->prev();
      }
      else
      {
        std::cout << "Interior point: " << current_segment.target() << std::endl;
        interior_points.push_back(current_segment.target());
      }

      circulator_current++;
      if (circulator_current == the_list.end())
      {
        circulator_current = the_list.begin();
      }
    } while (circulator_current != circulator_begin);
  }
}
//-----------------------------------------------------------------------------
// Check that all intersection segments are indeed edges
template<typename Polyhedron, int N>
bool check_splitting(const Polyhedron& a,
                     const std::list<std::list<boost::tuple<typename Polyhedron::Facet_handle,
                                                            typename Polyhedron::Facet_handle,
                                                            typename Polyhedron::Traits::Segment_3> > >& polylines,
                     const std::list<std::vector<typename Polyhedron::Halfedge_handle> >& intersection_list)
{
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Kernel::Segment_3 Segment;
  // typedef typename Kernel::Point_3 Point_3;
  typedef boost::tuple<Facet_handle, Facet_handle, Segment> IntersectionType;
  typedef std::list<IntersectionType> IntersectionList;


  std::cout << "Checking polyhedron" << std::endl;
  std::cout << "Number of intersection polylines: " << polylines.size() << std::endl;
  if (intersection_list.size() != polylines.size())
    std::cout << "Number of polylines doesn't match, " << intersection_list.size() << " != " << polylines.size() << std::endl;

  typename std::list<IntersectionList>::const_iterator polyline_it = polylines.begin();
  typename std::list<std::vector<Halfedge_handle> >::const_iterator edges_it = intersection_list.begin();
  while (polyline_it != polylines.end())
  {
    const IntersectionList &polyline = *polyline_it;
    const std::vector<typename Polyhedron::Halfedge_handle>& edges = *edges_it;

    std::cout << "Number of segments: " << polyline.size() << std::endl;
    if (polyline.size() != edges.size())
    {
      std::cout << "Number of segments doesn't match, " << polyline.size() << " != " << edges.size() << std::endl;
    }


    typename IntersectionList::const_iterator p_it = polyline.begin();
    typename std::vector<typename Polyhedron::Halfedge_handle>::const_iterator e_it = edges.begin();
    while (p_it != polyline.end())
    {
      std::cout << p_it->get<2>().source() << ", " << p_it->get<2>().target() << std::endl;
      std::cout << (*e_it)->vertex()->point() << ", " << (*e_it)->opposite()->vertex()->point() << std::endl << std::endl;
      // if (p_it->get<2>()

      p_it++;
      e_it++;
    }
    assert(e_it == edges.end());



    polyline_it++;
    edges_it++;
  }
  assert(edges_it == intersection_list.end());
  

  // Dummy return
  return true;
}
//-----------------------------------------------------------------------------
#endif 