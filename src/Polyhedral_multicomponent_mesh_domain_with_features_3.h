// Copyright (C) 2014 Benjamin Kehlet
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


#ifndef POLYHEDRAL_MULTICOMPONENT_MESH_DOMAIN_WITH_FEATURES_3_H
#define POLYHEDRAL_MULTICOMPONENT_MESH_DOMAIN_WITH_FEATURES_3_H

#include "Point3FuzzyStrictlyLess.h"
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

//-----------------------------------------------------------------------------
template<typename Polyhedron>
void recursive_remove(std::set<typename Polyhedron::Vertex_const_handle>& s,
                      typename Polyhedron::Vertex_const_handle h)
{
  typedef typename Polyhedron::Halfedge_around_vertex_const_circulator HV_const_circulator;
  typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
  typedef Polyhedron Polyhedron_type;

  const HV_const_circulator start = h->vertex_begin();
  HV_const_circulator current = start;
  do
  {
    Vertex_const_handle current_vertex = current->opposite()->vertex();
    assert(current_vertex != h);
    if (s.count(current_vertex))
    {
      s.erase(current_vertex);
      recursive_remove<Polyhedron_type>(s, current_vertex);
    }
    current++;
  } while (current != start);
}
//-----------------------------------------------------------------------------
// Scans the vertices of the polyhedron the polyhedron and returns a
// Polyhedron::Vertex_const_handle for each disconnected component.
template <typename Polyhedron, typename OutputIterator>
void get_disconnected_components(const Polyhedron& p, OutputIterator it)
{
  typedef Polyhedron Polyhedron_t;
  typedef typename Polyhedron_t::Vertex_const_handle Vertex_const_handle;

  // store all vertices in a set
  std::set<Vertex_const_handle> v;
  for (typename Polyhedron_t::Vertex_const_iterator vit = p.vertices_begin();
       vit != p.vertices_end(); vit++)
    v.insert(vit);

  while (!v.empty())
  {
    // Add the component to the output
    typename std::set<Vertex_const_handle>::iterator start_it = v.begin();
    Vertex_const_handle start = *start_it;
    v.erase(start_it);

    *it = start;
    it++;

    // Remove rest of component
    recursive_remove<Polyhedron_t>(v, start);
  }
}
//-----------------------------------------------------------------------------
// This class reimplements Construct_initial_points (from Polyhedral_mesh_domain)
// in order to make sure that all disconnected parts of the polyhedron are
// sufficiently covered. Otherwise the meshing algorithm may miss them entirely.
template< typename IGT_ >
class Polyhedral_multicomponent_mesh_domain_with_features_3
  : public CGAL::Polyhedral_mesh_domain_with_features_3< IGT_ >
{
 public:
  typedef typename CGAL::Polyhedral_mesh_domain_with_features_3< IGT_ > Base;
  typedef typename Base::Polyhedron Polyhedron;

// Passing the edge size with the constructor is a workaround. Ideally CGAL should pass it
// when calling construct_initial_points
Polyhedral_multicomponent_mesh_domain_with_features_3(const Polyhedron& p, double edge_size)
  : Base(p), edge_size(edge_size)
  {}

  ~Polyhedral_multicomponent_mesh_domain_with_features_3(){}

  struct Construct_initial_points
  {
    Construct_initial_points(const Polyhedral_multicomponent_mesh_domain_with_features_3& domain,
                             double edge_size)
     : r_domain_(domain), edge_size(edge_size) {}

    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int n = 8) const;

   private:
    const Polyhedral_multicomponent_mesh_domain_with_features_3& r_domain_;
    const double edge_size;
  };

  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this, edge_size);
  }

 private :
  const double edge_size;
};
//-----------------------------------------------------------------------------
template<typename Set, typename Polyhedron>
int recursive_insert(Set& set, typename Polyhedron::Vertex_const_handle v, int i, int n)
{
  std::pair<typename Set::iterator, bool> res = set.insert(v->point());
  if ( res.second )
    i++;

  typename Polyhedron::Halfedge_around_vertex_const_circulator start = v->vertex_begin(), current = start;
  do
  {
    if ( i ==n )
      break;

    i += recursive_insert<Set, Polyhedron>(set, current->opposite()->vertex(), i, n);

    current++;
  } while (current != start);

  return i;
}

//-----------------------------------------------------------------------------
template<typename IGT_>
template<class OutputIterator>
OutputIterator
Polyhedral_multicomponent_mesh_domain_with_features_3<IGT_>::
Construct_initial_points::operator()(OutputIterator pts, const int n) const
{
  std::cout << "Constructing initial points" << std::endl;

  typedef typename Polyhedral_multicomponent_mesh_domain_with_features_3::Polyhedron Polyhedron;
  typedef typename Polyhedron::Point_3 Point_3;
  typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;

  const Polyhedron& P  = r_domain_.polyhedron();
  std::list<Vertex_const_handle> components;
  get_disconnected_components(P, std::back_inserter(components));
  std::cout << "Number of components: " << components.size() << std::endl;

  // Store inserted points in a set with a fuzzy comparison operator
  // to ensure no points closer than the tolerance are inserted.
  typedef Point3FuzzyStrictlyLess<Point_3> CompareFunctor;
  typedef std::set<Point_3, CompareFunctor>  FuzzyPointSet;

  const CompareFunctor cf(edge_size);
  FuzzyPointSet inserted_points(cf);

  std::size_t current_index;
  {
    // get corners
    std::vector<std::pair<int, Point_3> > corners;
    r_domain_.get_corners(std::back_inserter(corners));
    current_index = corners.size();
    current_index++;
    for (typename std::vector<std::pair<int, Point_3> >::iterator it = corners.begin();
         it != corners.end(); it++)
    {
      inserted_points.insert(it->second);
    }
  }

  for (typename std::list<Vertex_const_handle>::iterator it = components.begin();
       it != components.end(); it++)
  {
    Vertex_const_handle current = *it;

    int i = recursive_insert<FuzzyPointSet, Polyhedron>(inserted_points, current, 0, n);
    std::cout << "Inserted " << i << " points from disconnected part" << std::endl;
  }

  return pts;
}
//-----------------------------------------------------------------------------
#endif
