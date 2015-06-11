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
#include "Polyhedron_utils.h"
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

#include <dolfin/log/log.h>
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
  void recursive_insert(Set& set,
                        std::set<typename Polyhedron::Vertex_const_handle>& visited,
                        typename Polyhedron::Vertex_const_handle v,
                        std::size_t n)
{
  std::pair<typename std::set<typename Polyhedron::Vertex_const_handle>::iterator, bool> v_insert = visited.insert(v);

  // If vertex is already visited, then return
  if (!v_insert.second)
  {
    return;
  }

  // Add point to set.
  set.insert(v->point());

  typename Polyhedron::Halfedge_around_vertex_const_circulator start = v->vertex_begin(), current = start;
  do
  {
    if ( set.size() >= n )
      break;

    recursive_insert<Set, Polyhedron>(set, visited, current->opposite()->vertex(), n);

    current++;
  } while (current != start);
}

//-----------------------------------------------------------------------------
template<typename IGT_>
template<class OutputIterator>
OutputIterator
Polyhedral_multicomponent_mesh_domain_with_features_3<IGT_>::
Construct_initial_points::operator()(OutputIterator pts, const int n) const
{
  typedef typename Polyhedral_multicomponent_mesh_domain_with_features_3::Polyhedron Polyhedron;
  typedef typename Polyhedron::Point_3 Point_3;
  typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;

  const Polyhedron& P  = r_domain_.polyhedron();
  std::list<Vertex_const_handle> components;
  mshr::PolyhedronUtils::get_disconnected_components(P, std::back_inserter(components));

  // Collect inserted points in a set with a fuzzy comparison operator
  // to ensure no points closer than the tolerance are inserted.
  typedef Point3FuzzyStrictlyLess<Point_3> CompareFunctor;
  typedef std::set<Point_3, CompareFunctor>  FuzzyPointSet;
  std::set<Vertex_const_handle> visited;

  //const CompareFunctor cf(edge_size);
  // TODO: Tune this parameter
  const double tolerance = edge_size*3;
  FuzzyPointSet inserted_points(CompareFunctor(tolerance*tolerance));

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

  // Insert n suraface points from each disconnected component
  for (typename std::list<Vertex_const_handle>::iterator it = components.begin();
       it != components.end(); it++)
  {
    Vertex_const_handle current = *it;
    recursive_insert<FuzzyPointSet, Polyhedron>(inserted_points, visited, current, n+inserted_points.size());
  }

  for (auto it = inserted_points.begin(); it != inserted_points.end(); it++)
  {
    *pts = std::make_pair(*it, current_index);
    pts++;
    current_index++;
  }

  return pts;
}
//-----------------------------------------------------------------------------
#endif
