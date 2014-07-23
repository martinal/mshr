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

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

template<typename Polyhedron>
void recursive_remove(std::set<typename Polyhedron::Vertex_const_handle>& s, typename Polyhedron::Vertex_const_handle h)
{
  typedef typename Polyhedron::Halfedge_around_vertex_const_circulator HV_const_circulator;
  typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
  typedef Polyhedron Polyhedron_type;

  const HV_const_circulator start = h->vertex_begin(), current = start;
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

template <typename Polyhedron, typename OutputIterator>
void get_connected_components(const Polyhedron& p, OutputIterator it)
{
  typedef Polyhedron Polyhedron_t;
  typedef typename Polyhedron_t::Vertex_const_handle Vertex_const_handle;

  // store all vertices in a set
  std::set<Vertex_const_handle> v;
  for (typename Polyhedron_t::Vertex_const_iterator vit = p.vertices_begin(); vit != p.vertices_end(); vit++)
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

template< typename IGT_ >
class Polyhedral_multicomponent_mesh_domain_with_features_3
  : public CGAL::Polyhedral_mesh_domain_with_features_3< IGT_ >
{
 public:
  typedef typename CGAL::Polyhedral_mesh_domain_with_features_3< IGT_ > Base;
  typedef typename Base::Polyhedron Polyhedron;

  Polyhedral_multicomponent_mesh_domain_with_features_3(const Polyhedron& p)
    : Base(p)
  {}

  ~Polyhedral_multicomponent_mesh_domain_with_features_3(){}

  struct Construct_initial_points
  {
    Construct_initial_points(const Polyhedral_multicomponent_mesh_domain_with_features_3& domain)
      : r_domain_(domain) {}

    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int n = 8) const
    {
      std::cout << "Constructing initial points" << std::endl;

      // Insert points to make sure points from all components are
      // represented initally
      typedef typename Polyhedral_multicomponent_mesh_domain_with_features_3::Polyhedron Polyhedron;
      typedef typename Polyhedron::Point_3 Point_3;
      typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;

      const Polyhedron& P  = r_domain_.polyhedron();
      std::list<Vertex_const_handle> components;
      get_connected_components(P, std::back_inserter(components));
      std::cout << "Number of components: " << components.size() << std::endl;

      std::size_t current_index;
      {
        // get number of corners
        std::vector<std::pair<int, Point_3> > corners;
        r_domain_.get_corners(std::back_inserter(corners));
        current_index = corners.size();
        current_index++;
      }

      for (typename std::list<Vertex_const_handle>::iterator it = components.begin();
           it != components.end(); it++)
      {
        const Vertex_const_handle v = *it;

        typename Polyhedron::Halfedge_around_vertex_const_circulator start = v->vertex_begin(), current = start;
        do
        {
          *pts++ = std::make_pair(current->opposite()->vertex()->point(), current_index);
          current_index++;
          std::cout << "Inserting point: " << current->opposite()->vertex()->point() << ", index: " << current_index << std::endl;
          
          current++;
        } while (current != start);
      }
    }

  private:
    const Polyhedral_multicomponent_mesh_domain_with_features_3& r_domain_;
  };

  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }


};


#endif
