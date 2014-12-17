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


#ifndef POLYHEDRON_UTILS_H__
#define POLYHEDRON_UTILS_H__

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

#endif
