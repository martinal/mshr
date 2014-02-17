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

#ifndef _H_TRIANGULATE_POLYHEDRON__
#define _H_TRIANGULATE_POLYHEDRON__

namespace 
{
template <class Polyhedron>
  inline void triangulate_facet(Polyhedron &p,
                                typename Polyhedron::Facet_handle f)
  {
    if (f->facet_degree() > 4)
      dolfin::dolfin_error("triangulate_polyhedron.h",
                   "Triangulate polyhedron",
                   "Triangulation of facets with degree > 4 not implemented");


    typename Polyhedron::Halfedge_handle h = f->halfedge();

    if (f->facet_degree() == 4)
    {
      p.split_facet(h, h->next()->next());
    }

  }
}

namespace mshr
{

template <class Polyhedron>
void triangulate_polyhedron(Polyhedron &p)
{
  for (typename Polyhedron::Facet_iterator it = p.facets_begin();
       it != p.facets_end(); it++)
  {
    triangulate_facet(p, it);

  }

}

}
#endif
