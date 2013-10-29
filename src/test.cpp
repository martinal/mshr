// Copyright (C) 2013 Benjamin Kehlet
//
// This file is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2013-10-03
// Last changed: 2013-10-03

#include <intersection_segments.h>
#include <PolyhedronFactory.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef typename Polyhedron::Facet_handle Facet_handle;
typedef typename Kernel::Segment_3 Segment;

int main(int argc, char** argv)
{
  Polyhedron a;

  make_tetrahedron(a, 
                   Point(1.0, 0.0, 0.0),
                   Point(2.0, 0.0, 0.0),
                   Point(1.5, 1.0, 0.0),
                   Point(1.5, .5, 10.0));

  Polyhedron b;
  make_tetrahedron(b,
                   Point(0.0, 0., .5),
                   Point(0.0, 0.0, 1.5),
                   Point(0.0, 1.0, 1.0),
                   Point(10.0, .5, 1.0));

  if (a.is_pure_triangle())
    std::cout << "a is pure triangle" << std::endl;

  if (b.is_pure_triangle())
    std::cout << "b is pure triangle" << std::endl;

  std::list<boost::tuple<Facet_handle, Facet_handle, Segment> > intersections;

  Polyhedron &biggest = a.size_of_facets() > b.size_of_facets() ? a : b;
  Polyhedron &smallest = a.size_of_facets() > b.size_of_facets() ? b : a;

  compute_intersections(biggest, smallest, std::back_inserter(intersections));
  split_facets(biggest, smallest, intersections);

  return 0;
}
