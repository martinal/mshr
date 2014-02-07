// Copyright (C) 2013 Benjamin Kehlet
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

#include <mshr/CSGCGALDomain3D.h>

#define CGAL_NO_DEPRECATED_CODE
#define CGAL_MESH_3_VERBOSE
//#define PROTECTION_DEBUG

#define CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
#define CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS

#include <CGAL/basic.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

namespace mshr
{
  namespace csg
  {

    // Exact polyhedron
    typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
    typedef Exact_Kernel::Triangle_3 Exact_Triangle_3;
    typedef CGAL::Nef_polyhedron_3<Exact_Kernel> Nef_polyhedron_3;
    typedef CGAL::Polyhedron_3<Exact_Kernel> Exact_Polyhedron_3;
    typedef Exact_Polyhedron_3::HalfedgeDS Exact_HalfedgeDS;
    typedef Nef_polyhedron_3::Point_3 Exact_Point_3;

  }
}



namespace mshr
{

struct CSGCGALDomain3DImpl
{
  csg::Nef_polyhedron_3 polyhedron;
};


CSGCGALDomain3D::CSGCGALDomain3D()
: impl(new CSGCGALDomain3DImpl)
{}
//-----------------------------------------------------------------------------
CSGCGALDomain3D::CSGCGALDomain3D(const mshr::CSGGeometry &csg)
: impl(new CSGCGALDomain3DImpl)
{
}
//-----------------------------------------------------------------------------
CSGCGALDomain3D::~CSGCGALDomain3D(){}
//-----------------------------------------------------------------------------
void CSGCGALDomain3D::remove_degenerated() {}

} // end namespace mshr
