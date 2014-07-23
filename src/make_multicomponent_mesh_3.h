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

#ifndef __MSHR_MAKE_MULTICOMPONENT_MESH_3_H
#define __MSHR_MAKE_MULTICOMPONENT_MESH_3_H

#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>


template<class C3T3, class MeshDomain, class MeshCriteria>
void make_multicomponent_mesh_3_impl(C3T3& c3t3,
                                     const MeshDomain&   domain,
                                     const MeshCriteria& criteria,
                                     const CGAL::parameters::internal::Exude_options& exude,
                                     const CGAL::parameters::internal::Perturb_options& perturb,
                                     const CGAL::parameters::internal::Odt_options& odt,
                                     const CGAL::parameters::internal::Lloyd_options& lloyd,
                                     const bool with_features,
                                     const CGAL::parameters::internal::Mesh_3_options& 
                                     mesh_options = CGAL::parameters::internal::Mesh_3_options())
{
  //std::cout << "Number of vertices initially: " << c3t3.triangulation().number_of_vertices() << std::endl;

  // Initialize c3t3 with points from the special features
  CGAL::internal::Mesh_3::C3t3_initializer< 
    C3T3,
    MeshDomain,
    MeshCriteria,
    CGAL::internal::Mesh_3::has_Has_features<MeshDomain>::value > () (c3t3,
                                                                      domain,
                                                                      criteria,
                                                                      with_features);

  // std::cout << "Number of vertices after features: " << c3t3.triangulation().number_of_vertices() << std::endl;
  
  // Inserts points from all connected components to the mesh
  CGAL::internal::Mesh_3::init_c3t3(c3t3, domain, criteria);
  // std::cout << "Number of vertices before meshing: " << c3t3.triangulation().number_of_vertices() << std::endl;
  
  // Build mesher and launch refinement process
  // Don't reset c3t3 as we just created it
  refine_mesh_3(c3t3, domain, criteria,
                exude, perturb, odt, lloyd, CGAL::parameters::no_reset_c3t3(), mesh_options);
}

#endif
