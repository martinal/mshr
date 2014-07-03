# Copyright (C) 2014 Benjamin Kehlet
#
# This file is part of mshr.
#
# mshr is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# mshr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with mshr. If not, see <http://www.gnu.org/licenses/>.
#

import dolfin
from mshr import *

dolfin.set_log_level(dolfin.TRACE)

# Define 3D geometry
sphere = Sphere(dolfin.Point(0, 0, 0), 0.5)
cone = Cylinder(dolfin.Point(0, 0, 0), dolfin.Point(0, 0, -1), .35, .1)

geometry = cone + sphere

# Geometry surfaces can be saved to off files
# which can be viewed by eg. MeshLab
meshing_domain = CSGCGALDomain3D(geometry)
meshing_domain.remove_degenerate_facets(1e-12)
meshing_domain.save_off("icecream.off")

# Test printing
dolfin.info("\nCompact output of 3D geometry:")
dolfin.info(geometry)
dolfin.info("\nVerbose output of 3D geometry:")
dolfin.info(geometry, True)

# Generate and plot mesh
m = generate_mesh(geometry, 16, "cgal")

dolfin.info(m)
dolfin.plot(m, "3D mesh")

dolfin.interactive()
