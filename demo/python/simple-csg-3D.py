# Copyright (C) 2012 Benjamin Kehlet
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
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#

import dolfin
from mshr import *

# Define 3D geometry
box = Box(0, 0, 0, 1, 1, 1)
sphere = Sphere(Point(0, 0, 0), 0.3)
cone = Cone(Point(0, 0, -1), Point(0, 0, 1), 1., .5)

domain = box + cone - sphere

# Test printing
info("\nCompact output of 3D geometry:")
info(domain)
info("\nVerbose output of 3D geometry:")
info(domain, True)

# Generate and plot mesh
m = dolfin.Mesh()
CSGMeshGenerator.generate(m, domain, 32)

mesh3d = Mesh(g3d, 32)
info(mesh3d)
plot(mesh3d, "3D mesh")

interactive()
