# Copyright (C) 2012-2014 Benjamin Kehlet
#
# This file is part of mshr.
#
# mshr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# mshr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with mshr.  If not, see <http://www.gnu.org/licenses/>.

import dolfin
from mshr import *

# Define 3D geometry
box = Box(dolfin.Point(0, 0, 0), dolfin.Point(1, 1, 1))
sphere = Sphere(dolfin.Point(0, 0, 0), 0.3)
cylinder = Cylinder(dolfin.Point(0, 0, -1), dolfin.Point(0, 0, 1), 1., .5)

domain = box + cylinder - sphere

# Test printing
dolfin.info("\nCompact output of 3D geometry:")
dolfin.info(domain)
dolfin.info("\nVerbose output of 3D geometry:")
dolfin.info(domain, True)

# Creating a mesh generator object gives access to parameters of the
# meshing backend
generator = CSGCGALMeshGenerator3D()
generator.parameters["edge_size"] = 0.025
generator.parameters["facet_angle"] = 25.0
generator.parameters["facet_size"] = 0.05

m = dolfin.Mesh()
generator.generate(domain, m)

dolfin.plot(m, "3D mesh")

dolfin.interactive()
