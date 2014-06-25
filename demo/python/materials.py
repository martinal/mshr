# Copyright (C) 2014 Benjamin Kehlet
#
# This file is part of mshr.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#

import dolfin
from mshr import *


dolfin.set_log_level(dolfin.TRACE)

# Define 2D geometry
domain = Rectangle(0., 0., 1., 1.) - Circle(0.0, 0.0, .35)
domain.set_subdomain(1, Circle(0, 0, .45))
domain.set_subdomain(2, Circle(0,0, .6))

# Generate and plot mesh
mesh2d = generate_mesh(domain, 45)
dolfin.plot(mesh2d, "2D mesh")

# Convert subdomains to mesh function for plotting
mf = dolfin.MeshFunction("size_t", mesh2d, 2, mesh2d.domains())
dolfin.plot(mf, "Subdomains")


dolfin.interactive()
