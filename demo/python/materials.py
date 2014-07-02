# Copyright (C) 2014 Benjamin Kehlet
#
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

dolfin.set_log_level(dolfin.TRACE)

# Define 2D geometry
domain = Rectangle(dolfin.Point(0., 0.), dolfin.Point(1., 1.)) - Circle(dolfin.Point(0.0, 0.0), .35)
domain.set_subdomain(1, Rectangle(dolfin.Point(.0, .0), dolfin.Point(.95, .95)))
domain.set_subdomain(2, Circle(dolfin.Point(0, 0), .45))
domain.set_subdomain(3, Circle(dolfin.Point(0,0), .6))

# Generate and plot mesh
mesh2d = generate_mesh(domain, 45)
dolfin.plot(mesh2d, "2D mesh")

# Convert subdomains to mesh function for plotting
mf = dolfin.MeshFunction("size_t", mesh2d, 2, mesh2d.domains())
dolfin.plot(mf, "Subdomains")

dolfin.interactive()
