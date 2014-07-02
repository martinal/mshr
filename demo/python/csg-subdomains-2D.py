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

dolfin.set_log_level(dolfin.TRACE)

# Define 2D geometry
domain =   Rectangle(dolfin.Point(0., 0.), dolfin.Point(5., 5.)) \
         - Rectangle(dolfin.Point(2., 1.25), dolfin.Point(3., 1.75)) \
         - Circle(dolfin.Point(1, 4), .25) \
         - Circle(dolfin.Point(4, 4), .25)
domain.set_subdomain(1, Rectangle(dolfin.Point(1., 1.), dolfin.Point(4., 3.)))
domain.set_subdomain(2, Rectangle(dolfin.Point(2., 2.), dolfin.Point(3., 4.)))

dolfin.info("\nVerbose output of 2D geometry:")
dolfin.info(domain, True)

# Generate and plot mesh
mesh2d = generate_mesh(domain, 45)
print mesh2d
dolfin.plot(mesh2d, "2D mesh")

# Convert subdomains to mesh function for plotting
mf = dolfin.MeshFunction("size_t", mesh2d, 2, mesh2d.domains())
dolfin.plot(mf, "Subdomains")

dolfin.interactive()
