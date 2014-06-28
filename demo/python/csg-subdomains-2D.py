# Copyright (C) 2012-2013 Benjamin Kehlet
#
# This file is part of DOLFIN.
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
# First added:  2012-11-12
# Last changed: 2013-03-15
# Begin demo

import dolfin
from mshr import *

dolfin.set_log_level(dolfin.TRACE)

# Define 2D geometry
domain = Rectangle(0., 0., 5., 5.) - Rectangle(2., 1.25, 3., 1.75) - Circle(1, 4, .25) - Circle(4, 4, .25)
domain.set_subdomain(1, Rectangle(1., 1., 4., 3.))
domain.set_subdomain(2, Rectangle(2., 2., 3., 4.))

dolfin.info("\nVerbose output of 2D geometry:")
dolfin.info(domain, True)

# # Plot geometry
# plot(domain, "2D Geometry (boundary)")

# Generate and plot mesh
mesh2d = generate_mesh(domain, 45)
print mesh2d
#dolfin.plot(mesh2d, "2D mesh")

# Convert subdomains to mesh function for plotting
mf = dolfin.MeshFunction("size_t", mesh2d, 2, mesh2d.domains())
dolfin.plot(mf, "Subdomains")


dolfin.interactive()
