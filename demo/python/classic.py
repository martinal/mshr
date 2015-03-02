# Copyright (C) 2015 Anders Logg
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

from dolfin import *
from mshr import *

# Create geometry
s1 = Sphere(Point(0, 0, 0), 1.4)
b1 = Box(Point(-1, -1, -1), Point(1, 1, 1))
c1 = Cylinder(Point(-2, 0, 0), Point(2, 0, 0), 0.8, 0.8)
c2 = Cylinder(Point(0, -2, 0), Point(0, 2, 0), 0.8, 0.8)
c3 = Cylinder(Point(0, 0, -2), Point(0, 0, 2), 0.8, 0.8)

geometry = s1*b1 - (c1 + c2 + c3)

# Create mesh
mesh = generate_mesh(geometry, 64)

# Save to file and plot
File("classic.pvd") << mesh
plot(mesh, interactive=True)
