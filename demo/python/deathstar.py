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
from math import pi, sin, cos, sqrt

# Parameters
R = 1.1
r = 0.4
t = 10
x = R*cos(float(t) / 180 * pi)
y = 0
z = R*sin(t)

# Create geometry
s1 = Sphere(Point(0, 0, 0), 1)
s2 = Sphere(Point(x, y, z), r)
b1 = Box(Point(-2, -2, -0.03), Point(2, 2, 0.03))
geometry = s1 - s2 - b1

# Create mesh
mesh = generate_mesh(geometry, 32)

# Save to file and plot
File("deathstar.pvd") << mesh
plot(mesh, interactive=True)
