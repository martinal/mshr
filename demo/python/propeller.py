# Copyright (C) 2014 Anders Logg
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

from mshr import *
from dolfin import *

# Set parameters
r = 0.125
R = 0.5
w = 0.3
h = 0.025
include_tip = False

# Define geometry
geometry = CSGGeometries.propeller(r, R, w, h, include_tip)

# Generate mesh
mesh = generate_mesh(geometry, 16)

# Rotate blades
print "Rotating blades..."
c = mesh.coordinates()
for i, (x, y, z) in enumerate(c):

    # Compute distance to axis
    _r = sqrt(x**2 + y**2)

    # Compute rotation angle
    v = -2*max(0, _r - r)

    # Rotate blades
    xx = x; yy = y; zz = z;
    if x > 0 and abs(y) < 2*h:
        yy = cos(v)*y - sin(v)*z
        zz = sin(v)*y + cos(v)*z
    elif x < 0 and abs(y) < 2*h:
        yy = cos(v)*y + sin(v)*z
        zz = -sin(v)*y + cos(v)*z
    elif y > 0 and abs(x) < 2*h:
        xx = cos(v)*x + sin(v)*z
        zz = -sin(v)*x + cos(v)*z
    elif y < 0 and abs(x) < 2*h:
        xx = cos(v)*x - sin(v)*z
        zz = sin(v)*x + cos(v)*z

    # Store coordinates
    c[i][0] = xx
    c[i][1] = yy
    c[i][2] = zz

# Plot mesh
plot(mesh, interactive=True)

# Save mesh to file
file = File("propeller.xml.gz")
file << mesh
