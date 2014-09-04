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
rotate_blades = True
include_tip = False
extra_rotation = True # only applied to inner mesh

# Define geometries
sphere = Sphere(Point(0, 0, 0), 2*R)
geometry_inside = CSGGeometries.propeller(r, R, w, h, rotate_blades, include_tip)
geometry_outside = sphere - geometry_inside

# Generate meshes
mesh_inside = generate_mesh(geometry_inside, 16)
mesh_outside = generate_mesh(geometry_outside, 16)

# Rotate blades
if extra_rotation:
    print "Rotating blades..."
    c = mesh_inside.coordinates()
    for i, (x, y, z) in enumerate(c):

        # Compute distance to axis
        _r = sqrt(x**2 + y**2)

        # Compute rotation angle
        v = -2*max(0, _r - r)

        # Rotate blades
        xx = x; yy = y; zz = z;
        if x > 0 and abs(y) < 5*h:
            yy = cos(v)*y - sin(v)*z
            zz = sin(v)*y + cos(v)*z
        elif x < 0 and abs(y) < 5*h:
            yy = cos(v)*y + sin(v)*z
            zz = -sin(v)*y + cos(v)*z
        elif y > 0 and abs(x) < 5*h:
            xx = cos(v)*x + sin(v)*z
            zz = -sin(v)*x + cos(v)*z
        elif y < 0 and abs(x) < 5*h:
            xx = cos(v)*x - sin(v)*z
            zz = sin(v)*x + cos(v)*z

        # Store coordinates
        c[i][0] = xx
        c[i][1] = yy
        c[i][2] = zz

# Save meshes to file
File("propeller_inside.xml.gz") << mesh_inside
File("propeller_outside.xml.gz") << mesh_outside

# Plot mesh
plot(mesh_inside)
plot(mesh_outside)
interactive()
