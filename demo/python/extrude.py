# Copyright (C) 2015 Benjamin Kehlet
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

# This demo illustrates how a 2D geometry can be extruded to 3D.

from dolfin import *
import mshr

g2d = mshr.Circle(Point(0,0), 1.2) + mshr.Circle(Point(0, 1.2), 1.2)

# Any simple 2D geometry can be extruded to 3D
g3d = mshr.Extrude2D(g2d, 
                     .2) # The z "thickness"

m = mshr.generate_mesh(g3d, 15)
plot(m, interactive=True)
