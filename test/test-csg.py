from mshr import *
from dolfin import *
import math

epsilon = 1e-11

def test_volume_exact(g, exact) :
    domain = CSGCGALDomain3D(g)
    vol = domain.volume()
    assert abs(vol-exact) < epsilon, "Computed volume of {} was {}, but expected {}.".format(g.str(False), vol, exact)

def test_volume_bounds(g, lower, upper) :
    domain = CSGCGALDomain3D(g)
    vol = domain.volume()
    assert vol <= upper and vol >= lower, "Computed volume of {} was {}, but expected to be in [{}, {}]".format(g.str(False), vol, lower, upper)

    
### Test volume of primitives
# Sphere
test_volume_bounds(Sphere(Point(1,1,1), 1), 3.8/3.*math.pi, 4./3.0*math.pi)

# TODO: Add primitives

# Test union
test_volume_exact(Box(Point(0,0,0), Point(2,2,2)) + Box(Point(-1,-1,-1), Point(1,1,1)), 15)

# TODO: Add operators

# Test translation
c = Cylinder(Point(-1, 1, 2), Point(3.5, -1, 2), 2, 3)
test_volume_exact(c, CSGCGALDomain3D(CSGTranslation(c, Point(2, 4, -5))).volume())

# Test rotation
e = Ellipsoid(Point(1,2,3), 2,4,6)
test_volume_exact(e, CSGCGALDomain3D(CSGRotation(e, Point(1,1,1), math.pi/2.)).volume())

# Test scaling
s = Sphere(Point(0,0,0), 2)
test_volume_exact(CSGScaling(s, 2), CSGCGALDomain3D(s).volume()*8)
