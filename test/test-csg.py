from mshr import *
from dolfin import *
import math

epsilon = 1e-15

def test_volume_exact(g, exact) :
    domain = CSGCGALDomain3D(g)
    vol = domain.volume()
    assert abs(vol-exact) < epsilon, "Computed volume of {} was {}, but expected {}".format(g.str(False), vol, exact)

def test_volume_bounds(g, lower, upper) :
    domain = CSGCGALDomain3D(g)
    vol = domain.volume()
    assert vol <= upper and vol >= lower, "Computed volume of {} was {}, but expected to be in [{}, {}]".format(g.str(False), vol, lower, upper)

    
test_volume_exact(Box(Point(0,0,0), Point(2,2,2)) + Box(Point(-1,-1,-1), Point(1,1,1)), 15)
test_volume_bounds(Sphere(Point(1,1,1), 1), 3.8/3.*math.pi, 4./3.0*math.pi)

