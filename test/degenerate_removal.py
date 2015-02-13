from mshr import *

# This union of spheres challenges the removal of degenerate facets
# since the intersection polyline of the two spheres matches lines in
# the triangulation of the spheres. Because of that roundoff errors
# introduces a lot of very small triangles when the union is carried
# out.

a = Sphere(dolfin.Point(0,0,0), .5)
b = Sphere(dolfin.Point(.5,0,0), .5)

domain = CSGCGALDomain3D(a+b)
domain.ensure_meshing_preconditions()

