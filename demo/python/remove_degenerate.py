import mshr
import dolfin

TOLERANCE = 1e-10

# This geometry generates a lot of degenerate facets
cone = mshr.Cylinder(dolfin.Point(-1.0, 1.0, 1.0), dolfin.Point(1.0, -1.0, -1.0), .5, .5)
cyl = mshr.Cone(dolfin.Point(1.0, -1.0, 1.0), dolfin.Point(-1.0, 1.0, -1.0), .5)
geometry = cone + cyl;

polyhedral_domain = mshr.CSGCGALDomain3D(geometry)

print "Degenerate facets after boolean operation: {0}".format(polyhedral_domain.num_degenerate_facets(TOLERANCE))
polyhedral_domain.remove_degenerate_facets(TOLERANCE)

dolfin.info(polyhedral_domain, True)

