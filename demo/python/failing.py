import mshr
import dolfin  

# Define 3D geometry  
cone = mshr.Cone(dolfin.Point(-1.0, 1.0, 1.0), dolfin.Point(1.0, -1.0, -1.0), .5, .5)
#csg::Cone cone2(Point(1.0, 1.0, 1.0), Point(-1.0, -1.0, -1.0), 0, .5)
cyl = mshr.Cylinder(dolfin.Point(1.0, -1.0, 1.0), dolfin.Point(-1.0, 1.0, -1.0), .5)
domain = cone + cyl;

# Test printing
dolfin.info("\nCompact output of 3D geometry:")
dolfin.info(domain)
dolfin.info("\nVerbose output of 3D geometry:")
dolfin.info(domain, True);

# Generate mesh
mesh3d = mshr.generate_mesh(domain, 10)

# Plot meshes
dolfin.plot(mesh3d, "3D mesh");
