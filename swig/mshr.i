
 %module mshr
 %{
    #include <numpy/arrayobject.h>

    #include <dolfin/common/types.h>
    #include <dolfin/mesh/MeshFunction.h>

    #include <mshr/CSGPrimitive.h>
    #include <mshr/CSGOperators.h>
    #include <mshr/CSGPrimitives2D.h>
    #include <mshr/CSGPrimitives3D.h>
 
    #include <mshr/CSGGeometry.h>
    #include <mshr/CSGGeometries3D.h>

    #include <mshr/CSGMeshGenerator.h>
    #include <mshr/CSGCGALMeshGenerator2D.h>
    #include <mshr/CSGCGALMeshGenerator3D.h>

 %}

// Global typemaps and forward declarations
%include "dolfin/swig/typemaps/includes.i"
%include "dolfin/swig/forwarddeclarations.i"

// Global exceptions
%include <exception.i>
%include "dolfin/swig/exceptions.i"

// STL SWIG string class
%include <std_string.i>

// Local shared_ptr declarations
%include <boost_shared_ptr.i>
%shared_ptr(dolfin::Variable)
%shared_ptr(dolfin::Hierarchical<dolfin::Mesh>)
%shared_ptr(dolfin::Mesh)


%include "dolfin/swig/common/pre.i"
%import(module="common") "dolfin/common/Variable.h"
%import(module="common") "dolfin/common/Hierarchical.h"

%include "dolfin/swig/mesh/pre.i"
%import(module="mesh") "dolfin/mesh/Mesh.h"

%shared_ptr(mshr::CSGGeometry)
%shared_ptr(mshr::CSGPrimitive)
%shared_ptr(mshr::CSGOperator)
%shared_ptr(mshr::CSGUnion)
%shared_ptr(mshr::CSGDifference)
%shared_ptr(mshr::CSGIntersection)
%shared_ptr(mshr::CSGPrimitive2D)
%shared_ptr(mshr::Circle)
%shared_ptr(mshr::Ellipse)
%shared_ptr(mshr::Rectangle)
%shared_ptr(mshr::Polygon)
%shared_ptr(mshr::CSGPrimitive3D)
%shared_ptr(mshr::Sphere)
%shared_ptr(mshr::Box)
%shared_ptr(mshr::Cone)
%shared_ptr(mshr::Cylinder)
%shared_ptr(mshr::Tetrahedron)
%shared_ptr(mshr::Surface3D)
%shared_ptr(mshr::CSGCGALMeshGenerator2D)
%shared_ptr(mshr::CSGCGALMeshGenerator3D)

%ignore CSGOperator;
%ignore CSGGeometry;
%ignore CSGPrimitive;

%include <mshr/CSGGeometry.h>
%include <mshr/CSGPrimitive.h>
%include <mshr/CSGOperators.h>
%include <mshr/CSGPrimitives2D.h>
%include <mshr/CSGPrimitives3D.h>
%include <mshr/CSGMeshGenerator.h>
%include <mshr/CSGCGALMeshGenerator2D.h>
%include <mshr/CSGCGALMeshGenerator3D.h>
%include <mshr/CSGGeometries3D.h>
