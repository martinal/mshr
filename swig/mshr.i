 %module mshr
 %{
    #include <numpy/arrayobject.h>

    #include <dolfin/common/types.h>
    #include <dolfin/mesh/MeshFunction.h>

    #include <mshr/CSGPrimitive.h>
    #include <mshr/CSGOperators.h>
    #include <mshr/CSGPrimitives2D.h>
    #include <mshr/CSGPrimitives3D.h>
    #include <mshr/CSGCGALDomain2D.h>
    #include <mshr/CSGCGALDomain3D.h>

    #include <mshr/CSGGeometry.h>
    #include <mshr/CSGGeometries3D.h>

    #include <mshr/MeshGenerator.h>
    #include <mshr/CSGCGALMeshGenerator2D.h>
    #include <mshr/CSGCGALMeshGenerator3D.h>
    #include <mshr/DolfinMeshUtils.h>

 %}

%init
%{
  import_array();
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
%include <std_shared_ptr.i>
%shared_ptr(dolfin::Variable)
%shared_ptr(dolfin::Hierarchical<dolfin::Mesh>)

%ignore dolfin::Variable::id;
%ignore dolfin::Variable::str;

%include "dolfin/swig/common/pre.i"
%import(module="dolfin") "dolfin/common/Variable.h"
%import(module="dolfin") "dolfin/common/Hierarchical.h"

%shared_ptr(dolfin::Mesh)

%ignore dolfin::Mesh::type;
%ignore dolfin::Mesh::hash;

%import "dolfin/swig/mesh/pre.i"
%import(module="dolfin") "dolfin/mesh/Mesh.h"

%ignore mshr::CSGGeometry::getType;
%ignore mshr::CSGUnion::_g0;
%ignore mshr::CSGUnion::_g1;
%ignore mshr::CSGUnion::getType;
%ignore mshr::CSGOperator::CSGOperator;

/* %ignore mshr::operator+(mshr::CSGGeometry& g0, mshr::CSGGeometry& g1); */
/* %ignore mshr::operator+(std::shared_ptr<mshr::CSGGeometry>, std::shared_ptr<mshr::CSGGeometry>); */
/* %ignore mshr::operator+(mshr::CSGGeometry&, std::shared_ptr<mshr::CSGGeometry>); */
/* %ignore mshr::operator+(std::shared_ptr<mshr::CSGGeometry>, mshr::CSGGeometry&); */

%ignore mshr::operator+;
%ignore mshr::operator-;
%ignore mshr::operator*;


%ignore mshr::CSGCGALDomain2D::operator=;
%ignore mshr::CSGCGALDomain2D::impl;
%ignore mshr::PSLG;

%shared_ptr(mshr::CSGGeometry)
%shared_ptr(mshr::CSGPrimitive)
%shared_ptr(mshr::CSGOperator)
%shared_ptr(mshr::CSGUnion)
%shared_ptr(mshr::CSGDifference)
%shared_ptr(mshr::CSGIntersection)
%shared_ptr(mshr::CSGTranslation)
%shared_ptr(mshr::CSGScaling)
%shared_ptr(mshr::CSGRotation)
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
%shared_ptr(mshr::CSGCGALDomain2D)
%shared_ptr(mshr::CSGCGALDomain3D)

%rename(_generate) mshr::generate;
%ignore mshr::get_boundary_mesh;

%include <mshr/CSGGeometry.h>
%include <mshr/CSGPrimitive.h>
%include <mshr/CSGOperators.h>
%include <mshr/CSGPrimitives2D.h>
%include <mshr/CSGPrimitives3D.h>
%include <mshr/CSGCGALDomain2D.h>
%include <mshr/CSGCGALDomain3D.h>
%include <mshr/MeshGenerator.h>
%include <mshr/CSGCGALMeshGenerator2D.h>
%include <mshr/CSGCGALMeshGenerator3D.h>
%include <mshr/CSGGeometries3D.h>
%include <mshr/DolfinMeshUtils.h>

%extend mshr::CSGGeometry {
  %pythoncode %{
     def __add__(self, other) :
         if isinstance(other, dolfin.Point) :
             return CSGTranslation(self, other)
         else :
             return CSGUnion(self, other)

     def __mul__(self, other) :
         from numbers import Number
         if isinstance(other, Number) :
             return CSGScaling(self, other)
         else:
             return CSGIntersection(self, other)

     def __sub__(self, other) :
         return CSGDifference(self, other)

  %}
 }

// Declare the generate function here. The signature must
// be in sync with mshr/MeshGenerator.h:generate()
%pythoncode {
      def generate_mesh(geometry, resolution, backend="cgal") :
          m = dolfin.Mesh()
          _generate(m, geometry, resolution, backend)
          return m
%}
