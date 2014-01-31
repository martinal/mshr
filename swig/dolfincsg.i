
 %module dolfincsg
 %{
     /* Put header files here or function declarations like below */
    #include <dolfincsg/CSGPrimitive.h>
    #include <dolfincsg/CSGOperators.h>
    #include <dolfincsg/CSGPrimitives2D.h>
    #include <dolfincsg/CSGPrimitives3D.h>
 
    #include <dolfincsg/CSGGeometry.h>
    #include <dolfincsg/CSGGeometries3D.h>

    #include <dolfincsg/CSGMeshGenerator.h>
    #include <dolfincsg/CSGCGALMeshGenerator2D.h>
    #include <dolfincsg/CSGCGALMeshGenerator3D.h>

 %}

%ignore CSGOperator;
%ignore CSGGeometry;
%ignore CSGPrimitive;
%include <dolfincsg/CSGGeometry.h>
%include <dolfincsg/CSGPrimitive.h>
%include <dolfincsg/CSGOperators.h>
%include <dolfincsg/CSGPrimitives2D.h>
%include <dolfincsg/CSGPrimitives3D.h>

%include <dolfincsg/CSGMeshGenerator.h>
%include <dolfincsg/CSGCGALMeshGenerator2D.h>
%include <dolfincsg/CSGCGALMeshGenerator3D.h>

%include <dolfincsg/CSGGeometries3D.h>
