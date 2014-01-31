
 %module mshr
 %{
     /* Put header files here or function declarations like below */
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
