# DTOADQ
GPU-Accelerated Light Transport on Sigend Distance Fields Video Emitter and Prototyping Tool.

# build
  Requires GLFW, OpenGL3.3, OpenCL1.2 with the cl\_khr\_gl\_sharing extension
    (use clinfo to check), OpenGL3.3, DUB and a D compiler.
  I've only tested this on NVIDIA hardware running Linux (Manjaro).

```
git clone https://github.com/AODQ/dtoadq
cd dtoadq && dub run
```

# usage
To create a project, create a new directory inside the projects directory.
  Scenes are .dtq, textures are .txt and generic functions are .cl. When DTOADQ
  generates the kernel, all generic functions from the directory the scene is
  located as well as the globals directory will be included. When any of the
  files in either directory change, all of the corresponding files and
  parameters are updated.

Generic Functions: All generic functions from the repository the scene is
  located as well as the globals directory are available. Here is an example
  generic function that, given a `float2`, will return a `float2` with the
fractional components of each element:
```
float2 fract2 ( float2 n ) {
  float2 t;
  return fract(n, &t);
}
```

Textures: Texture functions are not included, instead the textures they
  generate are included through the `\_\_read\_only image2d\_array\_t textures`
  parameter. To access them, inside the TEXTURESTART/END of your scene file,
  place the index followed by file location:
    `0 globals/textures/txCheckerboard.txt`
  These texture files generate a texture before the scene kernel call is made.
  The texture function call must match the filename, and any other function
    within will not be included for the scene call.
```
// from IQ: http://iquilezles.org/www/articles/voronoise/voronoise.htm
float voronoise ( float2 x, float u, float v ) {
  float2 p = floor(x),
         f = fract2(x);

  float k = 1.0f + 63.0f*pow(1.0f - v, 4.0f);
  float va = 0.0f, wt = 0.0f;
  for ( int j = -2; j != 2; ++ j )
  for ( int i = -2; i != 2; ++ i ) {
    float2 g = (float2)( (float)(i), (float)(j) );
    float3 o = (float3)(noise2to3(p+g)*(float3)(u, u, 1.0f));
    float2 r = g - f + o.xy;
    float d = dot(r, r);
    float w = pow(1.0f - smoothstep(0.0f, 1.414f, sqrt(d)), k);
    va += w*o.z;
    wt += w;
  }

  return va/wt;
}
float4 txCheckerboard ( float2 origin ) {
  return (float4)(voronoise(origin*20.0f, 0.2f, 2.2f));
}
```

Materials: The materials section of the DTQ scene file is parsed as a JSON
  string. Every value is a float from 0.0 to 1.0, but must be inputted as a
  string. The material describes a BSDF and these properties are allowed to be
  intermixed:
```
{"materials": [{
  "diffuse": "0.2", "specular": "0.0", "glossy": "0.2",
  "retroreflective": "0.0", "transmittive": "0.0"
}}
```
There will be support for more properties as well as microfacet via texture in the future.

Camera: The function `Update_Camera(Camera* camera, float time)` must exist
  within the CAMERASTART/END region. The Camera struct is defined as:
```
struct Camera {
  cl_float3 position, lookat, up;
  cl_int2 dimensions;
  cl_float fov;
  cl_int flags; // For imgui debug purposes
}
```


Emitter: The constant `__constant int EMITTER_AMT;` must be defined as a value
  greater than 0, and corresponds to the amount of emitters (light sources) in
  the given scene. The function `Emitter REmission ( int index, float3
  debug_values, float time )` must be also be defined, returning the
  corresponding emitter from the given index, defined by:
```
struct Emitter {
  float3 origin, emission;
  float radius;
}
```
Only area lights are supported. This is a function, and not parsed as JSON like
  materials, to allow the emitter's origin/emission/radius to change over time.

Map: The update map function must be defined
 `void Update_Map ( int avoid, float3 origin, SampledPt* pt_info, float time,
                   __read_only image2d_array_t textures, float3 debug_values )`
  When you are finished calculating the signed distance field to a primitive,
   you must call
   `MapUnionG(avoid, pt_info, sdf_distance, material_id, primitive_colour)`.
```
void Update_Map ( int avoid, float3 origin, SampledPt* pt_info, float time,
                  __read_only image2d_array_t textures, float3 debug_values ) {
  MapUnionG(avoid, pt_info, length(origin)-2.0f, 0, (float3)(0.9f, 0.2f, 0.3f));
}
```

Check the example scenes, functions and textures under `globals/`.
