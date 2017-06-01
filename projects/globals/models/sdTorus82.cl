!require projects/globals/functions/length8.cl
         projects/globals/models/sdSphere.cl;
float sdTorus82 ( float3 origin, float2 bounds ) {
  float2 q = (float2)(length(origin.xz) - bounds.x,  origin.y);
  return length8(q) - bounds.y;
}