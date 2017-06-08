
float sdScienceRotator ( float3 origin, float2 ringdim, float time ) {
  float torus = sdTorus82(origin, ringdim);
  float prism = sdSphere(
    origin + (float3)(cos(time)*ringdim.x, 0.0f, sin(time)*ringdim.x),
    ringdim.x*1.0f
  );
  return opSubtract(prism, torus);
}