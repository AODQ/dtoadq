
float3 fract3 ( float3 n ) {
  float3 t = (float3)(0.0f);
  return fract(n, &t);
}