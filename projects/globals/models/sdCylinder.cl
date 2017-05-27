
float sdCylinder ( float3 origin, float2 bounds ) {
  float2 d = fabs((float2)(length(origin.xz), origin.y)) - bounds;
  return fmin(fmax(d.x, d.y), 0.0f) + length(fmax(d, 0.0f));
}