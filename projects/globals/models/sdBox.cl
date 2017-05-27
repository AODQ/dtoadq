float sdBox ( float3 origin, float3 bounds ) {
  float3 d = fabs(origin) - bounds;
  return fmin(fmax(d.x, fmax(d.y, d.z)), 0.0f) + length(fmax(d, 0.0f));
}
