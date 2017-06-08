float sdBox2 ( float2 origin, float2 dist ) {
  float2 d = fabs(origin) - dist;
  return fmin(fmax(d.x, d.y), 0.0f) + length(fmax(d, 0.0f));
}