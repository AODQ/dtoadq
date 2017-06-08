
float opUnionRound ( float a, float b, float r ) {
  float2 u = max((float2)(r - a, r - b), (float2)(0.0f));
  return max(r, min(a, b)) - length(u);
}