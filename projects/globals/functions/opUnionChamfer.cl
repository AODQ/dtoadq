
float opUnionChamfer ( float a, float b, float r ) {
  return min(min(a, b), (a - r + b)*sqrt(0.5f));
}