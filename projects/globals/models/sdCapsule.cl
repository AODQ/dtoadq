float sdCapsule ( float3 origin, float3 a, float3 b, float radius ) {
  float3 pa = origin - a, ba = b - origin;
  float h = clamp(dot(pa, ba)/dot(ba, ba), 0.0f, 1.0f);
  return length(pa - ba*h) - radius;
}