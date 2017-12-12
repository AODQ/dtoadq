
float ffract ( float n ) {
  float t = 0.0f;
  return fract(n, &t);
}

float2 fract2 ( float2 n ) {
  float2 t = (float2)(0.0f);
  return fract(n, &t);
}

float lerp ( float a, float b, float c ) {
  return c < 0.0f ? a : c > 1.0f ? b :
         (mix(a, b, c));
}

float dabs ( float a ) {
  return a < 0.0f ? 0.0f : a;
}

float3 fract3 ( float3 n ) {
  float3 t = (float3)(0.0f);
  return fract(n, &t);
}

float sdTorus2 ( float3 p, float hole, float radius ) {
  return length2((float2)(length2(p.xz) - radius, p.y)) - hole;
}

float fmax2 ( float2 vec ) {return max(vec.x, vec.y);}
float fmin2 ( float2 vec ) {return min(vec.x, vec.y);}
float fmax3 ( float3 vec ) {return max(vec.x, max(vec.y, vec.z));}
float fmin3 ( float3 vec ) {return min(vec.x, min(vec.y, vec.z));}

float Shell ( float dist, float r ) {
  return fabs(dist) - r*0.5f;
}
