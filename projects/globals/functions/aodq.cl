
float2 fract2 ( float2 n ) {
  float2 t = (float2)(0.0f);
  return fract(n, &t);
}

float3 fract3 ( float3 n ) {
  float3 t = (float3)(0.0f);
  return fract(n, &t);
}

float fmax2 ( float2 vec ) {return max(vec.x, vec.y);}
float fmin2 ( float2 vec ) {return min(vec.x, vec.y);}
float fmax3 ( float3 vec ) {return max(vec.x, max(vec.y, vec.z));}
float fmin3 ( float3 vec ) {return min(vec.x, min(vec.y, vec.z));}

float Shell ( float dist, float r ) {
  return fabs(dist) - r*0.5f;
}