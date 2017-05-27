float sdHexPrism ( float3 origin, float2 bounds ) {
  float3 q = fabs(origin);
  float d1 = q.z - bounds.y,
        d2 = fmax((q.x*0.866025f+q.y*0.5f), q.y) - bounds.x;
  return length(fmax((float2)(d1, d2), 0.0f)) + fmin(fmax(d1, d2), 0.0f);
}