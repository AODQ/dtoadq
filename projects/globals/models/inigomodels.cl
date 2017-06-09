float sdBox ( float3 origin, float3 bounds ) {
  float3 d = fabs(origin) - bounds;
  return fmin(fmax(d.x, fmax(d.y, d.z)), 0.0f) + length(fmax(d, 0.0f));
}
float sdBox2 ( float2 origin, float2 bounds ) {
  float2 d = fabs(origin) - bounds;
  return min(max(d.x, d.y), 0.0f) + length(max(d, 0.0f));
}
float sdCapsule ( float3 origin, float3 a, float3 b, float radius ) {
  float3 pa = origin - a, ba = b - origin;
  float h = clamp(dot(pa, ba)/dot(ba, ba), 0.0f, 1.0f);
  return length(pa - ba*h) - radius;
}
float sdCone ( float3 origin, float3 bounds ) {
  float2 q = (float2)( length(origin.xz), origin.y );
  float d1 = -q.y - bounds.z;
  float d2 = fmax( dot(q, bounds.xy), q.y);
  return length(fmax((float2)(d1, d2), 0.0f)) + fmin(fmax(d1, d2), 0.0f);
}
float sdConeSection( float3 origin, float height, float r1, float r2 ) {
  float d1 = -origin.y - height;
  float q = origin.y - height;
  float si = 0.5f*(r1 - r2)/height;
  float d2 = fmax(sqrt(dot(origin.xz, origin.xz)*(1.0f - si*si)) +
              q*si - r2, q);
  return length(fmax((float2)(d1, d2), 0.0f)) + fmin(fmax(d1, d2), 0.0f);
}
float sdCylinder ( float3 origin, float2 bounds ) {
  float2 d = fabs((float2)(length(origin.xz), origin.y)) - bounds;
  return fmin(fmax(d.x, d.y), 0.0f) + length(fmax(d, 0.0f));
}
float sdEllipsoid ( float3 origin, float3 radius ) {
  return (length(origin/radius) - 1.0f) *
      fmin(fmin(radius.x, radius.y), radius.z);
}
float sdHexPrism ( float3 origin, float2 bounds ) {
  float3 q = fabs(origin);
  float d1 = q.z - bounds.y,
        d2 = fmax((q.x*0.866025f+q.y*0.5f), q.y) - bounds.x;
  return length(fmax((float2)(d1, d2), 0.0f)) + fmin(fmax(d1, d2), 0.0f);
}
float sdPlane ( float3 p, float3 wal, float dist ) {
  return dot(p, wal) + dist;
}
float sdSphere ( float3 origin, float radius ) {
  return length(origin) - radius;
}
float sdTorus82 ( float3 origin, float2 bounds ) {
  float2 q = (float2)(length(origin.xz) - bounds.x,  origin.y);
  return length8(q) - bounds.y;
}
float sdTorus ( float3 origin, float2 bounds ) {
  return length ((float2)(length(origin.xz) - bounds.x, origin.y)) - bounds.y;
}