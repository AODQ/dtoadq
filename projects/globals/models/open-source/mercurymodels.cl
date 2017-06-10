////////////////////////////////////////////////////////////////
//
//                           HG_SDF
//
//     GLSL LIBRARY FOR BUILDING SIGNED DISTANCE BOUNDS
//     (I ported to opencl)
//
//     version 2016-01-10
//
//     Check http://mercury.sexy/hg_sdf for updates
//     and usage examples. Send feedback to spheretracing@mercury.sexy.
//
//     Brought to you by MERCURY http://mercury.sexy
//
//
//
// Released as Creative Commons Attribution-NonCommercial (CC BY-NC)
//
////////////////////////////////////////////////////////////////


float sdSphere ( float3 p, float r ) { return length(p) - r; }

// Plane with normalized n at distance from origin
float sdPlane ( float3 p, float3 n, float dist ) {
  return dot(p, n) + dist;
}

// Cheap box (distance to corners is estimated
float sdCheapBox ( float3 p, float3 b ) { return fmax3(fabs(p) - b); }
float sdCheapBox2( float2 p, float2 b ) { return fmax2(fabs(p) - b); }

float sdBox ( float3 p, float3 b ) {
  float3 d = fabs(p) - b;
  return length(fmax(d, (float3)(0.0f))) + fmax3(fmin(d, (float3)(0.0f)));
}
float sdBox2( float2 p, float2 b ) {
  float2 d = fabs(p) - b;
  return length(fmax(d, (float2)(0.0f))) + fmax2(fmin(d, (float2)(0.0f)));
}

float sdCorner ( float2 p ) {
  return length(fmax(p, (float2)(0.0f))) + fmax2(fmin(p, (float2)(0.0f)));
}

// Cylinder standing upright on XZ plane
float sdCylinder ( float3 p, float r, float h ) {
  float d= length(p.xz) - r;
  d = fmax(d, fabs(p.y) - h);
  return d;
}

// A cylinder with round caps on both sides
float sdCapsule ( float3 p, float r, float c ) {
  return mix(length(p.zx) - r, length((float3)(p.x, fabs(p.y) - c, p.z)) - r,
                               step(c, fabs(p.y)));
}

// Line Segment between a and b
float sdLineSegment ( float3 p, float3 a, float3 b ) {
  float3 ab = b - a;
  float t = clamp(dot(p - a, ab)/dot(ab, ab), 0.0f, 1.0f);
  return length((ab*t + a) - p);
}

// Capsule between two end points, radius r
float sdCapsuleEnd ( float3 p, float3 a, float3 b, float r ) {
  return sdLineSegment ( p, a, b ) - r;
}

// Torus in XZ-plane
float sdTorus ( float3 p, float hole, float radius ) {
  return length((float2)(length(p.xz) - radius, p.y)) - hole;
}

// Circle line
float sdCircle ( float3 p, float r ) {
  float l = length(p.xz) - r;
  return length((float2)(p.y, l));
}

// Circular disc with no thickness, subtract to make rounded edge
float sdDisc ( float3 p, float r ) {
  float l = length(p.xz) - r;
  return l < 0.0f ? fabs(p.y) : length((float2)(p.y, l));
}

// Hexagonal prism, circumcircle variant
float sdHexagonCircumcircle ( float3 p, float2 h ) {
  float3 q = fabs(p);
  return fmax(q.y - h.y, fmax(q.x*sqrt(3.0f)*0.5f + q.z*0.5f, q.z) - h.x);
}

float sdHexagonIncircle ( float3 p, float2 h ) {
  return sdHexagonCircumcircle (p, (float2)(h.x*sqrt(3.0f)*0.5f, h.y));
}

// Cone with correct distance to tip and base circle
float sdCone ( float3 p, float radius, float height ) {
  float2 q = (float2)(length(p.xz), p.y);
  float2 tip = q - (float2)(0.0f, height);
  float2 mantle_dir = normalize((float2)(height, radius));
  float mantle = dot(tip, mantle_dir);
  float d = fmax(mantle, -q.y);
  float projected = dot(tip, (float2)(mantle_dir.y, -mantle_dir.x));

  if ( (q.y > height) && (projected < 0.0f) ) { d = fmax(d, length(tip)); }
  if ( (q.x > radius) && (projected > length((float2)(height, radius))) ) {
    d = fmax(d, length(q - (float2)(radius, 0.0f)));
  }
  return d;
}