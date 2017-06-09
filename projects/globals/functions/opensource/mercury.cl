
/// Mirror at an axis-aligned plane which is at a specified distance <dist>
/// from the origin
float opMirror ( float p, float dist, float* id ) {
  if ( id != 0 ) {
    *id = sign(p);
  }
  return fabs(p) - dist;
}

// Mirror in both dimensions and at the diagonal, yielding 1/8 of the space.
// translate by distance before mirror
float2 opMirrorOctant ( float2 p, float2 dist, float2* id ) {
  if ( id != 0 ) {
    *id = sign(p);
  }
  p.x = opMirror(p.x, dist.x, 0);
  p.y = opMirror(p.y, dist.y, 0);

  if ( p.y > p.x ) {
    p.xy = p.yx;
  }

  return p;
}

float opMod1 ( float p, float len, float* id ) {
  float hl = len*0.5f;
  if ( id != 0 ) {
    *id = floor((p + hl)/len);
  }
  return fmod(p + hl, len) - hl;
}

float2 opModMirror2 ( float2 p, float2 len, float2* id ) {
  float2 hl = len*0.5f;
  float2 c = floor((p + hl)/len);
  if ( id ) {
    *id = c;
  }
  p =  fmod(p + hl, len) - hl;
  p *= fmod(c, (float2)(2.0f))*2.0f - (float2)(1.0f);
  return p;
}

// Repeat around the origin by a fixed angle
float opModPolar ( float2 p, float repeats ) {
  float angle = 2*PI/repeats;
  float a = atan2(p.y, p.x) + angle/2.0f;
  float r = length(p);
  float c = floor(a/angle);
  a = fmod(a, angle) - angle/2.0f;
  p = (float2)(cos(a), sin(a))*r;
  // For odd number of repetitions, fix cell index of the cell in -x direction
  if ( fabs(c) >= (repeats/2.0f) ) c = fabs(c);
  return c;
}

float2 opRotate ( float2 p, float angle ) {
  return cos(angle)*p + sin(angle)*(float2)(p.y, -p.x);
}

float opUnionChamfer ( float a, float b, float r ) {
  return min(min(a, b), (a - r + b)*sqrt(0.5f));
}

float opUnionRound ( float a, float b, float r ) {
  float2 u = max((float2)(r - a, r - b), (float2)(0.0f));
  return max(r, min(a, b)) - length(u);
}