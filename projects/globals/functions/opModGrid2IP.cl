
// Repeat in 2 dimensions mirroring every second cell at diagonal
// in place
float2 opModGrid2IP ( float2 p, float2 len ) {
  float2 c = floor((p + len*0.5f)/len);
  p = fmod(p + len*0.5f, len) - len*0.5f;
  p *= fmod(c, (float2)(2.0f))*2.0f - (float2)(1.0f);
  p -= len/2.0f;
  if ( p.x > p.y ) p.xy = p.yx;
  return p;
}