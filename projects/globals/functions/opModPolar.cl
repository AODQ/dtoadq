// mercury hg sdf

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