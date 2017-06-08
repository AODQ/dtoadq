// mercury hg sdf

/// Mirror at an axis-aligned plane which is at a specified distance <dist>
/// from the origin
float opMirror ( float p, float dist, float* id ) {
  if ( id != 0 ) {
    *id = sign(p);
  }
  return fabs(p) - dist;
}