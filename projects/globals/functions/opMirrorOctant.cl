// mercury hg sdf

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