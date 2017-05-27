float length6 ( float2 p ) {
  p = p*p*p; p = p*p;
  return pow( p.x + p.y, 1.0f/6.0f );
}