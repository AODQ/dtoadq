float length8 ( float2 p ) {
  p = p*p; p = p*p; p = p*p;
  return pow ( p.x + p.y, 1.0f/8.0f );
}