
float2 opRotate ( float2 p, float angle ) {
  return cos(angle)*p + sin(angle)*(float2)(p.y, -p.x);
}