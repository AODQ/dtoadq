
float2 fract2 ( float2 n ) {
  float2 t = (float2)(0.0f);
  return fract(n, &t);
}