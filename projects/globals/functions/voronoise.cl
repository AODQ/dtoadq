
// from IQ: http://iquilezles.org/www/articles/voronoise/voronoise.htm
float voronoise ( float2 x, float u, float v ) {
  float2 p = floor(x),
         f = fract2(x);

  float k = 1.0f + 63.0f*pow(1.0f - v, 4.0f);
  float va = 0.0f, wt = 0.0f;
  for ( int j = -2; j != 2; ++ j )
  for ( int i = -2; i != 2; ++ i ) {
    float2 g = (float2)( (float)(i), (float)(j) );
    float3 o = (float3)(noise2to3(p+g)*(float3)(u, u, 1.0f));
    float2 r = g - f + o.xy;
    float d = dot(r, r);
    float w = pow(1.0f - smoothstep(0.0f, 1.414f, sqrt(d)), k);
    va += w*o.z;
    wt += w;
  }

  return va/wt;
}