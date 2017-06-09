
float length2 ( float2 p ) { return sqrt(p.x*p.x + p.y*p.y); }
float length6 ( float2 p ) {
  p = p*p*p; p = p*p;
  return pow( p.x + p.y, 1.0f/6.0f );
}
float length8 ( float2 p ) {
  p = p*p; p = p*p; p = p*p;
  return pow ( p.x + p.y, 1.0f/8.0f );
}


float noise1D ( float n ) {
  return noise2to1((float2)(n, n*1.3f));
}
float noise2to1 ( float2 n ) {
  float2 idk = (float2)(0.0f, 0.0f);
  float2 res = fract(n*1963.1844f, &idk);
  res += dot(res, res.yx+(float2)(19.19, 14.19));
  res = fract((res.xy + res.yx)*res.yx, &idk);
  return (res.y + res.x)/2.0f;
}
float3 noise2to3 ( float2 p ) {
  float3 q = (float3)( dot(p, (float2)(127.1f, 311.7f)),
                       dot(p, (float2)(269.5f, 183.3f)),
                       dot(p, (float2)(419.2f, 371.9f)) );
  return fract3(sin(q)*43758.5436f);
}


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