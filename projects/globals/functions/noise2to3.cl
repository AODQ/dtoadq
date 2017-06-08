
float3 noise2to3 ( float2 p ) {
  float3 q = (float3)( dot(p, (float2)(127.1f, 311.7f)),
                       dot(p, (float2)(269.5f, 183.3f)),
                       dot(p, (float2)(419.2f, 371.9f)) );
  return fract3(sin(q)*43758.5436f);
}