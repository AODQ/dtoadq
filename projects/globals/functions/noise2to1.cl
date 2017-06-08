
// Dave Hoskins hash
float noise2to1 ( float2 n ) {
  float2 idk = (float2)(0.0f, 0.0f);
  float2 res = fract(n*1963.1844f, &idk);
  res += dot(res, res.yx+(float2)(19.19, 14.19));
  res = fract((res.xy + res.yx)*res.yx, &idk);
  return (res.y + res.x)/2.0f;
}