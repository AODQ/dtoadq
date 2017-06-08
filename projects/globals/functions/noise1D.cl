
float noise1D ( float n ) {
  return noise2to1((float2)(n, n*1.3f));
}