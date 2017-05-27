float sdTorus ( float3 origin, float2 bounds ) {
  return length ((float2)(length(origin.xz) - bounds.x, origin.y)) - bounds.y;
}