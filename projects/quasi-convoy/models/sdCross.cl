RETURN float
PARAMS float3 origin float radius
REQUIRE sdBox2
COMMENT
BEGIN

float sdCross ( float3 origin, float dist ) {
  float da = sdBox2(origin.xy, (float2)(dist));
  float db = sdBox2(origin.yz, (float2)(dist));
  float dc = sdBox2(origin.zx, (float2)(dist));
  return fmin(da, fmin(db, dc));
}