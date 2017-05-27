RETURN float
PARAMS float2 origin float2 dist
REQUIRE
COMMENT
BEGIN

float sdBox2 ( float2 origin, float2 dist ) {
  float2 d = fabs(origin) - bounds;
  return fmin(fmax(d.x, d.y), 0.0f) + length(fmax(d, 0.0f));
}