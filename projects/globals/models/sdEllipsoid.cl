float sdEllipsoid ( float3 origin, float3 radius ) {
  return (length(origin/radius) - 1.0f) *
      fmin(fmin(radius.x, radius.y), radius.z);
}