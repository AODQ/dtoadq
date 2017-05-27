float sdPlane ( float3 p, float3 wal, float dist ) {
  return dot(p, wal) + dist;
}