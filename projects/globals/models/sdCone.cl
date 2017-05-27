float sdCone ( float3 origin, float3 bounds ) {
    float2 q = (float2)( length(origin.xz), origin.y );
    float d1 = -q.y - bounds.z;
    float d2 = fmax( dot(q, bounds.xy), q.y);
    return length(fmax((float2)(d1, d2), 0.0f)) + fmin(fmax(d1, d2), 0.0f);
}