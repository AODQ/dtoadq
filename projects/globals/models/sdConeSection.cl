float sdConeSection( float3 origin, float height, float r1, float r2 ) {
    float d1 = -origin.y - height;
    float q = origin.y - height;
    float si = 0.5f*(r1 - r2)/height;
    float d2 = fmax(sqrt(dot(origin.xz, origin.xz)*(1.0f - si*si)) +
               q*si - r2, q);
    return length(fmax((float2)(d1, d2), 0.0f)) + fmin(fmax(d1, d2), 0.0f);
}