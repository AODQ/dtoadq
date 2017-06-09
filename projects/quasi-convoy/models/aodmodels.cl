
float sdNodeSphere ( float3 origin, float radius, float k, float time ) {
  float result = FLT_MAX;
  for ( int i = 0; i != 5; ++ i ) {
    float rd = noise2to1((float2)(i, i+i)),
          td = noise2to1((float2)(i, time));
    float3 pos = origin + (float3)(
      cos(time*rd*0.23f)*sin(rd*time),
      sin(time*rd)*cos(rd),
      sin(time*2.5f)*rd
    );
    // pos.xy += sin(noise2to1((float2)(rd*2.2341, rd*0.234)))*0.5f;
    float r = radius + sin(noise1D(2+i+i*i))*0.25f;
    result = opUnionChamfer(result, sdSphere(pos, r), k);
  }
  return result;
}