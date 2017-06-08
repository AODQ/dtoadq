
float2 opModMirror2 ( float2* p, float2 len ) {
  float2 hl = len*0.5f;
  float2 c = floor((*p + hl)/len);
  *p =  fmod(*p + hl, len) - hl;
  *p *= fmod(c, (float2)(2.0f))*2.0f - (float2)(1.0f);
  return c;
}