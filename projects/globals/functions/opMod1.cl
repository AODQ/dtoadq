// mercury hg sdf

float opMod1 ( float* p, float len ) {
  float hl = len*0.5f;
  float c = floor((*p + hl)/len);
  *p = fmod(*p + hl, len) - hl;
  return c;
}