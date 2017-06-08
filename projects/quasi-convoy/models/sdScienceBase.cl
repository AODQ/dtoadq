
float sdScienceBase ( float3 origin, float scale, float time ) {
  origin.yx = opRotate(origin.yx, sin(time*2.0f)*0.686f);
  origin.yz = opRotate(origin.yz, cos(time*1.0f)*0.232f);
  float cone = sdConeSection ( origin, scale, scale*1.9f, scale*0.4f);
  return cone;
}