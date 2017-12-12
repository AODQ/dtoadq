TEXTURESSTART
TEXTURESEND
MATERIALSSTART
{
  
  "materials": [
    {
      "albedo":       "[0.7, 0.6, 0.4]",
      "diffuse":      "1.0",
      "specular":     "0.0",
      "glossy":       "0.0",
      "glossy_lobe":  "0.95",
      "transmittive": "0.0",
      "ior": "0.5",
      "roughness":    "0.6",
      "metallic":     "0.2",
      "fresnel":      "0.9",
      "subsurface":   "0.0",
      "anisotropic":  "0.8"
    },
    {
      "albedo":       "[0.2, 0.1, 0.4]",
      "diffuse":      "1.0",
      "specular":     "0.0",
      "glossy":       "0.0",
      "glossy_lobe":  "0.9",
      "transmittive": "0.0",
      "ior": "0.284",
      "roughness":    "0.8",
      "metallic":     "0.2",
      "fresnel":      "0.8",
      "subsurface":   "0.1",
      "anisotropic":  "0.3"
    },
    {
      "albedo":       "[0.65, 0.55, 0.42]",
      "diffuse":      "1.0",
      "specular":     "0.0",
      "glossy":       "0.0",
      "glossy_lobe":  "0.6",
      "transmittive": "0.0",
      "ior": "0.5",
      "roughness":    "0.912",
      "metallic":     "0.0",
      "fresnel":      "0.882",
      "subsurface":   "0.25",
      "anisotropic":  "0.706"
    },
    {
      "albedo":       "[1.0, 1.0, 1.0]",
      "diffuse":      "0.0",
      "specular":     "1.0",
      "glossy":       "0.0",
      "glossy_lobe":  "0.6",
      "transmittive": "0.0",
      "ior":          "0.5",
      "roughness":    "0.912",
      "metallic":     "0.0",
      "fresnel":      "0.882",
      "subsurface":   "0.25",
      "anisotropic":  "0.706"
    },
    {
      "albedo":       "[1.0, 1.0, 1.0]",
      "diffuse":      "0.0",
      "specular":     "0.0",
      "glossy":       "0.0",
      "glossy_lobe":  "0.6",
      "transmittive": "1.0",
      "ior":          "0.11",
      "roughness":    "0.912",
      "metallic":     "0.0",
      "fresnel":      "0.882",
      "subsurface":   "0.25",
      "anisotropic":  "0.706"
    }
  ]
}
MATERIALSEND

CAMERASTART
void Update_Camera ( Camera* camera, float time ) {
  // if ( camera->flags > 0 ) return;

  // camera->position = (float3)(0.0f);
  // camera->lookat   = (float3)(cos(time*0.5f), sin(time*0.5f), cos(time)*0.1f);
}
CAMERAEND

EMITTERSTART
__constant int EMITTER_AMT = 3;
Emitter REmission ( int index, float3 dval, float time ) {
  float f = (float)(index+1);
  float3 origin = (float3)(cos(time+f*234.0f)*2.0f-6.0f, 8.0f, sin(time+f*234.0f)*2.0f);
  if ( index == 0 ) // most light
    return (Emitter){origin, (float3)(0.1f)+dval.x, 1.0f, index};
  origin = (float3)(-0.7f, 1.32, -5.38f);
  origin += (float3)(sin(time)*2.5f, 0.0f, cos(time)*0.75f);
  if ( index == 1 ) // caustic light
    return (Emitter){origin, (float3)(0.5f)*dval.y, dval.z, index};
  origin = (float3)(0.5f, 0.0f, 1.18f);
  if ( index == 2 ) // glossy light
    return (Emitter){origin, (float3)(0.8f, 0.2f, 0.7f)*0.3f, 0.25f, index};
}
EMITTEREND

UPDATEMAPSTART
void Room ( int avoid, float3 origin, SampledPt* pt, float time, float3 dval ) {
  float dist;
  float3 colour = (float3)(1.0f);
  dist = sdPlane(origin, normalize((float3)(0.0f, 0.0f, 1.0f)), 10.5f);
  colour = (float3)(sin(origin.x*PI*5.0f)*sin(origin.y*PI*5.0f) >= 0.0f) + 0.5f;
  colour *= (float3)(1.0f, 0.0f, 0.0f);
  MapUnionG(avoid, pt, dist, 0, colour);
  dist = sdPlane(origin, normalize((float3)(0.0f, 0.0f, -1.0f)), 10.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  colour = (float3)(sin(origin.z*PI*5.0f)*sin(origin.y*PI*5.0f) >= 0.0f) + 0.5f;
  colour *= (float3)(0.0f, 0.5f, 1.0f);
  dist = sdPlane(origin, normalize((float3)(-1.0f, 0.0f, 0.0f)), 10.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  dist = sdPlane(origin, normalize((float3)(1.0f, 0.0f, 0.0f)), 10.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  colour = (float3)(sin(origin.x*PI*5.0f)*sin(origin.z*PI*5.0f) >= 0.0f) + 0.5f;
  colour *= (float3)(0.5f, 0.5f, 0.5f);
  dist = sdPlane(origin, normalize((float3)(0.0f, -1.0f, 0.0f)), 10.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  dist = sdPlane(origin, normalize((float3)(0.0f, 1.0f, 0.0f)), 1.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
}

void Panels ( int avoid, float3 origin, SampledPt* pt, float3 dval, float time){
  float dist = sdSphere(origin-(float3)(1.17f, 0.7f, 0.25f), 0.3f);
  MapUnionG(avoid, pt, dist, 1, (float3)(-1.0f));
}
void BSDFModel ( int avoid, float3 O, SampledPt* pt, float3 dval, float time){
  O += (float3)(0.3f);
  float id;
  O.x = opMod1(O.x, 0.5f, &id);
  float dist = sdBox(O, (float3)(0.1f))*0.5f;
  id *= 2349234.2349234f;
  MapUnionG(avoid, pt, dist, 2, noise2to3((float2)(id, noise1D(id))));
}

void Mirror ( int avoid, float3 O, SampledPt* pt, float3 dval, float time){
  O += (float3)(-1.825f, 0.504f, -0.467f);
  O.xz = opRotate(O.xz, -10.25f);
  float dist = sdBox(O, (float3)(0.78f, 0.75f, 0.002f));
  MapUnionG(avoid, pt, dist, 3, (float3)(-1.0f));
}

void Globe ( int avoid, float3 O, SampledPt* pt, float3 dval, float time){
  O -= (float3)(-0.99f,  0.0f, -4.83f);
  O.yz = opRotate(O.yz, 4.65f);
  O.xy = opRotate(O.xy, time*0.3f);
  float dist = sdAodMitsuba(O);
  dist = Shell(dist, 0.001f);
  MapUnionG(avoid, pt, dist, 4, (float3)(-1.0f));
}

void Update_Map ( int avoid, float3 origin, SampledPt* pt, float time,
                  __read_only image2d_array_t textures, float3 dval ) {
  Room(avoid, origin, pt, time, dval);
  Panels(avoid, origin, pt, dval, time);
  BSDFModel(avoid, origin, pt, dval, time);
  Mirror   (avoid, origin, pt, dval, time);
  Globe    (avoid, origin, pt, dval, time);
}
UPDATEMAPEND