TEXTURESSTART
TEXTURESEND
MATERIALSSTART
{
  
  "materials": [
       { "_c": "room 0",
      "albedo":       "[0.7, 0.6, 0.4]",
      "diffuse":      "1.0", "specular":     "0.0",
      "glossy":       "0.0", "glossy_lobe":  "0.95",
      "transmittive": "0.0", "ior": "0.5",
      "roughness":    "0.6", "metallic":     "0.2",
      "fresnel":      "0.9", "subsurface":   "0.0",
      "anisotropic":  "0.8"
    }, { "_c": "walkway 1",
      "albedo":       "[0.1, 0.1, 0.1]",
      "diffuse":      "0.75", "specular":     "0.0",
      "glossy":       "0.25", "glossy_lobe":  "0.9",
      "transmittive": "0.0", "ior": "0.284",
      "roughness":    "0.215", "metallic":     "0.268",
      "fresnel":      "0.5", "subsurface":   "0.364",
      "anisotropic":  "0.478"
    }, { "_c": "organism 2",
      "albedo":       "[0.67, 0.56, 0.79]",
      "diffuse":      "1.0", "specular":     "0.0",
      "glossy":       "0.0", "glossy_lobe":  "0.9",
      "transmittive": "0.0", "ior": "0.284",
      "roughness":    "0.127", "metallic":     "0.0",
      "fresnel":      "0.8", "subsurface":   "0.8",
      "anisotropic":  "0.7"
    }, { "_c": "interior 3",
      "albedo":       "[0.37, 0.28, 0.27]",
      "diffuse":      "1.0", "specular":     "0.0",
      "glossy":       "0.0", "glossy_lobe":  "0.9",
      "transmittive": "0.0", "ior": "0.284",
      "roughness":    "1.0", "metallic":     "0.5",
      "fresnel":      "0.5", "subsurface":   "0.0",
      "anisotropic":  "0.7"
    }, { "_c": "interior coat 4",
      "albedo":       "[0.05, 0.1, 0.2]",
      "diffuse":      "0.5", "specular":     "0.3",
      "glossy":       "0.2", "glossy_lobe":  "0.5",
      "transmittive": "0.0", "ior": "0.284",
      "roughness":    "0.8", "metallic":     "0.2",
      "fresnel":      "0.8", "subsurface":   "0.1",
      "anisotropic":  "0.3"
    }, { "_c": "pyramid extend 5",
      "albedo":       "[1.00, 0.81, 0.23]",
      "diffuse":      "0.2", "specular":     "0.0",
      "glossy":       "0.8", "glossy_lobe":  "0.9",
      "transmittive": "0.0", "ior": "0.284",
      "roughness":    "0.8", "metallic":     "0.2",
      "fresnel":      "0.8", "subsurface":   "0.1",
      "anisotropic":  "0.3"
    }, { "_c": "twisty thing 6",
      "albedo":       "[0.4, 0.36, 0.36]",
      "diffuse":      "1.0", "specular":     "0.0",
      "glossy":       "0.0", "glossy_lobe":  "0.9",
      "transmittive": "0.0", "ior": "0.284",
      "roughness":    "1.0", "metallic":     "0.6",
      "fresnel":      "0.2", "subsurface":   "0.2",
      "anisotropic":  "1.0"
    }, { "_c": "hour glass interior 7",
      "albedo":       "[1.0, 1.0, 1.0]",
      "diffuse":      " 0.0", "specular":     "0.0",
      "glossy":       "0.0", "glossy_lobe":  "0.9",
      "transmittive": "1.0", "ior": "0.07",
      "roughness":    "0.8", "metallic":     "0.2",
      "fresnel":      "0.8", "subsurface":   "0.1",
      "anisotropic":  "0.3"
    }, { "_c": "hour glass exterior 8",
      "albedo":       "[0.28, 0.06, 0.06]",
      "diffuse":      "1.0", "specular":     "0.0",
      "glossy":       "0.0", "glossy_lobe":  "0.9",
      "transmittive": "0.0", "ior": "0.07",
      "roughness":    "0.1", "metallic":     "0.2",
      "fresnel":      "0.8", "subsurface":   "0.1",
      "anisotropic":  "0.3"
    }, { "_c": "inner globe",
      "albedo":       "[0.13, 0.96, 0.73]",
      "diffuse":      "0.0", "specular":     "0.0",
      "glossy":       "0.0", "glossy_lobe":  "0.0",
      "transmittive": "1.0", "ior": "0.07",
      "roughness":    "0.8", "metallic":     "0.2",
      "fresnel":      "0.8", "subsurface":   "0.9",
      "anisotropic":  "0.8"
    }
  ]
}
MATERIALSEND

CAMERASTART
void Update_Camera ( Camera* camera, float time ) {
  // if ( camera->flags > 0 ) return;

  // camera->position = (float3)(0.0f, 0.37f, 4.24f);
  // camera->lookat = (float3)(0.0f, 0.51f, 0.0f);
  // camera->position = (float3)(0.0f);
  // camera->lookat   = (float3)(cos(time*0.5f), sin(time*0.5f), cos(time)*0.1f);
}
CAMERAEND

EMITTERSTART
__constant int EMITTER_AMT = 3;
Emitter REmission ( int index, float3 dval, float time ) {
  float f = (float)(index+1);
  float3 origin;
  origin = (float3)(0.0f, 18.0f, 2.3f);
  if ( index == 0 ) // most light
    return (Emitter){origin, (float3)(0.5f)*14.0f, 1.5f, index};
  origin = (float3)(-1.0f, 4.6f, 10.8f);
  if ( index == 1 ) // caustic light
    return (Emitter){origin, (float3)(0.5f, 0.5f, 0.7f)*11.0f, 0.5f, index};
  origin = (float3)(-0.5f, -0.6f, 1.3f);
  if ( index == 2 ) // glossy light
    return (Emitter){origin, (float3)(0.3f, 0.7f, 0.2f)*0.5f, 0.3f, index};
}
EMITTEREND

UPDATEMAPSTART
#define RCOL (float3)(0.2f, 0.2f, 0.2f)
void Room ( int avoid, float3 origin, SampledPt* pt, float time, float3 dval ) {
  float dist;
  float3 colour = (float3)(1.0f);
  dist = sdPlane(origin, normalize((float3)(0.0f, 0.0f, 1.0f)), 2.5f);
  colour = (float3)(sin(origin.x*PI*5.0f)*sin(origin.y*PI*5.0f) >= 0.0f) + 0.5f;
  colour *= RCOL;
  MapUnionG(avoid, pt, dist, 0, colour);
  dist = sdPlane(origin, normalize((float3)(0.0f, 0.0f, -1.0f)), 12.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  colour = (float3)(sin(origin.z*PI*5.0f)*sin(origin.y*PI*5.0f) >= 0.0f) + 0.5f;
  colour *= RCOL;
  dist = sdPlane(origin, normalize((float3)(-1.0f, 0.0f, 0.0f)), 2.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  dist = sdPlane(origin, normalize((float3)(1.0f, 0.0f, 0.0f)), 2.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  colour = (float3)(sin(origin.x*PI*5.0f)*sin(origin.z*PI*5.0f) >= 0.0f) + 0.5f;
  colour *= RCOL;
  dist = sdPlane(origin, normalize((float3)(0.0f, -1.0f, 0.0f)), 20.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  dist = sdPlane(origin, normalize((float3)(0.0f, 1.0f, 0.0f)), 1.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
}

void Globe ( int avoid, float3 O, SampledPt* pt, float3 dval, float time){
  float dist;
  float3 q = O;

  {//room
    O.xy = opMirrorOctant(O.xy, (float2)(0.4f), 0);
    O.yz = opMirrorOctant(O.yz, (float2)(0.4f), 0);
    // 0.0 -> 0.312
    float pyr = sdBox(O-(float3)(1.0f), (float3)(0.312f, 0.312f, 1.309f));
    dist = sdTorus(O-(float3)(1.0f), 0.115f, 0.874f);
    dist = fmin(pyr, dist);
    O = q;
    float3 a, b;
    a = (float3)(-1.717f, 0.697f, 1.702f);
    b = (float3)( 1.876f, 0.697f, 1.702f);
    a.y += sin(O.x*12.8f)*0.088f;
    b.y += sin(O.x*12.8f)*0.088f;
    a.z += cos(O.x*52.38f)*0.018f;
    b.z += cos(O.x*52.38f)*0.018f;
    pyr = sdCapsuleEnd(O, a, b, 0.029f);
    dist = opSubtractChamfer(dist, pyr, 0.022f)*0.5f;
  }
  MapUnionG(avoid, pt, dist, 3, (float3)(-1.0f));

  {//room armour
    O=q;
    O.xy = opMirrorOctant(O.xy, (float2)(0.4f), 0);
    O.yz = opMirrorOctant(O.yz, (float2)(0.4f), 0);
    float pyr = sdBox(O-(float3)(1.0f), (float3)(0.35f, 0.35f, 1.121f));
    dist = pyr;
  }
  MapUnionG(avoid, pt, dist, 4, (float3)(-1.0f));

  {
    O = q;
    float sp = sdSphere(O-(float3)(0.0f, 2.01f, 1.45f), 0.325f);
    dist = sdCapsuleEnd(O, (float3)(-1.3f, 1.74f, 1.53f),
                           (float3)(0.01f, 2.01f, 1.45f)+
                           (float3)(0.0f, 1.0f, 0.0f)*sin(O.y*1.513f)*0.1f*-4.0f,
                           0.072f);
    dist = opUnionRound(dist, sp, 0.144f);
    float3 a = (float3)(-1.44f, 0.16f, 1.75f);
    float3 b = (float3)(-1.3f, 1.74f, 1.53f);

    float3 ab = b - a;
    ab.z += sin(O.y*-3.843f)*0.345f;
    ab.x += cos(O.y*-10.695f)*0.09f;
    float t = clamp(dot(O - a, ab)/dot(ab, ab), 0.0f, 1.0f);
    sp = length((ab*t + a) - O) - 0.07f;
    dist = opUnionRound(dist, sp, 0.194f);

    a = (float3)(0.2f, 1.9f, 1.6f);
    b = (float3)(0.607f, 1.139f, 1.745f);
    a.z += sin(O.y*0.755f )*-0.098f;
    b.z += sin(O.y*11.978f)*-0.098f;
    sp = sdCapsuleEnd(O, a, b, 0.07f);
    dist = opUnionRound(dist, sp, 0.194f);
    a = b;
    b = (float3)(0.517f, 0.33f, 1.692f);
    sp = sdCapsuleEnd(O, a, b, 0.04f);
    dist = opUnionRound(dist, sp, 0.194f);
    dist *= 0.5f;
  }
  MapUnionG(avoid, pt, dist, 2, (float3)(-1.0f));

  {//walkway
    O = q;
    dist = sdBox(O-(float3)(0.00f, -0.16f, 6.05f),
                   (float3)(0.806f, 0.134f+sin(O.z*80.0f)*sin(O.x*80.0f)*0.01f,
                            4.587f));
    O -= (float3)(-20.35f, 0.296f, -27.11f);
    O.xz = opRotate(O.xz, PI/4.0f);
    O.xz = opMod2(O.xz, (float2)(0.15f), 0);
    float hol = sdBox(O, (float3)(0.05f));
    dist = opSubtractChamfer(dist, hol, 0.02f)*0.5f;
  }
  MapUnionG(avoid, pt, dist, 1, (float3)(-1.0f));
  {//cone meta ball holder
    O = q;
    O.y = opMirror(O.y, 1.417f, 0);
    O.xy = opRotate(O.xy, PI);
    dist = sdCone(O, 0.6f+sin(O.y)*0.2f, 0.78f);
  }
  MapUnionG(avoid, pt, dist, 5, (float3)(-1.0f));
  {//twisty thing
    O = q;
    O.zy = opMirrorOctant(O.zy, (float2)(0.828f, 0.843f), 0);
    O.xy = opRotate(O.xy, PI*0.5f);
    O.yz = opRotate(O.yz, 35.25f+time);
    O.x = opMod1(O.x, 0.825f, 0);
    float ring = sdTorus2(O, 0.070f, 0.807f+sin(time*2.0f)*0.1f);
    O = q;
    float ringdel = sdBox(O, (float3)(0.8f, 0.8f, 0.8f));
    ring = max(ring, -ringdel);
    dist = min(dist, ring);
  }
  MapUnionG(avoid, pt, dist, 6, (float3)(-1.0f));
  {//hour glass
    O = q;
    O.yz = opRotate(O.yz, 2.404f);
    O.xy = opRotate(O.xy, time*0.5f);
    float3 i = fabs(O);
    float2 h = (float2)(0.25f, 0.633f);
    float sp = fmax(i.y - h.y, fmax(i.x*sqrt(1.0f)*0.5f + i.z*0.5f, i.z) - h.x);
    float bnd = sdDisc(O, 2.237f) - -0.455f;
    dist = max(sp, -bnd);
    dist = opSubtractRound(sp, bnd, 0.941f)*0.5f;
    dist = Shell(dist, 0.0001f);
  }
  MapUnionG(avoid, pt, dist, 7, (float3)(-1.0f));
  {
    O.y = opMirror(O.y, 0.656f, 0);
    O.yz = opRotate(O.yz, 2.404f);
    O.xy = opRotate(O.xy, time*0.5f);
    O.xy = opRotate(O.xy, -0.002f + time*-0.5f);
    O.yz = opRotate(O.yz, 0.576f);


    float3 i = fabs(O);
    float2 h = (float2)(0.275f, 0.045f);
    float sp = fmax(i.y - h.y, fmax(i.x*sqrt(1.0f)*0.5f + i.z*0.5f, i.z) - h.x);
    dist = sp;
  }
  MapUnionG(avoid, pt, dist, 8, (float3)(-1.0f));
  float3 col = (float3)(0.0f);
  {//sand thing
    O = q;
    dist = sdSphere(O+(float3)(-0.605f, -4.493f, -11.015f), dval.x);
    dist = Shell(dist, 0.0001f);
  }
  MapUnionG(avoid, pt, dist, 9, (float3)(-1.0f));
}

void Update_Map ( int avoid, float3 origin, SampledPt* pt, float time,
                  __read_only image2d_array_t textures, float3 dval ) {
  Room(avoid, origin, pt, time, dval);
  Globe    (avoid, origin, pt, dval, time);
}
UPDATEMAPEND