TEXTURESSTART
TEXTURESEND
MATERIALSSTART
{
  
  "materials": [
       { "_c": "room / 0",
      "albedo":       "[0.7, 0.6, 0.4]",
      "diffuse":      "1.0", "specular":     "0.0",
      "glossy":       "0.0", "glossy_lobe":  "0.95",
      "transmittive": "0.0", "ior": "0.5",
      "roughness":    "0.6", "metallic":     "0.2",
      "fresnel":      "0.9", "subsurface":   "0.0",
      "anisotropic":  "0.8"
    }, { "_c": "bat / 1",
      "albedo":       "[0.7, 0.6, 0.4]",
      "diffuse":      "0.7", "specular":     "0.0",
      "glossy":       "0.3", "glossy_lobe":  "0.15",
      "transmittive": "0.0", "ior": "0.5",
      "roughness":    "0.4", "metallic":     "1.0",
      "fresnel":      "1.9", "subsurface":   "0.0",
      "anisotropic":  "1.0"
    }, { "_c": "box / 2",
      "albedo":       "[0.8, 0.2, 0.8]",
      "diffuse":      "1.0", "specular":     "0.0",
      "glossy":       "0.0", "glossy_lobe":  "0.95",
      "transmittive": "0.0", "ior": "0.5",
      "roughness":    "0.8", "metallic":     "0.4",
      "fresnel":      "0.5", "subsurface":   "1.0",
      "anisotropic":  "0.8"
    }, { "_c": "tra / 3",
      "albedo":       "[1.0, 1.0, 1.0]",
      "diffuse":      "0.0", "specular":     "0.0",
      "glossy":       "0.0", "glossy_lobe":  "0.95",
      "transmittive": "1.0", "ior": "0.131",
      "roughness":    "0.8", "metallic":     "0.4",
      "fresnel":      "0.5", "subsurface":   "1.0",
      "anisotropic":  "0.8"
    }
  ]
}
MATERIALSEND

CAMERASTART
void Update_Camera ( Camera* camera, float time ) {
  // if ( camera->flags > 0 ) return;
  camera->position = (float3)(-3.742067f, 2.68f, 6.910531f);
  camera->lookat = (float3)(0.72f, 0.61f, 0.0f);
  // camera->fov = 110.0f;
  camera->focal = 1.025f + 0.5f*sin(time);
  camera->radius = 0.00f;

  // camera->position = (float3)(22.467f, 1.815f, -2.7065f);
  // camera->lookat = (float3)(0.71f, 0.8f, 0.0f);
  // camera->lookat   = (float3)(cos(time*0.5f), sin(time*0.5f), cos(time)*0.1f);
}
CAMERAEND

EMITTERSTART
__constant int EMITTER_AMT = 2;
Emitter REmission ( int index, float3 dval, float time ) {
  float f = (float)(index+1);
  float3 origin;
  origin = (float3)(4.8f, 1.15f, 8.02f);
  if ( index == 0 ) // most light
    return (Emitter){origin, (float3)(0.5f)*13.47f, 0.756f, index};
  origin = -(float3)(13.47f, -15.288f, -9.859f);
  if ( index == 1 ) // caustic light
    return (Emitter){origin, (float3)(0.5f, 0.5f, 0.7f)*4.327f, 1.5f, index};
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
  origin.xy = opRotate(origin.xy, 0.923f);
  dist = sdPlane(origin, normalize((float3)(0.0f, 0.0f, 1.0f)), 3.5f);
  colour = (float3)(0.4f);
  MapUnionG(avoid, pt, dist, 0, colour);
  dist = sdPlane(origin, normalize((float3)(0.0f, 0.0f, -1.0f)), 12.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  dist = sdPlane(origin, normalize((float3)(-1.0f, 0.0f, 0.0f)), 14.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  dist = sdPlane(origin, normalize((float3)(1.0f, 0.0f, 0.0f)), 16.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  dist = sdPlane(origin, normalize((float3)(0.0f, -1.0f, 0.0f)), 20.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
  dist = sdPlane(origin, normalize((float3)(0.0f, 1.0f, 0.0f)), 11.5f);
  MapUnionG(avoid, pt, dist, 0, colour);
}

float sdTriPrism ( float3 p, float2 h ) {
  float3 q = fabs(p);
  float d1 = q.z-h.y;
  float d2 = fmax(q.x*0.866025f+p.y*0.5f,-p.y)-h.x*0.5f;
  return length(fmax((float2)(d1,d2),0.0f)) + fmin(fmax(d1,d2), 0.0f);
}


float Batman ( float3 p ) {
  float3 g = p;


  { // reflect space
    float3 pnormal = (float3)(0.0f, 0.0f, 1.0f);
    float t = dot(p, pnormal) + 2.82f;
    if ( t < 0.0f ) p = p - (2.0f*t)*pnormal;
  }
  opRotate(p.xy, 270.0f*PI/180.0f);

  float a, b, c, d;

  // ----- wing ----
  a = sdBox(p-(float3)( 1.68f, 0.0f,-1.65f),
              (float3)(1.7f, 0.1f, 1.1f));
  float h = p.x;
  // -- left half slice
  float3 q = p;
  q.x -= sin(p.y)*0.25f;
  q.z += sin(p.x)*0.45f;
  b = sdCylinder(q-(float3)( 1.2f, 0.0f,-0.1f), 1.2f, 0.2f);
  // -- bottom half slice
  q = p;
  q.x += sin(p.x)*0.85f;
  c = sdCylinder(q-(float3)(3.2f,0.0f,-1.7f), 1.2f, 0.2f);
  float body = fmax(fmax(a, -b), -c);

  q = p;
  q.x += sin(p.x)*0.3f;
  // q.x += cos(p.x)*0.1f;
  a = sdCylinder(q-(float3)( 0.2f, 0.0f,-2.5f), 0.3f, 0.2f);
  body = fmax(body, -a);

  // ------ head ------
  q = p;
  q -= (float3)(0.3f, 0.0f, -2.8f);
  opRotate(q.xy, 90.0f*PI/180.0f);
  opRotate(q.xz, 90.0f*PI/180.0f);
  opRotate(q.yx, 0.4f);
  a = sdTriPrism(q, (float2)(0.3f, 0.1f));
  body = fmin(body, a);
  // body = a;


  // --- bandings ---
  // ring
  q = p;
  a = sdTorus(q-(float3)(0.7f,0.0f,-2.9f),0.05f, 0.8f);
  body = fmax(fmax(body, -a), (body-a+0.1f)*sqrt(0.5f));
  // ridges
  q = p;
  q -= (float3)(-0.2f, 0.0f, -2.0f);
  float siz = 1.3f;
  q.z = fmod(q.z + 0.5f*siz, siz) - 0.5f*siz;
  // a = sdSphere(q, 0.3f);
  // body = fmax(fmax(body, -a), (body-a+0.1f)*sqrt(0.5f));


  return body*0.8f;
}

void Waves (int avoid, float3 O, SampledPt* pt, float time, float3 dval) {
  O += (float3)(6.105f, -4.041f, -6.865f);
  O.xy = opRotate(O.xy, -2.141f);
  O.yz = opRotate(O.yz, time);
  float d = Batman(O+(float3)(0.0f, 0.0f, -2.822f));
  MapUnionG(avoid, pt, d, 1, (float3)(-1.0f));
}

void Box ( int avoid, float3 O, SampledPt* pt, float time, float3 dval ) {
  O += (float3)(8.0f, -5.523f, -9.059f);
  O.xy = opRotate(O.xy, time);
  O.yz = opRotate(O.yz, time);
  float3 h = (float3)(0.5f);
  h.x += sin(O.y*(sin(time*2.5f)*0.17f))*7.0f;
  h.y += sin(O.z*45.0f + O.x*45.0f)*0.008f;
  h.z += cos(O.x*0.866f + O.y*1.988f)*0.562f;
  float d = sdBox(O, h)*0.5f;
  MapUnionG(avoid, pt, d, 2, (float3)(-1.0f));
}

void Trans(int avoid, float3 O, SampledPt* pt, float time, float3 dval ) {
  O += (float3)(4.665f, -2.466f, -7.023f);

  O.xy = opRotate(O.xy, time);
  O.yz = opRotate(O.yz, time);
  float d = sdCone(O, 0.3f, 1.0f);
  d = Shell(d, 0.0001f);
  MapUnionG(avoid, pt, d, 3, (float3)(-1.0f));
}

void Update_Map ( int avoid, float3 origin, SampledPt* pt, float time,
                  __read_only image2d_array_t textures, float3 dval ) {
  Room (avoid, origin, pt, time, dval);
  Waves(avoid, origin, pt, time, dval);
  Box  (avoid, origin, pt, time, dval);
  Trans(avoid, origin, pt, time, dval);
}
UPDATEMAPEND
