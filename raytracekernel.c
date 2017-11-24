#define MAX_DEPTH 16
#define MARCH_DIST //%MARCH_DIST.0f
#define MARCH_REPS //%MARCH_REPS
#define MAT_LEN //%MAT_LENGTH

#define TEXTURE_T __read_only image2d_array_t
#define SCENE_T(S, T) SceneInfo* S , TEXTURE_T T

#define DFLT2(V) ((V).x), ((V).y)
#define DFLT3(V) ((V).x), ((V).y), ((V).z)
#define DFLT4(V) ((V).x), ((V).y), ((V).z), ((V).w)

#define writeln(X)     if (Is_Debug()){printf(X "\n");                         }
#define writeint(X)    if (Is_Debug()){printf(#X " %d\n",                 (X));}
#define writefloat(X)  if (Is_Debug()){printf(#X " %f\n",                 (X));}
#define writefloat2(X) if (Is_Debug()){printf(#X " %f, %f\n",        DFLT2(X));}
#define writefloat3(X) if (Is_Debug()){printf(#X " %f, %f, %f\n",    DFLT3(X));}
#define writefloat4(X) if (Is_Debug()){printf(#X " %f, %f, %f, %f\n",DFLT4(X));}

__constant float MARCH_ACC = //%MARCH_ACC.0f/1000.0f;
// -----------------------------------------------------------------------------
// --------------- DEBUG -------------------------------------------------------
// Variadic functions not supported, so this is best you get :-(
bool Is_Debug ( ) {
  return get_global_id(0) == get_global_size(0)/2 &&
         get_global_id(1) == get_global_size(1)/2;
}
// -----------------------------------------------------------------------------
// --------------- GPU-CPU STRUCTS ---------------------------------------------
typedef struct T_Camera {
  float3 position, lookat, up;
  int2 dim;
  float fov;
  int flags;
} Camera;

typedef struct T_Material {
  // colour [set to (-1.0, -1.0, -1.0) to have map override it]
  float3 albedo;
  // sampling strategy
  float diffuse, specular, glossy, glossy_lobe;
  float transmittive;
  // PBR material
  float roughness, metallic, fresnel, subsurface, anisotropic;
} Material;

typedef struct T_Emitter {
  float3 origin, emission;
  float radius;
} Emitter;

typedef struct T_SceneInfo {
  float time;
  // __read_only image2d_array_t textures; IMAGES CANT BE USED AS FIELD TYPES:-(
  Material* materials;
  float3 debug_values;
} SceneInfo;
SceneInfo New_SceneInfo(float time, Material* materials,
                        float3 debug_values){
  SceneInfo si;
  si.time         = time;
  si.materials    = materials;
  si.debug_values = debug_values;
  return si;
}

__constant float PI   = 3.141592653589793f;
__constant float IPI  = 0.318309886183791f;
__constant float IPI2 = 0.159154943091895f;
__constant float TAU  = 6.283185307179586f;
__constant float ITAU = 0.159154943091895f;
// -----------------------------------------------------------------------------
// --------------- GENERAL STRUCTS ---------------------------------------------
typedef struct T_Ray {
  float3 origin, dir;
} Ray;

typedef struct T_SampledPt {
  float3 colour, origin, dir;
  float dist;
  int material_index;
} SampledPt;
// -----------------------------------------------------------------------------
// --------------- GENERAL UTILITIES      --------------------------------------
float Distance(float3 u, float3 v) {
  float x = u.x-v.x, y = u.y-v.y, z = u.z-v.z;
  return sqrt(x*x + y*y + z*z);
}
float sqr(float t) { return t*t; }
// -----------------------------------------------------------------------------
// --------------- MAP GEOMETRY FUNCTIONS --------------------------------------

//---MAP GEOMETRY INSERTION POINT---
//%MAPFUNCDECLARATIONS
//----------------------------------
//%MAPFUNCDEFINITIONS
//----------------------------------

void MapUnionG( int avoid, SampledPt* d1, float d, int mi, float3 c ) {
  if ( mi != avoid && d1->dist > d ) {
    d1->colour = c;
    d1->dist = d;
    d1->material_index = mi;
  }
}

// -----------------------------------------------------------------------------
// --------------- SCENE -------------------------------------------------------
//%SCENEINSERT
//------------------------
// -----------------------------------------------------------------------------
// --------------- MAP ---------------------------------------------------------
SampledPt Map ( int a, float3 origin, SCENE_T(si, Tx) ) {
  SampledPt res;
  res.dist = FLT_MAX;

  //---MAP INSERTION POINT---
  //%MAPINSERT
  //-------------------------

  // approx lighting with emissions
  for ( int i = 0; i != EMITTER_AMT; ++ i ) {
    Emitter e = REmission(i, si->debug_values, si->time);
    float dist = sdSphere(origin - e.origin, e.radius);
    MapUnionG(a, &res, dist, 0, (float3)(1.0f, 0.9f, 0.8f)*e.emission);
  }

  return res;
}

// -----------------------------------------------------------------------------
// --------------- GRAPHIC FUNCS -----------------------------------------------
float3 Normal ( float3 p, SCENE_T(si, Tx) ) {
  float2 e = (float2)(1.0f, -1.0f)*0.5773f*0.0005f;
  return normalize(
    e.xyy*Map(-1, p + e.xyy, si, Tx).dist +
    e.yyx*Map(-1, p + e.yyx, si, Tx).dist +
    e.yxy*Map(-1, p + e.yxy, si, Tx).dist +
    e.xxx*Map(-1, p + e.xxx, si, Tx).dist);
}

float3 reflect ( float3 V, float3 N ) {
  return V - 2.0f*dot(V, N)*N;
}

float3 refract(float3 V, float3 N, float refraction) {
  float cosI = -dot(N, V);
  float cosT = 1.0f - refraction*refraction*(1.0f - cosI*cosI);
  return (refraction*V) + (refraction*cosI - sqrt(cosT))*N;
}

// -----------------------------------------------------------------------------
// --------------- RAYTRACING/MARCH --------------------------------------------
SampledPt March ( int avoid, Ray ray, SCENE_T(si, Tx) ) {
  float distance = 0.0f;
  SampledPt t_info;
  for ( int i = 0; i < MARCH_REPS; ++ i ) {
    t_info = Map(avoid, ray.origin + ray.dir*distance, si, Tx);
    if ( t_info.dist < MARCH_ACC || t_info.dist > MARCH_DIST ) break;
    distance += t_info.dist;
  }
  if ( t_info.dist > MARCH_DIST ) {
    t_info.dist = -1.0f;
    return t_info;
  }
  t_info.dist = distance;
  return t_info;
}

//----POSTPROCESS---
//%POSTPROCESS

float Visibility_Ray(float3 orig, float3 other, float sphere_radius,
                     SCENE_T(si, Tx)) {
  float theoretical = Distance(orig, other) - sphere_radius*2.0f;
  float3 dir = normalize(other - orig);
  orig += dir*(MARCH_ACC*2.0f); // get out of our own intersection
  SampledPt ptinfo = March(-1, (Ray){orig, dir}, si, Tx);
  float actual = ptinfo.dist + MARCH_ACC*1.2f; // distance + leeway
  return 1.0f*(actual >= theoretical);
}

// --- some brdf utility functions .. ---
float3 Binormal ( float3 N ) {
  float3 axis = (fabs(N.x) < 1.0f ? (float3)(1.0f, 0.0f, 0.0f) :
                                    (float3)(0.0f, 1.0f, 0.0f));
  return normalize(cross(N, axis));
}

float Schlick_Fresnel ( float u ) {
  float f = clamp(1.0f - u, 0.0f, 1.0f);
  float f2 = f*f;
  return f2*f2*f; // f^5
}

float Smith_G_GGX_Correlated ( float L, float R, float a ) {
  return L * sqrt(R - a*sqr(R) + a);
}

float3 BRDF_F ( float3 wi, float3 N, float3 wo, Material* m ) {
  // get binormal, bitangent, half vec etc
  const float3 binormal  = Binormal(N),
               bitangent = cross(binormal, N),
               L         =  wo, V = -wi,
               H         = normalize(L+V);
  const float  cos_NV    = dot(N, V), cos_NL     = dot(N, L),
               cos_HV    = dot(H, V), cos_HL     = dot(H, L),
               Fresnel_L = Schlick_Fresnel(cos_NL),
               Fresnel_V = Schlick_Fresnel(cos_NV);
  // Diffusive component, just made it up
  float3 diffusive_albedo = m->albedo * pow(cos_HL * cos_NV, 0.5f) * IPI;

  float3 microfacet = (float3)(1.0f);

  /* if ( cos_NL < 0.0f || cos_NV < 0.0f ) */
  /*   return (float3)(0.0f); */

  { // ------- Fresnel
    // modified diffusive fresnel from disney, modified to use albedo & F0
    const float F0 = m->fresnel * m->metallic,
                Fresnel_diffuse_90 = F0 * sqr(cos_HL);
    microfacet *= (1.0f - F0) * diffusive_albedo +
                    mix(1.0f, Fresnel_diffuse_90, Fresnel_L) *
                    mix(1.0f, Fresnel_diffuse_90, Fresnel_V);
  }
  { // ------- Geometric
    // Heits 2014, SmithGGXCorrelated with half vec combined with anisotropic
    // term using disney's GTR2_aniso model
    const float Param  = sqr(0.5f + sqr(m->roughness)),
                Aspect = sqrt(1.0f - m->anisotropic*0.9f),
                Ax     = Param/Aspect, Ay = Param*Aspect,
                GGX_NV = Smith_G_GGX_Correlated(cos_HL, cos_NV, Ax),
                GGX_HL = Smith_G_GGX_Correlated(cos_NV, cos_HL, Ay);
    microfacet *= 0.5f / (GGX_NV*Ax + GGX_HL*Ay);
  }
  { // ------- Distribution
    // Hyper-Cauchy Distribution using roughness and metallic
    const float Param = 1.0f + m->roughness,
                Shape = (1.1f - sqrt(m->metallic)),
                tan_HL = length(cross(H, L))/cos_HL;
    const float Upper  = (Param - 1.0f)*pow(sqrt(2.0f), (2.0f*Param - 2.0f)),
                LowerL = (PI*sqr(Shape) * pow(cos_HL, 4.0f)),
                LowerR = pow(2.0f + sqr(tan_HL)/sqr(Shape), Param);
    microfacet *= (Upper / (LowerL * LowerR));
  }

  // Since microfacet is described using half vec, the following energy
  // conservation model may be used [Edwards et al. 2006]
  microfacet /= 4.0f * cos_HV * fmax(cos_NL, cos_NV);

  { // --------- Subsurface
    // modified disney retro reflection based off Hanrahan-Grueger BSSRDF
    //   approximation
    const float Rr_term = (2.0f * (0.5f + m->roughness) * sqr(dot(N, H)));
    const float3 Retro_reflection = diffusive_albedo * Rr_term *
      (Fresnel_L + Fresnel_V + (Fresnel_L*Fresnel_V*(Rr_term - 1.0f)));
    diffusive_albedo = mix(diffusive_albedo, Retro_reflection, m->subsurface);
  }

  return (microfacet + diffusive_albedo);
}

float3 Illuminate ( float3 O, float3 Wi, float3 col, Material* m,
                    SCENE_T(si, Tx) ) {
  float3 colour = (float3)(0.0f);
  for ( int i = 0; i != EMITTER_AMT; ++ i ) {
    Emitter emit = REmission(i, si->debug_values, si->time);

    float3 V = normalize(Wi),
           N = Normal(O, si, Tx),
           L = normalize(emit.origin - O);
    float3 pcolour = m->albedo.x < 0.0f ? col : m->albedo;
    float3 brdf = BRDF_F(V, N, L, m);
    float shadow = Visibility_Ray(O, emit.origin, emit.radius, si, Tx);

    colour += (brdf) * dot(N, L) * (emit.emission);
  }
  return colour / EMITTER_AMT;
}

SampledPt Colour_Pixel(Ray ray, SCENE_T(si, Tx)) {
  SampledPt result = March(-1, ray, si, Tx);
  float3 origin = ray.origin + ray.dir*result.dist;
  Material* m = si->materials + result.material_index;

  result.colour = Illuminate(origin, ray.dir, result.colour, m, si, Tx);

  return result;
}

// -----------------------------------------------------------------------------
// --------------- CAMERA ------------------------------------------------------
Ray Camera_Ray(Camera* camera) {
  float2 coord = (float2)((float)get_global_id(0), (float)get_global_id(1));
  float2 resolution = (float2)((float)camera->dim.x, (float)camera->dim.y);
  resolution.y *= 16.0f/9.0f;

  float2 mouse_pos = camera->lookat.xy;

  float2 puv = -1.0f + 2.0f * (coord/resolution);

  float input_angle = PI - 2.0f*PI*mouse_pos.x;
  float3 cam_pos = camera->position;
  float3 cam_target = cam_pos + (float3)(sin(input_angle),
                                         (3.0f * mouse_pos.y) - 1.0f,
                                         cos(input_angle));
  float3 cam_front = normalize(cam_target - cam_pos);
  float3 cam_right = normalize ( cross(cam_front, (float3)(0.0f, 1.0f, 0.0f)));

  float3 cam_up = normalize(cross(cam_right, cam_front));
  float3 ray_dir = normalize(puv.x*cam_right + puv.y*cam_up +
                             (180.0f - camera->fov)*PI/180.0f*cam_front);

  Ray ray;
  ray.origin = cam_pos;
  ray.dir = ray_dir;
  return ray;
}

// -----------------------------------------------------------------------------
// --------------- RAYTRACE KERNEL ---------------------------------------------
__kernel void DTOADQ_Kernel (
    __global unsigned char* img, // R G B ITER
    __write_only image2d_t output_image,
    __global unsigned char*     clear_img,
    __global Camera* camera_ptr,
    __global float* time_ptr,
    __read_only image2d_array_t textures,
    __global Material* g_materials,
    __global float* debug_val,
    __global int*   rng
  ) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  Camera camera = *camera_ptr;
  float3 dval = (float3)(debug_val[0], debug_val[1], debug_val[2]);
  float time = *time_ptr;
  // -- get old pixel, check if there are samples to be done
  //    (counter is stored in alpha channel)
  Update_Camera(&camera, time);

  Material materials[MAT_LEN];
  for ( int i = 0; i != MAT_LEN; ++ i ) {
    materials[i] = *(g_materials + i);
  }

  SceneInfo scene_info = New_SceneInfo(time, materials, dval);

  Ray ray = Camera_Ray(&camera);
  SampledPt result = Colour_Pixel(ray, &scene_info, textures);

  float3 colour;
  if ( result.dist >= 0.0f ) {
    colour = result.colour;
  } else {
    colour = (float3)(0.0f);
  }
  // convert to Y CB R, from matlab
  write_imagef(output_image, out, (float4)(colour, 1.0f));
  colour = (float3)(
      16.0f + (  65.481f*colour.x + 128.553f*colour.y + 24.966f*colour.z),
      128.0f + (- 37.797f*colour.x - 74.2030f*colour.y + 112.00f*colour.z),
      128.0f + ( 112.000f*colour.x - 93.7860f*colour.y - 18.214f*colour.z)
  );
  int irwx = (camera.dim.y - out.y)*camera.dim.x + out.x;
  int irwy = irwx + camera.dim.x*camera.dim.y;
  int irwz = irwx + camera.dim.x*camera.dim.y*2;
  img[irwx] = (unsigned char)(colour.x);
  img[irwy] = (unsigned char)(colour.y);
  img[irwz] = (unsigned char)(colour.z);
}
