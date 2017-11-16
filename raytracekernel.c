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
  // pdf
  float diffuse, specular, glossy, glossy_lobe;
  float transmittive, ior;
  // PBR material
  float pbr_roughness, pbr_metallic, pbr_fresnel;
  // disney material
  float d_subsurface, d_specular_tint, d_anisotropic, d_sheen, d_sheen_tint,
    d_clearcoat, d_clearcoat_gloss, d_specular;
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

__constant float PI  = 3.141592654f;
__constant float TAU = 6.283185307f;
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
  float m = clamp(1.0f - u, 0.0f, 1.0f);
  float m2 = m*m;
  return m2*m2*m;
}
float GTR1 ( float cosN_H, float a ) {
  if ( a >= 1.0f ) return 1.0f/PI;
  float a2 = a*a;
  float t = 1.0f + (a2-1.0f) * sqr(cosN_H);
  return (a2 - 1.0f)/(PI*log(a2)*t);
}
float GTR2 ( float cosN_H, float a ) {
  float a2 = a*a;
  float t = 1.0f + (a2 - 1.0f) * sqr(cosN_H);
  return a2/(PI*sqr(t));
}
float GTR2_Aniso (float cosH_N, float cosH_X, float cosH_Y, float ax, float ay){
  return 1.0f/(PI*ax*ay* sqr(sqr(cosH_X/ax) + sqr(cosH_Y/ay) + sqr(cosH_N)));
}
float Smith_G_GGX ( float cosN_V, float alpha_g ) {
  float a = sqr(alpha_g),
        b = sqr(cosN_V);
  return 1.0f/(cosN_V + sqrt(a + b - a*b));
}
float Smith_G_GGX_Aniso(float cosV_N, float cosV_X, float cosV_Y,
                        float ax, float ay) {
  return 1.0/(cosV_N + sqrt(sqr(cosV_X*ax) + sqr(cosV_Y*ay) + sqr(cosV_N)));
}
float3 Mon_To_Lin ( float3 x ) { return pow(x, (float3)(2.2)); }

float3 BRDF_F ( float3 wi, float3 wo, float3 N, Material* m, float3 pcolour ) {
  // get proper L, V, binormal (X) and bitangent (Y)
  float3 L = wo, V = -wi, X = Binormal(N), Y = cross(X, N);
  // get half vector (angle between L and V)
  float3 H = normalize(wi + wo);
  // useful normal angles
  float cosN_L = dot(N, L), cosN_V = dot(N, V),
        cosH_L = dot(H, L), cosH_N = dot(H, N);
  // check that wi -> wo is within hemisphere/not degenerate
  if ( cosN_L < 0.0f || cosN_V < 0.0f ) return (float3)(0.0f);

  // luminance approximation
  float3 albedo = Mon_To_Lin(pcolour);
  float3 lumination = albedo * (float3)(0.3f, 0.6f, 0.1f);

  // normalize lumination to isolate hue and saturation; get specular/sheen
  float3 tint = lumination > 0.0f ? albedo/lumination : (float3)(1.0f);
  float3 specular = mix(m->d_specular * 0.08f *
                        mix((float3)(1.0f), tint, (float3)(m->d_specular_tint)),
                        albedo, m->pbr_metallic);
  float3 sheen = mix((float3)(1.0f), tint, m->d_sheen_tint);

  // Diffuse Fresnel - go from 1 at normal incidence to 0.5 at grazing
  // and mix in diffuse retroreflection based on roughness
  // [retroreflection s when light reflects towards the light source]
  float F_L = Schlick_Fresnel(cosN_L),
        F_V = Schlick_Fresnel(cosN_V),
        F_diffuse = 0.5f + 2.0f * sqr(cosH_L) * m->pbr_roughness,
        F = mix(1.0f, F_diffuse, F_L) *
                  mix(1.0f, F_diffuse, F_V);

  // Subsurface-Scattering Approximation. Based on Hanrahan-Kruger BRDF
  // brdf approximation of isotropic BSSRDF. 1.25 scale is used to preserve
  // albedo. Fss90 is used to flatten retroreflection
  float F_ss90 = sqr(cosH_L)*m->pbr_roughness,
        F_ss   = mix(1.0f, F_ss90, F_L) *
                               mix(1.0f, F_ss90, F_V);
  float subsurface = 1.25f * (F_ss *
                              (1.0f/(cosN_L + cosN_V) - 0.5f) + 0.5f);

  // specular
  float aspect    = sqrt(1.0f - m->d_anisotropic * 0.9f),
        ax        = fmax(0.001f, sqr(m->pbr_roughness)/aspect),
        ay        = fmax(0.001f, sqr(m->pbr_roughness)*aspect),
        D_specular        = GTR2_Aniso(cosH_N, dot(H, X), dot(H, Y), ax, ay),
        F_H = Schlick_Fresnel(cosH_L);
  float3 F_specular = mix(specular, (float3)(1.0f), F_H);
  float G_specular = Smith_G_GGX_Aniso(cosN_L, dot(L, X), dot(L, Y), ax, ay) *
             Smith_G_GGX_Aniso(cosN_V, dot(V, X), dot(V, Y), ax, ay);

  // sheen
  float3 F_sheen = F_H * m->d_sheen * sheen;

  // clearcoat (ior = 1.5 -> f0 = 0.04)
  float D_clearcoat = GTR1(cosH_N, mix(0.1f, 0.001f, m->d_clearcoat_gloss)),
        F_clearcoat = mix(0.04f, 1.0f, F_H),
        G_clearcoat = Smith_G_GGX(cosN_L, 0.25f) * Smith_G_GGX(cosN_V, 0.25f);

  // final combinations
  float3 final_subsurface = mix(F, subsurface, m->d_subsurface) * albedo,
         final_albedo = (1.0f/PI) * final_subsurface + F_sheen,
         final_metallic = (1.0f - m->pbr_metallic),
         final_specular = G_specular * F_specular * D_specular,
         final_clearcoat = 0.25f * m->d_clearcoat * D_clearcoat *
                                    F_clearcoat * G_clearcoat;

  return final_albedo * final_metallic + final_specular + final_clearcoat;
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
    float3 brdf = BRDF_F(V, L, N, m, pcolour);
    float shadow = Visibility_Ray(O, emit.origin, emit.radius, si, Tx);

    colour += (brdf) * dot(N, L) * normalize(emit.emission);
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
