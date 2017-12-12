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

// math macros
#define SQR(T) ((T)*(T))
#define Spectrum float3

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
  float fov, focal, radius;
  int flags;
} Camera;


typedef struct T_Material {
  // colour [set to (-1.0, -1.0, -1.0) to have map override it]
  float3 albedo;
  // sampling strategy
  float diffuse, specular, glossy, glossy_lobe;
  float transmittive, ior;
  // PBR material
  float roughness, metallic, fresnel, subsurface, anisotropic;
} Material;

typedef struct T_Emitter {
  float3 origin, emission;
  float radius;
  int index;
} Emitter;

typedef struct T_SharedInfo {
  unsigned char clear_img;
  unsigned long finished_samples;
  unsigned char spp;
  uint2 rng_state;
} SharedInfo;

typedef struct T_SceneInfo {
  float time;
  // __read_only image2d_array_t textures; IMAGES CANT BE USED AS FIELD TYPES:-(
  __constant Material* materials;
  float3 debug_values;
  uint2 rng_state;
} SceneInfo;
SceneInfo New_SceneInfo(float time, __constant Material* materials,
                        float3 debug_values, uint2 rng){
  SceneInfo si;
  si.time         = time;
  si.materials    = materials;
  si.debug_values = debug_values;
  si.rng_state    = rng;
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

float Rand ( SceneInfo* si ) {
  enum { A=4294883355U };
  uint2 r = (*si).rng_state;
  uint res = r.x^r.y;
  uint hi = mul_hi(r.x, A);
  r.x = r.x*A + r.y;
  r.y = hi + (r.x<r.y);
  (*si).rng_state = r;
  return res/(float)(UINT_MAX);
}

float3 Float3_Max(float3 t, float g) {
  return (float3)(fmax(t.x, g), fmax(t.y, g), fmax(t.z, g));
}
float3 Float3_Min(float3 t, float g) {
  return (float3)(fmin(t.x, g), fmin(t.y, g), fmin(t.z, g));
}

float Sample_Uniform ( SceneInfo* si ) {
  return Rand(si);
}
float2 Sample_Uniform2 ( SceneInfo* si ) {
  return (float2)(Sample_Uniform(si), Sample_Uniform(si));
}
// -1 .. 1
float2 Sample_Uniform2_2 ( SceneInfo* si ) {
  return ((float2)(Sample_Uniform(si), Sample_Uniform(si))-(float2)(0.5f))*2.0f;
}
float3 Sample_Uniform3 ( SceneInfo* si ) {
  return (float3)(Sample_Uniform(si), Sample_Uniform2(si));
}
// -----------------------------------------------------------------------------
// --------------- GENERAL UTILITIES      --------------------------------------
float To_IOR ( float ior ) {
  return 0.0f + 2.5f*ior;
}
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

float3 refract(float3 I, float3 N, float ior) {
  float cos_NI = dot(N, I);
  float k = 1.0f - SQR(ior)*(1.0f - sqr(cos_NI));
  return k < 0.0f ? (float3)(0.0f) : ior*I - (ior*cos_NI + sqrt(k))*N;
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

float3 BRDF_F ( float3 wi, float3 N, float3 wo, __constant Material* m,
                float3 col ) {
  // get binormal, bitangent, half vec etc
  const float3 binormal  = Binormal(N),
               bitangent = cross(binormal, N),
               L         =  wo, V = -wi,
               H         = normalize(L+V);
  const float  cos_NV    = dot(N, V), cos_NL     = dot(N, L),
               cos_HV    = dot(H, V), cos_HL     = dot(H, L),
               Fresnel_L = Schlick_Fresnel(cos_NL),
               Fresnel_V = Schlick_Fresnel(cos_NV);
  // Diffusive component
  float3 diffusive_albedo = pow(col, 0.5f) * IPI;

  float3 microfacet = (float3)(1.0f);

  if ( isnan(cos_NL) || isnan(cos_NV) || cos_NL < 0.0f || cos_NV < 0.0f )
    return diffusive_albedo;

  { // ------- Fresnel
    // modified diffusive fresnel from disney, modified to use albedo & F0
    const float F0 = m->fresnel * m->metallic,
                Fresnel_diffuse_90 = F0 * SQR(cos_HL);
     microfacet *= (1.0f - F0) * diffusive_albedo +
                    mix(1.0f, Fresnel_diffuse_90, Fresnel_L) *
                    mix(1.0f, Fresnel_diffuse_90, Fresnel_V);
  }

  { // ------- Geometric
    // Heits 2014, SmithGGXCorrelated with half vec combined with anisotropic
    // term using GTR2_aniso model
    const float Param  = 0.5f + m->roughness,
                Aspect = sqrt(1.0f - m->anisotropic*0.9f),
                Ax     = Param/Aspect, Ay = Param*Aspect,
                GGX_NV = Smith_G_GGX_Correlated(cos_HL, cos_NV, Ax),
                GGX_HL = Smith_G_GGX_Correlated(cos_NV, cos_HL, Ay);
    microfacet *= 0.5f / (GGX_NV*Ax + GGX_HL*Ay);
  }
  { // ------- Distribution
    // Hyper-Cauchy Distribution using roughness and metallic
    const float Param = 1.2f + m->anisotropic,
                Shape = (1.1f - m->roughness*0.55f),
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
    // based off the Henyey-Greenstein equation
    const float R = 0.7f*(1.0 - m->roughness),
                M = 0.2f + m->subsurface;
    const float Rr_term = M * (1.0f - sqr(R))*(4.0f*IPI) *
                          (1.0f/(pow(1.0f + sqr(R) - 2.0f*R*cos_HL, 3.0f/2.0f)));
    const float3 Retro_reflection = diffusive_albedo * Rr_term *
                      (Fresnel_L + Fresnel_V + (Fresnel_L*Fresnel_V*(Rr_term)));
    diffusive_albedo = mix(diffusive_albedo, pow(Retro_reflection, (float3)(0.5f)),
                           m->subsurface*0.5f);
  }

  return mix((microfacet + diffusive_albedo),
             col,
             m->specular+m->transmittive);
}

float3 RColour ( float3 pt_colour, __constant Material* m ) {
  return (pt_colour.x >= 0.0f) ? pt_colour : m->albedo;
}

typedef struct T_IlluminateHorcrux {
  float3 colour;
  float3 N;
} IlluminateHorcrux;
IlluminateHorcrux Illuminate ( float3 O, float3 Wi, float3 col,
                               __constant Material* m, SCENE_T(si, Tx) ) {
  float3 colour = (float3)(0.0f);
  float3 N = Normal(O, si, Tx);
  for ( int i = 0; i != EMITTER_AMT; ++ i ) {
    Emitter emit = REmission(i, si->debug_values, si->time);

    float3 V = normalize(Wi),
           L = normalize(emit.origin - O);
    float3 brdf = BRDF_F(V, N, L, m, RColour(col, m));

    colour += (brdf) * emit.emission * dot(N, L);
  }
  return (IlluminateHorcrux){Float3_Min(colour / EMITTER_AMT, 1.0f), N};
}

typedef struct T_SpecularApproxHorcrux {
  float3 origin, colour, N;
  __constant Material* mat;
} SpecularApproxHorcrux;
SpecularApproxHorcrux Specular_Approx(Ray ray, SCENE_T(si, Tx)) {
  ray.origin += ray.dir*0.01f;
  SampledPt spec_result = March(-1, ray, si, Tx);
  if ( spec_result.dist < 0.0f )
    return (SpecularApproxHorcrux){ray.origin, ray.origin, NULL};
  __constant Material* sm = si->materials + spec_result.material_index;
  float3 colour = spec_result.colour;
  float3 origin = ray.origin + ray.dir*spec_result.dist;
  IlluminateHorcrux ih = Illuminate(origin, ray.dir, colour, sm, si, Tx);
  return (SpecularApproxHorcrux){origin, ih.colour, ih.N, sm};
}

SampledPt Colour_Pixel(Ray ray, SCENE_T(si, Tx)) {
  SampledPt result = March(-1, ray, si, Tx);
  if ( result.dist < 0.0f ) return result;
  float3 origin = ray.origin + ray.dir*result.dist;
  __constant Material* m = si->materials + result.material_index;

  IlluminateHorcrux res = Illuminate(origin, ray.dir, result.colour, m, si, Tx);
  result.colour = res.colour;
  // specular approximation
  if ( m->specular > 0.0f ) {
    ray = (Ray){origin, reflect(ray.dir, res.N)};
    SpecularApproxHorcrux sah = Specular_Approx(ray, si, Tx);
    if ( sah.mat )
      result.colour = mix(result.colour, sah.colour, m->specular);
  }
  // transmittive approximation
  if ( m->transmittive > 0.0f ) {
    ray = (Ray){origin, refract(ray.dir, res.N, To_IOR(m->ior))};
    SpecularApproxHorcrux first = Specular_Approx(ray, si, Tx);
    if ( first.mat ) {
      float3 realcol = mix(result.colour, first.colour, m->transmittive);
      if ( first.mat->transmittive > 0.0f ) {
        ray = (Ray){first.origin,
                    refract(ray.dir, first.N, To_IOR(first.mat->ior))};
        SpecularApproxHorcrux sec = Specular_Approx(ray, si, Tx);
        if ( sec.mat )
          realcol = mix(realcol, sec.colour, first.mat->transmittive);
      }
      result.colour = realcol;
    }
  }

  // glossy approximation [none now]

  return result;
}

// -----------------------------------------------------------------------------
// --------------- CAMERA ------------------------------------------------------
Ray Camera_Ray(Camera* camera, SceneInfo* si) {
  float2 coord = (float2)((float)get_global_id(0), (float)get_global_id(1));
  float2 resolution = (float2)((float)camera->dim.x, (float)camera->dim.y);
  resolution.y *= 16.0f/9.0f;

  float fov_r = (180.0f - camera->fov)*PI/180.0f;

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
  /* // ------ DOF & antialiasing ----- */
  // Adapted from https://www.shadertoy.com/view/lsX3DH
  float3 dof_puv = (float3)(puv, fov_r);
  float3 ray_dir = dof_puv.x*cam_right + dof_puv.y*cam_up + fov_r*cam_front;
  float3 dof_origin = camera->radius*(float3)(Sample_Uniform2_2(si), 0.0f);
  float3 dof_dir = normalize(dof_puv * camera->focal - dof_origin);
  cam_pos += dof_origin.x*cam_right + dof_origin.y*cam_up;
  ray_dir += dof_dir.x*cam_right + dof_dir.y*cam_up;
  ray_dir = normalize(ray_dir);
  return (Ray){cam_pos, ray_dir};
}

typedef struct _T_EvalPreviousOutputHorcrux {
  float4 old_colour;
  bool raycast_nav;
} _EvalPreviousOutputHorcrux;

_EvalPreviousOutputHorcrux Eval_Previous_Output(
    __global unsigned char* img, __write_only image2d_t output, int2 out,
    __global SharedInfo* shared_info, int pt) {
  _EvalPreviousOutputHorcrux hor;
  hor.old_colour = (float4)(0.0f);
  hor.raycast_nav = (shared_info->clear_img);
  if ( !hor.raycast_nav) {
    hor.old_colour = (float4)(img[pt+0]/255.0f, img[pt+1]/255.0f,
                              img[pt+2]/255.0f, (float)(img[pt+3]));
  } else {
    img[pt+0] = 0.0f; img[pt+1] = 0.0f; img[pt+2] = 0.0f; img[pt+3] = 0.0f;
    shared_info->finished_samples = 0;
  }
  return hor;
}

// -----------------------------------------------------------------------------
// --------------- RAYTRACE KERNEL ---------------------------------------------
__kernel void DTOADQ_Kernel (
    __global unsigned char*     img, // R G B ITER
    __write_only image2d_t      output_img,
    __global SharedInfo*        sinfo,
    __global Camera*            camera_ptr,
    __global float*             time_ptr,
    __read_only image2d_array_t textures,
    __constant Material*        g_materials,
    __global float*             debug_val_ptr,
    __global uint2*             rng_states
  ) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  Camera camera = *camera_ptr;
  float3 dval = (float3)(debug_val_ptr[0], debug_val_ptr[1], debug_val_ptr[2]);
  float time = *time_ptr;
  // -- get old pixel, check if there are samples to be done
  //    (counter is stored in alpha channel)
  Update_Camera(&camera, time);

  SceneInfo scene_info = New_SceneInfo(time, g_materials, dval,
                      rng_states[out.y*camera.dim.x + out.x]);

  Ray ray = Camera_Ray(&camera, &scene_info);
  SampledPt result = Colour_Pixel(ray, &scene_info, textures);

  float3 col;
  if ( result.dist >= 0.0f ) {
    col = result.colour;
  } else {
    col = (float3)(-1.0f);
  }

  float4 old_colour = (float4)(0.0f);
  int pix_pt = out.y*camera.dim.x*4 + out.x*4;
  {
    _EvalPreviousOutputHorcrux _hor =
      Eval_Previous_Output(img, output_img, out, sinfo, pix_pt);
    old_colour = _hor.old_colour;
    if ( old_colour.w >= (float)(64) ) return;
    if ( _hor.raycast_nav ) {
      write_imagef(output_img, out, (float4)(col, 1.0f));
      return;
    }
  }

  write_imagef(output_img, out, (float4)(col, 1.0f));
  if ( col.x >= 0.0f ) {
    old_colour = (float4)(mix(col, old_colour.xyz,
                              (old_colour.w/(old_colour.w+1.0f))),
                          old_colour.w+1.0f);
    float4 nold_colour;
    write_imagef(output_img, out, (float4)(old_colour.xyz, 1.0f));
    //
    img[pix_pt+0] = (unsigned char)(old_colour.x*255.0f);
    img[pix_pt+1] = (unsigned char)(old_colour.y*255.0f);
    img[pix_pt+2] = (unsigned char)(old_colour.z*255.0f);
    img[pix_pt+3] = (unsigned char)(old_colour.w);
  }
  rng_states[out.y*camera.dim.x + out.x] = scene_info.rng_state;
}
