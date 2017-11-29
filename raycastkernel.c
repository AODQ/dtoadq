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

bool Is_Debug ( ) {
  return get_global_id(0) == get_global_size(0)/2 &&
    get_global_id(1) == get_global_size(1)/2;
}

__constant float MARCH_ACC = //%MARCH_ACC.0f/1000.0f;
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
  Material* materials;
  float3 debug_values;
  uint2 rng_state;
} SceneInfo;
SceneInfo New_SceneInfo(float time, Material* materials,
                        float3 debug_values, uint2 rng){
  SceneInfo si;
  si.time         = time;
  si.materials    = materials;
  si.debug_values = debug_values;
  si.rng_state    = rng;
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
// --------------- GENERAL FUNCTIONS      --------------------------------------

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
    MapUnionG(a, &res, dist, -100 - i, e.emission);
  }

  return res;
}

// -----------------------------------------------------------------------------
// --------------- RAYTRACING/MARCH --------------------------------------------
SampledPt March ( int avoid, Ray ray, SCENE_T(si, Tx) ) {
  float distance = 0.0f;
  SampledPt t_info;
  for ( int i = 0; i < MARCH_REPS; ++ i ) {
    t_info = Map(avoid, ray.origin + ray.dir*distance, si, Tx);
    if ( t_info.dist <= MARCH_ACC || t_info.dist > MARCH_DIST ) break;
    distance += t_info.dist;
  }
  if ( t_info.dist > MARCH_DIST || t_info.dist < 0.0f ) {
    t_info.dist = -1.0f;
    return t_info;
  }
  t_info.dist = distance;
  return t_info;
}

float Distance(float3 u, float3 v) {
  float x = u.x-v.x, y = u.y-v.y, z = u.z-v.z;
  return sqrt(x*x + y*y + z*z);
}

float3 RColour ( float3 pt_colour, Material* m ) {
  return (pt_colour.x >= 0.0f) ? pt_colour : m->albedo;
}

//----POSTPROCESS---
//%POSTPROCESS

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

  Ray ray;

  /* // ------ DOF & antialiasing ----- */
  float3 dof_puv = (float3)(puv, fov_r);
  float3 ray_dir = dof_puv.x*cam_right + dof_puv.y*cam_up + fov_r*cam_front;
  float3 dof_origin = camera->radius*(float3)(Sample_Uniform(si), 0.0f, 0.0f);
  float3 dof_dir = normalize(dof_puv * camera->focal - dof_origin);
  cam_pos += dof_origin.x*cam_right + dof_origin.y*cam_up;
  ray_dir += dof_dir.x*cam_right + dof_dir.y*cam_up;
  ray_dir = normalize(ray_dir);
  ray.origin = cam_pos;
  ray.dir = ray_dir;
  return ray;
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
    __global Material*          g_materials,
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

  Material materials[MAT_LEN];
  for ( int i = 0; i != MAT_LEN; ++ i ) {
    Material tm = *(g_materials + i);
    materials[i] = tm;
  }
  SceneInfo scene_info = New_SceneInfo(time, materials, dval,
                                       rng_states[out.y*camera.dim.x + out.x]);

  Ray ray = Camera_Ray(&camera, &scene_info);
  SampledPt pt = March(-1, ray, &scene_info, textures);
  float3 col = (float3)(-1.0f);
  if ( pt.dist >= 0 ) {
    float3 O = ray.origin + pt.dist*ray.dir;
    for ( int i = 0; i != EMITTER_AMT; ++ i ) {
      Emitter emit = REmission(i, scene_info.debug_values, scene_info.time);
      if ( pt.material_index >= 0 && pt.material_index < MAT_LEN ) {
        Material m = materials[pt.material_index];
        float3 albedo = RColour(pt.colour, &m);
        // diffusive component i made up, no normal needed
        float3 V = -ray.dir, L = normalize(emit.origin-O);
        col = albedo*0.5f + dot(V, L)*0.1f;
      } else {
        col = (float3)(1.0f);
      }
    }
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
