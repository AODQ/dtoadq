module dtoadqkernel; immutable(string) DTOADQ_kernel = q{
#define MAX_DEPTH 16
#define MARCH_DIST //%MARCH_DIST.0f
#define MARCH_REPS //%MARCH_REPS

__constant float MARCH_ACC = //%MARCH_ACC.0f/1000.0f;
// -----------------------------------------------------------------------------
// --------------- DEBUG -------------------------------------------------------
// Variadic functions not supported, so this is best you get :-(
bool Is_Debug_Print ( ) {
  return get_global_id(0) == 200 &&
         get_global_id(1) == 200;
}
// -----------------------------------------------------------------------------
// --------------- GPU-CPU STRUCTS ---------------------------------------------
typedef struct T_Material {
  float3 base_colour;
  float metallic, subsurface, specular, roughness, specular_tint,
        anisotropic, sheen, sheen_tint, clearcoat, clearcoat_gloss,
        emission;
} Material;

typedef struct T_Camera {
  float3 position, lookat, up;
  int2 dim;
  float fov;
  int flags;
} Camera;

typedef struct RNG {
  ulong seed[16];
  ulong p;
} RNG;
// -----------------------------------------------------------------------------
// --------------- RANDOM ------------------------------------------------------
// --- random generation via xorshift1024star
ulong RNG_Next(RNG* rng) {
  const ulong s0 = (*rng).seed[(*rng).p];
  ulong       s1 = (*rng).seed[(*rng).p = ((*rng).p + 1)&15];
  s1 ^= s1 << 31; // a
  (*rng).seed[(*rng).p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b, c
  return (*rng).seed[(*rng).p] * 1181783497276652981L *
              (get_global_id(0) + 250) * (get_global_id(1) + 250);
}

float Uniform(RNG* rng, const float min, const float max) {
  return min + ((float)RNG_Next(rng) / (float)(ULONG_MAX/(max-min)));
}

float Uniform_Sample ( RNG* rng ) { return Uniform(rng, 0.0f, 1.0f); }

float3 Uniform_Float3(RNG* rng, const float min, const float max) {
  return (float3)(
    Uniform(rng, min, max),
    Uniform(rng, min, max),
    Uniform(rng, min, max)
  );
}
// -----------------------------------------------------------------------------
// --------------- NOISE -------------------------------------------------------
float2 fract2f ( float2 vec ) {float2 itptr; return fract(vec, &itptr);}
float  fract1f ( float  vec ) {float  itptr; return fract(vec, &itptr);}

float rand ( float2 n ) {
  return fract1f(sin(dot(n, (float2)(19.9898f, 4.1414f))) * 43758.5453f);
}
// -----------------------------------------------------------------------------
// --------------- MATRICES ----------------------------------------------------
typedef struct T_mat3 {
  float3 x, y, z;
} mat3;

float3 mat3_mul ( mat3 mat, float3 vec ) {
  return (float3)(
    mat.x.x*vec.x + mat.x.y*vec.y + mat.x.z*vec.z,
    mat.y.x*vec.x + mat.y.y*vec.y + mat.y.z*vec.z,
    mat.z.x*vec.x + mat.z.y*vec.y + mat.z.z*vec.z
  );
}

mat3 nmat3 ( float3 x_, float3 y_, float3 z_ ) {
  mat3 mat;
  mat.x = x_;
  mat.y = y_;
  mat.z = z_;
  return mat;
}

mat3 rotate_y ( float fi ) {
  return nmat3(
    (float3)(cos(fi), 0.0f, sin(fi)),
    (float3)(0.0f, 1.0f, 0.0f),
    (float3)(-sin(fi), 0.0f, cos(fi))
  );
}

mat3 rotate_x ( float fi ) {
  return nmat3(
    (float3)(1.0f, 0.0f, 0.0f),
    (float3)(0.0f, cos(fi), -sin(fi)),
    (float3)(0.0f, sin(fi),  cos(fi))
  );
}
// -----------------------------------------------------------------------------
// --------------- BASIC GRAPHIC FUNCTIONS -------------------------------------
__constant float PI = 3.1415926535f;

float3 reflect ( float3 V, float3 N ) {
  return V - 2.0f*dot(V, N)*N;
}

float sqr ( float t ) { return t*t; }

float3 refract(float3 V, float3 N, float refraction) {
  float cosI = -dot(N, V);
  float cosT = 1.0f - refraction*refraction*(1.0f - cosI*cosI);
  return (refraction*V) + (refraction*cosI - sqrt(cosT))*N;
}

// -----------------------------------------------------------------------------
// --------------- GENERAL STRUCTS ---------------------------------------------
typedef struct T_Ray {
  float3 origin, dir;
} Ray;

Ray New_Ray ( float3 o, float3 d ) {
  Ray ray;
  ray.origin = o;
  ray.dir    = d;
  return ray;
}

typedef struct T_IntersectionInfo {
  float dist;
  Material material;
  float3 origin, dir, normal;
} IntersectionInfo;


typedef struct T_MapInfo {
  float3 colour;
  float dist;
  int material_index;
} MapInfo;

MapInfo Create_Map_Info ( float d, int mi, float3 c ) {
  MapInfo info;
  info.dist = d;
  info.colour = c;
  info.material_index = mi;
  return info;
}
// -----------------------------------------------------------------------------
// --------------- MAP GEOMETRY FUNCTIONS --------------------------------------

//---MAP GEOMETRY INSERTION POINT---
//%MAPFUNCDECLARATIONS
//----------------------------------
//%MAPFUNCDEFINITIONS
//----------------------------------

void MapUnionG( int avoid, MapInfo* d1, MapInfo d2 ) {
  if ( d2.material_index != avoid && d1->dist > d2.dist )
    *d1 = d2;
}

// -----------------------------------------------------------------------------
// --------------- SCENE -------------------------------------------------------
//%SCENEINSERT
//------------------------
// -----------------------------------------------------------------------------
// --------------- MAP ---------------------------------------------------------
MapInfo Map ( int a, float3 origin, float time,
             __read_only image2d_array_t textures ) {
  MapInfo res = Create_Map_Info(FLT_MAX, -1, (float3)(0.3f));

  //---MAP INSERTION POINT---
  //%MAPINSERT
  //-------------------------

  return res;
}

// -----------------------------------------------------------------------------
// --------------- GRAPHIC FUNCS THAT NEED MAP ---------------------------------
float3 Normal ( float3 p, float time, __read_only image2d_array_t t ) {
  float2 eps = (float2)(0.001f,0.0);
	return normalize((float3)(
    (Map(-1, p+eps.xyy, time, t).dist - Map(-1, p-eps.xyy, time, t).dist),
    (Map(-1, p+eps.yxy, time, t).dist - Map(-1, p-eps.yxy, time, t).dist),
    (Map(-1, p+eps.yyx, time, t).dist - Map(-1, p-eps.yyx, time, t).dist)
  ));
}

/**
  FIXME !!
  Doesn't work! need to find hwo to properly get tangent!!
*/
float3 Tangent ( float3 normal ) {
  float3 t1 = cross(normal, (float3)(0.0f, 0.0f, 1.0f)),
         t2 = cross(normal, (float3)(0.0f, 1.0f, 0.0f));
  if ( length(t1) > length(t2) ) return t1;
  return t2;
}

// float3 BRDF ( float3 pos, float3 cam, float3 lpos, Material* material ) {
//   float3 N         = normalize(Normal(pos)),
//          L         = normalize(lpos - pos),
//          V         = normalize(cam - pos),
//          tangent   = normalize(cross(Tangent(N), N)),
//          bitangent = normalize(cross(N, tangent));
//   float3 result = Disney_BRDF(L, V, N, tangent, bitangent, material);

//   result *= dot(N, L);
//   return result;
// }

// -----------------------------------------------------------------------------
// --------------- RAYTRACING/MARCH --------------------------------------------
MapInfo March ( int avoid, Ray ray, float time,
                __read_only image2d_array_t textures ) {
  float distance = 0.0f;
  MapInfo t_info;
  for ( int i = 0; i < MARCH_REPS; ++ i ) {
    t_info = Map(avoid, ray.origin + ray.dir*distance, time, textures);
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

// uses cosine weight random hemisphere directions
// thanks to embree
float3 Hemisphere_Direction ( RNG* rng, float3 normal ) {
  const float phi = 2.0f*PI*Uniform_Sample(rng),
              vv  = 2.0f*(Uniform_Sample(rng) - 0.5f);
  const float cos_theta = sign(vv)*sqrt(fabs(vv)),
              sin_theta = sqrt(fmax(0.0f, 1.0f - (cos_theta*cos_theta)));
  return (float3)(cos(phi)*sin_theta, sin(phi)*sin_theta, cos_theta);
}

float3 Hemisphere_Weighted_Direction ( RNG* rng, float3 normal,
                                       float weight ) {
  float3 dir = Hemisphere_Direction(rng, normal);
  return normalize(dir + normal*weight);
}

float3 Orient_Normal ( float3 normal, float3 direction ) {
  return normal * (dot(normal, direction) < 0.0f ? -1.0f : 1.0f);
}

Ray TODO_BRDF_Reflect ( RNG* rng, const IntersectionInfo* info) {
  Ray rout;
  // if ( Uniform(rng, 0.0f, 1.0) > material.metallic ) {
    // diffuse
  if ( info->material.metallic < 0.2f )
    rout.dir = normalize(Hemisphere_Direction(rng, info->normal));
  else
    rout.dir = reflect(info->dir, info->normal);
  // } else {
  //   // specular/transmission
  //   float3 axis = Hemisphere_Weighted_Direction(rng, info.normal,
  //                         material.specular);
  //   if ( Uniform(rng, 0.0f, 1.0f) > material.anisotropic ) {
  //     //reflect
  //     rout.dir = normalize(info.dir - axis*dot(info.dir, axis)*2.0f);
  //   } else {
  //     //refract
  //     // TODO ! !
  //   }
  // }
  rout.origin = info->origin + rout.dir*0.2f;
  return rout;
}

float3 TODO_BRDF_Radiance ( float3 radiance, IntersectionInfo* info ){
  if ( info->material.metallic < 0.2f )
    return info->material.base_colour;
    // return (radiance+info->material.emission)*info->material.base_colour*PI;
    //         // cos(dot(info->normal, info->dir)/2.0f);
  else
    return radiance*0.8f;
}

typedef struct T_RayInfo {
  bool hit;
  float3 colour;
} RayInfo;

float3 Jitter ( float3 d, float phi, float sina, float cosa ) {
  float3 w = normalize(d), u = normalize(cross(w.yzx, w)), v = cross(w, u);
  return (u*cos(phi) + v*sin(phi))*sina + w*cosa;
}

RayInfo RPixel_Colour ( RNG* rng, const __global Material* material, Ray ray,
                        float time, __read_only image2d_array_t textures ) {
  //%KERNELTYPE
}

// -----------------------------------------------------------------------------
// --------------- CAMERA ------------------------------------------------------
Ray Camera_Ray(RNG* rng, Camera* camera) {
  float2 coord = (float2)((float)get_global_id(0), (float)get_global_id(1));
  coord += (float2)(Uniform(rng, -1.0f, 1.0f),
                    Uniform(rng, -1.0f, 1.0f));
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

  return New_Ray(cam_pos, ray_dir);
}

// -----------------------------------------------------------------------------
// --------------- KERNEL ------------------------------------------------------
__kernel void DTOADQ_Kernel (
        __read_only  image2d_t input_image,
        __write_only image2d_t output_image,
        __global bool* reset_image,
        __global RNG* rng_ptr,
        __global Camera* camera,
        __global float* time_ptr,
        __global Material* material, __global uint* material_size,
        __read_only image2d_array_t textures
      ) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  int image_rw_index = camera->dim.y*out.y + out.x;
  float time = *time_ptr;
  RNG rng = *rng_ptr;
  // -- get old pixel, check if there are samples to be done
  //    (counter is stored in alpha channel)
  float4 old_pixel = (float4)(0.0f);
  if ( !*reset_image ) {
    old_pixel = read_imagef(input_image, out);
  }

  Camera local_camera = *camera;
  Update_Camera(&local_camera, time);
  if ( old_pixel.w > 1024 ) return;
  // -- set up camera and stack

  if ( old_pixel.w < 1.0f ) {
    old_pixel *= 0.2f;
    old_pixel.w = 1.0f;
    // old_pixel *= 1.5f;
  }


  Ray ray = Camera_Ray(&rng, &local_camera);
  RayInfo result = RPixel_Colour(&rng, material, ray, time, textures);

  if ( result.hit ) {
    old_pixel =
      (float4)(
        mix(
          result.colour,
          old_pixel.xyz,
          (old_pixel.w/(old_pixel.w+1.0f))),
        old_pixel.w+1.0f);
  }

  write_imagef(output_image, out, old_pixel);

  if ( get_global_id(0) == 0 && get_global_id(1) == 0 ) {
    *rng_ptr = rng;
  }
}};
