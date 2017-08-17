#define MAX_DEPTH 8
#define MARCH_DIST //%MARCH_DIST.0f
#define MARCH_REPS //%MARCH_REPS

#define TEXTURE_T __read_only image2d_array_t
#define SCENE_T(S, T) SceneInfo* S, TEXTURE_T T
// EMIT_MAT = -120 - light_index
#define EMIT_MAT -120

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
  float diffuse, specular, glossy, retroreflective, transmittive;
} Material;

typedef struct T_Emitter {
  float3 origin, emission;
  float radius;
} Emitter;

typedef struct T_SharedInfo {
  unsigned char clear_img;
  unsigned long finished_samples;
  unsigned char spp;
} SharedInfo;

__constant float PI   = 3.141592654f;
__constant float IPI  = 0.318309886f;
__constant float TAU  = 6.283185307f;
__constant float ITAU = 0.159154943f;
// -----------------------------------------------------------------------------
// --------------- GENERAL STRUCTS ---------------------------------------------
typedef struct T_Ray {
  float3 origin, dir;
} Ray;

typedef struct T_BSDF_Sample_Info {
  float3 wi, wo;
} BSDF_Sample_Info;

typedef struct T_SampledPt {
  float3 colour, origin, dir, normal;
  float dist;
  BSDF_Sample_Info bsdf_info;
  int mat_index;
} SampledPt;
SampledPt SampledPt_From_Origin(float3 origin) {
  SampledPt pt;
  pt.origin = origin;
  pt.dist = 0.0f;
  return pt;
}

typedef struct T_SampledPath {
  SampledPt path [MAX_DEPTH+1];
  uint length;
  Emitter light_emitter;
} SampledPath;

typedef struct T_SceneInfo {
  float time;
  // __read_only image2d_array_t textures; IMAGES CANT BE USED AS FIELD TYPES:-(
  __global Material* materials;
  float3 debug_values;
  __global int* rng;
} SceneInfo;
SceneInfo New_SceneInfo(float time, __global Material* materials,
                        float3 debug_values, __global int* rng){
  SceneInfo si;
  si.time         = time;
  si.materials    = materials;
  si.debug_values = debug_values;
  si.rng          = rng;
  return si;
}
// -----------------------------------------------------------------------------
// --------------- GENERAL FUNCTIONS      --------------------------------------
// http://www0.cs.ucl.ac.uk/staff/ucacbbl/ftp/papers/langdon_2009_CIGPU.pdf
int Parker_Miller_Rand ( int seed ) {
  const int a = 16807,      // 7^5
            m = 2147483647; // 2^31 - 1
  return ((long)(seed * a))%m;
}

int Rand ( __global int* rng ) {
  int id = get_global_id(1)*get_global_size(0) + get_global_id(0);
  id %= 31;
  int rnum = Parker_Miller_Rand(rng[id]);
  rng[id] = rnum;
  return rnum;
}

float Uniform_Sample ( __global int* rng ) {
  return (float)((uint)(Rand(rng)))/(float)(UINT_MAX);
}
float2 Uniform_Sample2 ( __global int* rng ) {
  return (float2)(Uniform_Sample(rng), Uniform_Sample(rng));
}
float3 Uniform_Sample3 ( __global int* rng ) {
  return (float3)(Uniform_Sample(rng), Uniform_Sample2(rng));
}

float sqr(float t) { return t*t; }
float Distance(float3 u, float3 v) {
  float x = u.x+v.x, y = u.y+v.y, z = u.z+v.z;
  return sqrt(x*x + y*y + z*z);
}
float Distance_Sqr(float3 u, float3 v) {
  float x = u.x+v.x, y = u.y+v.y, z = u.z+v.z;
  return (x*x + y*y + z*z);
}

float Power_Heuristic ( float fn, float fpdf, float gn, float gpdf ) {
  float f = sqr(fn*fpdf),
        g = sqr(gn*gpdf);
  return f/(f + g);
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
    d1->mat_index = mi;
  }
}

// -----------------------------------------------------------------------------
// --------------- SCENE -------------------------------------------------------
//%SCENEINSERT
//------------------------
// -----------------------------------------------------------------------------
// --------------- MAP ---------------------------------------------------------
SampledPt Map ( int a, float3 origin, SCENE_T(si, Tx))  {
  SampledPt res;
  res.dist = FLT_MAX;

  //---MAP INSERTION POINT---
  //%MAPINSERT
  //-------------------------

  // lighting with emissions
  float3 light_emission = (float3)(1.0f, 0.9f, 0.8f);
  for ( int i = 0; i != EMITTER_AMT; ++ i ) {
    Emitter e = REmission(i, si->time);
    float dist = sdSphere(origin - e.origin, e.radius);
    int ta = a;
    if ( a == EMIT_MAT + 1 ) ta = EMIT_MAT - i;
    MapUnionG(a, &res, dist, EMIT_MAT - i, light_emission*e.emission);
  }

  return res;
}

// -----------------------------------------------------------------------------
// --------------- RAYMARCHING      --------------------------------------------
SampledPt March ( int avoid, Ray ray, SCENE_T(si, Tx)) {
  float distance = 0.0f;
  SampledPt t_info;
  for ( int i = 0; i < MARCH_REPS; ++ i ) {
    t_info = Map(avoid, ray.origin + ray.dir*distance, si, Tx);
    if ( t_info.dist < MARCH_ACC || t_info.dist > MARCH_DIST ) break;
    distance += t_info.dist;
    if ( t_info.mat_index != avoid ) avoid = -1;
  }
  t_info.dir    = ray.dir;
  if ( t_info.dist > MARCH_DIST ) {
    t_info.dist = -1.0f;
    return t_info;
  }
  t_info.dist = distance;
  t_info.origin = ray.origin + ray.dir*t_info.dist;
  return t_info;
}

// -----------------------------------------------------------------------------
// --------------- GRAPHIC FUNCS -----------------------------------------------
float3 Normal ( float3 p, SCENE_T(s, t)) {
  float2 e = (float2)(1.0f, -1.0f)*0.5773f*0.0005f;
  return normalize(
    e.xyy*Map(-1, p + e.xyy, s, t).dist + e.yyx*Map(-1, p + e.yyx, s, t).dist +
    e.yxy*Map(-1, p + e.yxy, s, t).dist + e.xxx*Map(-1, p + e.xxx, s, t).dist);
}

float3 reflect ( float3 V, float3 N ) {
  return V - 2.0f*dot(V, N)*N;
}

float3 refract(float3 V, float3 N, float refraction) {
  float cosI = -dot(N, V);
  float cosT = 1.0f - refraction*refraction*(1.0f - cosI*cosI);
  return (refraction*V) + (refraction*cosI - sqrt(cosT))*N;
}

float3 To_Cartesian ( float sin_theta, float cos_theta, float phi ) {
  return (float3)(cos(phi)*sin_theta, sin(phi)*sin_theta, cos_theta);
}
// Cosine weighted
float3 Hemisphere_Direction ( __global int* rng, float3 normal ) {
  const float phi = TAU*Uniform_Sample(rng),
              vv  = 2.0f*(Uniform_Sample(rng) - 0.5f);
  const float cos_theta = sign(vv)*sqrt(fabs(vv)),
              sin_theta = sqrt(fmax(0.0f, 1.0f - (cos_theta*cos_theta)));
  return To_Cartesian(sin_theta, cos_theta, phi);
}


float3 Cone_Direction( __global int* rng, float cos_theta_max,
                       float3 normal, float3 tangent, float3 binormal ) {
  float2 ur = (float2)(Uniform_Sample(rng), Uniform_Sample(rng));
  float cos_theta = (1.0f - ur.x) + ur.x*cos_theta_max,
        sin_theta = sqrt(1.0f - cos_theta*cos_theta);
  float phi = ur.y*TAU;
  float3 sample_v = To_Cartesian(sin_theta, cos_theta, phi);
  return sample_v.x*normal + sample_v.y*tangent + sample_v.z*binormal;
}

void Calculate_Binormals ( float3 N, float3* T, float3* B ) {
  *T = (fabs(N.x) > fabs(N.y)) ? normalize((float3)(-N.z, 0.0f,  N.x)) :
                                 normalize((float3)(0.0f, N.z,  -N.y));
  *B = cross(N, *T);
}
// -----------------------------------------------------------------------------
// --------------- BSDF    FUNCS -----------------------------------------------
float3 BSDF_Sample ( float3 wi, float3 normal, Material* m, float* out_pdf,
                     float3* out_bsdf, SCENE_T(si, Tx)) {
  float3 wo;
  wo = normalize(Hemisphere_Direction(si->rng, normal));
 
  *out_pdf = fabs(wo.z)*IPI;
  *out_bsdf = (float3)(IPI*m->diffuse);
  return wo;
}

float3 BSDF_Evaluate ( float3 wo, Material* m, float* out_pdf ) {
  *out_pdf = fabs(wo.z)*IPI;
  return (float3)(IPI*m->diffuse);
}

// -----------------------------------------------------------------------------
// --------------- LIGHT   FUNCS -----------------------------------------------
float Visibility_Ray(Ray ray, float max_dist, int mindex, SCENE_T(si, Tx)) {
  ray.origin += ray.dir*(MARCH_ACC+0.1f);
  SampledPt minfo = March(-1, ray, si, Tx);
  return step(max_dist, minfo.dist);
}

float3 Light_Emission(float distance) {
  return (float3)(1.0f, 0.98f, 0.95f)/SQR(distance);
}

float3 Light_Sample ( float3 wi, float3 P, float3 N, Emitter* m, float* pdf,
                      __global int* rng ) {
  float3 light_point = m->origin;
  float3 offset = normalize(Uniform_Sample3(rng));
  light_point += offset*m->radius*Uniform_Sample(rng);
  float3 light_dir = normalize(light_point - P);

  float cos_theta = fabs(dot(N, -light_dir));

  // if ( cos_theta < 0.001f ) return (float3)(0.0f);


  float area = 1.0f/(2*TAU*(m->radius*m->radius));
  *pdf = area*(sqr(length(light_dir)))/cos_theta;

  return light_dir;
}

// -----------------------------------------------------------------------------
// --------------- LIGHT TRANSPORT - BDPT --------------------------------------
int Generate_Subpath ( Ray ray, SampledPt* subpath, SCENE_T(si, Tx)) {
  int mindex = -1;
  int i;
  float roulette = 1.0f;
  for ( i = 0; i != MAX_DEPTH/2; ++ i ) {
    // // roulette
    // roulette -= 1.0f/(float)(MAX_DEPTH);
    // if ( Uniform_Sample(si->rng) > roulette ) break;
    // // march
    // SampledPt minfo = March(mindex, ray, si, Tx);
    // if ( minfo.dist < 0.0f ) break;
    // // importance sampling
    // mindex = minfo.mat_index;
    // Material material = si->materials[minfo.mat_index];
    // float3 hit = ray.origin + ray.dir*minfo.dist;
    // BSDF_Sample_Info bsdf_info = (BSDF_Sample_Info){
    //     ray.dir,
    //     BSDF_Sample(ray.dir, Normal(hit, si, Tx), &material, 0, si, Tx)
    // };
    // // set minfo/subpath element
    // minfo.origin = hit;
    // minfo.bsdf_info = bsdf_info;
    // subpath[i] = minfo;
    // // set ray
    // ray = (Ray){hit, bsdf_info.wo};
    // if ( minfo.mat_index <= EMIT_MAT ) break;
  }
  return i;
}

int Generate_BSDF_Path(Ray ray, SampledPt* subpath, SCENE_T(si, Tx)) {
  return Generate_Subpath(ray, subpath, si, Tx);
}

// does not do retro for you
int Generate_Light_Path (Emitter* light, SampledPt* subpath, SCENE_T(si, Tx)) {
  Ray ray = (Ray){light->origin, Uniform_Sample3(si->rng)};
  return Generate_Subpath(ray, subpath, si, Tx);
}

// TODO FIX
bool Check_Subpath(float3 xorigin, float3 yorigin, float xmat, float ymat,
                   SCENE_T(si, Tx)) {
  float3 angle = normalize(yorigin - xorigin);
  float xdist = March(xmat, (Ray){xorigin, normalize(yorigin - xorigin)},
                      si, Tx).dist;
  float ydist = March(ymat, (Ray){yorigin, normalize(xorigin - yorigin)},
                      si, Tx).dist;
  return true;
  // return ( fabs(ydist - xdist) < MARCH_ACC+0.01f );
}

bool Generate_Path ( Ray ray, SampledPath* p, SCENE_T(si, Tx)) {
  SampledPt light_path[MAX_DEPTH/2+1];
  // grab a random light
  int light_index = (int)(Uniform_Sample(si->rng)*EMITTER_AMT);
  Emitter light_emitter = REmission(light_index, si->time);
  p->light_emitter = light_emitter;
  // generate paths, use p->path as bsdf_path
  int bsdf_length  = Generate_BSDF_Path (ray,              p->path,  si, Tx);
  int light_length = Generate_Light_Path(&light_emitter, light_path, si, Tx);
  // --- configure last points on each path.. some special cases here as
  //    outlined in Veach's thesis
  if ( bsdf_length == 0 ) {
    if ( light_length == 0 ) { return false; }
    SampledPt* light = &light_path[light_length-1];
    if ( !Check_Subpath(ray.origin, light->origin, -1,light->mat_index, si, Tx))
      return false;
    light->bsdf_info.wo = normalize(ray.origin - light->origin);
  } else if ( light_length == 0 ) {
    SampledPt* bsdf = &p->path[bsdf_length-1];
    if ( !(Check_Subpath(bsdf->origin, light_emitter.origin,
            bsdf->mat_index, EMIT_MAT+1, si, Tx)))
      return false;
    bsdf->bsdf_info.wo = normalize(bsdf->origin - light_emitter.origin);
  } else {
    SampledPt* bsdf  = &p->path   [bsdf_length -1],
             * light = &light_path[light_length-1];
    // check connection made
    if ( !Check_Subpath(bsdf->origin, light->origin,
                bsdf->mat_index, light->mat_index, si, Tx) )
      return false;
    // set bsdf[$].wo = light[$].wo = |^ - $|
    bsdf->bsdf_info.wo = light->bsdf_info.wo =
          normalize(bsdf->origin - light->origin);
  }

  // retro light paths and concatenate
  for ( int i = 0; i != light_length; ++ i ) {
    float3 t_wi = light_path[i].bsdf_info.wi;
    light_path[i].bsdf_info.wi = -light_path[i].bsdf_info.wo;
    light_path[i].bsdf_info.wo = -t_wi;
    p->path[bsdf_length + light_length - i - 1] = light_path[i];
  }

  p->length = light_length + bsdf_length;
  return true;
}



// -----------------------------------------------------------------------------
// --------------- LIGHT TRANSPORT ---------------------------------------------
float MIS_Weight ( float pdf_i, float pdf_o ) {
  return pdf_i / ( pdf_i + pdf_o );
}
SampledPt Colour_Pixel ( Ray ray, SCENE_T(si, Tx) ) {
  float3 coeff = (float3)(1.0f);
  float3 colour = (float3)(0.0f);
  float bsdf_pdf = 1.0f;
  int mindex = -1;
  SampledPt result;
  result.dist = -1.0f;
  for ( int depth = 0; depth != 1; ++ depth ) {
    SampledPt minfo = March(mindex, ray, si, Tx);
    if ( minfo.dist < 0.0f ) break;
    // importance sampling
    mindex = minfo.mat_index;
    Material material = si->materials[mindex];
    float3 hit = ray.origin + ray.dir*minfo.dist;
    float3 normal = Normal(hit, si, Tx);
    if ( mindex <= EMIT_MAT ) {
      Emitter emitter = REmission(abs(mindex - EMIT_MAT), si->time);
      float light_pdf;
      Light_Sample(ray.dir, hit, normal, &emitter, &light_pdf, si->rng);
      float weight = MIS_Weight(bsdf_pdf, (1.0f/EMITTER_AMT)*light_pdf);
      // colour += coeff * emitter.emission * weight;
      // result.dist = 1.0f;
      // break;
    }

    { // light sample
      int light_index = (int)(Uniform_Sample(si->rng)*EMITTER_AMT);
      // if ( light_index == 0 ) light_index = 1;
      // light_index = 1;
      Emitter light = REmission(light_index, si->time);
      float light_pdf;
      float3 light_sample = Light_Sample(ray.dir, hit, normal, &light,
                                          &light_pdf, si->rng);
      SampledPt mt = March(mindex, (Ray){hit, light_sample}, si, Tx);
      if ( mt.mat_index <= EMIT_MAT ) {
        float3 factor = BSDF_Evaluate(light_sample, &material, &bsdf_pdf);
        float weight = MIS_Weight(light_pdf*(1.0f/EMITTER_AMT), bsdf_pdf);
        float cos_theta = light_sample.z;
        float3 contrib = (weight * cos_theta/((1.0f/EMITTER_AMT)*light_pdf)
                         ) * (factor*light.emission);
        colour += coeff * contrib;
        result.dist = 1.0f;
      }
    }

    BSDF_Sample_Info bsdf_info;
    {
      float3 bsdf;
      bsdf_info = (BSDF_Sample_Info){
          ray.dir,
          BSDF_Sample(ray.dir, normal, &material, &bsdf_pdf, &bsdf, si, Tx)
      };

      float dot_nl = fabs(dot(bsdf_info.wi, bsdf_info.wo));
      coeff *= minfo.colour/PI * (bsdf/bsdf_pdf) * dot_nl;
    }


    ray.dir = bsdf_info.wo;
    ray.origin = hit;
  }
  result.colour = colour;
  return result;
}
//   // Bidirectional Path Tracing . .
//   // Generate camera and light subpaths
//   SampledPath transport;
//   if ( !Generate_Path(ray, &transport, si, Tx) || transport.length == 0 ) {
//     SampledPt result;
//     result.dist = -1.0f;
//     result.colour = (float3)(0.0f);
//     return result;
//   }

//   // Evaluate path
//   SampledPt result;
//   result.dist = 1.0f;

//   /*

//     l₀ —————> l₁
//              /
//            /
//          /
//        /
//       b₁
//        \
//         \
//          \
//           b₀
// S^2
// S²
//     weight =
//       Lₑ(l₀, l₀->l₁)
//          G(l₀, l₁) V(l₀, l₁) Fᵣ(l₁, l₀->l₁, l₁->b₁)
//          G(l₁, b₁) V(l₁, b₁) Fᵣ(b₁, l₁->b₁, b₁->b₀)
//          ...
//       / probability_of_path
//   */

//   float3 colour = (float3)(0.0f),
//          coeff  = (float3)(1.0f);

//   for ( int t = 0; t < transport.length; ++ t ) {
//     SampledPt pt = transport.path[t];
//     Material mat = si->materials[pt.mat_index];
//     float3 albedo = pt.colour;
//     float3 normal = Normal(pt.origin, si, Tx);
//     float3 bsdf = BSDF(pt.bsdf_info, normal, &mat);
//     float  pdf  = BSDF_PDF(pt.origin, pt.bsdf_info, normal, &mat);
//     float3 hv = normalize(pt.bsdf_info.wi + pt.bsdf_info.wo);
//     hv *= hv;
//     coeff *= albedo*bsdf*dot(normal, pt.bsdf_info.wo)/pdf;
//   }

//   colour = coeff*transport.length*transport.length;

//   result.colour = colour;
//   return result;

//   // for ( int t = 1; t <= depth - 1; ++ t ) {
//   //   SampledPt pt = t < transport.bsdf_length ? transport.bsdf_path[t] :
//   //                       transport.light_path[t - transport.light_length];
//   //   Material mat = si->materials[pt.mat_index];
//   //   float3 albedo = pt.colour;
//   //   float3 normal = Normal(pt.origin, si, Tx);
//   //   if ( t < transport.bsdf_length ) {
//   //     coeff *= pt.bsdf_info.brdf * pt.bsdf_info.dir_cos/pt.bsdf_info.pdf;
//   //   } else {
//   //    // BSDF_Sample_Info bsdf = BSDF(pt.origin, pt.bsdf_info.wo,
//   //     //                                        pt.bsdf_info.wi,
//   //     //                                        mat, albedo, si, Tx);
//   //     // coeff *= pt.bsdf_info.brdf;
//   //   }
//   // }
//   // SampledPt result;
//   // result.dist = -1.0f;
//   // result.colour = (float3)(0.0f);
//   // float3 coeff = (float3)(1.0f);
//   // for ( int i = 0; i != 1; ++ i ) {
//   //   SampledPt pt = March(-1, ray, si, Tx);
//   //   if ( pt.mat_index == EMIT_MAT ) {
//   //     coeff = (float3)(1.0f);
//   //     break;
//   //   }
//   //   Material material = si->materials[pt.mat_index];
//   //   pt.normal = Normal(pt.origin, si, Tx);
//   //   coeff *= Integrate_Lighting(&pt, &material, -pt.dir, si, Tx);
//   // }
//   // result.colour = 10.0f * coeff;
//   // result.dist = 1.0f;
//   // return result;
// }

// -----------------------------------------------------------------------------
// --------------- CAMERA ------------------------------------------------------
Ray Camera_Ray(Camera* camera) {
  float2 coord = (float2)((float)get_global_id(0), (float)get_global_id(1));
  float2 resolution = (float2)((float)camera->dim.x, (float)camera->dim.y);
  resolution.y *= 16.0f/9.0f;

  float2 mouse_pos = camera->lookat.xy;

  float2 puv = -1.0f + 2.0f * (coord/resolution);

  float input_angle = PI - 2.0f*PI*mouse_pos.x;
  float3 cam_pos    = camera->position;
  float3 cam_target = cam_pos + (float3)(sin(input_angle),
            (3.0f * mouse_pos.y) - 1.0f, cos(input_angle));
  float3 cam_front = normalize(cam_target - cam_pos);
  float3 cam_right = normalize ( cross(cam_front, (float3)(0.0f, 1.0f, 0.0f)));

  float3 cam_up  = normalize(cross(cam_right, cam_front));
  float3 ray_dir = normalize(puv.x*cam_right + puv.y*cam_up +
                             (180.0f - camera->fov)*PI/180.0f*cam_front);

  Ray ray;
  ray.origin = cam_pos;
  ray.dir = ray_dir;
  return ray;
}

float4 Eval_Previous_Output(__global unsigned char* img,
            __write_only image2d_t output, int2 out,
            __global SharedInfo* shared_info, int pt) {
  float4 old_colour = (float4)(0.0f);
  if ( !(shared_info->clear_img) ) {
    old_colour = (float4)(img[pt+0]/255.0f, img[pt+1]/255.0f,
                          img[pt+2]/255.0f, (float)(img[pt+3]));
  } else {
    img[pt+0] = 0.0f; img[pt+1] = 0.0f; img[pt+2] = 0.0f; img[pt+3] = 0.0f;
    write_imagef(output, out, (float4)(old_colour.xyz, 1.0f));
    shared_info->finished_samples = 0;
  }
  return old_colour;
}

// -----------------------------------------------------------------------------
// --------------- KERNEL ------------------------------------------------------
__kernel void DTOADQ_Kernel (
    __global unsigned char*     img, // R G B ITER
    __write_only image2d_t      output_img,
    __global SharedInfo*        sinfo,
    __global Camera*            camera_ptr,
    __global float*             time_ptr,
    __read_only image2d_array_t textures,
    __global Material*          materials,
    __global float*             debug_val_ptr,
    __global int*               rng
  ) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  float spp = (float)(sinfo->spp);
  Camera camera = *camera_ptr;
  int pt = out.y*camera.dim.x*4 + out.x*4;
  float4 old_colour = Eval_Previous_Output(img, output_img, out, sinfo, pt);
  if ( old_colour.w >= spp ) return;
  float time = *time_ptr;
  float3 dval = (float3)(debug_val_ptr[0], debug_val_ptr[1],
                          debug_val_ptr[2]);
  SceneInfo scene_info = New_SceneInfo(time, materials, dval, rng);
  // -- get old pixel, check if there are samples to be done
  //    (counter is stored in alpha channel)
  Update_Camera(&camera, time);

  Ray ray = Camera_Ray(&camera);
  SampledPt result = Colour_Pixel(ray, &scene_info, textures);

  float3 colour;
  if ( result.dist >= 0.0f ) {
    colour = result.colour;
    old_colour = (float4)(mix(colour, old_colour.xyz,
                            (old_colour.w/(old_colour.w+1.0f))),
                        old_colour.w+1.0f);
    write_imagef(output_img, out, (float4)(old_colour.xyz, 1.0f));
    //
    img[pt+0] = (unsigned char)(old_colour.x*255.0f);
    img[pt+1] = (unsigned char)(old_colour.y*255.0f);
    img[pt+2] = (unsigned char)(old_colour.z*255.0f);
    img[pt+3] = (unsigned char)(old_colour.w);
    //
  }

  // convert to Y CB R, from matlab
  // if ( sinfo->finished_samples >= camera.dim.x*camera.dim.y ) {
  //   float3 colour;
  //   colour = (float3)(img[pt+0]/255.0f, img[pt+1]/255.0f, img[pt+2]/255.0f);
  //   colour = (float3)(
  //       16.0f + (  65.481f*colour.x + 128.553f*colour.y + 24.966f*colour.z),
  //       128.0f + (- 37.797f*colour.x - 74.2030f*colour.y + 112.00f*colour.z),
  //       128.0f + ( 112.000f*colour.x - 93.7860f*colour.y - 18.214f*colour.z)
  //   );
  //   int irwx = (camera.dim.y - out.y)*camera.dim.x + out.x;
  //   int irwy = irwx + camera.dim.x*camera.dim.y;
  //   int irwz = irwx + camera.dim.x*camera.dim.y*2;
  //   img[irwx] = (unsigned char)(colour.x);
  //   img[irwy] = (unsigned char)(colour.y);
  //   img[irwz] = (unsigned char)(colour.z);
  // }
}
// DTOADQ_kernel string
