module opencl_kernel; immutable(string) Test_raycast_string = q{
#define MAX_DEPTH 512
// -----------------------------------------------------------------------------
// --------------- DEBUG -------------------------------------------------------
bool Is_Debug_Print ( ) {
  return get_global_id(get_work_dim()/2) == 1 &&
         get_global_id(get_work_dim()/2) == 1;
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
} Camera;

typedef struct RNG {
  ulong seed[16];
  ulong p;
} RNG;
// -----------------------------------------------------------------------------
// --------------- RANDOM ------------------------------------------------------
// --- random generation via xorshift1024star
ulong RNG_Next(__global RNG* rng) {
  const ulong s0 = (*rng).seed[(*rng).p];
  ulong       s1 = (*rng).seed[(*rng).p = ((*rng).p + 1)&15];
  s1 ^= s1 << 31; // a
  (*rng).seed[(*rng).p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b, c
  return (*rng).seed[(*rng).p] * 1181783497276652981L *
              (get_global_id(0) + 250) * (get_global_id(1) + 250);
}

float Uniform(__global RNG* rng, const float min, const float max) {
  return min + ((float)RNG_Next(rng) / (float)(ULONG_MAX/(max-min)));
}

float3 Uniform_Float3(__global RNG* rng, const float min, const float max) {
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
  return fract1f(sin(dot(n, (float2)(19.9898, 4.1414))) * 43758.5453);
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

float length2(float2 p, int n) {
  return pow(pown(p.x, n) + pown(p.y, n), 1.0f/n);
}
float length3(float3 p, int n) {
  return pow(pown(p.x, n) + pown(p.y, n), 1.0f/n);
}
// -----------------------------------------------------------------------------
// --------------- RAY ---------------------------------------------------------
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
  int material_index;
  float3 origin, dir, normal;
} IntersectionInfo;
// -----------------------------------------------------------------------------
// --------------- BRDF/BSDF FUNCTIONS -----------------------------------------
float Schlick_Fresnel ( float u ) {
  return pown(clamp(1.0f - u, 0.0f, 1.0f), 5);
}

float GTR1 ( float NdotH, float a ) {
  if ( a >= 1.0f ) return 1.0f/PI;
  float a2 = a * a;
  float t = 1.0f + (a2 - 1.0f)*NdotH*NdotH;
  return (a2 - 1.0f)/(PI*log(a2)*t);
}

float GTR2 ( float NdotH, float a ) {
  float a2 = a*a;
  float t = 1.0f + (a2 - 1.0f)*NdotH*NdotH;
  return a2/(PI*t*t);
}

float GTR2_Aniso ( float NdotH, float HdotX, float HdotY, float ax, float ay ) {
  return 1.0f / ( PI * ax * ay * sqr(sqr(HdotX/ax) + sqr(HdotY/ay)
                                     + NdotH*NdotH));
}

float SmithG_GGX ( float NdotV, float alphaG ) {
  float a = alphaG*alphaG;
  float b = NdotV*NdotV;
  return 1.0f / (NdotV + sqrt(a + b - a*b));
}

float SmithG_GGX_Aniso ( float NdotV, float VdotX, float VdotY,
                         float ax, float ay ) {
  return 1.0f / (NdotV + sqrt(sqr(VdotX*ax) + sqr(VdotY*ay) + sqr(NdotV) ));
}


float3 Gamma_Correction ( float3 v ) {
  return pow(v, 2.2f);
}


float3 Disney_BRDF ( float3  L, float3 V, float3 N, float3 X, float3 Y,
                     Material* m ) {
  float NdotL = dot(N, L),
        NdotV = dot(N, V);
  if ( NdotL < 0.0f || NdotV < 0.0f ) return (float3)(0.0f);


  float3 H = normalize(L+V);
  float NdotH = dot(N, H),
        LdotH = dot(L, H);

  float3 Cdlin = Gamma_Correction(m->base_colour);
  float Cdlum = 0.3f*Cdlin.x + 0.6f*Cdlin.y + 0.1f*Cdlin.z; // emission approx
  // normalize lum to isolate hue+sat
  float3 Ctint = Cdlum > 0.0f ? Cdlin/Cdlum : (float3)(1.0f);
  float3 Cspec0 = mix(m->specular*0.08f*mix((float3)(1.0f), Ctint,
                                            m->specular_tint),
                      Cdlin, m->metallic);
  float3 Csheen = mix((float3)(1.0f), Ctint, m->sheen_tint);


  // Diffuse fresnel - 1 at normal incidence to 0.5 at grazing, mix in diffuse
  // retro reflection based on rougness
  float FL = Schlick_Fresnel(NdotL), FV = Schlick_Fresnel(NdotV);
  float Fd90 = 0.5f + 2.0f*LdotH*LdotH*m->roughness;
  float Fd = mix(1.0f, Fd90, FL) * mix(1.0f, Fd90, FV);

  // based on hanrahan krueger brdf approximation of isotropic bssrdf
  // 1.25 scale is to preserve albedo
  // fss90 used to flatten retroreflection based on roughness
  float Fss90 = LdotH*LdotH*m->roughness;
  float Fss = mix(1.0f, Fss90, FL) * mix(1.0f, Fss90, FV);
  float ss = 1.25f * (Fss * (1.0f / (NdotL + NdotV) - 0.5f) + 0.5f);

  // specular
  float aspect = sqrt(1.0f - m->anisotropic*0.9f);
  float ax = fmax(.001f, sqr(m->roughness)/aspect);
  float ay = fmax(.001f, sqr(m->roughness)*aspect);
  float Ds = GTR2_Aniso(NdotH, dot(H, X), dot(H, Y), ax, ay);
  float FH = Schlick_Fresnel(LdotH);
  float3 Fs = mix(Cspec0, (float3)(1.0f), FH);
  float Gs  = SmithG_GGX_Aniso(NdotL, dot(L, X), dot(L, Y), ax, ay);
        Gs *= SmithG_GGX_Aniso(NdotL, dot(V, X), dot(V, Y), ax, ay);

  // sheen
  float3 Fsheen = FH*m->sheen*Csheen;

  // clearcoat (ior = 1.5f -> F0 = 0.05)
  float Dr = GTR1(NdotH, mix(0.1f, 0.001f, m->clearcoat_gloss));
  float Fr = mix(0.04f, 1.0f, FH);
  float Gr = SmithG_GGX(NdotL, 0.25f) * SmithG_GGX(NdotV, 0.25f);

  return ((1.0f/PI) * mix(Fd, ss, m->subsurface)*Cdlin + Fsheen) *
         (1.0f-m->metallic) + Gs*Fs*Ds + 0.25f*m->clearcoat*Gr*Fr*Dr;
}

// -----------------------------------------------------------------------------
// --------------- MAP GEOMETRY FUNCTIONS --------------------------------------
float noise ( float2 n ) {
  const float2 d = (float2)(0.0f, 1.0f);
  float2 b = floor(n),
         f = smoothstep((float2)(0.0f), (float2)(1.0f), fract2f(n));
  return mix(mix(rand(b), rand(b + d.yx), f.x),
             mix(rand(b + d.xy), rand(b + d.yy), f.x),
             f.y);
}

float smin ( float a, float b ) {
  float k = -0.1f;
  float h = clamp(0.5f + 0.5f*(b - a)/k, 0.0f, 1.0f);
  return mix(b, a, h) - k*h*(1.0f - h);
}

float Boundary2 ( float3 p, float3 wal, float dist ) {
  return dot(p, wal) + dist;
}

float Boundary ( float3 p, float m) {
  return fmin(fmin(fmin(fmin(fmin(
    Boundary2(p, (float3)(1.0f, 0.0f, 0.0f),  m),
    Boundary2(p, (float3)(0.0f, 1.0f, 0.0f),  m)),
    Boundary2(p, (float3)(0.0f, 0.0f, 1.0f),  m)),
    Boundary2(p, (float3)(-1.0f, 0.0f, 0.0f), m)),
    Boundary2(p, (float3)(0.0f, -1.0f, 0.0f), m)),
    Boundary2(p, (float3)(0.0f, 0.0f, -1.0f), m));
}


float OP_Intersection ( float d1, float d2 ) {
  return max(-d2, d1);
}

float2 OP_Union ( float2 d1, float2 d2 ) {
  return mix(d2, d1, step(d1.x, d2.x));
}

float3 OP_Repeat ( float3 p, float3 c ) {
  return fmod(p, c) - 0.5f*c;
}

// -----------------------------------------------------------------------------
// --------------- MAP ---------------------------------------------------------
void Set_Map ( IntersectionInfo* pinfo, float dist, int mindex ) {
  if ( pinfo->dist > dist ) {
    pinfo->dist = dist;
    pinfo->material_index = mindex;
  }
}

float sdSphere ( float3 p, float r ) { return length(p) - r; }
float sdBox    ( float3 p, float3 b ) {
  float3 d = fabs(p) - b;
  return min(max(d.x, max(d.y, d.z)), 0.0f) + length(max(d, 0.0f));
}
float sdTorus  ( float3 p, float2 t ) {
  float2 q = (float2)(length(p.xz) - t.x, p.y);
  return length(q) - t.y;
}

float3 opRep ( float3 p, float3 c ) {
  return fmod(p, c) - 0.5f*c;
};

IntersectionInfo Map ( float3 p, float time ) {
  IntersectionInfo info;
  info.dist = FLT_MAX;

  Set_Map(&info, sdSphere(p + (float3)(-4.0f, 1.0f, 2.0f), 1.0f), 0);
  Set_Map(&info, Boundary(p, 10.0f), 1);
  Set_Map(&info, sdTorus(opRep(p, (float3)(8.0f, 8.0f, 8.0f)), (float2)(3.0f, 0.5f)), 2);

  return info;
}

// -----------------------------------------------------------------------------
// --------------- GRAPHIC FUNCS THAT NEED MAP ---------------------------------
float3 Normal ( float3 p, float time ) {
  const float Delta  = 0.001f,
              Delta2 = Delta*2.0f;

  const float3 X = (float3)(Delta, 0.0f, 0.0f),
               Y = (float3)(0.0f, Delta, 0.0f),
               Z = (float3)(0.0f, 0.0f, Delta);

  return (float3)(
    (Map(p + X, time).dist - Map(p - X, time).dist)/Delta2,
    (Map(p + Y, time).dist - Map(p - Y, time).dist)/Delta2,
    (Map(p + Z, time).dist - Map(p - Z, time).dist)/Delta2
  );
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

float3 BRDF ( float3 pos, float3 cam, float3 lpos, Material* material,
              float time ) {
  float3 N         = normalize(Normal(pos, time)),
         L         = normalize(lpos - pos),
         V         = normalize(cam - pos),
         tangent   = normalize(cross(Tangent(N), N)),
         bitangent = normalize(cross(N, tangent));
  float3 result = Disney_BRDF(L, V, N, tangent, bitangent, material);

  result *= dot(N, L);
  return result;
}

// -----------------------------------------------------------------------------
// --------------- RAYTRACING/MARCH --------------------------------------------
IntersectionInfo March ( Ray ray, float time ) {
  const float max_dist = 64.0f;
  float distance = 0.0f;
  IntersectionInfo t_info;
  for ( int i = 0; i < 64; ++ i ) {
    t_info = Map(ray.origin + ray.dir*distance, time);
    if ( t_info.dist < 0.05f || t_info.dist > max_dist ) break;
    distance += t_info.dist;
  }
  if ( t_info.dist > max_dist ) {
    t_info.dist = -1.0f;
    return t_info;
  }
  t_info.dist = distance;
  return t_info;
}

// uses cosine weight random hemisphere directions
float3 Hemisphere_Direction ( __global RNG* rng, float3 normal ) {
  float u1 = Uniform(rng,  0.0f, 1.0f),
        u2 = Uniform(rng,  0.0f, 1.0f);
  float theta = acos(sqrt(1.0f - u1));
  float phi   = 2.0f*PI*u2;

  float3 sdir = (fabs(normal.x) < 0.5f ? (float3)(1.0f, 0.0f, 0.0f) :
                                         (float3)(0.0f, 1.0f, 0.0f));
  float3 tdir = cross(normal, sdir);

  return theta*cos(phi)*sdir + theta*sin(phi)*tdir + theta*normal;
}

float3 Hemisphere_Weighted_Direction ( __global RNG* rng, float3 normal, float weight ) {
  float3 dir = Hemisphere_Direction(rng, normal);
  return normalize(dir + normal*weight);
}

float3 Orient_Normal ( float3 normal, float3 direction ) {
  return normal * (dot(normal, direction) < 0.0f ? -1.0f : 1.0f);
}

Ray TODO_BRDF_Reflect ( __global RNG* rng, Material material,
                           const IntersectionInfo info) {
  Ray rout;
  if ( Uniform(rng, 0.0f, 1.0) > material.metallic ) {
    // diffuse
    rout.dir = normalize(Hemisphere_Direction(rng, info.normal));
  } else {
    // specular/transmission
    float3 axis = Hemisphere_Weighted_Direction(rng, info.normal,
                          material.specular);
    if ( Uniform(rng, 0.0f, 1.0f) > material.anisotropic ) {
      //reflect
      rout.dir = normalize(info.dir - axis*dot(info.dir, axis)*2.0f);
    } else {
      //refract
      // TODO ! !
    }
  }
  rout.origin = info.origin + rout.dir*0.2f;
  return rout;
}

float3 TODO_BRDF_Radiance ( Material material, float3 radiance ) {
  return (radiance+material.emission)*material.base_colour;
}

float3 Raytrace ( __global RNG* rng, const __global Material* material,
                  Ray ray, int* hit_materials, float time ) {
  int depth = 0;
  Ray cray = ray;
  float terminate = 1.0f;
  // Trace a path from camera origin to some end
  while ( depth < MAX_DEPTH && Uniform(rng, 0.0f, 1.0f) <= terminate ) {
    IntersectionInfo info = March(ray, time);
    hit_materials[depth] = info.material_index;
    ++ depth;
    info.origin = ray.origin + ray.dir*info.dist;
    info.normal = Normal(info.origin, time);
    info.dir = ray.dir;
    Material mat = material[info.material_index];
    ray = TODO_BRDF_Reflect(rng, mat, info);
    terminate -= 0.1f + mat.emission; // TODO, fix me!
  }

  // Now we have a path, so we calculate the radiance from the light source
  float3 radiance = (float3)(0.0f, 0.0f, 0.0f);
  for ( int it = depth-1; it >= 0; -- it ) {
    radiance = TODO_BRDF_Radiance(material[hit_materials[it]], radiance);
    // [radiance *= ; TODO
  }
  return radiance;
}
// -----------------------------------------------------------------------------
// --------------- CAMERA ------------------------------------------------------
Ray Camera_Ray(__global RNG* rng, __global Camera* camera, float time) {
  float2 coord = (float2)((float)get_global_id(0), (float)get_global_id(1));
  float2 resolution = (float2)((float)camera->dim.x, (float)camera->dim.y);

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
  float3 ray_dir = normalize(puv.x*cam_right + puv.y*cam_up + 2.0f*cam_front);

  return New_Ray(cam_pos, ray_dir);
}

// -----------------------------------------------------------------------------
// --------------- KERNEL ------------------------------------------------------
__constant int SPP = 128;
__kernel void Kernel_Raycast(
        __write_only image2d_t output_image,
        __read_only  image2d_t input_image,
        __global RNG* rng,
        __global Camera* camera,
        __global float* time_ptr,
        __global Material* material, __global int* material_size
      ) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  // -- get old pixel, check if there are samples to be done
  //    (counter is stored in alpha channel)
  float4 old_pixel = read_imagef(input_image, out);
  uchar c = (uchar)old_pixel.w;
  if ( c >= SPP ) {
    return;
  }
  // -- set up camera and stack
  int hit_objects[MAX_DEPTH];

  // 50 samples . . .
  float3 result = (float3)(0.0f, 0.0f, 0.0f);
  float time = *time_ptr;
  Ray ray = Camera_Ray(rng, camera, *time_ptr);
  result += Raytrace(rng, material, ray, hit_objects, time);

  old_pixel += (float4)(result, 1.0f);

  write_imagef(output_image, out, old_pixel);


  // -- SUN GLARE :-) --
  // ray.dir += Uniform_Float3(rng, -2.5f, 2.5f);
  // float3 glare = Raytrace(ray, material, rng, time);
  // if ( glare.x > 1.5f ) {
  //   colour = glare/(float)SPP;
  //   write_imagef(output_image, out, (float4)(colour, 1.0f));
  // }
}};
