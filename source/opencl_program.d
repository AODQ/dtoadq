module opencl_program; immutable(string) Test_raycast_string = q{
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

IntersectionInfo Map ( float3 p, float time ) {
  IntersectionInfo info;
  info.dist = FLT_MAX;

  Set_Map(&info, sdSphere(p + (float3)(0.0f, 1.0f, 2.0f), 1.0f), 0);
  Set_Map(&info, Boundary(p, 10.0f), 1);

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
float3 Hemisphere_Direction ( float3 normal, __global RNG* rng ) {
  float u1 = Uniform(rng,  0.0f, 1.0f),
        u2 = Uniform(rng,  0.0f, 1.0f);
  float theta = acos(sqrt(1.0f - u1));
  float phi   = 2.0f*PI*u2;

  float x = sin(theta)*cos(phi),
        y = cos(theta),
        z = sin(theta)*sin(phi);

  float3 yvec = normal,
         xvec = normal,
         zvec;

  if (fabs(xvec.x)<=fabs(xvec.y) && fabs(xvec.x)<=fabs(xvec.z))
      xvec.x = 1.0f;
  else if (fabs(xvec.y)<=fabs(xvec.x) && fabs(xvec.y)<=fabs(xvec.z))
      xvec.y = 1.0f;
  else
      xvec.z = 1.0f;
  xvec = normalize(cross(xvec, y));
  zvec = normalize(cross(xvec, y));
  float3 dir = x*xvec + y*normal + z*zvec;
  return dir;
}

float3 Raytrace ( Ray ray, __global Material* material, __global RNG* rng,
                  float time ) {
  float3 colour = (float3)(0.0f), weight = (float3)(1.0f);
  int depth = 0;
  bool hit = false, material_hit = false;
  Material prev_material;
  float3 prev_pos;
  float3 cam_pos = ray.origin;
  while ( ++ depth ) {
    float depth_ratio = 1.0f-(1.0f/(float)(5 - depth));
    if ( Uniform(rng, 0.0f, 1.0f) >= depth_ratio ) {
      break;
    }
    IntersectionInfo info = March(ray, time);
    if ( info.dist < 0.0f ) {
      colour = (float3)(0.0f);
      hit = false;
      printf("MISSED?!\n");
      break;
    }

    Material m = *(material + info.material_index);

    if ( m.emission > FLT_MIN ) { // EMISSIVE
      float emittance = (m.emission/depth_ratio);
      if ( !material_hit )
        colour = (float3)(1.0f);
      else {
      colour = emittance *
        BRDF(
          prev_pos, ray.origin, ray.origin + ray.dir*info.dist,
          &prev_material, time
        )
      ;
      }
      hit = true;
      break;
    }

    material_hit = true;
    prev_pos = ray.origin + ray.dir*info.dist;
    prev_material = m;
    weight *= m.base_colour;
    ray.origin = prev_pos;
    ray.dir = Hemisphere_Direction(Normal(prev_pos, time), rng);
    ray.origin += ray.dir*0.2f;
  }

  if ( hit ) {
    return colour;
  }
  return (float3)(-1.0f);
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
__kernel void Kernel_Raycast(
        __write_only image2d_t output_image,
        __read_only  image2d_t input_image,
        __global RNG* rng,
        __global Camera* camera,
        __global float* time_ptr,
        __global ushort* converges, __global ushort* converges_size,
        __global Material* material, __global int* material_size
      ) {
  float time = *time_ptr;
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  Ray ray = Camera_Ray(rng, camera, time);

  float3 colour = Raytrace(ray, material, rng, time);

  // -- post effects --[B
  if ( colour.x > 0.0f ) {
    float gamma = 2.2f;
    colour = pow(colour, (float3)(2.0f/gamma));
    float4 old_colour = read_imagef(input_image, out) + (float4)(0.5f);
    ushort c = converges[out.y*camera->dim.x + out.x];
    if ( c < 5 ) {
      colour = mix(old_colour.xyz, colour, 1.0f/(float)(c));
      converges[out.y*camera->dim.x + out.x] += 1;
      if ( Is_Debug_Print() ) {
        printf("c: %f\n", 1.0f/(float)(c));
      }
      write_imagef(output_image, out, (float4)(colour, 1.0f));
    }
  }
}};
