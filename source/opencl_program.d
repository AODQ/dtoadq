module opencl_program; immutable(string) Test_raycast_string = q{
// --------------- DEBUG -------------------------------------------------------
bool Is_Debug_Print ( ) {
return get_global_id(0) == 1 && get_global_id(1) == 1;
}
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
// --------------- BASIC GRAPHIC FUNCTIONS -------------------------------------

float3 reflect ( float3 V, float3 N ) {
  return V - 2.0f*dot(V, N)*N;
}

__constant float PI = 3.1415926535f;
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
// --------------- MATERIAL/VOXEL/RAY ------------------------------------------

typedef struct T_Ray {
  float3 origin, dir, invdir;
  int3 sign;
} Ray;

Ray New_Ray ( float3 o, float3 d ) {
  Ray ray;
  ray.origin = o;
  ray.dir    = d;
  return ray;
}


// --------------- RANDOM ------------------------------------------------------
// // --- random generation via xorshift1024star
typedef struct RNG {
  ulong seed[16];
  ulong p;
} RNG;

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


float2 fract2f ( float2 vec ) {
  float2 itptr;
  return fract(vec, &itptr);
}

float fract1f ( float vec ) {
  float itptr;
  return fract(vec, &itptr);
}


float rand ( float2 n ) {
  return fract1f(sin(dot(n, (float2)(19.9898, 4.1414))) * 43758.5453);
}

float noise ( float2 n ) {
  const float2 d = (float2)(0.0f, 1.0f);
  float2 b = floor(n),
         f = smoothstep((float2)(0.0f), (float2)(1.0f), fract2f(n));
  return mix(mix(rand(b), rand(b + d.yx), f.x),
             mix(rand(b + d.xy), rand(b + d.yy), f.x),
             f.y);
}


// // https://pathtracing.wordpress.com/2011/03/03/cosine-weighted-hemisphere/
// float3 Random_Hemisphere_Direction(float3 normal, __global RNG* rng) {
//   int abort = 0;
//   float3 origin = (float3)(0.5f, 0.5f, 0.5f);
//   while ( true ) {
//     float3 v = normalize(Uniform(rng, 0.0f, 1.0f)) - origin;

//     if ( dot(normal, v) > 0.0f ) return v;
//     ++ abort;
//     if ( abort > 500 )
//       return v; // couldn't find a normal :/
//   }
// }

float smin ( float a, float b ) {
  float k = -0.1f;
  float h = clamp(0.5f + 0.5f*(b - a)/k, 0.0f, 1.0f);
  return mix(b, a, h) - k*h*(1.0f - h);
}

float PrimitiveA ( float3 point ) {
  return length(max(fabs(point) - (float3)(2.9f), 0.0f));
}

float Terrain ( float3 pos ) {
  return length(pos) - 2.0f;
}

float PrimitiveB ( float3 p ) {
  // float2 q = (float2)(length(point.xy) - 1.5f, point.z + cos(point.x));
  float r = smin(r, PrimitiveA(p));
  return smin(r, length(p + (float3)(0.0f, 0.0f,  1.0f)) - 1.0f);
}


float Boundary2 ( float3 p, float3 wal, float dist ) {
  return dot(p, wal) + dist;
}

float Boundary ( float3 p ) {
  return fmin(fmin(fmin(fmin(fmin(
    Boundary2(p, (float3)(1.0f, 0.0f, 0.0f), 5.0f),
    Boundary2(p, (float3)(0.0f, 1.0f, 0.0f), 5.0f)),
    Boundary2(p, (float3)(0.0f, 0.0f, 1.0f), 5.0f)),
    Boundary2(p, (float3)(-1.0f, 0.0f, 0.0f), 5.0f)),
    Boundary2(p, (float3)(0.0f, -1.0f, 0.0f), 5.0f)),
    Boundary2(p, (float3)(0.0f, 0.0f, -1.0f), 5.0f));
}

float Sun ( float3 p, float time ) {
  float c =  2.2325f + time*0.08f;
  float orig = p.z;
  p.y = fmod(p.y, c) - 0.5f*c;
  p.x = fmod(p.x, c) - 0.5f*c;
  return length(p + (float3)(sin(time)*0.02f,
                             cos(sqrt(time*2.5f))*0.51f,
                            -8.0f + cos(time*orig)*0.05f)) - 0.5f;
}


// --------------- INTERSECTION ------------------------------------------------
typedef struct T_IntersectionInfo {
  float dist;
  int material_index;
} IntersectionInfo;

void Set_Map ( IntersectionInfo* pinfo, float dist, int mindex ) {
  if ( pinfo->dist > dist ) {
    pinfo->dist = dist;
    pinfo->material_index = mindex;
  }
}

IntersectionInfo Map ( float3 p, float time ) {
  IntersectionInfo info;
  info.dist = FLT_MAX;
  Set_Map(&info, PrimitiveB(p), 0);

  Set_Map(&info, Sun(p, time), 1);
  return info;
}

IntersectionInfo Empty_IntersectionInfo ( ) {
  IntersectionInfo info;
  info.dist = 0.0f;
  info.material_index = -1;
  return info;
}

IntersectionInfo March ( Ray ray, float time ) {
  const float max_dist = 64.0f;
  float distance = 0.0f;
  IntersectionInfo t_info;
  for ( int i = 0; i < 64; ++ i ) {
    t_info = Map(ray.origin + ray.dir*distance, time);
    if ( t_info.dist < 0.05f || t_info.dist > max_dist ) break;
    distance += t_info.dist;
  }
  if ( t_info.dist > max_dist ) return Empty_IntersectionInfo();
  t_info.dist = distance;
  return t_info;
}

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

float smithG_GGX ( float NdotV, float alphaG ) {
  float a = alphaG*alphaG;
  float b = NdotV*NdotV;
  return 1.0f / (NdotV + sqrt(a + b - a*b));
}

float smithG_GGX_aniso ( float NdotV, float VdotX, float VdotY,
                         float ax, float ay ) {
  return 1.0f / (NdotV + sqrt(sqr(VdotX*ax) + sqr(VdotY*ay) + sqr(NdotV) ));
}


float3 mon2lin ( float3 v ) {
  return pow(v, 2.2f);
}


float3 disney_brdf ( float3  L, float3 V, float3 N, float3 X, float3 Y,
                     float3 base_colour,
                     float metallic,
                     float subsurface,
                     float specular,
                     float specular_tint,
                     float roughness,
                     float anisotropic,
                     float sheen,
                     float sheen_tint,
                     float clearcoat,
                     float clearcoat_gloss
) {
  float NdotL = dot(N, L),
        NdotV = dot(N, V);
  if ( NdotL < 0.0f || NdotV < 0.0f ) return (float3)(0.0f);


  float3 H = normalize(L+V);
  float NdotH = dot(N, H),
        LdotH = dot(L, H);

  float3 Cdlin = mon2lin(base_colour);
  float Cdlum = 0.3f*Cdlin.x + 0.6f*Cdlin.y + 0.1f*Cdlin.z; // luminance approx
  // normalize lum to isolate hue+sat
  float3 Ctint = Cdlum > 0.0f ? Cdlin/Cdlum : (float3)(1.0f);
  float3 Cspec0 = mix(specular*0.08f*mix((float3)(1.0f), Ctint, specular_tint),
                      Cdlin, metallic);
  float3 Csheen = mix((float3)(1.0f), Ctint, sheen_tint);


  // Diffuse fresnel - 1 at normal incidence to 0.5 at grazing, mix in diffuse
  // retro reflection based on rougness
  float FL = Schlick_Fresnel(NdotL), FV = Schlick_Fresnel(NdotV);
  float Fd90 = 0.5f + 2.0f*LdotH*LdotH*roughness;
  float Fd = mix(1.0f, Fd90, FL) * mix(1.0f, Fd90, FV);

  // based on hanrahan krueger brdf approximation of isotropic bssrdf
  // 1.25 scale is to preserve albedo
  // fss90 used to flatten retroreflection based on roughness
  float Fss90 = LdotH*LdotH*roughness;
  float Fss = mix(1.0f, Fss90, FL) * mix(1.0f, Fss90, FV);
  float ss = 1.25f * (Fss * (1.0f / (NdotL + NdotV) - 0.5f) + 0.5f);

  // specular
  float aspect = sqrt(1.0f - anisotropic*0.9f);
  float ax = fmax(.001f, sqr(roughness)/aspect);
  float ay = fmax(.001f, sqr(roughness)*aspect);
  float Ds = GTR2_Aniso(NdotH, dot(H, X), dot(H, Y), ax, ay);
  float FH = Schlick_Fresnel(LdotH);
  float3 Fs = mix(Cspec0, (float3)(1.0f), FH);
  float Gs  = smithG_GGX_aniso(NdotL, dot(L, X), dot(L, Y), ax, ay);
        Gs *= smithG_GGX_aniso(NdotL, dot(V, X), dot(V, Y), ax, ay);

  // sheen
  float3 Fsheen = FH*sheen*Csheen;

  // clearcoat (ior = 1.5f -> F0 = 0.05)
  float Dr = GTR1(NdotH, mix(0.1f, 0.001f, clearcoat_gloss));
  float Fr = mix(0.04f, 1.0f, FH);
  float Gr = smithG_GGX(NdotL, 0.25f) * smithG_GGX(NdotV, 0.25f);

  return ((1.0f/PI) * mix(Fd, ss, subsurface)*Cdlin + Fsheen) *
         (1.0f-metallic) + Gs*Fs*Ds + 0.25f*clearcoat*Gr*Fr*Dr;
}

float3 RGB_To_Float ( int r, int g, int b ) {
  return (float3)(r/255.0f, g/255.0f, b/255.0f);
}


typedef struct T_Material {
  float3 colour;
  float metallic, subsurface, specular, roughness, specular_tint,
        anisotropic, sheen, sheen_tint, clearcoat, clearcoat_gloss;
} Material;


// --------------- CAMERA ------------------------------------------------------
typedef struct T_Camera {
  float3 position, lookat, up;
  int2 dimensions;
  float fov;
} Camera;

Ray Camera_Ray(__global RNG* rng, __global Camera* cam, float time) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  float2 outf = (float2)((float)out.x, (float)out.y);
  float2 dim = (float2)((float)cam->dimensions.x, (float)cam->dimensions.y);
  float2 uv = ((2.0f * outf) - dim)/fmin(dim.x, dim.y);
  // uv.x += Uniform(rng, -0.001f, 0.001f);
  // uv.y += Uniform(rng, -0.001f, 0.001f);

  float3 eye_pos = (float3)(5.0f, 4.6f,  4.0f);
  float3 up = (float3)(0.0f, 1.0f, 0.0f);
  float3 forward = (float3)(0.0f, 0.1f, 1.0f);
  float3 right = cross(up, forward);

  return New_Ray(eye_pos,
                 normalize(up * uv.y + right * uv.x + forward)
  );
}


float3 BRDF ( float3 pos, float3 cam, float3 lpos, Material material, float time ) {
  float3 N = normalize(Normal(pos, time)),
         L = normalize(lpos - pos),
         V = normalize(cam - pos),
         tangent = normalize(cross((float3)(0.0f, 1.0f, 0.0f), N)),
         bitangent = normalize(cross(N, tangent));
  float3 result = disney_brdf(L, V, N, tangent, bitangent,
      material.colour,
      material.metallic,
      material.subsurface,
      material.specular,
      material.specular_tint,
      material.roughness,
      material.anisotropic,
      material.sheen,
      material.sheen_tint,
      material.clearcoat,
      material.clearcoat_gloss
  );

  result *= dot(N, L);
  return result;
}

// ------- kernel ---

__kernel void Kernel_Raycast(
        __write_only image2d_t output_image,
        __read_only  image2d_t input_image,
        __global RNG* rng,
        __global Camera* camera,
        __global float* time_ptr,
        __global Material* material, __global int* material_size
      ) {
  float time = *time_ptr;
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  Ray ray = Camera_Ray(rng, camera, time);

  IntersectionInfo info = March(ray, time);

  if ( info.dist > 0.0f ) {
    float3 point = ray.origin + ray.dir*info.dist;
    float3 colour;

    mat3 rotatey = rotate_y(sin(time*0.52f)*0.11f);
    colour = BRDF(
      point,
      ray.origin,
      mat3_mul(rotatey,
        (float3)(cos(time), -0.0f + sin(time)*2.0f, -2.0f)),
      material[info.material_index],
      time
    );

    colour += BRDF(
      point,
      ray.origin,
      mat3_mul(rotatey,
        (float3)(sin(time) - cos(time)*0.25, -1.0f + sin(time*0.25)*1.0f,
                 sin(time)+-2.0f)),
      material[info.material_index],
      time
    );

    float gamma = 2.2f;
    colour = pow(colour, (float3)(2.0f/gamma));
    write_imagef(output_image, out, (float4)(colour, 1.0f));
  } else {
    write_imagef(output_image, out, (float4)(0.0f, 0.0f, 0.0f, 0.0f));
  }
}};

