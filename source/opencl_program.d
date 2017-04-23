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
// typedef struct T_Material {
//   float3 colour;
//   float metallic, subsurface, specular, roughness, specular_tint,
//         anisotropic, sheen, sheen_tint, clearcoat, clearcoat_gloss,
//         emission;
// } Material;

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

float PrimitiveA ( float3 point ) {
  return length(max(fabs(point) - (float3)(0.9f), 0.0f));
}


float PrimitiveB ( float3 point ) {
  float2 q = (float2)(length(point.xy) - 1.5f, point.z + cos(point.x));
  return length(q) - 0.9f;
}

float Terrain ( float3 pos ) {
  return length(pos) - 2.0f;
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

float Sun ( float3 p ) {
  float c = 2.0f;
  p.z = fmod(p.z, c) - 0.5f*c;
  p -= (float3)(1.0f, 3.0f, 0.0f);
  float d = length(p) - 1.0f;
  float3 q = fabs(p);
  float prism = fmax(q.z - 0.5f,
                     fmax(q.z*0.866025f + q.y*0.5f, q.y)
                     - 0.5f);
  return fmax(d, prism);
}

float smin ( float a, float b ) {
  float k = -0.1f;
  float h = clamp(0.5f + 0.5f*(b - a)/k, 0.0f, 1.0f);
  return mix(b, a, h) - k*h*(1.0f - h);
}

// --------------- INTERSECTION ------------------------------------------------
float Map ( float3 p ) {
  float r = PrimitiveB(p);
  r = smin(r, PrimitiveA(p));
  r = smin(r, Terrain(p + (float3)(0.0f, 1.2f, 0.0f)));
  r = fmin(r, Boundary(p));
  r = fmin(r, Sun(p));
  return r;
}

float3 Normal ( float3 p ) {
  const float Delta = 0.01f;

  float3 X = (float3)(Delta, 0.0f, 0.0f),
         Y = (float3)(0.0f, Delta, 0.0f),
         Z = (float3)(0.0f, 0.0f, Delta);

  return normalize((float3)(
    Map(p + X) - Map(p - X),
    Map(p + Y) - Map(p - Y),
    Map(p + Z) - Map(p - Z)
  ));
}


float sqr ( float t ) { return t*t; }

float smithG_GGX_aniso ( float NdotV, float VdotX, float VdotY, float ax, float ay ) {
  return 1.0f / (NdotV + sqrt(sqr(VdotX*ax) + sqr(VdotY*ay) + sqr(NdotV) ));
}

float Schlick_Fresnel ( float u ) {
  float m = clamp(1.0f - u, 0.0f, 1.0f);
  float m2 = m*m;
  return m2*m2*m; // pow(m, 5)
}


float GTR2_Aniso ( float NdotH, float HdotX, float HdotY, float ax, float ay ) {
  return 1.0f / ( 3.141579 * ax * ay * sqr(sqr(HdotX/ax) + sqr(HdotY/ay) + NdotH*NdotH));
}

float GTR1 ( float NdotH, float a ) {
  if ( a >= 1.0f ) return 1/(3.14159);
  float a2 = a * a;
  float t = 1.0f + (a2 - 1.0f)*NdotH*NdotH;
  return (a2 - 1.0f)/(3.14159*log(a2)*t);
}

float smithG_GGX ( float NdotV, float alphaG ) {
  float a = alphaG*alphaG;
  float b = NdotV*NdotV;
  return 1.0f / (NdotV + sqrt(a + b - a*b));
}

float3 disney_brdf ( float3  L, float3 V, float3 N, float3 X, float3 Y,
                     float3 base_colour,
                     float specular_tint,
                     float specular,
                     float metallic,
                     float sheen_tint,
                     float sheen,
                     float roughness,
                     float anisotropic,
                     float clearcoat,
                     float clearcoat_gloss,
                     float subsurface

) {
  float NdotL = dot(N, L),
        NdotV = dot(N, V);
  float3 H = normalize(L+V);

  float NdotH = dot(N, H),
        LdotH = dot(L, H);

  float3 Cdlin = (float3)(pow(base_colour.x, 2.2f), pow(base_colour.y, 2.2f),
                          pow(base_colour.z, 2.2f));
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
  float ax = fmax(.001f, (roughness*roughness)/aspect);
  float ay = fmax(.001f, (roughness*roughness)*aspect);
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

  return ((1.0f/3.14157f) * mix(Fd, ss, subsurface)*Cdlin + Fsheen) *
         (1.0f-metallic) + Gs*Fs*Ds + 0.25f*clearcoat*Gr*Fr*Dr;
}

float3 BRDF ( float3 pos, float3 cam, float3 lpos ) {
  float3 N = normalize(Normal(pos)),
         L = normalize(lpos - pos),
         V = normalize(cam - pos),
         tangent = normalize(cross((float3)(1.0f, 0.0f, 0.0f), N)),
         bitangent = normalize(cross(N, tangent));
  return disney_brdf(L, V, N, tangent, bitangent,
              (float3)(0.6f, 0.0f, 1.0f), // base material
              0.2f, // specular_tint
              0.2f, // specular
              0.0f, // metallic
              0.9f, // sheen tint
              0.9f, // sheen
              0.2f,  // roughness
              0.8f,  // anisotrpoic
              0.2f,  // clearcoat
              0.2f,  // clearcoat_gloss
              0.8f  // subsurface
  );
}

float March ( Ray ray ) {
  const float max_dist = 8.0f;
  float distance = 0.0f;
  float interval;
  for ( int i = 0; i < 32; ++ i ) {
    interval = Map(ray.origin + ray.dir*distance);
    if ( interval < 0.00001f || interval > max_dist ) break;
    distance += interval;
  }
  if ( interval > max_dist ) return 0.0f;
  return distance;
}

typedef struct T_IntersectionInfo {
  bool intersection;
  float3 position;
  float distance;
  float3 normal, angle;
  float3 colour;
} IntersectionInfo;

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

  mat3 rotatey = rotate_y(sin(time*0.52f)*0.81f),
       rotatex = rotate_x(sin(time*0.59f)*0.81f);

  float3 eye_pos = mat3_mul(rotatex,
                   mat3_mul(rotatey, (float3)(0.0f, 0.0f, -4.0f)));
  float3 up = (float3)(0.0f, 1.0f, 0.0f);
  float3 forward = mat3_mul(rotatey, (float3)(0.0f, 0.0f, 1.0f));
  float3 right = cross(up, forward);

  return New_Ray(eye_pos,
                 normalize(up * uv.y + right * uv.x + forward)
  );
}

// ------- kernel ---

__kernel void Kernel_Raycast(
        __write_only image2d_t output_image,
        __read_only  image2d_t input_image,
        __global RNG* rng,
        __global Camera* camera,
        __global float* time
      ) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  Ray ray = Camera_Ray(rng, camera, *time);


  // ray = Camera_Ray(rng, cam, ray.dir);

  if ( Is_Debug_Print() ) {
    printf("ORIGIN: <%f, %f, %f> DIR: <%f, %f, %f>\n",
              ray.origin.x, ray.origin.y, ray.origin.z,
              ray.dir.x, ray.dir.y, ray.dir.z);
  }

  float dist = March(ray);

  if ( dist > 0.0f ) {
    float3 point = ray.origin + ray.dir*dist;

    float3 colour = BRDF(
      point,
      ray.origin,
      (float3)(40.0f, 20.0f, 200.0f + cos(*time)*2.0f)
    );

    write_imagef(output_image, out, (float4)(colour, 1.0f));
  } else {
    write_imagef(output_image, out, (float4)(0.0f, 0.0f, 0.0f, 1.0f));
  }
}};

