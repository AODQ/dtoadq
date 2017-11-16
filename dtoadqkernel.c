#define MAX_DEPTH 5
#define MARCH_DIST //%MARCH_DIST.0f
#define MARCH_REPS //%MARCH_REPS
#define MAT_LEN //%MAT_LENGTH

#define TEXTURE_T __read_only image2d_array_t
#define SCENE_T(S, T) SceneInfo* S, TEXTURE_T T
// EMIT_MAT = -120 - light_index
#define EMIT_MAT -120
#define REMIT(I) REmission(EMIT_MAT - I, si->debug_values, si->time)

#define DFLT2(V) ((V).x), ((V).y)
#define DFLT3(V) ((V).x), ((V).y), ((V).z)
#define DFLT4(V) ((V).x), ((V).y), ((V).z), ((V).w)

#define writeln(X)     if (Is_Debug()){printf(X "\n");                         }
#define writeint(X)    if (Is_Debug()){printf(#X " %d\n",                 (X));}
#define writeint2(X)   if (Is_Debug()){printf(#X " %d, %d\n",             (X));}
#define writeptr(X)    if (Is_Debug()){printf(#X " %p\n",                 (X));}
#define writefloat(X)  if (Is_Debug()){printf(#X " %f\n",                 (X));}
#define writefloat2(X) if (Is_Debug()){printf(#X " %f, %f\n",        DFLT2(X));}
#define writefloat3(X) if (Is_Debug()){printf(#X " %f, %f, %f\n",    DFLT3(X));}
#define writefloat4(X) if (Is_Debug()){printf(#X " %f, %f, %f, %f\n",DFLT4(X));}

// math macros
#define SQR(T) ((T)*(T))

//
#define Spectrum float3

__constant float MARCH_ACC = //%MARCH_ACC.0f/1000.0f;
__constant int   M_diffuse      = 1, M_glossy       = 2,
                 M_specular     = 3, M_transmittive = 4;
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

typedef struct T_SharedInfo {
  unsigned char clear_img;
  unsigned long finished_samples;
  unsigned char spp;
  uint2 rng_state;
} SharedInfo;

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
  float3 colour, origin, dir, normal;
  float dist;
  int mat_index;
} SampledPt;
SampledPt SampledPt_From_Origin(float3 origin) {
  SampledPt pt;
  pt.origin = origin;
  pt.dist = 0.0f;
  return pt;
}

typedef struct T_Vertex {
  float3 origin, albedo, normal, pdf_forward, pdf_backward;
  Material* material;
  int M_ID;
} Vertex;

typedef struct T_Subpath {
  Vertex vertices[MAX_DEPTH];
  uint length;
} Subpath;

typedef struct T_SceneInfo {
  float time;
  // __read_only image2d_array_t textures; IMAGES CANT BE USED AS FIELD TYPES:-(
  Material* materials;
  float3 debug_values;
  uint2 rng_state;
} SceneInfo;
SceneInfo New_SceneInfo(float time, Material* materials,
                        float3 debug_values, uint2 rng_state) {
  SceneInfo si;
  si.time         = time;
  si.materials    = materials;
  si.debug_values = debug_values;
  si.rng_state    = rng_state;
  return si;
}
// -----------------------------------------------------------------------------
// --------------- RANDOM FUNCTIONS       --------------------------------------
/*
  Using a high quality uniform RNG is very important for monte carlo, thus I
    use the MWC64X.
    http://cas.ee.ic.ac.uk/people/dt10/research/rngs-gpu-mwc64x.html

  The Warp Geneartor is, while superior, does not allow divergence, which makes
  it practically useless. For example, can't really know how much random numbers
  you'll need for a path you haven't generated yet of unknown length.
*/

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
// -----------------------------------------------------------------------------
// --------------- GENERAL FUNCTIONS      --------------------------------------

float Sample_Uniform ( SceneInfo* si ) {
  return Rand(si);
}
float2 Sample_Uniform2 ( SceneInfo* si ) {
  return (float2)(Sample_Uniform(si), Sample_Uniform(si));
}
float3 Sample_Uniform3 ( SceneInfo* si ) {
  return (float3)(Sample_Uniform(si), Sample_Uniform2(si));
}

float sqr(float t) { return t*t; }
float Distance(float3 u, float3 v) {
  float x = u.x-v.x, y = u.y-v.y, z = u.z-v.z;
  return sqrt(x*x + y*y + z*z);
}
float Distance_Sqr(float3 u, float3 v) {
  float x = u.x-v.x, y = u.y-v.y, z = u.z-v.z;
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
    Emitter e = REmission(i, si->debug_values, si->time);
    float dist = sdSphere(origin - e.origin, e.radius);
    MapUnionG(a, &res, dist, EMIT_MAT - i, light_emission*e.emission);
  }

  return res;
}

// -----------------------------------------------------------------------------
// --------------- RAYMARCHING      --------------------------------------------
SampledPt March ( int avoid, Ray ray, SCENE_T(si, Tx)) {
  float distance = 0.001f;
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

float3 RColour ( float3 pt_colour, Material* m ) {
  return (m->albedo.x < 0.0f) ? pt_colour : m->albedo;
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

// from PBRT 3rd ed
float2 Concentric_Sample_Disk(SceneInfo* si) {
  float2 u = Sample_Uniform2(si);
  // maps to [-1, 1]^2
  float2 offset = 2.0f * u - (float2)(1.0f, 1.0f);
  if ( offset.x == 0.0f && offset.y == 0.0f )
    return (float2)(0.0f);

  float theta, r;
  if ( fabs(offset.x) > fabs(offset.y) ) {
    r = offset.x;
    theta = (PI/4.0f) * (offset.y/offset.x);
  } else {
    r = offset.y;
    theta = (PI/2.0f) * (offset.x/offset.y);
  }
  return r * (float2)(cos(theta), sin(theta));
}

// from PBRT 3rd ed
float3 Sample_Cosine_Hemisphere ( SceneInfo* si, float* pdf ) {
  const float2 d = Concentric_Sample_Disk(si);
  float phi = sqrt(fmax(0.0f, 1.0f - d.x*d.x - d.y*d.y));
  pdf = phi/PI;
  return normalize((float3)(d.x, d.y, phi));
}
float3 Binormal ( float3 N ) {
  float3 axis = (fabs(N.x) < 1.0f ? (float3)(1.0f, 0.0f, 0.0f) :
                                    (float3)(0.0f, 1.0f, 0.0f));
  return normalize(cross(N, axis));
}
float Cosine_Hemisphere_PDF ( float3 N, float3 wi ) {
  return fmax(0.0f, dot(N, wi)) * IPI;
}

float3 Sample_Cosine_Sphere ( SceneInfo* si ) {
  // TODO: improve this, it's awful!
  float3 cos_hemi = Sample_Cosine_Hemisphere(si);
  cos_hemi *= Sample_Uniform(si) > 0.5f ? 1.0f : -1.0f;
  return cos_hemi;
}
float Cosine_Sphere_PDF ( float3 N, float3 wi ) {
  return fmax(0.0f, fabs(dot(N, wi))) * IPI2;
}

//// from PBRT
float3 Sample_Uniform_Cone ( float cos_theta_max, SceneInfo* si ) {
  float2 u = Sample_Uniform2(si);
  float cos_theta = (1.0f - u.x) + u.x*cos_theta_max,
        sin_theta = sqrt(1.0f - SQR(cos_theta));
  float phi = u.y*TAU;
  return To_Cartesian(cos(phi)*sin_theta, sin(phi)*sin_theta, cos_theta);
}
float Uniform_Cone_PDF ( float cos_theta_max ) {
  return 1.0f/(TAU * (1.0f - cos_theta_max));
}

// Probably only works on hemispheres
float3 Reorient_Angle ( float3 wi, float3 N ) {
  float3 binormal  = Binormal(N),
         bitangent = cross(binormal, N);
  return bitangent*wi.x + binormal*wi.y + wi.z*N;
}

// -----------------------------------------------------------------------------
// --------------- BSDF    FUNCS -----------------------------------------------
// --- some brdf utility functions .. ---
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

// -- diffuse --
float3 BRDF_Diffuse_Sample ( float3 wi, float3 N, SceneInfo* si ) {
  return Reorient_Angle(Sample_Cosine_Hemisphere(si), N);
}
float BRDF_Diffuse_PDF ( float3 wi, float3 wo, Material* mo ) {
  return fabs(wo.z)*IPI;
}

// -- specular --
float3 BRDF_Specular_Sample ( float3 wi, float3 N ) {
  return reflect(wi, N);
}
float BRDF_Specular_PDF ( float3 wi, float3 wo, Material* m ) {
  return 1.0f;
}

// -- glossy --
float3 BRDF_Glossy_Sample ( float3 wi, float3 N, float glossy_lobe,
                            SceneInfo* si) {
  return Reorient_Angle(Sample_Uniform_Cone(glossy_lobe, si), N);
}
float BRDF_Glossy_PDF ( float3 wi, float3 wo, Material* m ) {
  return Uniform_Cone_PDF(m->glossy_lobe);
}

// -- transmittive --
float3 BTDF_Sample ( float3 wi, float3 N, float ior ) {
  return refract(wi, N, ior);
}
float BTDF_PDF ( float3 wi, float3 wo, Material* m) {
  return 1.0;
}


/* The actual functions to call within the light transport algorithm. */
/*  Sample => Return a random sample from the material: */
/*              material id => if this hit is diffusive, glossy, specular or */
/*                             transmittive */
/*              bsdf        => the angle of the sample [wo] */
/*  PDF    => Probability Distribution Function of the material */
/*  F      => The 'luminosity' */

typedef struct T_BSDF_Sample_Results {
  int M_ID;
  float3 bsdf;
} BSDF_Sample_Results;
// Don't pass in a Vertex because it would be incomplete [while this
//    should compelte it, that's to the discretion of the user]
BSDF_Sample_Results BSDF_Sample ( float3 wi, float3 N, Material* m,
                                  SCENE_T(si, Tx)) {
  float transmittive_l = 0.0f,          transmittive_h = m->transmittive,
        glossy_l       = transmittive_h, glossy_h      = glossy_l  +m->glossy,
        specular_l     = glossy_h,      specular_h     = specular_l+m->specular;
        // diffuse takes up any remaining amount
  float pt = Sample_Uniform(si);
  float3 bsdf;
  int M_ID = 0;
  if       ( pt > transmittive_l && pt < transmittive_h ) {
    bsdf = BTDF_Sample(wi, N, m->ior);
    M_ID = M_transmittive;
  } else if ( pt > glossy_l && pt < glossy_h ) {
    bsdf = BRDF_Glossy_Sample(wi, N, m->glossy_lobe, si);
    M_ID = M_glossy;
  } else if ( pt > specular_l && pt < specular_h ) {
    writeln("specular");
    bsdf = BRDF_Specular_Sample(wi, N);
    M_ID = M_specular;
  } else {
    bsdf = BRDF_Diffuse_Sample(wi, N, si);
    M_ID = M_diffuse;
  }
  return (BSDF_Sample_Results){M_ID, bsdf};
}

float _BSDF_PDF ( float3 wi, float3 wo, float3 N, Material* m, int M ) {
  if( M==M_diffuse      ) return m->diffuse     *BRDF_Diffuse_PDF  (wi, wo, m );
  if( M==M_specular     ) return m->specular    *BRDF_Specular_PDF (wi, wo, m );
  if( M==M_glossy       ) return m->glossy      *BRDF_Glossy_PDF   (wi, wo, m );
  if( M==M_transmittive ) return m->transmittive*BTDF_PDF          (wi, wo, m );
  writeln("PDF???");
  return 0.0f;
}

float3 BRDF_F ( float3 wi, float3 wo, float3 N, Material* m, float3 pcolour ) {
  // unload material; I don't know why, but accessing material directly crashes
  // get proper L, V, binormal (X) and bitangent (Y)
  /* return pcolour; */
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
  return lumination;

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

float3 BTDF_F(float3 wi, float3 wo, float3 N, Material* m, float3 colour ) {
  return (1.0f/PI) * colour;
}

float3 _BSDF_F ( float3 wi, float3 wo, float3 N, Material* m, int M_ID,
                 float3 colour ) {
  if ( M_ID == M_transmittive )
    return BTDF_F(wi, wo, N, m, colour);
  return BRDF_F(wi, wo, N, m, colour);
}

float3 BSDF_F ( float3 wi, float3 wo, Vertex* V ) {
  return _BSDF_F ( wi, wo, V->normal, V->material, V->M_ID, V->colour );
}

float BSDF_PDF ( float3 wi, float3 wo, Vertex* V ) {
  return _BSDF_PDF ( wi, wo, V->normal, V->material, V->M_ID );
}

// -----------------------------------------------------------------------------
// --------------- LIGHT   FUNCS -----------------------------------------------
float Visibility_Ray(float3 orig, float3 other, SCENE_T(si, Tx)) {
  float theoretical = Distance(orig, other);
  float3 dir = normalize(other - orig);
  orig += dir*(0.001f);
  SampledPt ptinfo = March(-1, (Ray){orig, dir}, si, Tx);
  float actual = ptinfo.dist + MARCH_ACC + 0.1f;
  return 1.0f*(actual >= theoretical);
}

// -----------------------------------------------------------------------------
// --------------- LIGHT TRANSPORT ---------------------------------------------

/// Geometric term of e0 -> e1 no visibilty check, can be used to convert
///   a PDF w.r.t. SA to area by multiplying it by this factor
// ω -> A = PDF * |cosθ|/d²
// Adapted from PBRT 3d edition pg 1011
float Geometric_Term(Vertex* V0, Vertex* V1) {
  float3 vec = V1->origin - V0->origin;
  float3 dn = normalize(vec);
  return fabs(dot(V0->normal, dn)) * fabs(dot(V1->normal, dn));
}

float Geometric_Term_MIS(Vertex* V0, Vertex* V1) {
  float3 vec = V1->origin - V0->origin;
  float3 dn = normalize(vec);
  float  ds = sqr(length(vec));
  return fabs(dot(V0->normal, dn)) * fabs(dot(V1->normal, dn)) / ds;
}

Emitter Generate_Light_Subpath ( Subpath* path, float* light_weight,
                                 Spectrum* light_contrib, SCENE_T(si, Tx)) {
  // Grab a random light source to generate a subpath
  int light_index = Sample_Uniform(si)*EMITTER_AMT;
  Emitter light = REmission(light_index, si->debug_values, si->time);
  int mindex = EMIT_MAT - light_index;
  Ray ray;
  float pdf_forward;
  {// Generate ray position/angle and set vertex
    float3 N = Sample_Cosine_Sphere(si);
    float3 origin = light.origin + light.radius*N;
    float pdf_dir;
    float3 dir = Reorient_Angle(Sample_Cosine_Hemisphere(si, &pdf), N);

    ray = (Ray){origin, dir};
    pdf_forward = pdf_dir;

    Vertex* pv = path->vertices;
    pv->origin = origin;
    pv->albedo = light.emission;
    pv->normal = N;
    pv->pdf_forward = Cosine_Sphere_PDF(normal, wo);
    pv->material = NULL;

  }
  // generate path, Σ₀ᴿᴿ(Vᵢ + ωₒ*Rd)
  // where Rd is raymarch-distance, RR is russian-roulette,
  // ωₒ = i=0 ? light-ray-dir : fᵢ(ωᵢ, N)
  // fᵢ = BSDF importance-sampling (Alternatively, you could sample a random
  //         point on hemisphere Ω, but fᵢ ensures f(P, ωᵢ, ωₒ) > 0.0)
  // Also collect vertex data along the way
  path->length = 1;
  // setup V0, point on the surface of sphere
  Set_Vertex(path->vertices, light.emission, ray.origin, light_normal, NULL);

  { // --- set up V1 weight/contrib ---
    Vertex* V0 = path->vertices;
  }

  for ( int depth = 1; depth != 5; ++ depth ) {
    if ( Sample_Uniform(si) < 0.25f ) break;
    SampledPt ptinfo = March(mindex, ray, si, Tx);
    if ( ptinfo.dist < 0.0f )
      break;
    path->length += 1;
    // Vᵢ = Vᵢ₋₁ + ωᵢ*Rd
    float3 hit = ray.origin + ray.dir*ptinfo.dist;
    if ( ptinfo.mat_index < 0 ) break;
    float3 normal = Normal(hit, si, Tx);
    Material* lt_material = si->materials + ptinfo.mat_index;
    Set_Vertex(path->vertices + depth, RColour(ptinfo.colour, lt_material),
                hit, normal, lt_material);
    { // --- set up weight/contrib ---
      if ( depth == 1 ) { // s=1
        Vertex* V0 = path->vertices, * V1 = path->vertices + 1;
        float3 wo = normalize(V1->origin - V0->origin),
               wi = normalize(V0->origin - V1->origin);
        Spectrum bsdf_f = light.emission;
        float f_g   = Geometric_Term(V0, V1),
              b_g   = Geometric_Term(V1, V0),
              f_pdf = Cosine_Sphere_PDF(V0->normal, wo) * f_g,
              b_pdf = Cosine_Sphere_PDF(V0->normal, wi) * b_g;
        light_contrib[1] = bsdf_f * b_g * light_contrib[depth-1];
        light_weight[1] = (f_pdf/b_pdf) * (light_weight[depth-1]); /*XXX*/
        writefloat(light_weight[1]);
      } else { // s>1
        Vertex* V0 = path->vertices + depth - 2,
              * V1 = path->vertices + depth - 1,
              * V2 = path->vertices + depth;
        float3 wi = normalize(V1->origin - V0->origin),
               wo = normalize(V2->origin - V1->origin);
        Spectrum bsdf_f = BSDF_F(wi, wo, V1);
        float f_g   = Geometric_Term(V1, V2),
              b_g   = Geometric_Term(V1, V0),
              f_pdf = BSDF_PDF(wi, wo, V1) * f_g,
              b_pdf = BSDF_PDF(wo, wi, V1) * b_g;
        light_contrib[depth] = bsdf_f * b_g * light_contrib[depth-1];
        light_weight [depth] = (f_pdf*b_pdf) *  (light_weight[depth-1]);
      }
    }

    /* // W₀ = fᵢ(Wᵢ, N */
    BSDF_Sample_Results res = BSDF_Sample(ray.dir, normal, lt_material, si, Tx);
    mindex = ptinfo.mat_index;
    ray = (Ray){hit, res.bsdf};
    /* // M_ID is known b/c next vertex is known, though useless if path ends here */
    path->vertices[depth].M_ID = res.M_ID;
  }
  return light;
}

Spectrum BDPT_Integrate ( float3 pixel, float3 dir, SCENE_T(si, Tx)) {
  Subpath path;
  float    eye_weight,  light_weight [MAX_DEPTH];
  Spectrum eye_contrib, light_contrib[MAX_DEPTH];
  // α₀⁽ᴸᴱ⁾ = 1
  light_contrib[0] = eye_contrib = eye_weight = light_weight[0] = 1.0f;
  // instead of genearting a light and bsdf path, only the light path is
  // generated and the BSDF path is walked while evaluating the vertex behind it
  // to conserve GPU memory.
  Emitter light = Generate_Light_Subpath(&path, light_weight,
                                         light_contrib, si, Tx);
  Spectrum sample_colour = (Spectrum)(-1.0f, -1.0f, -1.0f);
  // --- generate eye path ---
  Ray ray = (Ray){pixel, dir};
  Vertex V0, V1, V2;
  {// set up V2 vertex [origin]
    Set_Vertex(&V2, (float3)(1.0f), ray.origin, (float3)(0.0f), NULL);
  }
  int mindex = -1;
  for ( int depth = 1; depth != 5; ++ depth ) {
    /* if ( Sample_Uniform(si) < 0.25f ) break; */
    SampledPt ptinfo = March(mindex, ray, si, Tx);
    if ( ptinfo.dist < 0.0f )
      break;

    // Even though evaluation is on V1, in the case that an emitter is hit,
    // there is enough information to compute this value. Then the loop must be
    // broken
    if ( ptinfo.mat_index <= EMIT_MAT ) { // S=0, T=n strategy
      if ( depth == 1 ) { // directly hit light from camera
        sample_colour = normalize(REMIT(ptinfo.mat_index).emission);
        break;
      }
      // just throw it out, it's such a small contribution that it's
      // not even worth computing
      break; // light sources don't reflect [though, they could]
    }

    // Vᵢ = Vᵢ₋₁ + ωᵢ*Rd
    float3 hit = ray.origin + ray.dir*ptinfo.dist;
    float3 normal = Normal(hit, si, Tx);

    // V0 –> V1 –> V2
    // V2 is current, but we're evaluating V1
    V0 = V1; V1 = V2;
    {// set up vertex
      Material* m = si->materials + ptinfo.mat_index;
      Set_Vertex(&V2, RColour(ptinfo.colour, m), hit, normal, m);
    }
    // previous vertices
    {// --- calculate contribution and weights ---
      if ( depth == 1 ) {
        // what's defined; V1 -> V2
        // This is ultimately meaningless. The biggest thing to understand
        // as that we're evaluating the lens here, which doesn't exist in the
        // scene. Even though we know the colour of V2, we do not have enough
        // information for its BSDF, so even though it seems that V2 is being
        // skipped, it's not.
        float3 wo = normalize(V2.origin - V1.origin),
               wi = normalize(V0.origin - V1.origin);
        Spectrum bsdf_f = 1.0f; // this is for lens, don't have one yet
        float f_g   = Geometric_Term(&V0, &V1), // CANT work w/o lens
              b_g   = Geometric_Term(&V1, &V0), // neither can this
              f_pdf = 1.0f * f_g, // no lens... pdf must be 1.0f
              b_pdf = 1.0f * b_g;
        eye_contrib = 1.0f * eye_contrib;
        eye_weight  = 1.0f;
      } else {
        // what's defined; V0 -> V1 -> V2
        // Even though it might not make senst to eval the colour of V1 in terms
        // of a connection [if V1 goes to the light source, V2 becomes
        // meaningless], this is handled in the connection code. This
        // information becomes relavent in cases such as V0 -> V1 -> V2 -> Y1
        float3 wi = normalize(V1.origin - V0.origin),
               wo = normalize(V2.origin - V1.origin);
        Spectrum bsdf_f = BSDF_F(wi, wo, &V1);
        float f_g   = Geometric_Term(&V1, &V2),
              b_g   = Geometric_Term(&V1, &V0),
              f_pdf = BSDF_PDF(wi, wo, &V1) * f_g,
              b_pdf = BSDF_PDF(wo, wi, &V1) * b_g;
        eye_contrib = bsdf_f * f_g * eye_contrib;
        eye_weight  = (b_pdf*f_pdf) * (eye_weight + 1.0f);
      }
    }

    // setup next ray and mindex
    mindex = ptinfo.mat_index;
    ray.origin = hit;
    {// W₀ = fᵢ(Wᵢ, N)
      BSDF_Sample_Results sres = BSDF_Sample(ray.dir, normal, V2.material,
                                             si, Tx);
      ray.dir = sres.bsdf;
      V2.M_ID = sres.M_ID  ;
    }

    // --- calculate connection term ---
    // Depth == 1 is the camera lens, which doesn't exist. The evaluations
    // would be occuring at the "middle of the screen" and therefore nearly
    // completely useless, not being worth the additional computation
    if ( depth == 1 ) continue;
    for ( int light_it = 0; light_it < path.length; ++ light_it ) {
      Vertex* light_vertex = (path.vertices + light_it);
      Spectrum contribution;
      float vis = 0.0f;
      Vertex* Z2, * Z1, * Y1, * Y2;
      if ( light_it == 0 ) { // S == 1 connection strategy
        // Y1 <--> Z1 -> Z2 [Y2 doesn't exist]
        Z2 = &V0; Z1 = &V1; Y1 = light_vertex; Y2 = NULL;
        float3 wi = normalize(Z1->origin - Y1->origin),
               wo = normalize(Z2->origin - Z1->origin);
        Spectrum bsdf_f = BSDF_F(wi, wo, Z1);
        float con_g = Geometric_Term(Y1, Z1);
        contribution = eye_contrib * light_contrib[light_it] * con_g * bsdf_f;
      } else {/*  // s >= 1, t >= 1 connection contrib */
        // Y2 -> Y1 <--> Z1 -> Z2
        Y2 = light_vertex - 1; Y1 = light_vertex;
        Z2 = &V0;              Z1 = &V1;
        float3 wi = normalize(Y1->origin - Y2->origin),
               wo = normalize(Z1->origin - Y1->origin);
        Spectrum light_bsdf = BSDF_F(wi, wo, Y1);
        wi = wo;
        wo = normalize(Z2->origin - Z1->origin);
        Spectrum eye_bsdf   = BSDF_F(wi, wo, Z1);
        float con_g = Geometric_Term_MIS(Y1, Z1);
        contribution = eye_contrib * light_contrib[light_it] * con_g *
                       light_bsdf * eye_bsdf;
      }
      // visibility ray
      vis += Visibility_Ray(Z1->origin, Y1->origin, si, Tx);
      // calculate MIS
      float MIS = 1.0f/(1.0f + eye_weight + light_weight[light_it]);
      // add contrib
      if ( sample_colour.x < 0.0f && vis > 0.0f )
        sample_colour = (float3)(0.0f);
      contribution = (float3)(
        fmax(0.0f, contribution.x),
        fmax(0.0f, contribution.y),
        fmax(0.0f, contribution.z)
      ) ;
      sample_colour += (contribution * vis)*MIS;
    }
  }
  // throw away results that go over 1
  if ( sample_colour.x > 1.0f || sample_colour.y > 1.0f || sample_colour.z > 1.0f )
    sample_colour = (float3)(-1.0f);
  sample_colour = (float3)(
    sample_colour.x,
    sample_colour.y,
    sample_colour.z
  );
  return sample_colour;
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
    // perform cheap raycast for sake of fast navigation.
    /* write_imagef(output, out, (float4)(hor.old_colour.xyz, 1.0f)); */
    shared_info->finished_samples = 0;
  }
  return hor;
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
    __global Material*          g_materials,
    __global float*             debug_val_ptr,
    __global uint2*             rng_states
  ) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  float spp = (float)(sinfo->spp);
  // --- construct camera --
  Camera camera = *camera_ptr;
  Ray ray = Camera_Ray(&camera);
  // pixel pos
  int pix_pt = out.y*camera.dim.x*4 + out.x*4;
  // -- set up scene info ---
  Material materials[MAT_LEN];
  for ( int i = 0; i != MAT_LEN; ++ i ) {
    materials[i] = *(g_materials + i);
    materials[i] = *(g_materials + i);
  }
  float time = *time_ptr;
  float3 dval = (float3)(debug_val_ptr[0], debug_val_ptr[1],
                          debug_val_ptr[2]);
  SceneInfo scene_info = New_SceneInfo(time, materials, dval,
                                       rng_states[out.y*camera.dim.x + out.x]);

  Update_Camera(&camera, time);

  // grav previous results (has to come after all initialization for fast nav
  //  raycast)
  float4 old_colour;
  {
    _EvalPreviousOutputHorcrux _hor =
          Eval_Previous_Output(img, output_img, out, sinfo, pix_pt);
    old_colour = _hor.old_colour;
    if ( old_colour.w >= spp ) return;
    if ( _hor.raycast_nav ) {
      // Laughably bad raycast to make it easy to navigate a scene
      SampledPt pt = March(-1, ray, &scene_info, textures);
      if ( pt.mat_index >= 0 )
        pt.colour = RColour(pt.colour, materials+pt.mat_index);
      pt.colour += pt.dist/MARCH_DIST;
      write_imagef(output_img, out, (float4)(pt.colour, 1.0f));
      return;
    }
  }


  // --- integrate ---
  Spectrum colour = BDPT_Integrate(ray.origin, ray.dir, &scene_info, textures);

  // --- store results ---
  // random
  rng_states[out.y*camera.dim.x + out.x] = scene_info.rng_state;

  // colour
  if ( colour.x >= 0.0f && colour.y >= 0.0f && colour.z >= 0.0f ) {
    old_colour = (float4)(mix(colour, old_colour.xyz,
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
  //
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
