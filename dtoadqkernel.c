#define MAX_EYE_DEPTH 5
#define MAX_EMI_DEPTH 2
#define MARCH_DIST //%MARCH_DIST.0f
#define MARCH_REPS //%MARCH_REPS
#define MAT_LEN //%MAT_LENGTH

#define TEXTURE_T __read_only image2d_array_t
#define SCENE_T(S, T) SceneInfo* S, TEXTURE_T T
// EMIT_MAT = -120 - light_index
#define EMIT_MAT -120
#define REMITM(I) (EMIT_MAT - (I))
#define REMIT(I) REmission(REMITM(I), si->debug_values, si->time)

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
#define Spectrum float3

__constant float MARCH_ACC = //%MARCH_ACC.0f/1000.0f;
__constant int   DO_NAVIGATION = 1;
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
  float3 origin, irradiance, normal, mat_col;
  float pdf_fwd, pdf_bwd;
  __constant Material* material;
  int M_ID;
} Vertex;

float Calc_Prob(Vertex* t) {
  float pdf_bwd = t->pdf_bwd,
        pdf_fwd = t->pdf_fwd;
  if ( pdf_bwd == 0.0f ) pdf_bwd = 1.0f;
  if ( pdf_fwd == 0.0f ) pdf_fwd = 1.0f;
  return pdf_bwd/pdf_fwd;
}

typedef struct T_Subpath {
  Vertex vertices[MAX_EMI_DEPTH];
  uint length;
} Subpath;

typedef struct T_SceneInfo {
  float time;
  // __read_only image2d_array_t textures; IMAGES CANT BE USED AS FIELD TYPES:-(
  __constant Material* materials;
  float3 debug_values;
  uint2 rng_state;
} SceneInfo;
SceneInfo New_SceneInfo(float time, __constant Material* materials,
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
// -1 .. 1
float2 Sample_Uniform2_2 ( SceneInfo* si ) {
  return ((float2)(Sample_Uniform(si), Sample_Uniform(si))-(float2)(0.5f))*2.0f;
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

float3 Float3_Max(float3 t, float g) {
  return (float3)(fmax(t.x, g), fmax(t.y, g), fmax(t.z, g));
}
float3 Float3_Min(float3 t, float g) {
  return (float3)(fmin(t.x, g), fmin(t.y, g), fmin(t.z, g));
}
float3 Float3_Clamp(float3 t, float a, float b) {
  return Float3_Max(Float3_Min(t, b), a);
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

float3 RColour ( float3 pt_colour, __constant Material* m ) {
  return (pt_colour.x >= 0.0f) ? pt_colour : m->albedo;
}

bool Is_Delta ( Vertex* vtx ) {
  return vtx->material == NULL ? 0 : ( vtx->material->specular     != 0.0f ||
                                       vtx->material->transmittive != 0.0f);
}

bool Is_Black ( Vertex* vtx ) {
  return vtx->irradiance.x == 0.0f &&
         vtx->irradiance.y == 0.0f &&
         vtx->irradiance.z == 0.0f;
}

float Map_PDF(float pdf) {
  return pdf == 0.0f ? 1.0f : pdf;
}

float3 reflect ( float3 V, float3 N ) {
  return V - 2.0f*dot(V, N)*N;
}

float3 refract(float3 I, float3 N, float ior) {
  float cos_NI = dot(N, I);
  float k = 1.0f - SQR(ior)*(1.0f - sqr(cos_NI));
  return k < 0.0f ? (float3)(0.0f) : ior*I - (ior*cos_NI + sqrt(k))*N;
}

float3 To_Cartesian ( float cos_theta, float phi ) {
  float sin_theta = sqrt(fmax(0.0f, 1.0f - cos_theta));
  return (float3)(cos(phi)*sin_theta, sin(phi)*sin_theta, cos_theta);
}

float3 Binormal ( float3 N ) {
  float3 axis = (fabs(N.x) < 1.0f ? (float3)(1.0f, 0.0f, 0.0f) :
                                    (float3)(0.0f, 1.0f, 0.0f));
  return normalize(cross(N, axis));
}

// Most sample functions from embree
float3 Sample_Cosine_Hemisphere ( SceneInfo* si, float* pdf ) {
  const float2 u = Sample_Uniform2(si);
  const float cos_theta = sqrt(u.y);
  *pdf = cos_theta * IPI;
  return To_Cartesian(cos_theta, TAU * u.x);
}
float PDF_Cosine_Hemisphere ( float3 wi, float3 N ) {
  return fabs(dot(wi, N)) * IPI;
}

// from Embree renderer (though, unlike PBRT's rather unique method of
// calculating a cosine hemisphere, I've seen this method in many places)
float3 Sample_Cosine_Sphere ( SceneInfo* si, float* pdf ) {
  float2 u = Sample_Uniform2(si);
  float phi = TAU * u.x,
        vv  = 2.0f*(u.y - 0.5f),
        cos_theta = sign(vv) * sqrt(fabs(vv));
  *pdf = 2.0f * fabs(cos_theta) * IPI;
  return To_Cartesian(sign(vv) * sqrt(fabs(vv)), phi);
}
float PDF_Cosine_Sphere ( float3 wi ) {
  return 2.0f * fabs(wi.z) * IPI;
}

float _PDF_Uniform_Cone ( float cos_theta_max ) {
  return 1.0f/(4.0f*PI*sqr(sin(0.5f*cos_theta_max)));
}
float3 Sample_Uniform_Cone ( float cos_theta_max, float* pdf, SceneInfo* si ) {
  float2 u = Sample_Uniform2(si);
  float phi = TAU*u.x,
        cos_theta = 1.0f - u.y*(1.0f - cos(cos_theta_max));
  *pdf = _PDF_Uniform_Cone(cos_theta_max);
  return To_Cartesian(cos_theta, phi);
}
float PDF_Uniform_Cone ( float3 wi, float3 N, float cos_theta_max ) {
  return _PDF_Uniform_Cone(cos_theta_max);
}

float Light_PDF ( Emitter* m ) {
  float pdf_amt = 1.0f/(EMITTER_AMT),
        pdf_pos = 1.0f/(2.0f*TAU*SQR(m->radius));
  return pdf_amt*pdf_pos;
}

float3 Reorient_Angle ( float3 wi, float3 N ) {
  float3 binormal  = Binormal(N),
         bitangent = cross(binormal, N);
  return bitangent*wi.x + binormal*wi.y + wi.z*N;
}

// -----------------------------------------------------------------------------
// --------------- BSDF    FUNCS -----------------------------------------------
// --- some brdf utility functions .. ---
float Schlick_Fresnel ( float u ) {
  float f = clamp(1.0f - u, 0.0f, 1.0f);
  float f2 = f*f;
  return f2*f2*f; // f^5
}

float Smith_G_GGX_Correlated ( float L, float R, float a ) {
  return L * sqrt(R - a*sqr(R) + a);
}

// -------- samples
float3 BRDF_Diffuse_Sample ( float3 wi, float3 N, float* pdf, SceneInfo* si ) {
  return Reorient_Angle(Sample_Cosine_Hemisphere(si, pdf), N);
}
float BRDF_Diffuse_PDF  ( float3 wi, float3 wo, float3 N ) {
  return PDF_Cosine_Hemisphere(wo, N);
}
float3 BRDF_Glossy_Sample ( float3 wi, float3 N, float* pdf, float glossy_lobe,
                            SceneInfo* si) {
  float3 NN = reflect(wi, N);
  return Reorient_Angle(Sample_Uniform_Cone(glossy_lobe, pdf, si), NN);
}
float BRDF_Glossy_PDF   ( float3 wi, float3 wo, float3 N, float cos_theta_max ) {
  if ( cos_theta_max < 0.001f ) return 1.0f;
  return PDF_Uniform_Cone(wo, N, cos_theta_max);
}
float3 BRDF_Specular_Sample ( float3 wi, float3 N, float* pdf ) {
  *pdf = 0.0f;
  return reflect(wi, N);
}
float BRDF_Specular_PDF ( float3 wi, float3 wo, float3 N ) {
  return 0.0f;
}
// Misnomer, it should be either BSDF or BTDF ; TODO make all pdf/sample/etc bsdf
float3 BRDF_Transmittive_Sample ( float3 wi, float3 N, float* pdf, float ior ) {
  *pdf = 0.0f;
  return refract(wi, N, ior*5.0f);
}
float BRDF_Transmittive_PDF ( float3 wi, float3 wo, float3 N, float ior ) {
  return 0.0f;
}


// Actual BRDF function that returns the albedo of the surface
Spectrum _BRDF_F ( float3 wi, float3 N, float3 wo, __constant Material* m,
                   float3 col ) {
  // get binormal, bitangent, half vec etc
  const float3 binormal  = Binormal(N),
               bitangent = cross(binormal, N),
               L         =  wo, V = -wi,
               H         = normalize(L+V);
  const float  cos_NV    = dot(N, V), cos_NL     = dot(N, L),
               cos_HV    = clamp(dot(H, V), 0.0f, 1.0f),
               cos_HL    = clamp(dot(H, L), 0.0f, 1.0f),
               Fresnel_L = Schlick_Fresnel(cos_NL),
               Fresnel_V = Schlick_Fresnel(cos_NV);
  // Diffusive component
  float3 diffusive_albedo = col*IPI;

  float3 microfacet = (float3)(1.0f);


  // probably transmittive
  if ( cos_NL <= 0.0f || cos_NV <= 0.0f ) {
    diffusive_albedo *= PI*m->transmittive;
    return diffusive_albedo;
  }

  { // ------- Fresnel
    // modified diffusive fresnel from disney, modified to use albedo & F0
    const float F0 = m->fresnel * m->metallic,
                Fresnel_diffuse_90 = F0 * SQR(cos_HL);
    float3 F = (1.0f - F0) * diffusive_albedo +
            mix(1.0f, Fresnel_diffuse_90, Fresnel_L) *
            mix(1.0f, Fresnel_diffuse_90, Fresnel_V);
    microfacet *= Float3_Clamp(F, 0.0f, 1.0f);
  }

  { // ------- Geometric
    // Heits 2014, SmithGGXCorrelated with half vec combined with anisotropic
    // term using GTR2_aniso model
    const float Param  = 0.5f + m->roughness,
                Aspect = sqrt(1.0f - m->anisotropic*0.9f),
                Ax     = Param/Aspect, Ay = Param*Aspect,
                GGX_NV = Smith_G_GGX_Correlated(cos_HL, cos_NV, Ax),
                GGX_HL = Smith_G_GGX_Correlated(cos_NV, cos_HL, Ay);
    float G = 0.5f / (GGX_NV*Ax + GGX_HL*Ay);
    microfacet *= Float3_Clamp(G, 0.0f, 1.0f);
  }
  { // ------- Distribution
    // Hyper-Cauchy Distribution using roughness and metallic
    const float Param = 1.2f + m->anisotropic,
                Shape = (1.1f - m->roughness*0.55f),
                tan_HL = length(cross(H, L))/cos_HL;
    const float Upper  = (Param - 1.0f)*pow(sqrt(2.0f), (2.0f*Param - 2.0f)),
                LowerL = (PI*sqr(Shape) * pow(cos_HL, 4.0f)),
                LowerR = pow(2.0f + sqr(tan_HL)/sqr(Shape), Param);
    float D = (Upper/(LowerL*LowerR));
    microfacet *= Float3_Clamp(D, 0.0f, 1.0f);
  }

  // Since microfacet is described using half vec, the following energy
  // conservation model may be used [Edwards et al. 2006]
  microfacet /= 4.00f * cos_HV * fmax(cos_NL, cos_NV);

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

  float3 result = Float3_Clamp((diffusive_albedo + microfacet), 0.0f, 1.0f);
  result = mix(result, col, m->specular+m->transmittive);
  return isnan(result.x) ? (float3)(0.0f) : result;
}

int Weighted_Uniform ( float4 weights, float rand, int amt ) {
  if ( rand < weights.x ) return 0;
  rand -= weights.x;
  if ( rand < weights.y ) return 1;
  rand -= weights.y;
  if ( rand < weights.z ) return 2;
  rand -= weights.z;
  if ( rand < weights.w ) return 3;

  writeln("WEIGHTED UNIFORM ERROR");
  writefloat(rand);
  writefloat4(weights);
  return 0;
}

// BSDF Sample helper function that sets sample dir, pdf and f
void _BSDF_Sample ( float3 wi, float3 N, Vertex* vtx,
                    float3* sample_dir, float* sample_pdf, float3* sample_f,
                    SceneInfo* si ) {
  // get sample_dir and pdf which is based off material brdf samples
  __constant Material* m = vtx->material;
  float4 weights = (float4){m->diffuse, m->glossy, m->specular, m->transmittive};
  int indx = Weighted_Uniform(weights, Sample_Uniform(si), 4);

  switch ( indx ) {
    case 0:
      *sample_dir = BRDF_Diffuse_Sample(wi, N, sample_pdf, si);
      *sample_pdf *= m->diffuse;
    break;
    case 1:
      *sample_dir = BRDF_Glossy_Sample(wi, N, sample_pdf, m->glossy_lobe, si);
      *sample_pdf *= m->glossy;
    break;
    case 2:
      *sample_dir = BRDF_Specular_Sample(wi, N, sample_pdf);
      *sample_pdf *= m->specular;
    break;
    case 3:
      *sample_dir = BRDF_Transmittive_Sample(wi, N, sample_pdf, m->ior);
      *sample_pdf *= m->transmittive;
    break;
  }

  // get sample_f
  *sample_f = _BRDF_F(wi, N, *sample_dir, m, vtx->mat_col);
}

// generates PDF, wo doesn't really do anything
float BSDF_PDF ( Vertex* vtx, float3 wi, float3 wo) {
  float pdf = 0.0f;
  __constant Material* m = vtx->material;
  float3 N = vtx->normal;

  pdf += m->diffuse      * BRDF_Diffuse_PDF     (wi, wo, N);
  pdf += m->glossy       * BRDF_Glossy_PDF      (wi, wo, N, m->glossy_lobe);
  pdf += m->specular     * BRDF_Specular_PDF    (wi, wo, N);
  pdf += m->transmittive * BRDF_Transmittive_PDF(wi, wo, N, m->ior);

  return pdf;
}

// Connect subpaths together
Spectrum Subpath_Connection ( Vertex* V0, Vertex* V1, Vertex* V2 ) {
  float gterm = fabs(dot(V1->normal, normalize(V2->origin - V1->origin)));
  float3 wi = normalize(V1->origin - V0->origin),
         wo = normalize(V2->origin - V1->origin);
  float pdf = BSDF_PDF(V1, wi, wo);
  Spectrum brdf = _BRDF_F(wi, V1->normal, wo, V1->material, V1->mat_col);
  if ( pdf == 0.0f )
    pdf = 1.0f;
  return brdf*gterm/pdf;
}

/* sets vtx with random sample from material
*/
void BSDF_Sample ( float3 O, float3 wi, float3 N, __constant Material* m,
                   Ray* ray, float3* irradiance, float* pdf,
                   Vertex* prev_vtx, Vertex* vtx, SceneInfo* si){
  // vertex info
  vtx->origin   = O;
  vtx->normal   = N;
  vtx->material = m;
  // sample info
  float3 sample_dir, sample_f;
  float sample_pdf;
  _BSDF_Sample(wi, N, vtx, &sample_dir, &sample_pdf, &sample_f, si);
  // edge info
  float3 edge = O - prev_vtx->origin;
  float edge_dist = length(edge);
  edge /= edge_dist; // normalize
  // vertex info
  vtx->pdf_fwd  = (*pdf) * dot(edge, prev_vtx->normal);
  vtx->irradiance = *irradiance;
  // prev vertex pdf bwd
  prev_vtx->pdf_bwd = BSDF_PDF(vtx, sample_dir, wi) * fabs(dot(-edge, N));
  // delta check [specular/transmittive]
  if ( Is_Delta(vtx) ) {
    prev_vtx->pdf_bwd = vtx->pdf_fwd = 0.0f;
  }
  // misc scope info
  sample_pdf = sample_pdf == 0.0f ? 1.0f : sample_pdf;
  *irradiance = sample_f * (*irradiance) * fabs(dot(N, sample_dir))/sample_pdf;
  *pdf    = sample_pdf;
  *ray    = (Ray){O, sample_dir};
}


// -----------------------------------------------------------------------------
// --------------- LIGHT   FUNCS -----------------------------------------------
SampledPt March_Visibility ( int avoid, Ray ray, SCENE_T(si, Tx)) {
  float distance = 0.005f;
  SampledPt t_info;
  for ( int i = 0; i < MARCH_REPS; ++ i ) {
    t_info = Map(avoid, ray.origin + ray.dir*distance, si, Tx);
    if ( t_info.dist < 0.00001f || t_info.dist > MARCH_DIST ) break;
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
float Visibility_Ray(float3 orig, float3 other, SCENE_T(si, Tx)) {
  float theoretical = Distance(orig, other);
  float3 dir = normalize(other - orig);
  orig += dir*(0.1f);
  SampledPt ptinfo = March_Visibility(-1, (Ray){orig, dir}, si, Tx);
  float actual = ptinfo.dist + MARCH_ACC + 0.1f;
  return 1.0f*(actual >= theoretical*0.9f);
}

// -----------------------------------------------------------------------------
// --------------- LIGHT TRANSPORT ---------------------------------------------
float _Geometric_Term_PDF_Glossy ( Vertex* V0, Vertex* V1, Vertex* V2 ) {
  if ( !V0 || !V2 ) return 0.0f;
  float3 wi = normalize(V1->origin - V0->origin),
         wo = normalize(V2->origin - V1->origin);
  float3 NN = reflect(wi, V1->normal);
  return (dot(wo, NN) > cos(V1->material->glossy_lobe))*1.0f;
}
float Geometric_Term_PDF ( Vertex* E0, Vertex* E1, Vertex* L1, Vertex* L0 ) {
  if ( E1->material != NULL && E1->material->glossy > 0.0f ) {
    return _Geometric_Term_PDF_Glossy(E0, E1, L1);
  }
  /* if ( L1->material != NULL && L1->material->glossy > 0.0f ) { */
  /*   return _Geometric_Term_PDF_Glossy(E1, L1, L0); */
  /* } */
  if ( Is_Delta(E1) || Is_Delta(L1) ) {
    return 0.0f;
  }
  return 1.0f;
}
/// Geometric term of e0 -> e1 w/ visibilty check, can be used to convert
///   a PDF w.r.t. SA to area by multiplying it by this factor
// ω -> A = PDF * |cosθ|/d²
// Adapted from PBRT 3d edition pg 1011
float Geometric_Term_Novis(Vertex* V0, Vertex* V1, SCENE_T(si, Tx)) {
  float3 d = V0->origin - V1->origin;
  float g = 1.0f/(SQR(length(d)));
  d *= sqrt(g);
  return g * fabs(dot(V0->normal, d)) * fabs(dot(V1->normal, d));
}
float Geometric_Term(Vertex* E0, Vertex* E1, Vertex* L1, Vertex* L0,
                     SCENE_T(si, Tx)) {
  float visible = Visibility_Ray(E1->origin, L1->origin, si, Tx);
  visible *= Geometric_Term_PDF(E0, E1, L1, L0);
  return Geometric_Term_Novis(E1, L1, si, Tx) * visible;
}

float PDF_Combine ( float pdf_bwd, float pdf_fwd ) {
  // dirac delta distribution handler
  return Map_PDF(pdf_bwd)/Map_PDF(pdf_fwd);
}

__constant int Eval_Vertex_Enum_Hit   = 0,
               Eval_Vertex_Enum_Miss  = 1,
               Eval_Vertex_Enum_Light = 2;
int Eval_Vertex ( Ray* ray, float3* irradiance, float* pdf,
                  Vertex* V0, Vertex* V1, SCENE_T(si, Tx) ) {
  // raymarch
  if ( V1->material && V1->material->transmittive > 0.0 )
    (*ray).origin += 0.1f*(*ray).dir;
  SampledPt ptinfo = March(-1, *ray, si, Tx);
  if ( ptinfo.dist < 0.0f ) return Eval_Vertex_Enum_Miss;
  // grab info
  float3 O = ray->origin + ray->dir * ptinfo.dist,
         N = Normal(O, si, Tx);
  V1->origin = O;
  V1->normal = N;
  // exit if hit light
  int mindex = ptinfo.mat_index;
  if ( mindex <= EMIT_MAT ) return mindex;
  __constant Material* mat = si->materials + ptinfo.mat_index;
  // set col
  float3 tcol = ptinfo.colour;
  V1->mat_col = RColour(tcol, mat);
  // sample
  BSDF_Sample(O, ray->dir, N, mat, ray, irradiance, pdf, V0, V1, si);
  // no point in continuing anymore if energy lost
  // though it might be wise to continue anyway to avoid branching
  /* if ( radiance->x == 0.0f && radiance->y == 0.0f && radiance->z == 0.0f ) */
  /*   return Eval_Vertex_Enum_Miss; */
  return Eval_Vertex_Enum_Hit;
}

float PDF3(Vertex* V0, Vertex* V1, Vertex* V2) {
  float3 wi, wo;
  wo = V2->origin - V1->origin; // use as edge temporary
  float edge_dist = length(wo);
  wo /= edge_dist;
  if ( V0 ) wi = normalize(V0->origin - V1->origin);

  float pdf;
  // null -> L1 [pt on surface] -> E1
  // cosΩ(|E₁-L₁|, L₁N) * |E₁ · ωₒ| / ||L₁ - E₁||²
  if ( !V0 ) pdf = PDF_Cosine_Hemisphere(-wo, V1->normal);
  else       pdf = BSDF_PDF(V1, wi, wo);


  return pdf * fabs(dot(V1->normal, wo))/(SQR(edge_dist));
}

float Calculate_MIS ( Vertex* L0, Vertex* L1, Vertex* E0, Vertex* E1,
                      float light_sum, float light_prob,
                      float eye_sum,   float eye_prob  ) {
  // calculate pdf bwd
  if ( E1 ) E1->pdf_bwd = PDF3(L0, L1, E1);
  if ( E0 ) E0->pdf_bwd = PDF3(L1, E1, E0);
  if ( L1 ) L1->pdf_bwd = PDF3(E0, E1, L1);
  if ( L0 ) L0->pdf_bwd = PDF3(E1, L1, L0);

  // calculate eye probability/mis with connection
  if ( E0 ) {
    eye_prob *= Calc_Prob(E0); eye_sum += SQR(eye_prob);
  }
  eye_prob *= Calc_Prob(E1); eye_sum += SQR(eye_prob);

  if ( L0 ) {
    light_prob *= Calc_Prob(L0); light_sum += SQR(light_prob);
  }
  light_prob *= Calc_Prob(L1); light_sum += SQR(light_prob);


  return 1.0f/(1.0f + eye_sum + light_sum);
}

Emitter Generate_Light_Subpath ( Subpath* path, float* light_sum,
                                 float* light_prob, SCENE_T(si, Tx)) {
  // Grab a random light source to generate a subpath
  int light_index = Sample_Uniform(si)*EMITTER_AMT;
  Emitter light = REmission(light_index, si->debug_values, si->time);
  int mindex = EMIT_MAT - light_index;
  Ray ray;
  float3 irradiance;
  float pdf_fwd;
  {// Generate ray position/angle, set initial vertex and set pdf_fwd
    float pdf_pos, pdf_dir;
    float3 N      = Sample_Cosine_Sphere(si, &pdf_pos),
           origin = light.origin + light.radius*N,
           dir    = Reorient_Angle(Sample_Cosine_Hemisphere(si, &pdf_dir), N);

    ray = (Ray){origin, dir};
    pdf_fwd = pdf_dir;

    float pdf_amt = 1.0f/EMITTER_AMT;

    Vertex* pv = path->vertices;
    pv->origin = origin;
    irradiance = pv->irradiance = light.emission/(pdf_amt*pdf_pos);
    pv->normal = N;
    pv->pdf_fwd = pdf_dir;
    pv->material = NULL;
    path->length = 1;
  }

  float t_prob = 1.0f, t_sum = 0.0f;

  light_prob[0] = 1.0f;

  // Generate path
  for ( int depth = 1; depth != MAX_EMI_DEPTH; ++ depth ) {
    Vertex* L0 = path->vertices+depth-1,
          * L1 = path->vertices+depth;
    int res = Eval_Vertex(&ray, &irradiance, &pdf_fwd, L0, L1, si, Tx);
    if ( res != Eval_Vertex_Enum_Hit ) break;
    // weights for L0 now that its bdf pwd is defined
    t_prob = light_prob[depth] = t_prob*Calc_Prob(L0);
    // delta handling
    if ( !Is_Delta(L0) && !Is_Delta(L1) )
      t_sum  = light_prob[depth] = t_sum+t_prob;
    path->length += 1;
  }
  return light;
}

Spectrum Connect_Light_Surface (Vertex* E0, Vertex* E1, Vertex* L1,
                                float t_sum, float t_prob, SCENE_T(si, Tx)) {
  Spectrum contribution = Subpath_Connection(E0, E1, L1);
  contribution = E1->irradiance * L1->irradiance;
  contribution *= Geometric_Term_Novis(E1, L1, si, Tx);
  // calculate MIS
  Vertex TE0 = *E0, TE1 = *E1;

  float l_prob = Calc_Prob(L1);
  float real_mis = Calculate_MIS(NULL, L1, &TE0, &TE1,
                                 l_prob, l_prob, t_sum, t_prob);
  return contribution*real_mis;
}

void Set_Sample_Colour(float3* sample_col, float3 con) {
  if ( sample_col->x < 0.0f && (con.x > 0.0f || con.y > 0.0f || con.z > 0.0f)) {
    *sample_col = (float3)(0.0f);
  }
  con = Float3_Max(Float3_Min(con, 1.0f), 0.0f);
  *sample_col += con;
}

Spectrum BDPT_Integrate ( float3 pixel, float3 dir, SCENE_T(si, Tx)) {
  Subpath path;
  float light_sum_arr[MAX_EMI_DEPTH];
  float light_prob_arr[MAX_EMI_DEPTH];
  // α₀⁽ᴸᴱ⁾ = 1
  // instead of generating a light and bsdf path, only the light path is
  // generated and the BSDF path is walked while evaluating the vertex behind it
  // to conserve GPU memory.
  Emitter light = Generate_Light_Subpath(&path, light_sum_arr,
                                         light_prob_arr, si, Tx);
  Spectrum sample_colour = (Spectrum)(-1.0f, -1.0f, -1.0f);
  // --- generate eye path ---
  Ray ray = (Ray){pixel, dir};
  Vertex E0, E1, E2;
  float3 irradiance = (float3)(1.0f);
  float pdf_fwd = 1.0f;
  { // Generate initial vertex and pdf_fwd
    E2.origin = ray.origin;
    E2.irradiance = (float3)(1.0f);
    E2.normal = ray.dir;
    E2.pdf_fwd = 1.0f;
    E2.material = NULL;
  }

  int mindex = -1;
  // Generate path, in non-gpu bdpt pseudo/D-code
  /* Connect ( eye_verts[2..$] , light_verts[1..$] ) */
  /* One immediate problem is that, since the eye path is connected during
    generation, and the next vertex needs to be known, the results are stored in
    E0, E1. E2 is the current evaluation vertex.
    where &#R = ready to render

       Eval Shift Eval Shift Eval Rendr       Shift  Eval   Rendr
   E0  NUL  NUL   NUL  &0R   &0R  &0R   . . . &N-2R  &N-1R  &N-1R
   E1  NUL  &0    &0R  &1    &1R  &1R         &N-1   &N-1R  &N-1R
   E2  &0   NUL   &1   NUL   &2   &2          NUL    &N     &N
  */
  float t_prob = 1.0f, t_sum = 0.0f;
  for ( int eye_depth = 1; eye_depth != MAX_EYE_DEPTH; ++ eye_depth ) {
    // shift edges
    E0 = E1;
    E1 = E2;
    {
    Vertex CL1;
    int CL1_connected = 0;
    {// eval path
        int result = Eval_Vertex(&ray, &irradiance, &pdf_fwd,
                                &E1, &E2, si, Tx);
        /* if ( pdf_fwd == 0.0f ) break; */
        if ( result == Eval_Vertex_Enum_Miss )
          break;
        if ( result <= EMIT_MAT ) {
          // direct light hit
          Emitter tlight = REMIT(result);
          if ( eye_depth == 1 ) {
            return (float3)(1.0f);//tlight.emission;
          }
          // E0 gaurunteed at this point
          CL1.origin = E2.origin;
          // cheat to get caustics to work better :-)
          CL1.irradiance = tlight.emission*Light_PDF(&tlight);
          CL1.normal = Normal(CL1.origin, si, Tx);
          CL1_connected = 2; // 2 = direct
        }
      }
      // Skip the first iteration because, although there is no camera lens,
      // E1 is still a point on the eye (just pixel). The actual first vertex
      // hasn't been evaluated yet, it needs to be set to E1!
      // In S=1 case, CL1 shouldn't even be defined
      if ( eye_depth == 1 ) continue;
      // connect to random point on light
      if ( !CL1_connected ) {
        float3 lorig = normalize((Sample_Uniform3(si)-0.5f)*2.0f);
        lorig = (lorig*light.radius*1.00f) + light.origin;
        float3 ldir = normalize(light.origin-E1.origin);
        SampledPt lres = March(-1, (Ray){E1.origin+ldir*0.005f, ldir},
                              si, Tx);
        CL1.origin = E1.origin + ldir*lres.dist;
        CL1.irradiance = light.emission/Light_PDF(&light);
        CL1.normal = Normal(CL1.origin, si, Tx);
        float pdf = Geometric_Term_PDF(&E0, &E1, &CL1, NULL);
        CL1_connected = (lres.mat_index != REMITM(light.index) || pdf <= 0.0f)?
                        -1:1;
      }
      if ( CL1_connected > 0 ) {
        Spectrum c = Connect_Light_Surface(&E0, &E1, &CL1,
                                           t_prob, t_sum, si,Tx);
        Set_Sample_Colour(&sample_colour, c);
        if ( CL1_connected == 2 ) return sample_colour; // direct light hit
      }
    }
    /* Perform connection, in pseudo/D-code: */
    /*    light_path.Each!( L => L.Connect(E1) ) */
    for ( int light_depth = 1; light_depth < path.length; ++ light_depth ) {
      Vertex* L1 = (path.vertices + light_depth),
            * L0 = light_depth == 0 ? NULL : (path.vertices + light_depth-1);
      // if delta continue

      Spectrum contribution = E1.irradiance * L1->irradiance;

      // As there is no T=1 case, L1 -> E1 -> E0 [Le in s=1 case] must happen
      contribution *= Subpath_Connection(L1, &E1, &E0);

      // For s == 1, L1 = L0, a point on the surface of emitter, there is no
      // L0 -> L1 -> E1 connection
      if ( light_depth > 0 ) {
        // Connect L0 -> L1 -> E1
        contribution *= Subpath_Connection( L0,  L1, &E1);
      }

      // Geometric connection term [includes visibility check]
      contribution *= Geometric_Term(&E0, &E1, L1, L0, si, Tx);

      // Connection failed (ei, Specular->Diffuse path didn't match or
      //    dot(L, V) < 0.0f etc)
      if ( contribution.x == 0.0f && contribution.y == 0.0f &&
           contribution.z == 0.0f ) {
        continue;
      }

      // ------------ Breaking the connectioning ------------
      // Create temporary as their bwd pdfs are overwriten in connection
      // strategy
      // Not necessary for TE1 but still done for consistency
      Vertex TL1, TL0, TE0, TE1;
      TE0 = E0; TE1 = E1;
      TL1 = *L1;
      if ( L0 ) TL0 = *L0;

      float real_mis = 0.0f;

      real_mis =
        Calculate_MIS(L0?&TL0:NULL, &TL1, &TE0, &TE1,
                      light_sum_arr[light_depth], light_prob_arr[light_depth],
                      t_sum, t_prob);

      // add contrib
      if ( sample_colour.x < 0.0f &&
              (contribution.x > 0.0f || contribution.y > 0.0f ||
               contribution.z > 0.0f) ) {
        sample_colour = (float3)(0.0f);
      }
      sample_colour += Float3_Min(contribution*real_mis, 1.0f);
    }
    // Calculate eye prob/mis after light calculation to avoid using a buffer
    /*
      depth| 0    1    2    3   4  ... N
      -----|----------------
      prob | X    X   1.0f  V0  V0    V0
           | X    X    X    X   V1    ..
           | X    X    X    X   X     VN

    */
    t_prob *= Calc_Prob(&E0);
    // Add MIS only for non-delta distribution
    if ( !Is_Delta(&E0) && !Is_Delta(&E1) )
      t_sum  += t_prob;
  }
  return Float3_Clamp(sample_colour, 0.0f, 1.0f);
}

// -----------------------------------------------------------------------------
// --------------- CAMERA ------------------------------------------------------
// from gllookat
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
  /* float3 dof_puv = (float3)(puv, fov_r); */
  /* float3 ray_dir = dof_puv.x*cam_right + dof_puv.y*cam_up + fov_r*cam_front; */
  /* float3 dof_origin = camera->radius*(float3)(Sample_Uniform2_2(si), 0.0f); */
  /* float3 dof_dir = normalize(dof_puv * camera->focal - dof_origin); */
  /* cam_pos += dof_origin.x*cam_right + dof_origin.y*cam_up; */
  /* ray_dir += dof_dir.x*cam_right + dof_dir.y*cam_up; */
  /* ray_dir = normalize(ray_dir); */
  /* return (Ray){cam_pos, ray_dir}; */
  return (Ray){cam_pos, normalize(puv.x*cam_right + puv.y*cam_up +
                                  fov_r*cam_front)};
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
    __constant Material*        g_materials,
    __global float*             debug_val_ptr,
    __global uint2*             rng_states
  ) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  float spp = (float)(sinfo->spp);
  // --- construct camera --
  Camera camera = *camera_ptr;
  float time = *time_ptr;
  // pixel pos
  int pix_pt = out.y*camera.dim.x*4 + out.x*4;
  // -- set up scene info ---
  float3 dval = (float3)(debug_val_ptr[0], debug_val_ptr[1],
                          debug_val_ptr[2]);
  SceneInfo scene_info = New_SceneInfo(time, g_materials, dval,
                                       rng_states[out.y*camera.dim.x + out.x]);
  if ( Sample_Uniform(&scene_info) < 0.3f ) {
    rng_states[out.y*camera.dim.x + out.x] = scene_info.rng_state;
    return;
  }
  // -- update camera / generate ray --
  Update_Camera(&camera, time);
  Ray ray = Camera_Ray(&camera, &scene_info);

  // grav previous results (has to come after all initialization for fast nav
  //  raycast)
  float4 old_colour;
  {
    _EvalPreviousOutputHorcrux _hor =
          Eval_Previous_Output(img, output_img, out, sinfo, pix_pt);
    old_colour = _hor.old_colour;
    if ( (int)(old_colour.w) == spp ) {
      sinfo->finished_samples += 1;
      old_colour.w += 1;
    }
    if ( old_colour.w >= spp ) return; /*  */
    if ( _hor.raycast_nav && DO_NAVIGATION ) {
      // Laughably bad raycast to make it easy to navigate a scene
      SampledPt pt = March(-1, ray, &scene_info, textures);
      if ( pt.mat_index >= 0 )
        pt.colour = RColour(pt.colour, g_materials+pt.mat_index);
      pt.colour = (float3)((pt.colour.x+pt.colour.y+pt.colour.z)/3.0f);
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
}
