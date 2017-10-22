#define MAX_DEPTH 5
#define MARCH_DIST //%MARCH_DIST.0f
#define MARCH_REPS //%MARCH_REPS
#define MAT_LEN 4

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
  float3 colour, origin, normal;
  Material* material;
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

float Uniform_Sample ( SceneInfo* si ) {
  return Rand(si);
}
float2 Uniform_Sample2 ( SceneInfo* si ) {
  return (float2)(Uniform_Sample(si), Uniform_Sample(si));
}
float3 Uniform_Sample3 ( SceneInfo* si ) {
  return (float3)(Uniform_Sample(si), Uniform_Sample2(si));
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

// from PBRT 3rd ed
float2 Concentric_Sample_Disk(SceneInfo* si) {
  float2 u = Uniform_Sample2(si);
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
float3 Sample_Cosine_Hemisphere ( SceneInfo* si ) {
  const float2 d = Concentric_Sample_Disk(si);
  float phi = sqrt(fmax(0.0f, 1.0f - d.x*d.x - d.y*d.y));
  return (float3)(d.x, d.y, phi);
}
float3 Sample_Cosine_Hemisphere_N ( float3 N, SceneInfo* si ) {
  float3 wo = normalize(Sample_Cosine_Hemisphere(si));
  float3 binormal = (fabs(N.x) < 1.0f ? (float3)(1.0f, 0.0f, 0.0f) :
                                        (float3)(0.0f, 1.0f, 0.0f));
  binormal = normalize(cross(N, binormal));
  float3 bitangent = cross(binormal, N);
  return bitangent*wo.x + binormal*wo.y + wo.z*N;
}
float Cosine_Hemisphere_PDF ( float3 N, float3 wi ) {
  return fmax(0.0f, dot(N, wi)) * IPI;
}

float3 Sample_Cosine_Sphere ( SceneInfo* si ) {
  // TODO: improve this, it's awful!
  float3 cos_hemi = Sample_Cosine_Hemisphere(si);
  cos_hemi *= Uniform_Sample(si) > 0.5f ? 1.0f : -1.0f;
  return cos_hemi;
}
float Cosine_Sphere_PDF ( float3 N, float3 wi ) {
  return fmax(0.0f, fabs(dot(N, wi))) * IPI2;
}

//// from PBRT
float3 Uniform_Sample_Cone ( float2 u, float cos_theta_max ) {
  float cos_theta = (1.0f - u.x) + u.x*cos_theta_max,
        sin_theta = sqrt(1.0f - SQR(cos_theta));
  float phi = u.y*TAU;
  return To_Cartesian(cos(phi)*sin_theta, sin(phi)*sin_theta, cos_theta);
}

float Uniform_Cone_PDF ( float cos_theta_max ) {
  return 1.0f/(TAU * (1.0f - cos_theta_max));
}

void Calculate_Binormals ( float3 N, float3* T, float3* B ) {
  *T = (fabs(N.x) > fabs(N.y)) ? normalize((float3)(-N.z, 0.0f,  N.x)) :
                                 normalize((float3)(0.0f, N.z,  -N.y));
  *B = cross(N, *T);
}
// -----------------------------------------------------------------------------
// --------------- BSDF    FUNCS -----------------------------------------------
float3 BSDF_Sample ( float3 wi, float3 N, Material* m, SCENE_T(si, Tx)) {
  float3 diffuse = Sample_Cosine_Hemisphere_N(N, si);
  return diffuse;
}


float3 BSDF_F ( float3 wi, float3 wo, Material* m ) {
  return (float3)(IPI);
}

float BSDF_PDF ( float3 P, float3 wo, Material* m ) {
  return fabs(wo.z)*IPI;
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
  float  ds = sqr(length(vec));
  return fabs(dot(V0->normal, dn)) * fabs(dot(V1->normal, dn)) / ds;
}

Emitter Generate_Light_Subpath ( Subpath* path, float* light_weight,
                                 Spectrum* light_contrib, SCENE_T(si, Tx)) {
  // Grab a random light source to generate a subpath
  int light_index = Uniform_Sample(si)*EMITTER_AMT;
  Emitter light = REmission(light_index, si->debug_values, si->time);
  int mindex = EMIT_MAT - light_index;
  Ray ray;
  float3 light_normal;
  {// Generate ray position/angle
    float3 N = Sample_Cosine_Sphere(si);
    float3 origin = light.origin + light.radius*N;
    float3 dir = Sample_Cosine_Hemisphere_N(N, si);
    ray = (Ray){origin, dir};
    light_normal = N;
  }
  // generate path, Σ₀iᴿᴿ(Vᵢ + ωₒ*Rd)
  // where Rd is raymarch-distance, RR is russian-roulette, and
  // ωₒ = i=0 ? light-ray-dir : fᵢ(ωᵢ, N)
  // fᵢ = BSDF importance-sampling (Alternatively, you could sample a random
  //         point on hemisphere Ω, but fᵢ ensures f(P, ωᵢ, ωₒ) > 0.0)
  // Also collect vertex data along the way
  path->length = 1;
  // setup V0, point on the surface of sphere
  path->vertices[0].colour = light.emission;
  path->vertices[0].origin = ray.origin;
  path->vertices[0].normal = light_normal;
  path->vertices[0].material = NULL;

  { // --- set up V1 weight/contrib ---
    Vertex* V0 = path->vertices;
  }

  for ( int depth = 1; depth != 4; ++ depth ) {
    if ( Uniform_Sample(si) < 0.25f ) break;
    SampledPt ptinfo = March(mindex, ray, si, Tx);
    if ( ptinfo.dist < 0.0f )
      break;
    path->length += 1;
    // Vᵢ = Vᵢ₋₁ + ωᵢ*Rd
    float3 hit = ray.origin + ray.dir*ptinfo.dist;
    float3 normal = Normal(hit, si, Tx);
    Material* lt_material = si->materials + ptinfo.mat_index;
    {// set up vertex
      path->vertices[depth].colour = ptinfo.colour;
      path->vertices[depth].origin = hit;
      path->vertices[depth].normal = Normal(hit, si, Tx);
      path->vertices[depth].material = lt_material;
    }
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
      } else { // s>1
        Vertex* V0 = path->vertices + depth - 2,
              * V1 = path->vertices + depth - 1,
              * V2 = path->vertices + depth;
        float3 wi = normalize(V1->origin - V0->origin),
               wo = normalize(V2->origin - V1->origin);
        Spectrum bsdf_f = BSDF_F(wo, wi, V1->material)*V1->colour;
        float f_g   = Geometric_Term(V1, V2),
              b_g   = Geometric_Term(V1, V0),
              f_pdf = BSDF_PDF(wi, wo, V1->material) * f_g,
              b_pdf = BSDF_PDF(wo, wi, V1->material) * b_g;
        light_contrib[depth] = bsdf_f * b_g * light_contrib[depth-1];
        light_weight [depth] = (f_pdf/b_pdf) *  (light_weight[depth-1]+1.0f);
      }
    }

    path->length = 1;
    // W₀ = fᵢ(Wᵢ, N)
    float3 wo = BSDF_Sample(ray.dir, normal, lt_material, si, Tx);
    mindex = ptinfo.mat_index;
    ray.origin = hit;
    ray.dir = wo;
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
    V2.colour = (float3)(0.0f);
    V2.origin = ray.origin;
    V2.normal = ray.dir;
  }
  int mindex = -1;
  for ( int depth = 1; depth != 2; ++ depth ) {
    /* if ( Uniform_Sample(si) < 0.25f ) break; */
    SampledPt ptinfo = March(mindex, ray, si, Tx);
    if ( ptinfo.dist < 0.0f )
      break;

    // Vᵢ = Vᵢ₋₁ + ωᵢ*Rd
    float3 hit = ray.origin + ray.dir*ptinfo.dist;
    float3 normal = Normal(hit, si, Tx);

    // V0 –> V1 –> V2
    // V2 is current, but we're evaluating V1
    V0 = V1; V1 = V2;
    {// set up vertex
      V2.colour = ptinfo.colour;
      V2.origin = hit;
      V2.normal = Normal(hit, si, Tx);
    }
    // previous vertices
    {// --- calculate contribution and weights ---
      if ( depth == 1 ) {
        // what's defined; V1 -> V2
        // but there's no physical camera lens, so you can't do normal stuff
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
        float3 wi = normalize(V1.origin - V0.origin),
               wo = normalize(V2.origin - V1.origin);
        Spectrum bsdf_f = BSDF_F(wi, wo, V1.material)*V1.colour;
        float f_g = Geometric_Term(&V1, &V2),
              b_g = Geometric_Term(&V1, &V0),
              f_pdf = BSDF_PDF(wi, wo, V1.material) * f_g,
              b_pdf = BSDF_PDF(wo, wi, V1.material) * b_g;
        eye_contrib = bsdf_f * f_g * eye_contrib;
        eye_weight = (b_pdf/f_pdf) * (eye_weight + 1.0f);
      }
    }

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

    // setup next ray and mindex
    mindex = ptinfo.mat_index;
    ray.origin = hit;
    // W₀ = fᵢ(Wᵢ, N)
    ray.dir = BSDF_Sample(ray.dir, normal, V2.material, si, Tx);

    // --- calculate connection term ---
    /* if ( depth == 1 ) */
    /*   continue; // we're evaluating the "delta camera" which is pointless */

    for ( int light_it = 0; light_it < path.length; ++ light_it ) {
      Vertex* light_vertex = (path.vertices + light_it);
      Spectrum contribution;
      float vis = 0.0f;
      Vertex* Z2, * Z1, * Y1, * Y2;
      if ( light_it == 0 ) { // S == 1 connection strategy [vertex on light]
        // Y1 -> Z1 -> Z2 [Y2 doesn't exist]
        Z2 = &V1; Z1 = &V2; Y1 = light_vertex; Y2 = NULL;
        float3 wi = normalize(Z1->origin - Y1->origin),
               wo = normalize(Z2->origin - Z1->origin);
        Spectrum bsdf_f = BSDF_F(wi, wo, Z1->material) * Z1->colour;
        float con_g = Geometric_Term(Y1, Z1);
        contribution = eye_contrib * light_contrib[light_it+1] * con_g * bsdf_f;
      } else { // s >= 1, t >= 1 connection contrib
        // Y2 -> Y1 <--> Z1 -> Z2
        Y2 = (path.vertices + light_it - 1); Y1 = light_vertex;
        Z2 = &V1;                            Z1 = &V2;
        float3 wi = normalize(Y1->origin - Y2->origin),
               wo = normalize(Z1->origin - Y1->origin);
        Spectrum light_bsdf = BSDF_F(wi, wo, Y1->material) * Y1->colour;
        /*  */ wi = wo;
               wo = normalize(Z2->origin - Z1->origin);
        Spectrum eye_bsdf   = BSDF_F(wi, wo, Z1->material) * Z1->colour;
        float con_g = Geometric_Term(Y1, Z1);
        contribution = eye_contrib * light_contrib[light_it+1] * con_g *
                        light_bsdf * eye_bsdf;
      }
      // visibility ray
      vis = Visibility_Ray(Z1->origin, Y1->origin, si, Tx);
      // calculate MIS
      float MIS = 1.0f/( 1.0f + eye_weight + light_weight[light_it+1]);
      writefloat(light_weight[light_it+1]);
      // add contrib
      if ( sample_colour.x < 0.0f && vis > 0.0f )
        sample_colour = (float3)(0.0f);
      contribution = (float3)(
        fmax(0.0f, contribution.x),
        fmax(0.0f, contribution.y),
        fmax(0.0f, contribution.z)
      );
      sample_colour += MIS * (contribution * vis);
    }
  }
  sample_colour = (float3)(
    fmin(1.0f, sample_colour.x),
    fmin(1.0f, sample_colour.y),
    fmin(1.0f, sample_colour.z)
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
      // Laughably bad "raytracing" just to make it easy to navigate a scene
      SampledPt pt = March(-1, ray, &scene_info, textures);
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
