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
} SharedInfo;

__constant float PI   = 3.141592653589793f;
__constant float IPI  = 0.318309886183791f;
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
  float rng;
} SceneInfo;
SceneInfo New_SceneInfo(float time, Material* materials,
                        float3 debug_values, float rng){
  SceneInfo si;
  si.time         = time;
  si.materials    = materials;
  si.debug_values = debug_values;
  si.rng          = rng;
  return si;
}
// -----------------------------------------------------------------------------
// --------------- GENERAL FUNCTIONS      --------------------------------------
float Rand ( SceneInfo* si ) {
  float t;
  float val = fract(sin(dot((float2)(get_global_id(0), get_global_id(1)),
                            (float2)(si->rng*100.0f+2.8f,
                                     si->rng*100.0f+4.3f))
                       )*48561.4982f, &t);
  si->rng = val;
  return val;
}

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

float2 Sample_Concentric_Disk ( float2 u ) {
  float2 offset = 2.0f*u - (float2)(1.0f);
  if ( offset.x == 0.0f && offset.y == 0.0f ) return (float2)(0.0f);

  float theta, r;
  if ( fabs(offset.x) > fabs(offset.y) ) {
    r = offset.x;
    theta = (PI/4.0f)*(offset.y/offset.x);
  } else {
    r = offset.y;
    theta = (PI/2.0f) - (PI/4.0f)*(offset.x/offset.y);
  }
  return r*(float2)(cos(theta), sin(theta));
}

float3 Sample_Cosine_Hemisphere ( float2 u ) {
  float2 d = Sample_Concentric_Disk(u);
  float z = sqrt(fmax(0.0f, 1.0f - SQR(d.x) - SQR(d.y)));
  return (float3)(d.x, d.y, z);
}

float Cosine_Hemisphere_PDF ( float3 N, float3 wi ) {
  return fmax(0.0f, dot(N, wi)) * IPI;
}
// Cosine weighted
float3 Hemisphere_Direction ( SceneInfo* si ) {
  const float phi = TAU*Uniform_Sample(si),
              vv  = 2.0f*(Uniform_Sample(si) - 0.5f);
  const float cos_theta = sign(vv)*sqrt(fabs(vv)),
              sin_theta = sqrt(fmax(0.0f, 1.0f - (cos_theta*cos_theta)));
  return To_Cartesian(sin_theta, cos_theta, phi);
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
float3 BSDF_Sample ( float3 wi, float3 normal, Material* m, SCENE_T(si, Tx)) {
  float3 wo;
  wo = normalize(Hemisphere_Direction(si));
  // Now I have to orient it with the normal . . .
  // whatever
  return wo;
}


float3 BSDF_F ( float3 wi, float3 wo, Material* m ) {
  return (float3)(1.0f);
}

float BSDF_PDF ( float3 P, float3 wo, Material* m ) {
  return fabs(wo.z)*IPI;
}

// -----------------------------------------------------------------------------
// --------------- LIGHT   FUNCS -----------------------------------------------
float Visibility_Ray(float3 orig, float3 other, SCENE_T(si, Tx)) {
  float theoretical = Distance(orig, other);
  float3 dir = normalize(other - orig);
  orig += dir*(MARCH_ACC*2.0f + 0.2f);
  SampledPt ptinfo = March(-1, (Ray){orig, dir}, si, Tx);
  float actual = ptinfo.dist + MARCH_ACC + 0.2f;
  return 1.0f*(actual >= theoretical);
}

/*
  The PDF from a light source on the surface of a light emitter. Since, as of
  now, there is only support for one light source, the PDF is a constant.
  Consider the surface area of a sphere is (2.0f*tau*r^2). Thus Pr(X <= x) =
  x/(SA). p(x) = 1/(SA). Then to check, the integral from 0 to SA of p(x)dx = 1.
*/
float Light_PDF ( float3 P, Emitter* m ) {
  float rad_sqr = m->radius*m->radius;
  return 1.0f/(2.0f*TAU*rad_sqr);
}

float3 Light_Sample_Li ( float3 P, float3* wi, float* pdf, Emitter* m,
                         SCENE_T(si, Tx)) {
  float3 O = m->origin;
  float theta = 2.0f * Uniform_Sample(si) * PI,
        phi   = PI * Uniform_Sample(si);
  O += (float3)(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi))
       * (2.0f * (Uniform_Sample(si)-0.5f)) * m->radius;
  *wi = O - P;
  if ( length(*wi) == 0.0f ) *pdf = 0.0f;
  else {
    *wi = normalize(*wi);
    float sin_max = m->radius*m->radius/Distance_Sqr(O, P);
    *pdf = Light_PDF(P, m);
    if ( *pdf == INFINITY ) *pdf = 0.0f;
  }
  return O;
}


// -----------------------------------------------------------------------------
// --------------- LIGHT TRANSPORT ---------------------------------------------

/*
  typedef struct T_Vertex {
  float3 colour, origin, normal;
  int mat_index;
  } Vertex;

  typedef struct T_Subpath {
  Vertex vertices[MAX_DEPTH];
  uint length;
  } Subpath;
 */
Emitter Generate_Light_Subpath ( Subpath* path, SCENE_T(si, Tx)) {
  // Grab a random light source to generate a subpath
  int light_index = Uniform_Sample(si)*EMITTER_AMT;
  Emitter light = REmission(light_index, si->debug_values, si->time);
  int mindex = EMIT_MAT - light_index;
  Ray ray;
  {// Generate ray position/angle
    float3 origin =
      (normalize(Uniform_Sample3(si)) - (float3)(0.5f)) * (float3)(2.0f);
    origin = light.origin + origin*light.radius;
    ray = (Ray){origin, Uniform_Sample3(si)};
  }
  // generate path, Σ₀iᴿᴿ(Vᵢ + ωₒ*Rd)
  // where Rd is raymarch-distance, RR is russian-roulette, and
  // ωₒ = i=0 ? light-ray-dir : fᵢ(ωᵢ, N)
  // fᵢ = BSDF importance-sampling (Alternatively, you could sample a random
  //         point on hemisphere Ω, but fᵢ ensures f(P, ωᵢ, ωₒ) > 0.0)
  // Also collect vertex data along the way
  path->length = 1;
  // setup V0
  path->vertices[0].colour = light.emission;
  path->vertices[0].origin = ray.origin;
  path->vertices[0].normal = Normal(ray.origin, si, Tx);
  path->vertices[0].material = NULL;
  for ( int depth = 1; depth != 2; ++ depth ) {
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

    // W₀ = fᵢ(Wᵢ, N)
    float3 wo = BSDF_Sample(ray.dir, normal, lt_material, si, Tx);
    mindex = ptinfo.mat_index;
    ray.origin = hit;
    ray.dir = wo;
  }
  path->length = 0;
  return light;
}

/// Geometric term of e0 -> e1 no visibilty check, can be used to convert
///   a PDF w.r.t. SA to area
float Geometric_Term(Vertex* V0, Vertex* V1) {
  float3 w = V1->origin - V0->origin;
  float inv_dist_sqr = 1.0f / sqr(length(w));
  return inv_dist_sqr * fabs(dot(V1->normal, w*sqrt(inv_dist_sqr)));
}

Spectrum BDPT_Integrate ( float3 pixel, float3 dir, SCENE_T(si, Tx)) {
  Subpath light_path;
  // instead of genearting a light and bsdf path, only the light path is
  // generated and the BSDF path is walked while evaluating the vertex behind it
  // to conserve GPU memory.
  Emitter light = Generate_Light_Subpath(&light_path, si, Tx);
  // generate weights and PDFs
  float    eye_weight,  light_weight [MAX_DEPTH];
  Spectrum eye_contrib, light_contrib[MAX_DEPTH];
  // α₀⁽ᴸᴱ⁾ = 1
  light_contrib[0] = eye_contrib = eye_weight = light_weight[0] = 1.0f;

  if ( light_path.length > 1 ) { // α₁
    // -------- unweighted contribution --------
    // αᴸ₁ = Lₑ⁰(y₀ –> y₁) / Pₐ(y₀)
    // Pₐ(y₀) = pᴸ₁
    Vertex* V0 = light_path.vertices, * V1 = light_path.vertices + 1;
    // Pₐ(y₀) = Pσ(y₀) * G(y₀ –> y₁)
    float sigma_pdf = Light_PDF(V0->origin, &light);
    float forward_g = Geometric_Term(V0, V1);
    float forward_pdf = sigma_pdf * forward_g;
    // consider an area light: Lₑ(P, ωₒ) ≡ Lₑ(P); Lₑ(V₀ –> V₁) = Lₑ⁰(V₀)
    float3 Le = light.emission;
    /* writefloat(forward_pdf); */
    light_contrib[1] = Le * forward_pdf;
    // ---------- MIS weight -------------------
    // <-PA / ->PA
    float backward_g = Geometric_Term(V1, V0);
    float backward_pdf = sigma_pdf * backward_g;
    /* writefloat(backward_g); */
    /* writefloat(backward_pdf); */
    light_weight[1] = backward_pdf * (1.0f/forward_pdf);
    /* writefloat(light_weight[1]); */
  }
  for ( int i = 2; i < light_path.length; ++ i ) { // α₂ –> αₜ
    // -------- unweighted contribution --------
    Vertex* V0 = light_path.vertices + i - 2, * V1 = light_path.vertices + i - 1,
          * V2 = light_path.vertices + i;
    //       fₛ(V2 –> V1 –> V0)
    // αᴸᵢ = —————————————————————– * αᴸᵢ₋₁
    //         Pσ(V1 –> V0)
    float3 wi = normalize(V1->origin - V0->origin),
           wo = normalize(V2->origin - V1->origin);
    Spectrum bsdf_f = BSDF_F(wo, wi, V1->material);
    float backward_g = Geometric_Term(V1, V0);
    float backward_pdf = BSDF_PDF(wo, wi, V1->material) * backward_g;
    light_contrib[i] = fabs(backward_pdf) < 0.01f ? 0.0f :
                       bsdf_f/backward_pdf * light_contrib[i-1];
    // ---------- MIS weight -------------------
    float forward_g  = Geometric_Term(V1, V2);
    /* backward_pdf *= backward_g; //XXX not part of unweighted contrib? */
    float forward_pdf = BSDF_PDF(wi, wo, V1->material) * forward_g;
    light_weight[i] = backward_pdf/forward_pdf * light_weight[i-1];
  }

  Spectrum sample_colour = (Spectrum)(-1.0f);
  // Pretty much same exact code as generate light subpath, except with bdpt,
  //  and I leave out the near-identical equations from L path equations
  Ray ray = (Ray){pixel, dir};
  Vertex V0, V1, V2;
  int mindex = -1;
  for ( int depth = 0; depth != 2; ++ depth ) {
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
    {// calculate contribution and weights
      if ( depth == 0 ) {
        // What's defined: V2 [we're not evaluating anything]
        // αᴱ₀ = 1
        eye_contrib = eye_weight = 1.0f; // unnecessary, but unfortunately the
                                         // if statement is
      } else if ( depth == 1 ) {
        // What's defined: V1 -> V2 [we're evaluating V1] This is for a lens.
        // So, it's ultimately completely useless until an
        // actual lens exists in the SDF scene. Note that this isn't an actual
        // point on a surface of the scene (that's V2, next up for eval), this
        // is just an infinitesimal point based off where the camera currently
        // is, just like a delta light.
        // -------- unweighted contribution --------
        // αᴱ₁ = Wₑ⁰(z₀)/Pₐ(z₀)
        Spectrum We = 1.0f; // XXX This is for lens, we don't have one.
        // Pₐ = Pσ(z̆ᵢ₋₁ –> z̆ᵢ) * G(zᵢ₋₁ <–> zᵢ)
        float sigma_pdf = 1.0f; // Still mathematically correct until a camera
                                // is implemented [if ever]
        float forward_g = 1.0f; // based off camera too
        eye_contrib = We/(sigma_pdf * forward_g);
        // ---------- MIS weight -------------------
        eye_weight = 1.0f;
      } else {
        // What's defined: V0 -> V1 -> V2 [we're evaluating V1]
        // -------- unweighted contribution --------
        //       fₛ(V0 -> V1 -> V2)
        // αᴱᵢ = —————————————————— αᴱᵢ₋₁
        //        Pα(V1 -> V2)
        float3 wi = normalize(V1.origin - V0.origin),
               wo = normalize(V2.origin - V1.origin);
        Spectrum forward_bsdf = BSDF_F(wi, wo, V1.material);
        float forward_g    = Geometric_Term(&V1, &V2);
        float forward_pdf = BSDF_PDF(V1.origin, wo, V1.material) * forward_g;
        eye_contrib = (forward_bsdf/forward_pdf) * eye_contrib;
        // ---------- MIS weight -------------------
        float backward_g   = Geometric_Term(&V1, &V0);
        float backward_pdf = BSDF_PDF(V1.origin, wi, V1.material) * backward_g;
        /* forward_pdf       *= forward_g; */
        eye_weight = (backward_pdf/forward_pdf) * eye_weight;
      }
    }

    // Even though evaluation is on V1, in the case that an emitter is hit,
    // there is enough information to compute this value. Then the loop must be
    // broken
    if ( ptinfo.mat_index <= EMIT_MAT ) { // S=0, T=n strategy
      writeint(depth);
      if ( depth == 0 ) { // directly hit light from camera
        sample_colour = (float3)(1.0f);
        break;
      }
      // This is a very weird situation, we basically have to calculate the PDF
      // of the light source? I believe?
      // Well, we know αᴸ₀ is 1.0f, and we just computed αᴱₜ, and our connection
      // contribution is Lₑ(zₜ₋₁ → zₜ₋₂), but there is no directional component.
      // So in other words, it is the simple C*(0, t) = αᴱₜ * Lₑ(zₜ₋₁)
      // This makes sense intuitively, this is just simple path tracing.
      Emitter path_emitter = REMIT(ptinfo.mat_index);
      Spectrum contribution = eye_contrib * path_emitter.emission;
      // The difficult component is calculating the weight, not even PBRT gives
      // a very good explanation. TODO I believe it's similar to above with just
      // light path weight = 1.0, and using existing eye weight
      if ( sample_colour.x < 0.0f )
        sample_colour = (float3)(0.0f);
      float MIS = 1.0f/(1.0f + light_weight[0] + eye_weight);
      sample_colour += (contribution) * MIS;
      break;
    }

    if ( depth > 0 ) {
      // iterate over every light index & add to contribution sample_colour
      for ( int light_it = 0; light_it != light_path.length; ++ light_it ) {
        Vertex* light_vertex = (light_path.vertices + light_it);
        // calculate unweighted contribution alphaL/E
        Spectrum contribution = eye_contrib * light_contrib[light_it];
        // calculate unweighted contribution connection edge
        float vis = Visibility_Ray(ray.origin, light_vertex->origin, si, Tx);
        if ( light_it == 0 ) { // S == 1 connection strategy TODO
          contribution = (float3)(1.0f);
        } else if ( depth == 0 ) { // T == 1 connection strategy TODO
        } else { // s >= 1, t >= 1 connection contrib
          Vertex* Y2 = (light_path.vertices + light_it - 1), *Y1 = light_vertex,
                * Z2 = &V1, * Z1 = &V2;
          // Y2 -> Y1 <--> Z1 -> Z2
          float3 wi = normalize(Y1->origin - Y2->origin),
                 wo = normalize(Z1->origin - Y1->origin);
          Spectrum light_path_bsdf = BSDF_F(wi, wo, Y1->material);
          /*  */ wi = wo;
                 wo = normalize(Z2->origin - Z1->origin);
          Spectrum eye_path_bsdf   = BSDF_F(wi, wo, Z1->material);
          float G = Geometric_Term(Y1, Z1);
          contribution *= light_path_bsdf * eye_path_bsdf * G;
        }
        // calculate MIS
        float MIS = 1.0f/( 1.0f + eye_weight + light_weight[light_it] );
        // add contrib
        if ( sample_colour.x < 0.0f )
          sample_colour = (float3)(0.0f);
        sample_colour += (contribution * vis) * MIS;
      }
    }

    // setup next ray and mindex
    mindex = ptinfo.mat_index;
    ray.origin = hit;
    // W₀ = fᵢ(Wᵢ, N)
    ray.dir = BSDF_Sample(ray.dir, normal, V2.material, si, Tx);
  }
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
    __global float*             rng
  ) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  float spp = (float)(sinfo->spp);
  // --- construct camera --
  Camera camera = *camera_ptr;
  Ray ray = Camera_Ray(&camera);
  // -- get old pixel, check if there are samples to be done
  //    (counter is stored in alpha channel)
  int pt = out.y*camera.dim.x*4 + out.x*4;
  // -- set up scene info ---
  Material materials[MAT_LEN];
  for ( int i = 0; i != MAT_LEN; ++ i ) {
    materials[i] = *(g_materials + i);
  }
  float time = *time_ptr;
  float3 dval = (float3)(debug_val_ptr[0], debug_val_ptr[1],
                          debug_val_ptr[2]);
  SceneInfo scene_info = New_SceneInfo(time, materials, dval, *rng);
  Uniform_Sample(&scene_info); // throw away first random result

  Update_Camera(&camera, time);

  // grav previous results (has to come after all initialization for fast nav
  //  raycast)
  float4 old_colour;
  {
    _EvalPreviousOutputHorcrux _hor =
          Eval_Previous_Output(img, output_img, out, sinfo, pt);
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
  if ( colour.x >= 0.0f ) {
    old_colour = (float4)(mix(colour, old_colour.xyz,
                            (old_colour.w/(old_colour.w+1.0f))),
                        old_colour.w+1.0f);
    float4 nold_colour;
    writefloat4(nold_colour);
    write_imagef(output_img, out, (float4)(old_colour.xyz, 1.0f));
    //
    img[pt+0] = (unsigned char)(old_colour.x*255.0f);
    img[pt+1] = (unsigned char)(old_colour.y*255.0f);
    img[pt+2] = (unsigned char)(old_colour.z*255.0f);
    img[pt+3] = (unsigned char)(old_colour.w);
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
