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
  int mat_index;
} Vertex;

typedef struct T_Edge {
  float3 wi, wo;
} Edge;

typedef struct T_Subpath {
  Vertex vertices[MAX_DEPTH];
  Edge   edges   [MAX_DEPTH];
  uint length;
} Subpath;

typedef struct T_SceneInfo {
  float time;
  // __read_only image2d_array_t textures; IMAGES CANT BE USED AS FIELD TYPES:-(
  __global Material* materials;
  float3 debug_values;
  float rng;
} SceneInfo;
SceneInfo New_SceneInfo(float time, __global Material* materials,
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
  /* In order to do a normal, using SDFs, you must perform a gradient
       approximation: the central difference on the SDF at P, for example the
        difference quotient: (f(x+h) - f(x))/h
     The most obvious way to calculate this is the six-point gradient; where
       for each axis you compute f(x+h) - f(x-h), and normalize the result,
       this returns a vector pointing in the direction where the SDF map
       changes the most, so the gradient is the same as a normal.
     A 4 point gradient is possible, you just have to multiply each point
       of the gradient by epsilon before normalizing. Computing these gradients
       is expensive, because you have to map multiple times, so the speedup
       is well worth the (rather unnoticeable) accuracy.
     It's worth noting that with an additional normal, you can also compute a
       bumpmap, and that with simple object you can calculate the gradient
       by hand rather than this general-purpoe solution. But that's not an
       option with dtoadq.
  */
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
  return (float3)(IPI*m->diffuse, IPI*m->diffuse, IPI*m->diffuse);
}

float BSDF_PDF ( float3 wi, float3 wo, Material* m ) {
  return fabs(wo.z)*IPI;
}

// -----------------------------------------------------------------------------
// --------------- LIGHT   FUNCS -----------------------------------------------
float Visibility_Ray(float3 orig, float3 other, int goal_index,
                     SCENE_T(si, Tx)) {
  float max_dist = Distance(orig, other);
  float3 dir = normalize(other - orig);
  orig -= dir*(MARCH_ACC+0.1f);
  /* writefloat3(ray.origin); */
  SampledPt minfo = March(-1, (Ray){orig, dir}, si, Tx);
  return minfo.mat_index == goal_index;
}

float Light_PDF ( float3 P, Emitter* m ) {
  float sin_theta = (m->radius*m->radius/Distance_Sqr(m->origin, P));
  return 1.0f/(TAU * (1.0f - sqrt(max(0.0f, 1.0f - sin_theta))));
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

SampledPt Colour_Pixel ( Ray ray, SCENE_T(si, Tx) ) {
  float3 coeff = (float3)(0.0f);
  float3 colour = (float3)(0.0f);
  /* float bsdf_pdf = 1.0f; */
  int mindex = -1;
  SampledPt result;
  result.dist = -1.0f;
  float emitter_div = 1.0f/(float)(EMITTER_AMT);
  for ( int depth = 0; depth != 4; ++ depth ) {
    SampledPt minfo = March(mindex, ray, si, Tx);
    if ( minfo.dist < 0.0f ) break;
    // importance sampling
    mindex = minfo.mat_index;
    Material material = si->materials[mindex];
    float3 hit = ray.origin + ray.dir*minfo.dist;
    float3 normal = Normal(hit, si, Tx);
    if ( mindex <= EMIT_MAT ) {
      Emitter emitter = REmission(abs(mindex - EMIT_MAT), si->debug_values,
                                  si->time);
      /* float light_pdf; */
      /* Light_Sample(ray.dir, hit, normal, &emitter, &light_pdf, si); */
      /* float weight = MIS_Weight(bsdf_pdf, emitter_div*light_pdf); */
      float3 emit = emitter.emission * (1.0f/(float)(depth+1.0f));
      colour += coeff * emit;
      result.dist = 1.0f;
      break;
    }


    int light_index = (int)(Uniform_Sample(si)*EMITTER_AMT);
    Emitter light = REmission(light_index, si->debug_values, si->time);

    { // light sample
      float light_pdf, bsdf_pdf;
      float3 wi;
      float3 Li = Light_Sample_Li(hit, &wi, &light_pdf, &light, si, Tx);
      /* if ( light_pdf > 0.0f ) { */
        float3 f = BSDF_F(wi, ray.dir, &material) * fabs(dot(wi, normal))
                   * (minfo.colour/PI);
        bsdf_pdf = BSDF_PDF(wi, ray.dir, &material);
        bool visible = Visibility_Ray(hit, Li, light_index, si, Tx);
        /* visible &= light_pdf > 0.0f; */
        visible = true;
        float weight = Power_Heuristic(light_pdf, 1.0f,  bsdf_pdf, 1.0f);
        colour += f * light.emission * visible * 100.0f *
                  (weight/(emitter_div*light_pdf));

        colour = fabs(colour);
        result.dist += (float)(visible)*2.0f;
      /* } */
    }

    { // bsdf sample
      float3 wo       = BSDF_Sample(ray.dir, normal, &material, si, Tx);
      float3 f        = BSDF_F(ray.dir, wo, &material)*fabs(dot(ray.dir, wo))
                        * (minfo.colour/PI);
      float bsdf_pdf  = BSDF_PDF(ray.dir, wo, &material);
      float light_pdf = Light_PDF(hit, &light);
      float weight    = Power_Heuristic(1.0f/bsdf_pdf,          2.0f,
                                       (light_pdf*emitter_div), 2.0f);

      coeff += f * (weight/bsdf_pdf);

      ray.dir = wo;
      ray.origin = hit;
    }
  }
  result.colour = min((float3)(1.0f), max((float3)(0.0f), fabs(colour)));
  return result;
}


#if 0
/*
  typedef struct T_Vertex {
  float3 colour, origin, normal;
  int mat_index;
  } Vertex;

  typedef struct T_Edge {
  float3 wi, wo;
  } Edge;

  typedef struct T_Subpath {
  Vertex vertices[MAX_DEPTH];
  Edge   edges   [MAX_DEPTH];
  uint length;
  } Subpath;
 */
Emitter Generate_Light_Subpath ( Subpath* path, Scene_T(si, Tx)) {
  // Grab a random light source to generate a subpath
  int light_index = Uniform_Sample(si)*EMITTER_AMT;
  Emitter light = REmission(light_index, si->debug_values, si->time);
  int mindex = EMIT_MIT - light_index;
  Ray ray;
  {// Generate ray position/angle
    float3 origin = (normalize(Uniform_Sample(si)) - (float3)(0.5f)) *
                  (float2)(2.0f);
      origin = light.origin + origin*light.radius;
    ray = (Ray){origin, Uniform_Sample3(si)};
  }
  // generate path, Σ₀iᴿᴿ(Vᵢ + ωₒ*Rd)
  // where Rd is raymarch-distance, RR is russian-roulette, and
  // ωₒ = i=0 ? light-ray-dir : fᵢ(ωᵢ, N)
  // fᵢ = BSDF importance-sampling (Alternatively, you could sample a random
  //         point on hemisphere Ω, but fᵢ ensures f(P, ωᵢ, ωₒ) > 0.0)
  // Also collect vertex and edge data along the way
  for ( int depth = 0; depth != 3; ++ depth ) {
    SampledPt ptinfo = March(mindex, ray, si, Tx);
    if ( minfo.dist < 0.0f )
      break;

    // Vᵢ = Vᵢ₋₁ + ωᵢ*Rd
    float3 hit = ray.origin + ray.dir*minfo.dist;
    float3 normal = Normal(hit, si, Tx);
    {// set up vertex
      path.vertices[depth].colour = ptinfo.colour;
      path.vertices[depth].origin = hit;
      path.vertices[depth].normal = Normal(hit);
      path.vertices[depth].mat_index = ptinfo.mat_index;
    }
    {// set up edge
      path.edges[depth].wi = ray.dir;
      if ( depth > 0 ) {
        path.edges[depth-1].wo = ray.dir;
      }
    }

    // W₀ = fᵢ(Wᵢ, N)
    float3 wo = BSDF_Sample(ray.dir, normal, &material, si, Tx);
    mindex = minfo.mat_index;
    ray.origin = hit;
    ray.dir = wo;
  }
  return light;
}
Spectrum Bidirectional_Pathtrace ( float3 pixel, float3 dir, SCENE_T(si, Tx)) {
  /**
    Remember the rendering equation, that uses the solid angle
      Lₒ(P, ωₒ) = Lₑ(P, ωₒ) + ∫ Fₛ(P, ωᵢ, ωₒ) Lᵢ(P, ωᵢ) cosθᵢ dωᵢ
                              Ω

    To get from solid angle to path, consider:
      Fₛ(v –> v` –> v``) = Fₛ(P, ωᵢ, ωₒ)
      Lₑ(v –> v`)        = Lₑ(P, ωₒ)
      cosθ  , is transformed to area with G(v <–> v`)
    Lₒ(v –> v`) = Lₑ(v –> v`) * ∫M Lₒ(v₀–>v₁)G(v₀<–>v₁)fₛ(v₍ₖ₋₁₎–>vₖ) dA
      where A is area
      M is all set of possible paths
    This is known as the three-point form or light transport equation, we
      recursively expand it to get this:

        ∞                              k-1
        Σ  ∫( Lₑ(V₀ –> V₁) G(V₀ <–> V₁) Π( fₛ(Vᵢ₋₁–>Vᵢ–>Vᵢ₊₁)G(Vᵢ<–>Vᵢ₊₁)
       k=1 ₼ᵏ⁺¹                        i=1
                     * Wₑʲ(Vₖ₋₁ –>Vₖ) dA(v₀) ̇̇̇ dA(vₖ)
           ))

    This becomes, = Iⱼ = ∫Ω ( pⱼ(v̆) dμ(v̆) )
      where v̆ is a path z₀–>zₛ–>yₜ–>y₀ (bsdf -> light) [or v₀–>vₖ, k=s+t-1]
      Ω is the set of paths of any length
      μ is the union of a set of paths
    since μ is the area-product measure, we have to use the area-integral
      of the transport equation. Well, this is actually the recursively
      expanded form of the transport equation:

    Of course, to compute this, we need to apply monte carlo sampling:
      Iⱼ ≈ fⱼ(v̆)/p(v̆)
       where p is the probability distribution function. It should be noted
       that this PDF must be with respects to the area, so we must convert the
       PDF w.r.t. solid angle.

    The amount that each path contributes to Iⱼ differ greatly, so even naively
      weighing each contribution by something like 1.0/(s+t+1), will produce very
      noisy images, with results similar to naive path tracing.
    In order to get good contributions, you must combine them in an optimal
      manner using multiple importance sampling, as described in Veach's thesis:
        F = ΣₛΣₜWₛₜ(v̆ₛₜ)/(pₛₜ(v̆ₛₜ))
        But this can be rewritten to
            ΣₛΣₜCₛ,ₜ = (Wₛₜ(v̆₍ₛₜ₎))*(αᴸₛ*cₛₜ*αᴱₜ)
        where (αᴸₛ*cₛₜ*αᴱₜ) is the unweighted contribution (Cᵘₛₜ)
              α is the contribution from light L or bsdf B
              c is the contribution from connecting L and B, zₛ–>yₜ
              W is the weight
      For clarity,
        Cᵘₛₜ = (αᴸₛ * cₛₜ * αᴱₜ) ≡ fⱼ(V̆ₛₜ)/p(x̆ₛₜ)
      And if you were to expand this, you would get the light transport equation,
        divided by a probablity distribution function.
      Thus F = ΣₛΣₜWₛₜ*Cᵘₛₜ
  **/
  Subpath light_path;
  // Instead of generating a light and bsdf path, I generate the light path and
  // then walk the bsdf path, in order to conserve GPU memory
  Emitter light = Generate_Light_Subpath(&light_path, si, Tx);
  // -- generate weights and PDFS --
  Spectrum bsdf_weight, light_weights[MAX_DEPTH];
  Spectrum bsdf_contrib,  light_contrib[MAX_DEPTH];
  // Calculate weights/pdfs for light path... I calculated pdfs as described in
  //   veach's paper.
  // α₀⁽ᴸᴱ⁾ = 1
  light_contrib[0] = light.colour;
  bsdf_contrib = 1.0f;
  if ( light_path.length > 1 ) { // α₁
    // αᴸ₁ = Lₑ⁰(y₀ –> y₁) / Pₐ(y₁)
    // Calculate Pₐ then Lₑ⁰
    // Pₐ(y₀) = pᴸ₁
    // where Pₐ(i) [PDF w.r.t. area] = Pσ(v̆₍ᵢ₋₁₎ –> v̆ᵢ) * G(v̆₍ᵢ₋₁₎ <–> v̆ᵢ)
    //       Pσ = PDF w.r.t. solid angle
    //       G(v <–> v') = V(v –> v')*(cosθₒ * cosθ'ᵢ)/||v - v'||²
    //       V = visibility test, no need to calculate it here; check the
    //             expansion of the light transport eq., they all cancel out
    //             except for cₛₜ – makes sense intuitively as the path has
    //             already been constructed: v –> v' must be visible
    Vertex* V1 = light_path.vertices, * V0 = light_path.vertices - 1;
    Edge  * E1 = light_path.edges,    * E0 = light_path.edges    - 1;
    // Pσ * G(v <–> v')
    float sigma_pdf = Light_PDF(V1.origin, &light),
          g = (dot(E1.wi, V0.normal)*dot(E1.wi, V1.normal) / // cosθₒ * cosθ'ᵢ
               sqr(length(V0.origin - V1.origin)));          // sqr magnitude
    float light_pdf = sigma_pdf * g;
    // Lₑ⁰(y₀ –> y₁) [spatial and directional components of Le]
    // Specifically, the Lₑ component of the rendering eq. is split into two:
    //   Lₑ(V₀–>V₁) = Lₑ⁰(V₀)* Lₑ¹(V₀–>V₁)
    //   XXX This is mostly meaningless, because I only deal with area lights.
    //         consider Lₑ(P, ωₒ) ≡ Lₑ(P):
    //            Lₑ(V₀ –> V₁) = Lₑ⁰(V₀)
    float Le = light.colour;
    // αᴸ₁ = Lₑ⁰(y₀ –> y₁) / Pₐ(y₁)
    light_contrib[1] = Le / light_pdf;
  }
  for ( int i = 2; i < light_path.length(); ++ i ) { // α₂ –> αₜ
    Vertex* V0 = light_path.vertices - 1, * V1 = light_path.vertices,
            V2 = light_path.vertices + 1;
    Edge  * E0 = light_path.edges    - 1, * E1 = light_path.edges,
            E2 = light_path.vertices + 1;
    //       fₛ(yᵢ₊₁ –> yᵢ –> yᵢ₋₁)
    // αᴸᵢ = —————————————————————– * αᴸᵢ₋₁
    //         Pσ(yᵢ –> yᵢ₋₁)
    Spectrum bsdf_f = BSDF_F(E1.wo, E1.wi, &material);
    float sigma_pdf = BSDF_PDF(E1.wo, E1.wi, &material);
    light_contrib[i] = bsdf_f/sigma_pdf * (light_contrib[i-1]);
  }

  Spectrum sample_colour = (Spectrum)(0.0f);
  // Pretty much same exact code as generate light subpath, except with bdpt,
  //  and I leave out the near-identical equations from L path equations
  Ray ray = Ray{pixel, dir};
  Vertex V0, V1, V2;
  Edge E0, E1, E2;
  Spectrum prev_contrib;
  for ( int depth = 0; depth != 3; ++ depth ) {
    SampledPt ptinfo = March(mindex, ray, si, Tx);
    if ( minfo.dist < 0.0f )
      break;

    // V0 –> V1 –> V2
    // V2 is current, but we're evaluating V1
    V0 = V1; V1 = V2;
    {// set up vertex
      V2.colour = ptinfo.colour;
      V2.origin = hit;
      V2.normal = Normal(hit);
      V2.mat_index = ptinfo.mat_index;
    }
    {// set up edge
      E2.wi = ray.dir;
    }
    // Vᵢ = Vᵢ₋₁ + ωᵢ*Rd
    float3 hit = ray.origin + ray.dir*minfo.dist;
    float3 normal = Normal(hit, si, Tx);
    // W₀ = fᵢ(Wᵢ, N)
    float3 wo = BSDF_Sample(ray.dir, normal, &material, si, Tx);
    // previous vertices/edges
    {// calculate contribution and weights
      if ( depth == 0 ) {
        // αᴱ₀ = 1
        bsdf_contrib = 1.0f;
      } else if ( depth == 1 ) {
        // ———————————— contribution ————————————
        // αᴱ₁ = Wₑ⁰(z₀)/Pₐ(z₀)
        float We = V1.colour; // XXX correct ?
        // Pₐ = Pσ(z̆ᵢ₋₁ –> z̆ᵢ) * G(zᵢ₋₁ <–> zᵢ)
        float sigma_pdf = 1.0f/PI, // TODO, based off camera ??
        g = (dot(E1.wi, V1.normal)*dot(E0.wo, V1.normal) / // cosθₒ * cosθ'ᵢ
            sqr(length(V0.origin - V1.origin)));          // sqr magnitude
        bsdf_contrib = We/(sigma_pdf * g);
        prev_contrib = bsdf_contrib;
        // ———————————— weight       ————————————
      } else {
        //       fₛ(zᵢ₋₁ –> zᵢ₋₂ –> zᵢ₋₃)
        // αᴱᵢ = ———————————————————————— αᴱᵢ₋₁
        //        Pα(zᵢ₋₂ –> zᵢ₋₁)
        float bsdf = BSDF_F(E1.wi, E2.wi, &si.materials);
        float sigma_pdf = BSDF_PDF(E1.wo, E1.wi, &material);
        bsdf_contrib = bsdf/sigma_pdf * prev_contrib;
        prev_contrib = bsdf_contrib;
      }
      // Calculate pᵢ(y̆)
      // calculate weight
      /*
                    pₛ(x̆)
        ωₛₜ(x̆) = —————————–
                 Σᵢ Pᵢ(x̆)
        As described by Veach's thesis; there's a few problems with this:
          1) PDFs can under/overflow floats; a vertex's area density is
             inversely proportional to the square of the scene dimensions.
          2) Has a complexity of O(n²)
          3) Produces very lengthy code that is hard to debug

        PBR 3rd Ed, chapter 16 pg 1014-1016 explains in detail how to solve this
          while simplifying the equation:

                 ( 1.0f,                i = s
                 |
                 | p←(xᵢ)
                 | ——————– rᵢ₊₁(x̆)      i < s
         rᵢ(x̆) = { p→(xᵢ)
                 |
                 | p→(xᵢ₋₁)
                 | ————————– rᵢ₋₁(x̆)    i > s
                 ( p←(xᵢ₋₁)

                                              1.0f
        Thus, something like Wₛₜ(x̆) = ——————————–——————————–
                                       rₛ(x̆) + rₜ(x̆) + 1.0f

      Well, I've found there's a bit more simplification possible. The thing is,
        there isn't a single path x̆, there's a light and eye path, the latter of
        which is an input range. So thus you can reduce:

                                1.0f
          Wₛₜ(y̆, z̆) = –––––––––––––––——————————–
                       Rₛ₋₁(y̆) + Rₜ₋₁(z̆) + 1.0f

              ( 1.0f               i=0
              |
      Rᵢ(x̆) = { Pᵢ←(x̆)
              | —————– * Rᵢ₋₁(x̆)  i>0
              ( pᵢ→(x̆)

      There's no need to swap the Pᵢ directions, because they are in respects to
        the light path and not the entire path (swapping them to form to the
        camera path). Expanding the recursion we get

               i    Pᵢ←(x̆)
      Rᵢ(x̆) =  Π    —————–
              n=0   Pᵢ→(x̆)

      pᵢ = the probability of a path to have been created. Specifically, pₛ is
        how the path was actually generated, while p₀ ̣̣̣ pₛ₋₁ , p₀ ̣̣̣ pₜ₋₁
        represents the ways the path could have been generated
      pᵢ = pᵢ,₍ₛ₊ₜ₎₋ᵢ(x̆ₛₜ) for i 0 .. s+t
      */
    }
    {
    }
    // ωₛₜ ≡ ωₛₜ(x̆ₛₜ) ; value depends on pdfs with which x̆ is generated by s+t+1
    //                  possibilities
    // pᵢ = pᵢ₍ₛ₊ₜ₎₋ᵢ(x̆ₛₜ) for i = 0 .. s+t
    // ωₛₜ = p²ₛ
    //       ——–
    //      ΣᵢP²ᵢ

    mindex = minfo.mat_index;
    ray.origin = hit;
    ray.dir = wo;
  }
  // Iterate over
  for ( int lindex = light_subpath.length-1; lindex >= 0; -- lindex ) {
    int min_bindex = max(2-lindex, 2),
        max_bindex = min(bsdf_subpath.length - 1, MAX_DEPTH+1-lindex);

    for ( int bindex = max_bindex; bindex >= min_bindex; -- bindex ) {
      bool was_sample_direct = false;
      int remaining = MAX_DEPTH - lindex - bindex + 1;

      Spectrum path_colour = (Spectrum)(0.0f);

      // connect paths
    }
  }
}
#endif
// Bidirectional Path Tracing . .
// Generate camera and light subpaths
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
    __global float*             rng
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
  SceneInfo scene_info = New_SceneInfo(time, materials, dval, *rng);
  Uniform_Sample(&scene_info);
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
