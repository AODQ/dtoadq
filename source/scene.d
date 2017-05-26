module scene;
import opencl;
import globals;

struct Material {
  cl_float3 base_colour;
  float metallic, subsurface, specular, roughness, specular_tint,
        anisotropic, sheen, sheen_tint, clearcoat, clearcoat_gloss,
        emission;
  float padding;
}

struct RNG {
  cl_ulong[16] seed;
  cl_ulong     p;

  static RNG Generate_New ( ) {
    RNG rng;
    import std.random;
    foreach ( ref i; rng.seed ) {
      i = uniform(0, cl_ulong.max).to!cl_ulong;
    }
    return rng;
  }
}

import std.random : uniform;
auto Default_Material ( ) {
  return Material(
    To_CLFloat3([uniform(0.0f, 1.0f), uniform(0.0f, 1.0f),
                 uniform(0.0f, 1.0f)]),
    uniform(0.0f, 1.0f), uniform(0.0f, 1.0f), uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f), uniform(0.0f, 1.0f), uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f), uniform(0.0f, 1.0f), uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f),
    0.0f
  );
}

auto Create_Material ( float[] vals ) {
  return Material(
    To_CLFloat3(cast(float[3])vals[0..3]),
    vals[ 3], vals[ 4], vals[ 5], vals[ 6], vals[ 7], vals[ 8], vals[ 9],
    vals[10], vals[11], vals[12], vals[13],
  );
}

