module ocl.structs;
import ocl.opencl;

struct Camera {
  cl_float3 position, lookat, up;
  cl_int2 dimensions;
  cl_float fov;
  cl_int flags;

  this ( float[3] pos, float[3] dir, int[2] dim ) {
    static import stl;
    position = To_CLFloat3(pos);
    stl.writeln(position);
    lookat = To_CLFloat3(dir);
    stl.writeln(lookat);
    up = To_CLFloat3([0.0f, 1.0f, 0.0f]);
    dimensions.x = dim[0]; dimensions.y = dim[1];
    fov = 110.0f;
    flags = 0;
  }
}

struct ImageMetaData {
  ubyte clear_img = true;
  ulong finished_samples;
  ubyte spp = 255;
}

struct RNG {
  cl_ulong[16] seed;
  cl_ulong     p;

  static RNG New ( ) {
    RNG rng;
    import functional, std.random, std.conv : to;
    rng.seed = iota(0, 16).map!(n => uniform(0, cl_ulong.max).to!cl_ulong)
                          .array;
    return rng;
  }
}

struct OCLMaterial {
  // colour [set to (-1.0, -1.0, -1.0) to have map override it]
  cl_float3 albedo;
  // sampling strategy
  @("ModifyFloat") {
    @("Normalize") {
      cl_float diffuse, specular, glossy;
    }
    cl_float glossy_lobe;
    cl_float transmittive;
    // PBR material
    cl_float roughness, metallic, fresnel, subsurface, anisostropic;
  }
  cl_float padding_1, padding_2;
}

auto Set_OCLMaterial ( float[3] colour, float[] vals ) {
  return OCLMaterial(To_CLFloat3(colour),
          vals [0 ], vals [1 ], vals [2 ], vals [3 ], vals [4 ], vals [5 ],
          vals [6 ], vals [7 ], vals [8 ], vals [9 ],
          0.0f, 0.0f);
}

// CLFloat doesn't play nicely with cimgui so have to create this
struct Material {
  // albedo [set to (-1.0, -1.0, -1.0) to have map override it]
  float[3] albedo;
}

import std.traits, std.string : format;
private template TypeMap(string name) {
  mixin(`alias field = %s.%s;`.format(fullyQualifiedName!OCLMaterial, name));
  enum TypeMap = typeof(field).stringof;
}
