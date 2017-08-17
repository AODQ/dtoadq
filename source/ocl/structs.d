module ocl.structs;
import ocl.opencl;

struct Camera {
  cl_float3 position, lookat, up;
  cl_int2 dimensions;
  cl_float fov;
  cl_int flags;

  this ( float[3] pos, float[3] dir, int[2] dim ) {
    position = To_CLFloat3(pos);
    lookat = To_CLFloat3(dir);
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
    import functional;
    rng.seed = iota(0, 16).map!(n => n = n.uniform(0, cl_ulong.max).to!cl_ulong)
                          .array;
    return rng;
  }
}

struct Material {
  float diffuse, specular, glossy, retroreflective, transmittive;
}

auto Default_Material ( ) {
  return Material(1.0f);
}

auto Create_Material ( float[] vals ) {
  return Material(
    vals[ 0], vals[ 1], vals[ 2], vals[ 3]
  );
}
