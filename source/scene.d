module scene;
import opencl;
import globals;

struct Material {
  cl_float3 base_colour;
  float metallic, subsurface, specular, roughness, specular_tint,
        anisotropic, sheen, sheen_tint, clearcoat, clearcoat_gloss;
  float padding, padding2;
}


import std.random : uniform;
auto Default_Material ( ) {
  return Material(
    To_CLFloat3([uniform(0.0f, 1.0f), uniform(0.0f, 1.0f),
                 uniform(0.0f, 1.0f)]),
    uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f),
    uniform(0.0f, 1.0f)
  );
}



struct Camera {
  cl_float3 position, lookat, up;
  cl_int2 dimensions;
  cl_float fov;
}

auto Construct_Camera(float[3] pos, float[3] dir, int[2] dim) {
  cl_int2 dimensions;
  dimensions.x = dim[0]; dimensions.y = dim[1];
  return Camera(
    To_CLFloat3(pos), To_CLFloat3(dir), To_CLFloat3([0.0f, 1.0f, 0.0]),
    dimensions,
    60.0f
  );
}
