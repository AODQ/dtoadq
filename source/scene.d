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

struct Camera {
  cl_float3 position, lookat, up;
  cl_int2 dimensions;
  cl_float fov;
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

auto Construct_Camera(float[3] pos, float[3] dir, int[2] dim) {
  cl_int2 dimensions;
  dimensions.x = dim[0]; dimensions.y = dim[1];
  return Camera(
    To_CLFloat3(pos), To_CLFloat3(dir), To_CLFloat3([0.0f, 1.0f, 0.0]),
    dimensions,
    60.0f
  );
}

bool Update_Camera ( ref Camera camera ) {
  import input, std.math;
  bool cam_update;
  auto lookvec = camera_buffer.lookat;
  auto pos = camera_buffer.position;

  if ( RMouse_X1() ) {
    cam_update = true;
    lookvec[0] += (RMouse_X-RMouse_X_Stick)/Win_dim;
    lookvec[1] -= (RMouse_Y-RMouse_Y_Stick)/Win_dim;
    camera_buffer.lookat = lookvec;
    Unstick();
  }

  float vel = 0.1f + RKey_Input(82)*0.3f - RKey_Input(70)*0.09f;
  import gln;

  if ( RKey_Input(65) || RKey_Input(68) ) { // AD
    cam_update = true;
    float angle = PI - lookvec[0]*2.0*PI;
    float vel_k = ((RKey_Input(65)?1:-1)*vel);
    pos[0] +=  angle.cos*vel_k;
    pos[2] += -angle.sin*vel_k;
  }
  if ( RKey_Input(87) || RKey_Input(83) ) { // WS
    float angle = PI - lookvec[0]*2.0*PI;
    cam_update = true;
    float vel_k = ((RKey_Input(87)?1:-1)*vel);
    pos[0] +=  angle.sin*vel_k;
    pos[2] +=  angle.cos*vel_k;
  }
  if ( RKey_Input(81) || RKey_Input(69) ) { // QE
    cam_update = true;
    pos[1] += (RKey_Input(81)?1:-1)*vel;
  }

  camera_buffer.position = pos;
  return cam_update;
}
