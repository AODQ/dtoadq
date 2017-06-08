module camera;
import globals, opencl;

struct Camera {
  cl_float3 position, lookat, up;
  cl_int2 dimensions;
  cl_float fov;
  cl_int flags;
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
  auto lookvec = camera.lookat;
  auto pos = camera.position;

  if ( RMouse_X1() ) {
    cam_update = true;
    lookvec[0] += (RMouse_X-RMouse_X_Stick)/camera.dimensions.x;
    lookvec[1] -= (RMouse_Y-RMouse_Y_Stick)/camera.dimensions.y;
    Unstick();
  }

  float vel = 0.1f + RKey_Input(82)*0.3f - RKey_Input(70)*0.09f;

  if ( RKey_Input(65) || RKey_Input(68) ) { // AD
    cam_update = true;
    float angle = PI - lookvec[0]*2.0f*PI;
    float vel_k = ((RKey_Input(65)?1:-1)*vel);
    pos[0] +=  angle.cos*vel_k;
    pos[2] += -angle.sin*vel_k;
  }
  if ( RKey_Input(87) || RKey_Input(83) ) { // WS
    float angle = PI - lookvec[0]*2.0f*PI;
    cam_update = true;
    float vel_k = ((RKey_Input(87)?1:-1)*vel);
    pos[0] +=  angle.sin*vel_k;
    pos[2] +=  angle.cos*vel_k;
  }
  if ( RKey_Input(81) || RKey_Input(69) ) { // QE
    cam_update = true;
    pos[1] += (RKey_Input(81)?1:-1)*vel;
  }

  camera.lookat = lookvec;
  camera.position = pos;
  return cam_update;
}
