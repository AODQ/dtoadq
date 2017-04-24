module raytracer;
import opencl;
import opencl_program : Test_raycast_string;
import globals;
import core.time : MonoTime;
import scene;

struct RNG {
  cl_ulong[16] seed;
  cl_ulong p;
};

RNG Generate_New_RNG() {
  RNG rng;
  import std.random;
  foreach ( ref i; rng.seed ) {
    i = uniform(0, cl_ulong.max).to!(cl_ulong);
  }
  return rng;
}

private OpenCLProgram program;
private OpenCLImage img_buffer_write, img_buffer_read;
private OpenCLSingleton!RNG rng_buffer;
private OpenCLSingleton!Camera camera_buffer;
private OpenCLSingleton!float timer_buffer;
private OpenCLBuffer!Material material_buffer;
private Material[] material;
private MonoTime timer_start;
private bool reset_camera;
private uint texture;

void Initialize ( ){
  // --- opencl        ---
  opencl.Initialize();
  // --- variable init ---
  auto rng = Generate_New_RNG();
  material = [ Default_Material(), Default_Material() ];
  auto camera = Construct_Camera([1.0f, 0.0f,  0.0f], [ 0.0f, 0.0f, -1.0f],
                                                [Img_dim, Img_dim]);
  // --- kernel init ---
  program = Compile(Test_raycast_string);
  program.Set_Kernel("Kernel_Raycast");

  // --- parameter buffer init ---
  auto RO = BufferType.read_only, WO = BufferType.write_only;
  img_buffer_write = program.Set_Image_Buffer(WO, Img_dim);
  img_buffer_read  = program.Set_Image_Buffer(RO, Img_dim);
  rng_buffer       = program.Set_Singleton!RNG(RO, rng);
  camera_buffer    = program.Set_Singleton!Camera(RO, camera);
  timer_buffer     = program.Set_Singleton!float(RO, 0.0f);
  material_buffer  = program.Set_Buffer!Material(RO, material);

  // --- opengl ---
  import derelict.opengl3.gl3;
  glGenTextures(1, &texture);
  import gl_renderer;
  GL_Renderer_Initialize();
}

void Remove ( ) {
  program.Clean_Up();
}

bool Should_Update ( float timer ) {
  static prev_timer = 0.0f;
  // -- only update at ~30 frames --
  bool status = timer - prev_timer > 32.0f/1000.0f;
  if ( status )
    prev_timer = timer;
  return status;
}

bool Update_Camera ( ) {
  import input, std.math;
  bool cam_update;
  auto lookvec = camera_buffer.data[0].lookat;
  auto pos = camera_buffer.data[0].position;
  if ( RMouse_Right() ) {
    cam_update = true;
    lookvec[0] += (RMouse_X-RMouse_X_Stick)/Win_dim/10.0f;
    lookvec[1] -= (RMouse_Y-RMouse_Y_Stick)/Win_dim/10.0f;
    lookvec[1]  = fmax(0.0f, fmin(1.0f, lookvec[1]));
    camera_buffer.data[0].lookat = lookvec;
  }

  if ( RMouse_Middle() ) {
    cam_update = true;
    lookvec[0] += (RMouse_X-RMouse_X_Stick)/Win_dim;
    lookvec[1] -= (RMouse_Y-RMouse_Y_Stick)/Win_dim;
    camera_buffer.data[0].lookat = lookvec;
    Unstick();
  }

  float vel = 0.1f + RKey_Input(32)*0.3f;
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

  camera_buffer.data[0].position = pos;
  return cam_update;
}

void Update ( float timer ) {
  static bool updated_last_frame;
  static int previous_img_dim = 0;
  if ( updated_last_frame )
    program.Read_Image(img_buffer_write);
  updated_last_frame = false;
  Render();
  import gui;
  bool material_changed = Imgui_Render ( material, camera_buffer.data[0] );
  if ( !timer.Should_Update ) return;

  if ( previous_img_dim != Img_dim && previous_img_dim != 0 ) {
    program.Modify_Image_Size(img_buffer_write, Img_dim, Img_dim);
    // program.Write(img_buffer_write);
    camera_buffer.data[0].dimensions = [Img_dim, Img_dim];
    program.Write(camera_buffer);
  }
  previous_img_dim = Img_dim;
  timer_buffer.data[0] = timer;
  if ( material_changed ) {
    material_buffer.data = material;
    program.Write(material_buffer);
  }
  if ( Update_Camera() ) {
    program.Write(camera_buffer);
  }
  program.Write(timer_buffer);
  program.Run([Img_dim, Img_dim, 1], [1, 1, 1]);
  updated_last_frame = true;
}

void Render ( ) {
  CLImage_To_GLImage(CLImage(img_buffer_write.data,
                             cast(int)img_buffer_write.width,
                             cast(int)img_buffer_write.height));
  import gl_renderer;
  GL_Render(texture, Img_dim, Img_dim);
}

private void CLImage_To_GLImage ( CLImage image ) {
    import derelict.opengl3.gl3;
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width, image.height, 0,
               GL_RGBA, GL_FLOAT, cast(void*)image.buffer.ptr);
  glBindTexture(GL_TEXTURE_2D, 0);
}
