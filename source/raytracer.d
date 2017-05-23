module raytracer;
import opencl;
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
private OpenCLGLImage[2] img_buffer; // Ping-pong ogl/ocl interop as described:
      // https://github.com/nvpro-samples/gl_cl_interop_pingpong_st
private OpenCLSingleton!RNG rng_buffer;
private OpenCLSingleton!bool img_reset_buffer;
private OpenCLSingleton!Camera camera_buffer;
private OpenCLSingleton!float timer_buffer;
private OpenCLBuffer!Material material_buffer;
private Material[] material;
private MonoTime timer_start;
private bool reset_camera;
private uint texture;

void Load_Default_File ( ) {
  static import std.file;
  auto tstr = std.file.read("testmap.cl").to!string;
  import functional;
  string mapfunc;
  string[] material_str_list;
  bool hit;
  foreach ( i; tstr.splitLines() ) {
    if ( i == "-MATERIALS" ) {
      hit = true;
      continue;
    }
    if ( hit ) material_str_list ~= i;
    else       mapfunc           ~= i;
  }
  foreach ( m; material_str_list ) {
    float[] fvals = m.split[1 .. $].map!(to!float).array;
    material ~= Create_Material(fvals);
  }
}

void Initialize ( ){
  Load_Default_File();
  opencl.Initialize();
  // --- opengl ---
  import gl_renderer;
  GL_Renderer_Initialize();
  // --- kernel init ---
  Recompile(true); // Just compiles but same as recompiling anyways
}

enum Kernel_Type { Raycast, Raytrace, MLT, }
enum Kernel_Flag { Show_Normals,           }
enum Kernel_Var  { March_Dist, March_Reps }
struct Kernel_Information {
  Kernel_Type type;
  bool[Kernel_Flag] flags;
  int [Kernel_Var ] vars;
}
private Kernel_Information kernel_info;
static this ( ) {
  alias KF = Kernel_Flag, KV = Kernel_Var;
  kernel_info = Kernel_Information(
    Kernel_Type.Raycast,
    [ KF.Show_Normals : false ],
    [ KV.March_Dist   : 64,
      KV.March_Reps   : 128 ]
  );
}
private string map_function;
/** Does not recompile for you */
auto Set_Kernel_Type ( Kernel_Type type ) {
  kernel_info.type = Kernel_Type.Raycast;
}
/** Does not recompile for you */
auto Set_Kernel_Flag ( Kernel_Flag flag, bool status ) {
  kernel_info.flags[flag] = status;
}
/** Does not recompile for you */
auto Set_Kernel_Var   ( Kernel_Var var, int value ) {
  kernel_info.vars[var] = value;
}
/** Does not recompile for you */
void Set_Map_Function ( string val ) { map_function = val; }

private string RPixel_Colour_String ( ) {
  import opencl_kernel;
  final switch ( kernel_info.type ) with ( Kernel_Type ) {
    case Raycast:  return Raycast_pixel_colour_function;
    case Raytrace: return Raytrace_pixel_colour_function;
    case MLT:      return MLT_pixel_colour_function;
  }
}

void Recompile ( bool reset_all = false ) {
  import opencl_kernel, std.string : replace;
  if ( program !is null ) program.Clean_Up();
  string kernel = DTOADQ_kernel;
  {
    alias KV = Kernel_Var;
    kernel = kernel
      .replace("//%MAP", map_function)
      .replace("//%PIXELCOLOUR", RPixel_Colour_String)
      .replace("//%MARCH_DIST", kernel_info.vars[KV.March_Dist].to!string)
      .replace("//%MARCH_REPS", kernel_info.vars[KV.March_Reps].to!string);
  }
  if ( kernel_info.flags[Kernel_Flag.Show_Normals] ) {
    kernel = kernel.replace("//#define SHOW_NORMALS", "#define SHOW_NORMALS");
  }
  program = Compile(kernel);
  program.Set_Kernel("Kernel_Pathtrace");

  auto rng = Generate_New_RNG();
  auto camera = Construct_Camera([1.0f, 0.0f,  0.0f], [ 0.0f, 0.0f, -1.0f],
                                                [Img_dim, Img_dim]);
  if ( !reset_all ) {
    camera   = camera_buffer.data[0];
  } else {
  }

  // --- parameter buffer init ---
  auto RO = BufferType.read_only, WO = BufferType.write_only;
  writeln("DOING GL THING");
  img_buffer[0]    = program.Set_Image_GL_Buffer(RO, Img_dim);
  writeln("DOING GL THING DONESO");
  img_buffer[1]    = program.Set_Image_GL_Buffer(WO, Img_dim);
  img_reset_buffer = program.Set_Singleton!bool(RO, true);
  rng_buffer       = program.Set_Singleton!RNG(RO, rng);
  camera_buffer    = program.Set_Singleton!Camera(RO, camera);
  timer_buffer     = program.Set_Singleton!float(RO, 0.0f);
  material_buffer  = program.Set_Buffer!Material(RO, material);
  writeln("Recompiled!");
}


void Remove ( ) {
  if ( program ) program.Clean_Up();
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

  if ( RMouse_X1() ) {
    cam_update = true;
    lookvec[0] += (RMouse_X-RMouse_X_Stick)/Win_dim;
    lookvec[1] -= (RMouse_Y-RMouse_Y_Stick)/Win_dim;
    camera_buffer.data[0].lookat = lookvec;
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

  camera_buffer.data[0].position = pos;
  return cam_update;
}


static int flip_img = 0;
void Update ( float timer ) {
  program.Acquire_Ownership(img_buffer[0]);
  program.Acquire_Ownership(img_buffer[1]);
  static bool updated_last_frame;
  static int previous_img_dim = 0;
  updated_last_frame = false;
  Render();
  bool material_changed;
  {
    import gui.gui;
    material_changed = Imgui_Render ( material, camera_buffer.data[0] );
  }
  bool reset_img_buffers = false;
  if ( !timer.Should_Update ) return;

  // --- check if dimensions of image changed ---
  if ( previous_img_dim != Img_dim && previous_img_dim != 0 ) {
    Recompile();
    import functional;
    // -- update size of buffers --
    program.Resize_Image_GL_Buffer(img_buffer[0], Img_dim, Img_dim);
    program.Resize_Image_GL_Buffer(img_buffer[1], Img_dim, Img_dim);
    // -- update camera --
    camera_buffer.data[0].dimensions = [Img_dim, Img_dim];
    reset_img_buffers = true;
    // -- write results --
    program.Write(camera_buffer);
  }
  previous_img_dim = Img_dim;

  // --- check if user has changed material ---
  if ( material_changed ) {
    reset_img_buffers = true;
    material_buffer.data = material;
    program.Write(material_buffer);
  }

  // --- check if user has moved camera ---
  if ( Update_Camera() ) {
    reset_img_buffers = true;
    program.Write(camera_buffer);
  }

  // --- check if the image needs to be reset ---
  img_reset_buffer.data[0] = reset_img_buffers;
  program.Write(img_reset_buffer);

  // timer
  timer_buffer.data[0] = timer;
  program.Write(timer_buffer);

  // -- update kernel arguments --
  program.Change_Kernel_Arg(img_buffer[flip_img  ], 0);
  program.Change_Kernel_Arg(img_buffer[flip_img^1], 1);
  flip_img ^= 1;
  program.Change_Kernel_Arg(img_reset_buffer, 2);
  program.Change_Kernel_Arg(rng_buffer,  3);
  program.Change_Kernel_Arg(camera_buffer, 4);
  program.Change_Kernel_Arg(timer_buffer, 5);
  program.Change_Kernel_Arg(material_buffer, 6);
  // -- run --
  program.Run([Img_dim, Img_dim, 1], [1, 1, 1]);
  updated_last_frame = true ;
  program.Release_Ownership(img_buffer[0]);
  program.Release_Ownership(img_buffer[1]);
}

void Render ( ) {
  import gl_renderer;
  GL_Render(img_buffer[flip_img].gl_texture, Img_dim, Img_dim);
}

/** Still need this in case the GPU doesn't have
      GL_CL interoperability for some reason . . .
*/
private void CLImage_To_GLImage ( CLImage image ) {
    import derelict.opengl3.gl3;

  // glBindTexture(GL_TEXTURE_2D, texture);
  // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width, image.height, 0,
  //              GL_RGBA, GL_FLOAT, cast(void*)image.buffer.ptr);
  // glBindTexture(GL_TEXTURE_2D, 0);
}
