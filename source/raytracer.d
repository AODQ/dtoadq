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
  // --- kernel init ---
  static import std.file;
  writeln("FORMAT STR");
  auto tstr = std.file.read("testmap.cl").to!string;
  import functional;
  string mapfunc;
  string[] material_str_list;
  bool hit;
  foreach ( i; tstr.splitLines() ) {
    writeln("`", i, "` == `-MATERIALS` ?");
    if ( i == "-MATERIALS" ) {
      hit = true;
      continue;
    }
    if ( hit ) material_str_list ~= i;
    else       mapfunc           ~= i;
  }
  writeln("MAPFUNC: ", mapfunc);
  writeln("\n\nmaterial list: ", material_str_list);
  foreach ( m; material_str_list ) {
    float[] fvals = m.split[1 .. $].map!(to!float).array;
    material ~= Create_Material(fvals);
  }
  Recompile(true); // Just compiles but same as recompiling anyways
  // --- opengl ---
  import derelict.opengl3.gl3;
  glGenTextures(1, &texture);
  import gl_renderer;
  GL_Renderer_Initialize();
  writeln("done with renderer");
}

enum Kernel_Type { Raycast, Raytrace, MLT, }
enum Kernel_Flag { Render_Normals,         }
struct Kernel_Information {
  Kernel_Type type;
  bool[Kernel_Flag] flags;
}
private Kernel_Information kernel_information;
static this ( ) {
  with ( Kernel_Flag ) {
    kernel_information = Kernel_Information(
      Kernel_Type.Raycast,
      [ Render_Normals : false ]
  );}
}
private string map_function;
/** Does not recompile for you */
auto Set_Kernel_Type ( Kernel_Type type ) {
  kernel_information.type = Kernel_Type.Raycast;
}
/** Does not recompile for you */
auto Set_Kernel_Flag ( Kernel_Flag flag, bool status ) {
  kernel_information.flags[flag] = status;
}
/** Does not recompile for you */
void Set_Map_Function ( string val ) { map_function = val; }

private string RPixel_Colour_String ( ) {
  import opencl_kernel;
  final switch ( kernel_information.type ) with ( Kernel_Type ) {
    case Raycast:  return Raycast_pixel_colour_function;
    case Raytrace: return Raytrace_pixel_colour_function;
    case MLT:      return MLT_pixel_colour_function;
  }
}

void Recompile ( bool reset_all = false ) {
  import opencl_kernel, std.string : replace;
  if ( program !is null ) program.Clean_Up();
  string kernel = DTOADQ_kernel;
  kernel = kernel.replace("//%MAP%//", map_function)
                 .replace("//%PIXELCOLOUR%//", RPixel_Colour_String);
  program = Compile(kernel);
  program.Set_Kernel("Kernel_Pathtrace");

  auto rng = Generate_New_RNG();
  auto camera = Construct_Camera([1.0f, 0.0f,  0.0f], [ 0.0f, 0.0f, -1.0f],
                                                [Img_dim, Img_dim]);
  if ( !reset_all ) {
    camera = camera_buffer.data[0];
  }

  // --- parameter buffer init ---
  auto RO = BufferType.read_only, WO = BufferType.write_only;
  img_buffer_write = program.Set_Image_Buffer(WO, Img_dim);
  img_buffer_read  = program.Set_Image_Buffer(RO, Img_dim);
  rng_buffer       = program.Set_Singleton!RNG(RO, rng);
  writeln("CAMERA BUFFER: ", camera_buffer);
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
  bool material_changed;
  {
    import gui.gui;
    material_changed = Imgui_Render ( material, camera_buffer.data[0] );
  }
  bool reset_img_buffers = false;
  if ( !timer.Should_Update ) return;

  // --- check if dimensions of image chanegd ---
  if ( previous_img_dim != Img_dim && previous_img_dim != 0 ) {
    import functional;
    // -- update size of buffers --
    program.Modify_Image_Size(img_buffer_write, Img_dim, Img_dim);
    program.Modify_Image_Size(img_buffer_read, Img_dim, Img_dim);
    // -- update camera --
    camera_buffer.data[0].dimensions = [Img_dim, Img_dim];
    // -- update image buffers --
    reset_img_buffers = true;
    float[] buffer;
    buffer.length = 4*Img_dim*Img_dim;
    buffer = buffer.map!(n => 0.0f).array;
    img_buffer_read.data = img_buffer_write.data = buffer;
    // -- write results --
    program.Write(camera_buffer);
    program.Write(img_buffer_write);
    // (no need to write read as it's done every frame anyways)
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
  if ( reset_img_buffers ) {
    import functional;
    // {import functional; img_buffer_write.data.each!((ref n) => n = 0.0f);}
    for ( int i = 3; i < img_buffer_write.data.length; i += 4 )
      img_buffer_write.data[i] = 0.0f;
    program.Write(img_buffer_write);
    // again, no need to write to read as it's done directly after
  }

  // -- update image read and timer buffers --
  img_buffer_read.data = img_buffer_write.data;
  program.Write(img_buffer_read);
  timer_buffer.data[0] = timer;
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
