module dtoadq;
import globals, scene, kernelinfo;
static import DIMG = dtoadqimage;
static import OCL  = opencl;
static import GL   = gl_renderer;
import camera;
import derelict.opencl.cl : cl_event;

// -- kernel vars --
private DIMG.Image image_buffer;
private RNG        rng_buffer;
private bool       img_reset_buffer;
private Camera     camera_buffer;
private Material[] material_buffer;


void Initialize ( ) {
  // -- GL -> OCL --
  GL.Renderer_Initialize();
  OCL.Initialize("DTOADQ_Kernel");
  // -- initialize variables --
  DIMG.Initialize();
  rng_buffer = RNG.Generate_New();
  img_reset_buffer = true;
  camera_buffer = Construct_Camera([1.0f, 0.0f, 0.0f], [0.0f, 0.0f, -1.0f],
                                   [image_buffer.x, image_buffer.y]);
  material_buffer = [ Default_Material() ];
  Set_Image_Buffer(DIMG.Resolution.r640_360, true);
}

void Compile ( ) {
  import opencl_kernel, std.string : replace;
  string kernel = DTOADQ_kernel;
  {
    alias KV = Kernel_Info.Var;
    kernel = kernel
      .replace("//%MAP", RKernel_Info.map_function)
      .replace("//%KERNELTYPE", RKernel_Type_String)
      .replace("//%MARCH_DIST", RKernel_Info.vars[KV.March_Dist].to!string)
      .replace("//%MARCH_REPS", RKernel_Info.vars[KV.March_Reps].to!string);
  }
  if ( RKernel_Info.flags[Kernel_Info.Flag.Show_Normals] ) {
    kernel = kernel.replace("//#define SHOW_NORMALS", "#define SHOW_NORMALS");
  }
  OCL.Compile(kernel);
}

void Set_Image_Buffer ( DIMG.Resolution resolution, bool force = false ) {
  if ( force || image_buffer.resolution != resolution ) {
    image_buffer = DIMG.RImage(resolution);
    camera_buffer.dimensions.x = image_buffer.x;
    camera_buffer.dimensions.y = image_buffer.y;
    writeln("IMG X: ", image_buffer.x, " Y: ", image_buffer.y);
    writeln("IMAGE BUFFER: ", image_buffer);
    import glfw;
  }
}

cl_event image_unlock_event;

void Update ( float timer ) {
  // -- camera/gui/etc --
  Update_Camera(camera_buffer);
  // -- kernel --
  if ( Should_Recompile(false) ) {
    writeln("recompilin");
    Compile();
  }
  image_buffer.Lock();
  OCL.Run(
    image_buffer.RRead(),
    image_buffer.RWrite(),
    img_reset_buffer,
    rng_buffer,
    camera_buffer,
    timer,
    material_buffer, material_buffer.length,
    image_buffer.x, image_buffer.y
  );
  auto image_unlock_event = image_buffer.Unlock();
  img_reset_buffer ^= 1;
}

void Render ( ) {
  import gui.gui;
  Imgui_Render(material_buffer, camera_buffer);
  GL.Render(image_buffer.RRender, image_unlock_event);
}

void Clean_Up ( ) {
  OCL.Clean_Up();
}
