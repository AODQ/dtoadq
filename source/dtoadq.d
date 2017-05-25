module dtoadq;
import globals, scene, kernelinfo;
static import DIMG = dtoadqimage;
static import OCL  = opencl;
static import GL   = gl_renderer;

// -- kernel vars --
private DIMG.Image image_buffer;
private RNG        rng_buffer;
private bool       img_reset_buffer;
private Camera     camera_buffer;
private Material[] material_buffer;


void Initialize ( ) {
  // -- GL -> OCL --
  GL.Renderer_Initialize();
  OCL.Initialize("DTOAD_Kernel");
  // -- initialize variables --
  DIMG.Initialize();
  image_buffer = DIMG.RImage(DIMG.Resolution.r640_360);
  rng_buffer = RNG.Generate_New();
  img_reset_buffer = true;
  camera_buffer = Construct_Camera([1.0f, 0.0f, 0.0f], [0.0f, 0.0f, -1.0f],
                                   [image_buffer.x, image_buffer.y]);
  material_buffer = [ Default_Material() ];
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

void Set_Image_Buffer ( DIMG.Resolution resolution ) {
  if ( image_buffer.resolution != resolution ) {
    image_buffer = DIMG.RImage(resolution);
  }
}

void Update ( float timer ) {
  if ( Should_Recompile(false) ) {
    Compile();
  }
  image_buffer.Lock();
  OCL.Run(
    image_buffer.RRead(),
    image_buffer.RWrite(),
    img_reset_buffer,
    rng_buffer,
    camera_buffer,
    material_buffer, material_buffer.length,
    image_buffer.x, image_buffer.y
  );
  image_buffer.Unlock();
}

void Render ( ) {
  GL.Render(image_buffer.RRender);
}

void Clean_Up ( ) {
  OCL.Clean_Up();
}
