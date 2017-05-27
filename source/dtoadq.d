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
  camera_buffer = Construct_Camera([2.5f, 0.0f, 6.0f], [-0.06f, 0.4f, -1.0f],
                                   [image_buffer.x, image_buffer.y]);
  material_buffer = [ Default_Material() ];
  Set_Image_Buffer(DIMG.Resolution.r640_360, true);
}

void Compile ( ) {
  OCL.Compile(Parse_Kernel());
}

void Set_Image_Buffer ( DIMG.Resolution resolution, bool force = false ) {
  if ( force || image_buffer.resolution != resolution ) {
    img_reset_buffer = true;
    image_buffer = DIMG.RImage(resolution);
    camera_buffer.dimensions.x = image_buffer.x;
    camera_buffer.dimensions.y = image_buffer.y;
    import glfw;
  }
}

cl_event image_unlock_event;

void Update ( float timer ) {
  // -- camera/gui/etc --
  img_reset_buffer |= Update_Camera(camera_buffer);
  img_reset_buffer |= Should_Recompile();
  // -- kernel --
  if ( Should_Recompile(false) ) {
    Compile();
  }
  image_buffer.Lock();
  OCL.Run(
    image_buffer.RRead(),
    image_buffer.RWrite(),
    img_reset_buffer,
    OCL.CLStoreMem(rng_buffer),
    camera_buffer,
    timer,
    material_buffer, material_buffer.length,
    image_buffer.x, image_buffer.y
  );
  auto image_unlock_event = image_buffer.Unlock();
  img_reset_buffer = false;
  // -- read rng --
}

void Render ( ) {
  import gui.gui;
  img_reset_buffer |= Imgui_Render(material_buffer, camera_buffer);
  GL.Render(image_buffer.RRender, image_unlock_event);
}

void Clean_Up ( ) {
  OCL.Clean_Up();
}
