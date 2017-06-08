module dtoadq;
import globals, scene;
static import KI   = kernelinfo;
static import DIMG = dtoadqimage;
static import OCL  = opencl;
static import GL   = gl_renderer;
import camera;
import derelict.opencl.cl : cl_event;
import parser : Reparse_Files;

// -- kernel vars --
private DIMG.Image image_buffer;
private RNG        rng_buffer;
private bool       img_reset_buffer;
private Camera     camera_buffer;
private Material[] material_buffer;


void Initialize ( bool init_gl = true ) {
  // -- GL -> OCL --
  if ( init_gl ) GL.Renderer_Initialize();
  OCL.Initialize();
  // -- initialize variables --
  DIMG.Initialize();
  rng_buffer = RNG.Generate_New();
  img_reset_buffer = true;
  camera_buffer = Construct_Camera([2.5f, 0.0f, 6.0f], [-0.06f, 0.4f, -1.0f],
                                   [image_buffer.x, image_buffer.y]);
  material_buffer = [ Default_Material() ];
  Set_Image_Buffer(DIMG.Resolution.r640_360, true);
  KI.Set_Map_Function(KI.ProceduralType.Model,
                      KI.FileType.CL,
                      "projects/globals/models/sdSphere.cl");
}

void Compile ( string kernel_name ) {
  import parser;
  OCL.Compile(Parse_Kernel(), kernel_name);
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

static float timer = 0.0f;

void Set_Time ( float t ) {
  timer = t;
}

bool allow_time_change = true;

void Add_Time ( float t ) {
  if ( allow_time_change )
    timer += t;
}

float RTime ( ) { return timer; }

void Update_DTOADQ ( bool video_update, ref float[] data ) {
  auto imgs = DIMG.RImages();
  // -- buffer --
  bool should_reparse = Reparse_Files();
  // -- camera/gui/etc --
  img_reset_buffer |= Update_Camera(camera_buffer);
  img_reset_buffer |= KI.Should_Recompile();
  img_reset_buffer |= should_reparse;
  // -- kernel --
  if ( KI.Should_Recompile(false) || should_reparse ) {
    writeln("Recompiling");
    Compile("DTOADQ_Kernel");
    return;
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
    OCL.CLPredefinedMem(imgs),
    // kernel size
    image_buffer.x, image_buffer.y
  );
  image_buffer.Unlock();
  img_reset_buffer = false;
  if ( video_update ) {
    return image_buffer.RFloatData(data);
  }
}

DIMG.GLImage gl_image;

void Update_Texture ( ) {
  bool should_reparse = Reparse_Files();
  if ( gl_image is null ) {
    gl_image = new DIMG.GLImage(1024, 1024);
  }
  // we only want to run once per change
  if ( KI.Should_Recompile(false) || should_reparse ) {
    float[] data = DIMG.Kernel_Run_Texture(1024, KI.RFilename);
    gl_image.Update(data);
  }
}

void Update ( bool video_update, ref float[] data ) {
  if ( KI.RFile_Type == KI.FileType.TXT ) {
    Update_Texture();
  } else {
    Update_DTOADQ(video_update, data);
  }
}

void Render ( ) {
  import gui.gui;
  img_reset_buffer |= Imgui_Render(material_buffer, camera_buffer);
  if ( KI.RFile_Type == KI.FileType.TXT ) {
    if ( gl_image !is null ) {
      GL.Render(gl_image.gl_texture, image_unlock_event);
    }
  } else {
    GL.Render(image_buffer.RRender, image_unlock_event);
  }
}

void Clean_Up ( ) {
  OCL.Clean_Up();
}
