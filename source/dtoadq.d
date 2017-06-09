module dtoadq;
import globals;
import oclstructs;
static import KI   = kernelinfo;
static import DIMG = dtoadqimage;
static import OCL  = opencl;
static import GL   = gl_renderer;
import derelict.opencl.cl : cl_event;
import parser : Reparse_Files;

private class DTOADQ {
  // -- kernel arguments
  DIMG.Image image;
  DIMG.GLImage gl_image; // only for texture
  ubyte[] rw_image;
  RNG rng;
  bool image_reset = true;
  Camera camera;
  Material[] material;
  float timer = 0.0f;
  // -- opengl/opencl events
  cl_event image_unlock_event;
  // -- kernel debug
  float[] debug_vals;
  bool allow_time_change = true;
  bool running = true;

  // -- funcs
  this ( ) {
    debug_vals = [ 0.0f, 0.0f, 0.0f ];
    rng = RNG.Generate_New();
    camera = Construct_Camera([2.5f, 0.0f, 6.0f], [-0.06f, 0.4f, -1.0f],
                                    [image.x, image.y]);
    material = [ Default_Material() ];
    Set_Image_Buffer(DIMG.Resolution.r640_360, true);
  }

  void Set_Image_Buffer ( DIMG.Resolution resolution, bool force = false ) {
    if ( force || image.resolution != resolution ) {
      image_reset = true;
      image = DIMG.RImage(resolution);
      rw_image.length = image.x*image.y*4;
      camera.dimensions.x = image.x;
      camera.dimensions.y = image.y;
    }
  }

  void Compile ( string kernel_name ) {
    import parser;
    OCL.Compile(Parse_Kernel(), kernel_name);
  }

  void Update_DTOADQ ( ) {
    auto imgs = DIMG.RImages();
    // -- buffer --
    bool should_reparse = Reparse_Files();
    // -- camera/gui/etc --
    image_reset |= camera.Update|KI.Should_Recompile()|should_reparse;
    // -- kernel --
    if ( KI.Should_Recompile(false) || should_reparse ) {
      Compile("DTOADQ_Kernel");
      writeln("Recompiled");
      return;
    }
    image.Lock();
    OCL.Run(
      OCL.CLStoreMem(rw_image),
      image.RWrite(),
      camera,
      timer,
      OCL.CLPredefinedMem(imgs),
      debug_vals,
      // kernel size
      image.x, image.y
    );
    image.Unlock();
    image_reset = false;
  }
  void Update_Render ( ) {
    Update_DTOADQ();
    static import VI = videorender;
    running = !VI.Update(rw_image[0..1920*1080*3]);
  }
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
  void Render ( ) {
    import gui.gui : Imgui_Render;
    image_reset |= Imgui_Render(material, camera);
    if ( KI.RProcedural_Type == KI.ProceduralType.Texture ) {
      if ( gl_image !is null )
        GL.Render(gl_image.gl_texture, image_unlock_event);
    } else
      GL.Render(image.RRender, image_unlock_event);
  }
}

private DTOADQ dtoadq;

void Initialize () {
  // initialize opengl, opencl, images
  GL.Renderer_Initialize();
  OCL.Initialize();
  DIMG.Initialize();
  // initialize DTOADQ
  dtoadq = new DTOADQ();
}

void Set_Time ( float t ) {dtoadq.timer = t;}
auto RAllow_Time_Change_Ptr ( ) { return &dtoadq.allow_time_change; }
float RTime ( ) { return dtoadq.timer; }
void Add_Time ( float t ) {
  if ( dtoadq.allow_time_change )
    dtoadq.timer += t;
}

void Update () {
  if ( !dtoadq.running ) return;
  if ( KI.RProcedural_Type == KI.ProceduralType.Texture )
    return dtoadq.Update_Texture();
  if ( KI.RKernel_Type == KI.KernelType.VideoRender )
    return dtoadq.Update_Render();
  dtoadq.Update_DTOADQ();
}
void Render ( ) { dtoadq.Render(); }
void Clean_Up ( ) { OCL.Clean_Up(); }


void Set_Image_Buffer ( DIMG.Resolution resolution ) {
  dtoadq.Set_Image_Buffer(resolution);
}

bool RRunning ( ) { return dtoadq.running; }

auto RDebug_Vals_Ptr ( ) { return dtoadq.debug_vals.ptr; }
