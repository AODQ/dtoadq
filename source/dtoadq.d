module dtoadq;
import globals;
import oclstructs;
static import KI   = kernelinfo;
static import DIMG = dtoadqimage;
static import OCL  = opencl;
static import GL   = gl_renderer;
static import stl;
import derelict.opencl.cl : cl_event;
import parser : Reparse_Files;

private class DTOADQ {
  // -- kernel arguments
  DIMG.Image image;
  DIMG.GLImage gl_image; // only for texture
  ubyte[] rw_image;
  int[] rng;
  SharedInfo shared_info;
  Camera camera;
  Material[] material;
  float timer = 0.0f;
  // -- opengl/opencl events
  cl_event image_unlock_event;
  // -- kernel debug
  float[] debug_vals;
  bool update_timer = false,
       update_kernel = true;
  bool running = true;

  // -- funcs
  this ( ) {
    debug_vals = [ 0.0f, 0.0f, 0.0f ];
    camera = Construct_Camera([2.5f, 0.0f, 6.0f], [-0.06f, 0.4f, -1.0f],
                                    [image.x, image.y]);
    shared_info.finished_samples = 0;
    material = [
      Material(0.0f, 1.0f, 0.0f, 0.0f, 0.0f),
      Material(0.0f, 1.0f, 0.0f, 0.0f, 0.0f),
      Material(0.0f, 1.0f, 0.0f, 0.0f, 0.0f),
      Material(1.0f, 0.0f, 0.0f, 0.0f, 0.0f),
      Material(0.0f, 1.0f, 0.0f, 0.0f, 0.0f),
      Material(0.0f, 1.0f, 0.0f, 0.0f, 0.0f),
      Material(0.0f, 1.0f, 0.0f, 0.0f, 0.0f),
      Material(0.0f, 1.0f, 0.0f, 0.0f, 0.0f),
      Material(0.0f, 1.0f, 0.0f, 0.0f, 0.0f)
    ];
    Set_Image_Buffer(DIMG.Resolution.r640_360, true);
  }

  void Set_Image_Buffer ( DIMG.Resolution resolution, bool force = false ) {
    if ( force || image.resolution != resolution ) {
      shared_info.clear_img = true;
      image = DIMG.RImage(resolution);
      rw_image.length = image.x*image.y*4;
      camera.dimensions.x = image.x;
      camera.dimensions.y = image.y;
    }
  }

  bool Compile ( string kernel_name ) {
    import parser;
    return OCL.Compile(Parse_Kernel(), kernel_name);
  }

  void Update_DTOADQ ( ) {
    auto imgs = DIMG.RImages();
    // -- buffer --
    bool should_reparse = Reparse_Files();
    // -- camera/gui/etc --
    static float previous_timer = 0.0f;
    shared_info.clear_img |= camera.Update|KI.Should_Recompile()|should_reparse;
    shared_info.clear_img |= stl.abs(previous_timer - timer) > 0.0f;
    previous_timer = timer;
    // -- kernel --
    if ( KI.Should_Recompile(false) || should_reparse ) {
      import std.datetime;
      if ( Compile("DTOADQ_Kernel") ) {
        writeln("Recompiled @ ", Clock.currTime);
      } else {
        writeln("Failed to compiled @ ", Clock.currTime);
      }
      return;
    }
    // -- random numbers --
    {
      import stl;
      rng = iota(0, 32).map!(n => n = uniform(-int.max, int.max)).array;
    }
    //
    image.Lock();
    OCL.Run(
      OCL.CLStoreMem(rw_image),
      image.RWrite(),
      OCL.CLStoreMem(shared_info),
      camera,
      timer,
      OCL.CLPredefinedMem(imgs),
      material,
      debug_vals,
      rng,
      // kernel size
      image.x, image.y
    );
    image.Unlock();
    shared_info.clear_img = false;
  }
  void Update_Render ( ) {
    Update_DTOADQ();
    static import VI = videorender;
    auto len = DIMG.RResolution_Length(VI.RRender_Resolution);
    if ( shared_info.finished_samples >= len ) {
      running = !VI.Update(rw_image[0..len*3]);
      shared_info.finished_samples = 0;
      shared_info.clear_img = true;
    }
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
    shared_info.clear_img |= Imgui_Render(material, camera);
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
auto RUpdate_Timer_Ptr ( ) { return &dtoadq.update_timer; }
auto RUpdate_Kernel_Ptr     ( ) { return &dtoadq.update_kernel; }
float RTime ( ) { return dtoadq.timer; }
void Add_Time ( float t ) {
  if ( dtoadq.update_timer )
    dtoadq.timer += t;
}

void Update () {
  if ( !dtoadq.running ) return;
  if ( !dtoadq.update_kernel ) return;
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
