module dtoadq;
import derelict.opencl.cl : cl_event;
import parser : Reparse_Files;
static import core.shared_info;

private cl_event image_unlock_event;
private float[] debug_values;
private bool update_timer  = false,
             update_kernel = true;
private bool running = true;

private void Initialize ( ) {
  debug_vals = [ 0.0f, 0.0f, 0.0f ];
  camera = opencl.Construct_Camera([2.5f, 0.0f, 6.0f], [-0.06f, 0.4f, -1.0f],
                                  [image.x, image.y]);
  shared_info.finished_samples = 0;
  material = [];
  Set_Image_Buffer(DIMG.Resolution.r160_140, true);
}

  // -- funcs
  this ( ) {
  }

  void Set_Image_Buffer ( DIMG.Resolution resolution, bool force = false ) {
    core.shared_info.Set_Image_Buffer(resolution, force);
  }

  bool Compile ( string kernel_name ) {
    import parser;
    return OCL.Compile(Parse_Kernel(), kernel_name);
  }

  void Update_DTOADQ ( ) {
    auto imgs = core.image.RImages();
    bool should_reparse = Reparse_Files() | core.shared_info.Update_Buffer();
    // -- kernel --
    if ( kernel.info.Should_Recompile(false) || should_reparse ) {
      import std.datetime;
      stl.writeln("Compile " ~ (Compile("DTOADQ_Kernel")?"Success":"Failure") ~
                  ": " ~ Clock.currTime);
      return;
    }
    { // Run kernel
      static import info = core.shared_info;
      info.image.Lock();
      OCL.Run(
        OCL.CLStoreMem(info.rw_image),
        info.image.RWrite(),
        OCL.CLStoreMem(info.image_meta_data),
        info.camera,
        info.timer,
        OCL.CLPredefinedMem(info.imgs),
        info.material,
        info.debug_vals,
        info.rng,
        // kernel size
        info.image.x, info.image.y
      );
      info.image.Unlock();
      info.image_meta_data.clear_img = false;
    }
  }

  void Update_Render ( ) {
    Update_DTOADQ();
    // static import VI = videorender;
    // auto len = DIMG.RResolution_Length(VI.RRender_Resolution);
    // if ( shared_info.finished_samples >= len ) {
    //   running = !VI.Update(rw_image[0..len*3]);
    //   shared_info.finished_samples = 0;
    //   shared_info.clear_img = true;
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
  glfw.Renderer_Initialize();
  ocl.Initialize();
  core.image.Initialize();
  // kernelinfo configure
  kernel.info.Configure();
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

void Set_Material(Material[] material) {
  dtoadq.material = material.dup;
}

static import JSON = std.json;
void Set_Material(string json_value) {
  Material[] materials;
  JSON.JSONValue val = JSON.parseJSON(json_value);
  foreach ( json; val["materials"].array ) {
    materials ~= Material(json["diffuse"].str.to!float,
                          json["specular"].str.to!float,
                          json["glossy"].str.to!float,
                          json["retroreflective"].str.to!float,
                          json["transmittive"].str.to!float);
  }
  Set_Material(materials);
}
