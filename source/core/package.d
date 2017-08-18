module core;
public static import kernel = core.kernel,
                     renderer = core.renderer, config = core.config,
                     info = core.info;
public import core.image;
static import ocl, stl, glfw, gui;
static import core.image;
import core.info : RType, KernelType;
static import shared_info = core.shared_info;

void Initialize () {
  // initialize config and glfw
  static import configuration, gui;
  glfw.Initialize();
  // initialize opengl, opencl, images
  ocl.Initialize();
  core.image.Initialize();
  // gui
  // gui.Initialize();
  configuration.Configure();
}

private void Update_Kernel ( ) {
  if ( !kernel.update_kernel ) return;
  if ( RType == KernelType.Texture     )
    return kernel.Update_Texture();
  if ( RType == KernelType.VideoRender )
    return kernel.Update_Video_Render();
  // both raytrace and dtq
  kernel.Update_DTOADQ;
}

void Update ( ) {
  // update imgui/glfw3
  glfw.Update();
  // update timer
  static float curr_time = 0.0f, prev_time = 0.0f;
  curr_time = glfw.RTime();
  if ( kernel.update_timer )
    Set_Time(RTime() + curr_time - prev_time);
  prev_time = curr_time;
  // update kernel
  Update_Kernel();
}

void Render ( ) {
  static import derelict.opengl3.gl3, derelict.imgui.imgui;
  // clear screen -> process gui/render gl -> render imgui -> swap buffer
  glfw.Start_Buffer();
  renderer.Render();
  glfw.Finish_Buffer();
}

void Clean_Up ( ) {
  glfw.Clean_Up();
  ocl.Clean_Up();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

///
void Set_Time ( float t ) { shared_info.timer = t;    }
///
float RTime   (         ) { return shared_info.timer; }
///

///
void Set_Material(ocl.Material[] material) {
  shared_info.material = material.dup;
}
///
void Set_Material(stl.json.JSONValue json_value) {
  import stl : to;
  ocl.Material[] materials;
  foreach ( json; json_value["materials"].array ) {
    materials ~= ocl.Material(json["diffuse"         ].str.to!float,
                              json["specular"        ].str.to!float,
                              json["glossy"          ].str.to!float,
                              json["retroreflective" ].str.to!float,
                              json["transmittive"    ].str.to!float);
  }
  Set_Material(materials);
}

/// Sets image buffer (camera, gl images, etc) based off resolution.
void Set_Image_Buffer ( core.image.Resolution resolution, bool force = false ) {
  shared_info.Set_Image_Buffer(resolution, force);
}

/// For use with imgui
bool* RUpdate_Timer_Ptr ( ) { return &kernel.update_timer; }
/// For use with imgui
bool* RUpdate_Kernel_Ptr( ) { return &kernel.update_kernel; }
/// For use with imgui
float* RDebug_Vals_Ptr ( ) { return kernel.debug_values.ptr; }

bool RRunning ( ) { return kernel.running; }

int RKernel_Var(info.KernelVar var) {
  return info.RVar(var);
}
