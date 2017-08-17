module gui;
import globals;
import derelict.imgui.imgui;
import oclstructs;
static import Kernel = kernelinfo;
static import DTOADQ = dtoadq;
static import KI = kernelinfo;

private void Render_Materials ( ref Material[] materials, ref bool change ) {
  string RMaterial_Mixin ( ) {
    import std.traits, std.string : format;
    auto member_names = [ __traits(derivedMembers, Material) ];
    string mix = `foreach ( i; 0 .. materials.length ) {`;
    mix ~= `if ( igTreeNode(igAccum("Material ", i)) ) {`;
      mix ~= `auto m = materials.ptr + i;`;
      foreach ( member; member_names ) {
        mix ~= q{change |= igSliderFloat(
          igAccum("(", i.to!string, "): %s"),
                            &m.%s, 0.0f, 1.0f);} .format(member, member);
      }
    mix ~= `igTreePop(); }`;
    return mix ~ `}`;
  }

  igBegin("Material Properties");
    import std.conv;
    mixin(RMaterial_Mixin());
  igEnd();
}

bool Imgui_Render ( ref Material[] materials, ref Camera camera ) {
  bool change;
  igBegin("Project Details");
    igText("FPS: %.3f ms/frame (%.1f FPS)", 1000.0f / igGetIO().Framerate,
                                                      igGetIO().Framerate);
    import functional;
    if ( igCollapsingHeader("Statistics") ) {
      igText("CAMERA POSITION --");
      foreach ( info; zip(iota(0, 3), ["X", "Y", "Z"]) )
        change |= igInputFloat(info[1].ptr, &camera.position[info[0]]);
      igText(igAccum("Camera Angle",
        camera.lookat[0..3].map!(n => cast(int)(n*100.0f)/100.0f)));
      change |= igSliderFloat("FOV", &camera.fov, 50.0f, 140.0f);
    }
    // -- render options --
    if ( igCollapsingHeader("Render Options") ) {
      static bool Show_Normals = false;
      static int kernel_type = -1, resolution, pkernel_type,
                 march_dist, march_reps;
      static float march_acc;
      if ( kernel_type == -1 ) {
        kernel_type = KI.RKernel_Type();
        march_dist = KI.RVar(KI.KernelVar.March_Dist);
        march_reps = KI.RVar(KI.KernelVar.March_Reps);
        march_acc  = KI.RVar(KI.KernelVar.March_Acc)/1000.0f;
      }
      pkernel_type = kernel_type;
      igRadioButton("Raytrace", &kernel_type, 0); igSameLine();
      igRadioButton("MLT",      &kernel_type, 1);
      alias Kernel_Type = Kernel.KernelType,
            Kernel_Var  = Kernel.KernelVar;
      if ( kernel_type != pkernel_type ) {
        Kernel.Set_Kernel_Type(cast(Kernel_Type)kernel_type);
      }

      if ( igSliderInt("March Distance", &march_dist, 1, 1024) ) {
        Kernel.Set_Kernel_Var(Kernel_Var.March_Dist, march_dist);
      }
      if ( igSliderInt("March Repetitions", &march_reps, 1, 1024) ) {
        Kernel.Set_Kernel_Var(Kernel_Var.March_Reps, march_reps);
      }
      if ( igSliderFloat("March Accuracy", &march_acc, 0.00f, 0.1f) ) {
        int acc_var = cast(int)(march_acc*1000);
        Kernel.Set_Kernel_Var(Kernel_Var.March_Acc, acc_var);
      }

      static import DIMG = dtoadqimage;
      foreach ( it, str; DIMG.RResolution_Strings ) {
        igRadioButton(str.toStringz, &resolution, cast(int)it);
      }

      DTOADQ.Set_Image_Buffer(cast(DIMG.Resolution)resolution);
    }

    static bool open_file_browser = false;
    igCheckbox("File Browser", &open_file_browser);
    if ( open_file_browser ) {
      static import Files = gui.files;
      Files.Update(open_file_browser);
    }

    static bool open_materials = false;
    igCheckbox("Materials", &open_materials);
    if ( open_materials )
      Render_Materials(materials, change);

    static bool open_editor = false;
    igCheckbox("Editor", &open_editor);
    if ( open_editor ) {
      Editor(change);
    }
    igSeparator();

    static bool small_stats_open = false;
    igCheckbox("SmallStats", &small_stats_open);
    if ( small_stats_open ) {SmallStats(camera);}
    igSeparator();
  igEnd();

  return change;
}

private void SmallStats ( ref Camera camera ) {
  string TS ( float r ) {
    import std.math : round;
    return (round(r*100.0f)/100.0f).to!string;
  }
  igBegin("SmallStats");
    auto p = camera.position,
         a = camera.lookat;
    igText(igAccum("P %s, %s, %s".format(TS(p[0]), TS(p[1]), TS(p[2]))));
    igText(igAccum("A %s, %s".format(TS(a[0]), TS(a[1]))));
  igEnd();
}

private void Editor ( ref bool change ) {
  igBegin("Editor");
    float timer = 0.0f;
    timer = DTOADQ.RTime();
    if ( igDragFloat("Timer", &timer, 0.01f, 0.0f, 120.0f,
                    timer.to!string.toStringz, 1.0f) ) {
      DTOADQ.Set_Time(timer);
    }
    auto update_timer  = DTOADQ.RUpdate_Timer_Ptr,
         update_kernel = DTOADQ.RUpdate_Kernel_Ptr;
    igCheckbox("Update Time", update_timer);
    igCheckbox("Update kernel", update_kernel);
    change |= *update_timer;
    if ( *update_timer ) {
      static float start_timer = 0.0f, end_timer = 120.0f;
      igInputFloat("START", &start_timer);
      igInputFloat("END",   &end_timer);
      if ( timer > end_timer || timer < start_timer )
        DTOADQ.Set_Time(start_timer);
    }
    foreach ( dbg; 0 .. 3 ) {
      auto dbgstr = ("DBG " ~ (cast(char)('X'+dbg)).to!string).toStringz;
      change |=
        igDragFloat(dbgstr, &DTOADQ.RDebug_Vals_Ptr[dbg], 0.01f, -100.0,100.0,
                  DTOADQ.RDebug_Vals_Ptr[dbg].to!string.toStringz, 2.0f);
    }
    igEnd();
}

// -- allow variadic calls to igFN, ei, igText(igAccum("lbl: ", l), ..);
auto igAccum(T...)(T t) {
  string res = "";
  foreach ( i; t ) {
    res ~= to!string(i);
  }
  return res.toStringz;
}

import opencl : cl_float4, To_CLFloat3;
bool CLFloat3_Colour_Edit(T...)(T t, ref cl_float4 vec ) {
  float[3] arr = [ vec[0], vec[1], vec[2] ];
  bool change = igColorEdit3(t.Accum.toStringz, arr);
  vec = To_CLFloat3(arr);
  return change;
}
