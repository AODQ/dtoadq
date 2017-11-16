module gui.render;
static import ocl, core.info, gui, core;
import derelict.imgui.imgui;
import gui.modules;

bool Imgui_Render ( ref ocl.Camera camera ) {
  bool change;
  igBegin("Project Details");
    igText("FPS: %.3f ms/frame (%.1f FPS)", 1000.0f / igGetIO().Framerate,
                                                      igGetIO().Framerate);
    import functional;
    if ( igCollapsingHeader("Statistics") ) {
      igText("CAMERA POSITION --");
      foreach ( info; zip(iota(0, 3), ["X", "Y", "Z"]) )
        change |= igInputFloat(info[1].ptr, &camera.position[info[0]]);
      igText(gui.Accum("Camera Angle",
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
        kernel_type = core.info.RType();
        march_dist = core.info.RVar(core.info.KernelVar.March_Dist);
        march_reps = core.info.RVar(core.info.KernelVar.March_Reps);
        march_acc  = core.info.RVar(core.info.KernelVar.March_Acc)/1000.0f;
      }
      pkernel_type = kernel_type;
      { // render type
        alias KT = core.info.KernelType;
        igRadioButton("DTQ", &kernel_type, cast(int)KT.DTQ);      igSameLine();
        igRadioButton("RT",  &kernel_type, cast(int)KT.Raytrace); igSameLine();
        igRadioButton("RC",  &kernel_type, cast(int)KT.Raycast);
      }
      if ( kernel_type != pkernel_type ) {
        core.info.Set_Kernel_Type(cast(core.info.KernelType)kernel_type);
      }

      if ( igSliderInt("March Distance", &march_dist, 1, 1024) ) {
        core.info.Set_Kernel_Var(core.info.KernelVar.March_Dist, march_dist);
      }
      if ( igSliderInt("March Repetitions", &march_reps, 1, 1024) ) {
        core.info.Set_Kernel_Var(core.info.KernelVar.March_Reps, march_reps);
      }
        if ( igSliderFloat("March Accuracy", &march_acc, 0.00f, 0.1f) ) {
        int acc_var = cast(int)(march_acc*1000);
        core.info.Set_Kernel_Var(core.info.KernelVar.March_Acc, acc_var);
      }

      static import core.image;
      foreach ( it, str; core.image.Resolution_Strings ) {
        igRadioButton(str.toStringz, &resolution, cast(int)it);
      }

      core.Set_Image_Buffer(cast(core.image.Resolution)resolution);
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
      Render_Materials(change);

    static bool open_editor = false;
    igCheckbox("Editor", &open_editor);
    if ( open_editor ) {
      Render_Editor(change);
    }
    igSeparator();

    static bool small_stats_open = false;
    igCheckbox("SmallStats", &small_stats_open);
    if ( small_stats_open ) {Render_SmallStats(camera);}
    igSeparator();
  igEnd();

  return change;
}
