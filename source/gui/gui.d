module gui.gui;
import globals;
import derelict.imgui.imgui;
import scene, camera : Camera;
import gui.nodegraphrenderer : Update_Node_Graph;
static import Kernel = kernelinfo;
static import DTOADQ = dtoadq;

private void Render_Materials ( ref Material[] materials, ref bool change ) {
    igBegin("Materials");
      foreach ( i; 0 .. materials.length ) {
        auto m = &materials[i];
        import std.conv : to;
        auto istr = i.to!string;
        change |= CLFloat3_Colour_Edit("Base Colour", istr, m.base_colour);
        change |= gdSliderNorm("Metallic",        istr, m.metallic);
        change |= gdSliderNorm("Subsurface",      istr, m.subsurface);
        change |= gdSliderNorm("Specular",        istr, m.specular);
        change |= gdSliderNorm("Specular_tint",   istr, m.specular_tint);
        change |= gdSliderNorm("Roughness",       istr, m.roughness);
        change |= gdSliderNorm("Anisotropic",     istr, m.anisotropic);
        change |= gdSliderNorm("Sheen",           istr, m.sheen);
        change |= gdSliderNorm("Sheen_tint",      istr, m.sheen_tint);
        change |= gdSliderNorm("Clearcoat",       istr, m.clearcoat);
        change |= gdSliderNorm("Clearcoat_gloss", istr, m.clearcoat_gloss);
        change |= gdSliderNorm("Emission",        istr, m.emission);
      }
    igEnd();
}

bool Imgui_Render ( ref Material[] materials, ref Camera camera ) {
  bool change;
  igBegin("Project Details");
    import functional;
    if ( igCollapsingHeader("Statistics") ) {
      igText("FPS: %.3f ms/frame (%.1f FPS)", 1000.0f / igGetIO().Framerate,
                                                        igGetIO().Framerate);
      gdText("CAMERA POSITION --");
      change |= gdInputFloat("X", camera.position[0]);
      change |= gdInputFloat("Y", camera.position[1]);
      change |= gdInputFloat("Z", camera.position[2]);
      gdText("Camera Angle    ", camera.lookat[0..3]
                                .map!(n => cast(int)(n*100.0f)/100.0f));
      change |= gdSlider("FOV", camera.fov, 50.0f, 140.0f);
    }
    // -- render options --
    if ( igCollapsingHeader("Render Options") ) {
      static bool Show_Normals = false;
      static int kernel_type = 0, resolution = 0, pkernel_type,
                march_dist = 64,
                march_reps = 128;

      pkernel_type = kernel_type;
      igRadioButton("Raycast",  &kernel_type, 0); igSameLine();
      igRadioButton("Raytrace", &kernel_type, 1); igSameLine();
      igRadioButton("MLT",      &kernel_type, 2);
      alias Kernel_Type = Kernel.KernelInfo.Type,
            Kernel_Flag = Kernel.KernelInfo.Flag,
            Kernel_Var  = Kernel.KernelInfo.Var;
      if ( kernel_type != pkernel_type ) {
        Kernel.Set_Kernel_Type(cast(Kernel_Type)kernel_type);
      }

      if ( igCheckbox("Show Normals", &Show_Normals) ) {
        Kernel.Set_Kernel_Flag(Kernel_Flag.Show_Normals, Show_Normals);
      }

      if ( igSliderInt("March Distance", &march_dist, 1, 2048) ) {
        Kernel.Set_Kernel_Var(Kernel_Var.March_Dist, march_dist);
      }
      if ( igSliderInt("March Repetitions", &march_reps, 1, 256) ) {
        Kernel.Set_Kernel_Var(Kernel_Var.March_Reps, march_reps);
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
      final switch ( KI.RProcedural_Type ) {
        case KI.ProceduralType.Function: break;
        case KI.ProceduralType.Texture: break;
        case KI.ProceduralType.Model:
          import gui.modeleditor;
          change |= Update_Model(open_editor);
        break;
        case KI.ProceduralType.Scene: break;
      }
    }

    static bool open_node_graph = false;
    igCheckbox("Node Graph", &open_node_graph);
    if ( open_node_graph ) {
      Update_Node_Graph();
    }
  igEnd();

  return change;
}

private string Accumulator(T...)(T t) {
  string res = "";
  foreach ( i; t ) {
    res ~= to!string(i);
  }
  return res;
}
void gdText(T...)(T t) {
  igText(t.Accumulator.toStringz);
}

void gdTextColour(T...)(T t, float r, float g, float b) {
  igTextColored(ImVec4(r, g, b, 1.0f), t.Accumulator.toStringz);
}


bool gdButton(T...)(T t) {
  return igButton(t.Accumulator.toStringz);
}

bool gdSlider(T...)(T t, ref float f, float low, float hi) {
  return igSliderFloat(t.Accumulator.toStringz, &f, low, hi);
}
bool gdSliderNorm(T...)(T t, ref float f) {
  return igSliderFloat(t.Accumulator.toStringz, &f, 0.0f, 1.0f);
}

bool gdInputFloat(T...)(T t, ref float f) {
  return igInputFloat(t.Accumulator.toStringz, &f);
}

bool gdInputInt(T...)(T t, ref int f) {
  return igInputInt(t.Accumulator.toStringz, &f);
}

bool gdNewWindow(T...)(T id) {
  bool t = true;
  igBegin(Accumulator(id).toStringz, &t);
  return t;
}

bool gdInputText(T...)(T t, ref string input) {
  char[] iptr = input.dup ~ '\0';
  iptr.length = 20;
  auto change = igInputText(t.Accumulator.toStringz, iptr.ptr, iptr.length, 0);
  input = iptr.ptr.fromStringz.idup;
  return change;
}

auto gdCalcTextSize ( string str ) {
  // -- The igCalcTextSize function call wasn't working so i just implemented
  //    it myself rather boring with cimgui again ... --
  return ImVec2(str.length*13, 14);
}

auto gdRMousePos ( ) {
  ImVec2 vec;
  igGetMousePos(&vec);
  return vec;
}

auto gdRMouseDelta ( ) {
  ImVec2 vec;
  igGetMouseDragDelta(&vec, 0, 1.0f);
  return vec;
}

auto gdMenuItem(string name, string shortcut = "", bool selected = false,
                bool enabled = true) {
  return igMenuItem(name.toStringz, shortcut.toStringz, selected, enabled);
}

import opencl : cl_float4, To_CLFloat3;
bool CLFloat3_Colour_Edit(T...)(T t, ref cl_float4 vec ) {
  float[3] arr = [ vec[0], vec[1], vec[2] ];
  bool change = igColorEdit3(t.Accumulator.toStringz, arr);
  vec = To_CLFloat3(arr);
  return change;
}
