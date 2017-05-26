module gui.gui;
import globals;
import derelict.imgui.imgui;
import scene, camera : Camera;
import gui.nodegraphrenderer : Update_Node_Graph;
static import Kernel = kernelinfo;
static import DTOADQ = dtoadq;

bool Imgui_Render ( ref Material[] materials, ref Camera camera ) @trusted {
  bool change;
  bool close;
  igBegin("Material Properties", &close);
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
  bool closecam;
  import functional;
  if ( igCollapsingHeader("Statistics") ) {
    igText("FPS: %.3f ms/frame (%.1f FPS)", 1000.0f / igGetIO().Framerate,
                                                      igGetIO().Framerate);
    gdText("Camera Position ", camera.position[0..3]);
    gdText("Camera Angle    ", camera.lookat[0..3]
                              .map!(n => cast(int)(n*100.0f)/100.0f));
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
    alias Kernel_Type = Kernel.Kernel_Info.Type,
          Kernel_Flag = Kernel.Kernel_Info.Flag,
          Kernel_Var  = Kernel.Kernel_Info.Var;
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


  bool menuthingasdf = true;
  Update_Node_Graph();

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


bool gdButton(T...)(T t) {
  return igButton(t.Accumulator.toStringz);
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

bool gdInputText(T...)(T t, ref string f) {
  return igInputText(t.Accumulator.toStringz, f.ptr);
}

bool gdNewWindow(T...)(T id) {
  bool t = true;
  igBegin(Accumulator(id).toStringz, &t);
  return t;
}

void gdInputText(T...)(T t, char* input_ptr, size_t input_length) {
  igInputText(t.Accumulator.toStringz, input_ptr, input_length, 0);
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
