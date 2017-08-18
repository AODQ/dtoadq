module gui;
static import ocl, core.info;
import derelict.imgui.imgui, std.string : toStringz;

bool Imgui_Render ( ref ocl.Material[] materials, ref ocl.Camera camera ) {
  static import gui.render;
  return gui.render.Imgui_Render ( materials, camera );
}

// -- allow variadic calls to igFN, ei, igText(gui.Accum("lbl: ", l), ..);
auto Accum(T...)(T t) {
  import std.conv : to;
  string res = "";
  foreach ( i; t ) {
    res ~= to!string(i);
  }
  return res.toStringz;
}

static import ocl;
bool CLFloat3_Colour_Edit(T...)(T t, ref ocl.cl_float4 vec ) {
  float[3] arr = [ vec[0], vec[1], vec[2] ];
  bool change = igColorEdit3(t.Accum.toStringz, arr);
  vec = ocl.To_CLFloat3(arr);
  return change;
}
