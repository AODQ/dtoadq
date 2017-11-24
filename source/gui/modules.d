module gui.modules;
static import stl, ocl, functional, core.info, core;
import derelict.imgui.imgui;
import std.conv : to; import std.string : format;
import std.traits;

void Mat_Normalize(T)(T* important, T*[] x) {
  import std.math : abs;
  float rval = *important;
  // TODO: i'm sure there's a better way, you have to play with it a little at
  //       1.0 unfortunately
  foreach ( cnt; 0 .. 250 ) {
    float len = 0.0f;
    *important = rval;
    // accum
    foreach ( i; x )
      len += *i;
    foreach ( i; x ) {
      *i /= len;
      *i = cast(int)((*i) * 100.0f)/100.0f;
    }
    if ( abs(*important - rval) >= 0.01f ) break;
  }
}

void Render_Materials ( ref bool change ) {
  // return a GUI renderer for materials, with their corresponding
  // attributes
  static import sinfo = core.shared_info;
  import std.typecons;
  import std.traits, std.string : format;
  import functional : zip;

  igBegin("Material Properties");
    import std.conv;
    alias OCLM = ocl.OCLMaterial;
    float* importance_member;
    OCLM*  importance_mat;
    // ---- gui render
    foreach ( i; 0 .. sinfo.material.length ) {
      if ( !igTreeNode(gui.Accum("Material ", i)) ) continue;
      auto m  = &sinfo.material[i],
          om = &sinfo.ocl_material[i];
      float* importance = null;
      // colour
      auto col3_acc = gui.Accum("(%s): albedo".format(i.to!string));
      change |= igColorEdit3(col3_acc, m.albedo);
      // floatz
      float* fmem;
      string acc;
      static foreach(member; stl.AllFilteredAttributes!(OCLM, `"ModifyFloat"`)){
        mixin(q{fmem = &om.%s;}.format(member));
        acc = "(%s): %s".format(i.to!string, member);
        if ( igSliderFloat(gui.Accum(acc), fmem, 0.0f, 1.0f) ) {
          change = true;
          static if (stl.HasAttribute!(OCLM, member, `"Normalize"`)){
            importance_member = fmem;
            importance_mat = om;
          }
        }
      }
      igTreePop();
    }

    if ( importance_mat && importance_member ) {
      float*[] mats;
      static foreach ( member; __traits(allMembers, OCLM) ) {
        static if ( stl.HasAttribute!(ocl.OCLMaterial, member, `"Normalize"`) )
          mixin(q{mats ~= &importance_mat.%s;}.format(member));
      }
      Mat_Normalize(importance_member, mats);
    }
  igEnd();
}

void Render_SmallStats ( ref ocl.Camera camera ) {
  string TS ( float r ) {
    import std.math : round;
    return stl.to!string(round(r*100.0f)/100.0f);
  }
  igBegin("SmallStats");
    auto p = camera.position,
         a = camera.lookat;
    igText(gui.Accum("P %s, %s, %s".format(TS(p[0]), TS(p[1]), TS(p[2]))));
    igText(gui.Accum("A %s, %s".format(TS(a[0]), TS(a[1]))));
  igEnd();
}


void Render_Editor ( ref bool change ) {
  igBegin("Editor");
    float timer = 0.0f;
    timer = core.RTime();
    if ( igDragFloat("Timer", &timer, 0.01f, 0.0f, 120.0f,
                    stl.toStringz(to!string(timer)), 1.0f) ) {
                    core.Set_Time(timer);
    }
    auto update_timer  = core.RUpdate_Timer_Ptr,
         update_kernel = core.RUpdate_Kernel_Ptr;
    igCheckbox("Update Time", update_timer);
    igCheckbox("Update kernel", update_kernel);
    change |= *update_timer;
    if ( *update_timer ) {
      static float start_timer = 0.0f, end_timer = 120.0f;
      igInputFloat("START", &start_timer);
      igInputFloat("END",   &end_timer);
      if ( timer > end_timer || timer < start_timer )
        core.Set_Time(start_timer);
    }
    foreach ( dbg; 0 .. 3 ) {
      auto dbgstr = stl.toStringz("DBG " ~ (cast(char)('X'+dbg)).to!string);
      change |=
        igDragFloat(dbgstr, &core.RDebug_Vals_Ptr[dbg], 0.01f, -100.0,100.0,
                  stl.toStringz(core.RDebug_Vals_Ptr[dbg].to!string), 2.0f);
    }
    igEnd();
}
