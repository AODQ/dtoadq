module gui.modules;
static import stl, ocl, core.info, core;
import derelict.imgui.imgui;
import std.conv : to; import std.string : format;


void Render_Materials ( ref ocl.Material[] materials, ref bool change ) {
  string RMaterial_Mixin ( ) {
    import std.traits, std.string : format;
    auto member_names = [ __traits(derivedMembers, ocl.Material) ];
    string mix = `foreach ( i; 0 .. materials.length ) {`;
    mix ~= `if ( igTreeNode(gui.Accum("Material ", i)) ) {`;
      mix ~= `auto m = materials.ptr + i;`;
      foreach ( member; member_names ) {
        mix ~= q{change |= igSliderFloat(
          gui.Accum("(", i.to!string, "): %s"),
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
