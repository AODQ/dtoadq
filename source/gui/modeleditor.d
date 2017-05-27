module gui.modeleditor;
import globals, derelict.imgui.imgui, gui.gui;

private float[3][200] low, high;

static this ( ) {
  foreach ( ref l; low  ) foreach ( ref t; l ) t = 0.0f;
  foreach ( ref l; high ) foreach ( ref t; l ) t = 1.0f;
}

bool Update_Model ( ref bool open ) {
  auto model = KI.RModel_Info_Ptr;
  bool change = false;
  igBegin("Model Editor", &open);
    igText(model.name.toStringz);
    igSeparator();
    foreach ( it, ref dat; model.params ) {
      gdText(" --- ", dat.label);
      switch ( dat.type ) {
        default: assert(0);
        case KI.ModelInfo.Data.TFloat:
          writeln("IT: ", it);
          gdInputFloat("low", it, 0, low[it][0]); igSameLine();
          gdInputFloat("hi", it, 0, high[it][0]);
          change |= gdSlider(dat.label, dat.tfloat, low[it][0], high[it][0]);
        break;
        case KI.ModelInfo.Data.TFloat2:
        case KI.ModelInfo.Data.TFloat3:
          foreach ( fit, ref val; dat.tfloatarr ) {
            gdInputFloat("low", it, fit, low[it][fit]);
            gdInputFloat("hi", it, fit, high[it][fit]);
            string str = dat.label;
            switch ( fit ) {
              default: str ~= " " ~ fit.to!string; break;
              case 0: str ~= " x"; break;
              case 1: str ~= " y"; break;
              case 2: str ~= " z"; break;
              case 3: str ~= " w"; break;
            }
            change |= gdSlider(str, val, low[it][fit], high[it][fit]);
          }
        break;
        case KI.ModelInfo.Data.TInt:
        break;
      }
      igSeparator();
    }
  igEnd();
  if ( change ) KI.Set_Recompile();
  return change;
}
