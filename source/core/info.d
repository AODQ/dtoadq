module core.info;
static import stl, ocl, parser;

enum KernelType { DTQ, Raytrace, Raycast, Texture, VideoRender };
enum KernelVar  { March_Dist, March_Reps, March_Acc }
private struct KernelInfo {
  string filename;

  KernelType type;
  int [KernelVar] vars;
  bool recompile;
}

private KernelInfo kernel_info;

void Set_Kernel_Type ( KernelType type ) {
  kernel_info.recompile = true;
  kernel_info.type = type;
  parser.Change_Kernel(type);
}

void Set_Kernel_Var  ( KernelVar var, int value ) {
  kernel_info.recompile = true;
  kernel_info.vars[var] = value;
}

void Set_Recompile ( ) { kernel_info.recompile = true; }
void Clear_Recompile ( ) { kernel_info.recompile = false; }

void Set_Map_Function ( string filename ){
  Set_Kernel_Type(RType);
  kernel_info.recompile = true;
  kernel_info.filename  = filename;
  parser.Change_Scene(filename);
}

void Set_Texture_Function ( string filename ){
  Set_Kernel_Type(KernelType.Texture);
  kernel_info.recompile = true;
  kernel_info.filename  = filename;
  parser.Change_Scene(filename);
}

void Set_Any_Function ( string filename ) {
  switch ( stl.RMatching_File_Extensions(filename, "dtq", "txt") ) {
    default: assert(false, "Unknown file type");
    case "dtq": Set_Map_Function(filename); break;
    case "txt": Set_Texture_Function(filename); break;
  }
}

int        RVar ( KernelVar var ) { return kernel_info.vars [var ]; }
KernelType RType    ()            { return kernel_info.type;        }
string     RFilename()            { return kernel_info.filename;    }

bool Should_Recompile ( bool silent = true ) {
  bool ret = kernel_info.recompile;
  kernel_info.recompile = ret&silent;
  return ret;
}
