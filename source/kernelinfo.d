module kernelinfo;
import globals;
static import OCL = opencl;

enum ProceduralType { Texture, Scene };
enum KernelType { Raytrace, MLT, VideoRender };
enum KernelVar  { March_Dist, March_Reps, March_Acc }
private struct KernelInfo {
  string filename;

  KernelType type;
  int [KernelVar] vars;
  ProceduralType procedural_type;
  bool recompile;
}

private KernelInfo kernel_info;

ProceduralType String_To_ProceduralType ( string str ) {
  auto ind = str.indexOf(".");
  if ( ind != -1 ) str = str[ind+1 .. $];
  switch ( str ) {
    default:
      writeln("'", str, "' filetype not supported");
    return ProceduralType.Scene;
    case "dtq": return ProceduralType.Scene;
    case "txt": return ProceduralType.Texture;
  }
}

auto Set_Kernel_Type ( KernelType type ) {
  kernel_info.recompile = true;
  kernel_info.type = type;
}

auto Set_Kernel_Var  ( KernelVar var, int value ) {
  kernel_info.recompile = true;
  kernel_info.vars[var] = value;
}

void Set_Recompile ( ) { kernel_info.recompile = true; }
void Clear_Recompile ( ) { kernel_info.recompile = false; }

void Set_Map_Function ( string filename ){
  kernel_info.recompile = true;
  kernel_info.filename  = filename;
  kernel_info.procedural_type = String_To_ProceduralType(filename);
}

auto RVar ( KernelVar var ) { return kernel_info.vars [var ]; }
auto RType            () { return kernel_info.type;            }
auto RFilename        () { return kernel_info.filename;        }
auto RProcedural_Type () { return kernel_info.procedural_type; }
auto RKernel_Type     () { return kernel_info.type;            }

auto Should_Recompile ( bool silent = true ) {
  auto ret = kernel_info.recompile;
  kernel_info.recompile = ret&silent;
  return ret;
}

auto Configure ( ) {
  static import config;
  kernel_info.filename = config.RDefault_Project();
}

static this ( ) {
  with ( KernelInfo ) {
    kernel_info = KernelInfo(
      "projects/globals/defaultscene.dtq",
      KernelType.MLT,
      [ KernelVar.March_Dist : 355, KernelVar.March_Reps : 256,
       KernelVar.March_Acc : 10 ],
      ProceduralType.Scene
    );
    Set_Recompile();
  }
}
