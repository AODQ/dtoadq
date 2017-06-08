module kernelinfo;
import globals;
static import OCL = opencl;

enum ProceduralType { Function, Texture, Model, Scene }
enum FileType       { CL, DTQ, TXT }
FileType String_To_FileType ( string str ) {
  switch ( str ) {
    default: assert(0, str ~ " filetype not supported");
    case "cl":  return FileType.CL;
    case "dtq": return FileType.DTQ;
    case "txt": return FileType.TXT;
  }
}

// -- kernel compile shit --
struct KernelInfo {
  enum Type { Raycast, Raytrace, MLT, }
  enum Flag { Show_Normals,           }
  enum Var  { March_Dist, March_Reps, March_Acc }

  string filename;

  Type type;
  bool[Flag] flags;
  int [Var ] vars;
  ProceduralType procedural_type;
  FileType       file_type;

  this ( Type type_, bool[Flag] flags_, int[Var] vars_ ) {
    flags = flags_; vars = vars_; type = type_;
  }
}

private KernelInfo kernel_info;
private bool       recompile     = true,
                   reparse_files = true;

auto Set_Kernel_Type ( KernelInfo.Type type ) {
  recompile = true;
  kernel_info.type = type;
}

auto Set_Kernel_Flag ( KernelInfo.Flag flag, bool status ) {
  recompile = true;
  kernel_info.flags[flag] = status;
}

auto Set_Kernel_Var  ( KernelInfo.Var var, int value ) {
  recompile = true;
  kernel_info.vars[var] = value;
}

void Set_Map_Function ( ProceduralType ptype, FileType ftype, string filename ){
  recompile                   = true;
  reparse_files               = true;
  kernel_info.filename        = filename;
  kernel_info.procedural_type = ptype;
  kernel_info.file_type       = ftype;
}

auto RType ( ) { return kernel_info.type; }
auto RFlag ( KernelInfo.Flag flag ) { return kernel_info.flags[flag]; }
auto RVar  ( KernelInfo.Var  var  ) { return kernel_info.vars [var ]; }
auto RFilename ( ) { return kernel_info.filename; }
auto RProcedural_Type ( ) { return kernel_info.procedural_type; }
auto RFile_Type ( )  { return kernel_info.file_type; }

auto Should_Recompile ( bool silent = true ) {
  auto ret = recompile;
  recompile = recompile&silent;
  return ret;
}

auto Should_Reparse_Files ( bool silent = true ) {
  auto ret = reparse_files;
  reparse_files = reparse_files&silent;
  return ret;
}

immutable(KernelInfo) RKernel_Info ( ) { return cast(immutable)kernel_info; }

auto RKernel_Type_String ( ) {
  import opencl_kernel;
  final switch ( kernel_info.type ) with ( KernelInfo.Type ) {
    case Raycast:  return Raycast_kernel_function;
    case Raytrace: return Raytrace_kernel_function;
    case MLT:      return MLT_kernel_function;
  }
}

static this ( ) {
  with ( KernelInfo ) {
    kernel_info = KernelInfo(
      Type.Raycast,
      [ Flag.Show_Normals : false ],
      [ Var.March_Dist : 64, Var.March_Reps : 128, Var.March_Acc : 1 ]
    );
  }
}

struct ModelInfo {
  struct Data {
    enum { TFloat, TFloat2, TFloat3, TInt }
    union {
      float   tfloat;
      float[] tfloatarr;
      int     tint;
    }
    string label, rlabel;
    int type;
    string str_type; // TODO replace type with str_type
    bool label_special, override_special;

    auto To_String ( bool no_override = false ) {
      import std.conv : to;
      if ( !no_override && label_special && !override_special ) return rlabel;
      switch ( type ) {
        default: assert(0, "Unknown type");
        case TInt:     return tint   .to!string;
        case TFloat:   return OCL.To_OpenCL_Float(tfloat);
        case TFloat2: case TFloat3:
          return OCL.To_OpenCL_Float_Array(tfloatarr);
      }
    }
    this ( string type_, string label_ ) {
      rlabel = label = label_;
      import std.stdio;
      switch ( label ) {
        default: break;
        case "otime":
          rlabel = "time";
        goto case;
        case "time":
          label_special    = true;
          override_special = true;
        break;
      }
      str_type = type_;
      switch ( type_ ) {
        default: assert(0, "Unsupported type: " ~ type_);
        case "float":  type = TFloat;  tfloat = 0.0f;        break;
        case "float2": type = TFloat2; tfloatarr = [0.0f, 0.0f]; break;
        case "float3": type = TFloat3; tfloatarr = [0.0f, 0.0f, 0.0f]; break;
        case "int":    type = TInt;    tint = 0;             break;
      }
    }
  }
  Data[] params;
  string name;
}
ModelInfo model_info;

auto RModel_Info_Ptr ( ) { return &model_info; }

void Set_Recompile ( ) { recompile = true; }
