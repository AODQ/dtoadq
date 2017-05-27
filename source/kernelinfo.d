module kernelinfo;
enum ProceduralType { Function, Texture, Model, Scene }
enum FileType       { CL, JSON }
FileType String_To_FileType ( string str ) {
  switch ( str ) {
    default: assert(0, str ~ " filetype not supported");
    case "cl":   return FileType.CL;
    case "json": return FileType.JSON;
  }
}

// -- kernel compile shit --
struct KernelInfo {
  enum Type { Raycast, Raytrace, MLT, }
  enum Flag { Show_Normals,           }
  enum Var  { March_Dist, March_Reps, }

  Type type;
  bool[Flag] flags;
  int [Var ] vars;
  string filename;
  ProceduralType procedural_type;
  FileType       file_type;
}

private KernelInfo kernel_info;
private bool        recompile = true;

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
  recompile = true;
  kernel_info.filename        = filename;
  kernel_info.procedural_type = ptype;
  kernel_info.file_type       = ftype;
}

auto RProcedural_Type ( ) { return kernel_info.procedural_type; }
auto RFile_Type ( )  { return kernel_info.file_type; }

auto Should_Recompile ( bool silent = true ) {
  auto ret = recompile;
  recompile = recompile&silent;
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
      [ Var.March_Dist : 64, Var.March_Reps : 128 ]
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
    string label;
    int type;
    auto To_String ( ) {
      import std.conv : to;
      import opencl : To_OpenCL_Float, To_OpenCL_Float_Array;
      switch ( type ) {
        default: assert(0, "Unknown type");
        case TInt:    return tint   .to!string;
        case TFloat:  return tfloat.To_OpenCL_Float;
        case TFloat2: case TFloat3:
          return tfloatarr.To_OpenCL_Float_Array;
      }
    }
    this ( string type_, string label_ ) {
      label = label_;
      import std.stdio;
      writeln("TYPE: ", type_);
      writeln("LABEL: ", label_);
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

private string Parse_Model_Func ( string filename, inout string data ) {
  // -- first create model info --
  import globals, std.string, std.regex;
  filename = filename[filename.lastIndexOf('/')+1 .. $];
  filename = filename[0 .. filename.indexOf('.')];
  if ( model_info.name == filename ) goto PLSDELETELATER;
  model_info.name = filename;
  model_info.params.length = 0;
  // we have to find the declaration of the function ..
  // %s filename ( %s, ... ) {
  {
  string extract_func_sig = // 1 = return, 2 = params
    r"(\w+)\s*" ~ filename ~ r"\s*\(\s*([a-zA-Z0-9,\s]+)\)";
  string extract_params = // 1 = type, 2 = name (have to insert , after)
    r"\s*(\w+)\s*(\w+)\s*,";
  auto func_sig = matchFirst(data, extract_func_sig);
  if ( func_sig.empty() ) {
    writeln("could not parse function signature for: ", filename);
    return "";
  }
  auto matches = matchAll(func_sig[2]~",", extract_params);
  while ( !matches.empty ) {
    string type  = matches.front[1];
    string label = matches.front[2];
    matches.popFront();
    model_info.params ~= ModelInfo.Data(type, label);
    writeln("MODEL INFO: ", model_info.params[$-1]);
  }
  }

  PLSDELETELATER:
  // -- then set up the model for rendering --
  import functional;
  string output = "res = opU(a, res, (float2)(%s(origin + %s), 0.0f));";
  string args = model_info.params.map!(n => n.To_String() ~ ",")
                          .tee!(n => writeln("RESULTS: ", n))
                          .joiner.array[0..$-1].to!string;
  return output.format(filename, args);
}
private string Parse_Map_Func ( string filename, inout string data ) {
  switch ( kernel_info.procedural_type ) {
    case ProceduralType.Model: return Parse_Model_Func(filename, data);
    default: assert(0, "Not supported yet...");
  }
}

void Set_Recompile ( ) { recompile = true; }

string Parse_Kernel ( ) {
   import globals, opencl_kernel, std.string : replace;
   static import std.file;
   string kernel = DTOADQ_kernel;
   string map_func, model_func;
   if ( std.file.exists(kernel_info.filename) ) {
     model_func = std.file.read(kernel_info.filename).to!string;
     map_func = Parse_Map_Func(kernel_info.filename, model_func);
   } else
     map_func = q{
       res = opU(a, res, (float2)(length(origin)-1.25f, 0.0f));
     };
   {
     alias KV = KernelInfo.Var;
     kernel = kernel
       .replace("//%MODEL", model_func)
       .replace("//%MAP", map_func)
       .replace("//%KERNELTYPE", RKernel_Type_String)
       .replace("//%MARCH_DIST", kernel_info.vars[KV.March_Dist].to!string)
       .replace("//%MARCH_REPS", kernel_info.vars[KV.March_Reps].to!string);
     writeln("KERNEL: ", kernel);
   }
   if ( RKernel_Info.flags[KernelInfo.Flag.Show_Normals] ) {
     kernel = kernel.replace("//#define SHOW_NORMALS", "#define SHOW_NORMALS");
   }
   return kernel;
 }
