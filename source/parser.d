module parser;
import globals, kernelinfo;

private bool[string] parsed_files;

private string Parse_Model_Func ( string filename, inout string data ) {
  // -- first boolcreate model info --
  import globals, std.string, std.regex;
  parsed_files[filename] = true;
  filename = filename[filename.lastIndexOf('/')+1 .. $];
  filename = filename[0 .. filename.indexOf('.')];
  model_info.name = filename;
  model_info.params.length = 0;
  // we have to find the declaration of the function ..
  // %s filename ( %s, ... ) {
  string extract_require_sig = r"!require\s+([\w\.\s]+);";
  string extract_func_sig = // 1 = return, 2 = params
    r"(\w+)\s*" ~ filename ~ r"\s*\(\s*([a-zA-Z0-9,\s]+)\)";
  string extract_params = // 1 = type, 2 = name (have to insert , after)
    r"\s*(\w+)\s*(\w+)\s*,";
  // -- require --
  {
    import functional;
    auto ma = data.matchAll(extract_require_sig);
    if ( ma.exists ) {
      string[] req_files = ma.front[1].split(" ").filter!(n => n != "").array;
      foreach ( fil; req_files ) {
        if ( fil !in parsed_files ) {
          static import std.file;
          Parse_Model_Func(fil, std.file.read(file).to!string);
        }
      }
    }
  }

  // -- signature --
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


private string Parse_Map_Func ( string filename, inout string data ) {
  switch ( RProcedural_Type ) {
    case ProceduralType.Model: return Parse_Model_Func(filename, data);
    default: assert(0, "Not supported yet...");
  }
}

string Parse_Kernel ( ) {
   import globals, opencl_kernel, std.string : replace;
   static import std.file;
   string kernel = DTOADQ_kernel;
   string map_func, model_func;
   parsed_files = ["":false];
   auto filename = RFilename;
   if ( std.file.exists(filename) ) {
     model_func = std.file.read(filename).to!string;
     map_func = Parse_Map_Func(filename, model_func);
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
       .replace("//%MARCH_DIST", RVar(KV.March_Dist).to!string)
       .replace("//%MARCH_REPS", RVar(KV.March_Reps).to!string);
   }
   if ( RKernel_Info.flags[KernelInfo.Flag.Show_Normals] ) {
     kernel = kernel.replace("//#define SHOW_NORMALS", "#define SHOW_NORMALS");
   }
   return kernel;
 }
