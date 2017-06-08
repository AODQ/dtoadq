module parser;
import globals;
static import KI = kernelinfo;
static import DIMG = dtoadqimage;
static import std.file;
import std.datetime : SysTime;

MapParser[string] parsed_files;
SceneParser       parsed_scene;
bool source_file_structure_updated,
     scene_file_structure_updated;


private class MapParser {
  string filename, funcname, declaration, definition, data;
  SysTime last_modification;

  this ( string filename_ ) {
    assert(std.file.exists(filename_), "FILE " ~ filename_ ~ " NOT FOUND");
    filename = filename_;
    funcname = To_Funcname(filename);
    Reprocess_File(true);
  }

  bool Reprocess_File ( bool force = false ) {
    if ( !std.file.exists(filename) ) return true;
    auto tlast_modification = RModification_Time(filename);
    if ( !force && tlast_modification == last_modification )
      return false;
    last_modification = tlast_modification;
    data = std.file.read(filename).to!string;
    definition = Parse_Definition(data);

    auto old_declaration = declaration;
    declaration = Parse_Declaration(data, funcname);

    source_file_structure_updated = true;
    return true;
  }
}

private class SceneParser {
  string filename, data, var_match;
  SysTime last_modification;

  this ( string filename_ ) {
    static import std.file;
    assert(std.file.exists(filename_), "FILE " ~ filename_ ~ " NOT FOUND");
    filename = filename_;
    data = std.file.read(filename).to!string;
    Reprocess_File(true);
  }

  bool Reprocess_File ( bool force = false ) {
    if ( !std.file.exists(filename) ) return true;
    auto tlast_modification = RModification_Time(filename);
    if ( !force && tlast_modification == last_modification )
      return false;
    last_modification = tlast_modification;
    data = std.file.read(filename).to!string;

    if ( force || var_match != RScene_Map_Match("MATERIALS") )
      scene_file_structure_updated = true;

    return true;
  }

  string Mixin_Debug_Vars() {
    string debug_vars;
    foreach ( data; KI.model_info.params ) {
      debug_vars ~= data.str_type ~ " " ~ data.label ~ " = " ~
                    data.To_String ~ ";\n";
    }
    return debug_vars;
  }
}

auto RScene_Map_Match ( string label ) in {
  assert(parsed_scene !is null);
} body {
  import std.regex;
  try {
    auto r = parsed_scene.data.matchFirst(label~r"START([\s\S]*)"~label~"END");
    r.popFront();
    return r.front();
  } catch ( Exception e ) {
    writeln("Label ", label, " was not found in scene");
    return "";
  }
}

private auto RModification_Time ( string filename ) in {
  assert(std.file.exists(filename), "FILE " ~ filename ~ " NOT FOUND");
} body {
  SysTime access, modif;
  std.file.getTimes(filename, access, modif);
  return modif;
}

private void Create_Map_Parser ( string filename ) {
  filename = filename.strip();
  if ( filename !in parsed_files )
    parsed_files[filename] = null;
  parsed_files[filename] = new MapParser(filename);
}

/** public access to map parser */
public auto Parse_File ( string filename ) {
  return new MapParser(filename);
}

private string To_Funcname ( string filename ) {
  string funcname = filename[filename.lastIndexOf('/')+1 .. $];
  funcname        = funcname[0 .. funcname.indexOf('.')];
  return funcname;
}
immutable static string
  Extract_params  = // 1 = type, 2 = name (have to insert , after)
                    r"\s*(\w+)\s*(\w+)\s*[,\)]";

private void Recompile_Files ( string[] files ) {
  foreach ( fil; files ) {
    auto location = fil in parsed_files;
    if ( !location ) { Create_Map_Parser(fil); }
    else
      location.Reprocess_File();
  }
}

/// returns the regex match of the declaration "ret fn(..)"
private string Parse_Declaration ( inout string data, inout string funcname ){
  import functional, std.regex;
  string Extract_declaration = r"(\w+)\s*" ~ funcname ~
                               r"\s*\(\s*([_a-zA-Z0-9,\s\*]+)\)";
  auto dec = matchFirst(data, Extract_declaration);
  if ( dec.empty ) {
    writeln("Could not parse declaration from " ~ data);
    writeln("REGEX: ", Extract_declaration);
    return "INVALID";
  }
  return dec.hit ~ ";\n";
}

private string Parse_Definition ( inout string data ) {
  import std.regex, std.string;
  return data;
}

private string Parse_Map ( inout string data, inout string funcname ) {
  import functional;
  string output = "MapUnionG(a, &res, Create_Map_Info(%s(origin + %s), 0,
                             (float3)(0.2f)));";
  string args = KI.model_info.params.map!(n => n.To_String() ~ ",")
                             .joiner.array[0..$-1].to!string;
  return output.format(funcname, args);
}

public auto Create_Model_Info ( MapParser model ) {
  import functional, std.regex;
  KI.ModelInfo model_info;
  model_info.name = model.funcname;
  auto params = matchAll(model.declaration, Extract_params);
  while ( !params.empty ) {
    model_info.params ~= KI.ModelInfo.Data(params.front[1], params.front[2]);
    params.popFront();
  }
  return model_info;


}

private string[2] RParsed_File_Info ( ) {
  string[2] info;
  foreach ( mapparser; parsed_files ) {
    info[0] ~= mapparser.declaration;
    info[1] ~= mapparser.definition  ~ "\n";
  }
  return info;
}

string Ext (T...) ( string str, T extensions ) {
  foreach ( ext; extensions ) {
    if ( str.length <= ext.length ) continue;
    if (  str[$-ext.length .. $] == ext ) return ext;
  }
  return "";
}

bool Valid_File_To_Parse(T)(T file) {
  import std.regex;
  if ( !file.isFile ) return false;
  if ( matchFirst(cast(string)(file)[0..2], r".*\.#")) return false;
  if ( Ext(file, "dtq") != "" ) return false;
  return true;
}

void Reparse_Directory ( string directory ) {
  import functional, std.regex;
  std.file.dirEntries(directory, std.file.SpanMode.breadth)
    .filter!(n => n.Valid_File_To_Parse)
    .map!(n => n.name).array.Recompile_Files();
}

bool Recheck_New_Files ( string directory ) {
  import functional, std.regex;
  return std.file.dirEntries(directory, std.file.SpanMode.breadth)
    .filter!(n => n.Valid_File_To_Parse)
    .filter!(n => n.name !in parsed_files).array.length > 0;
}

bool Reparse_Files ( ) {
  /// -- reparse model/functions --
  static int counter = 0;
  if ( ++counter > 250 ) {
    counter = 0;
    if ( "projects/".Recheck_New_Files() ) parsed_files.clear();
  }
  import functional;
  if ( parsed_files.length == 0 ) {
    "projects/".Reparse_Directory;
    return true;
  }
  bool status = false;
  parsed_files.each!(n => status |= n.Reprocess_File);
  /// -- reparse map ---
  if ( KI.RFile_Type == KI.FileType.DTQ ) {
    // if one does exist we check filename in case scene was swapped
    if ( parsed_scene is null || parsed_scene.filename != KI.RFilename ) {
      parsed_scene = new SceneParser(KI.RFilename);
    }
    status |= parsed_scene.Reprocess_File;
  }
  return status;
}

private string Parse_DTOADQ_Kernel ( ) {
   import globals, opencl_kernel, std.string : replace;
   // -- grab default kernel
   string kernel = DTOADQ_kernel;
   auto filename = KI.RFilename;
   if ( KI.Should_Reparse_Files(false) ) {
     parsed_files.clear();
   }

   Reparse_Files();

  if ( std.file.exists(filename) ) {
    switch ( KI.RFile_Type() ) {
      default: assert(0);
      case KI.FileType.CL: //  // MODEL VIEW
      Recompile_Files([filename]);
      auto mapper = parsed_files[filename];
      if ( source_file_structure_updated ) {
        KI.model_info = Create_Model_Info(mapper);
      }
      source_file_structure_updated = false;
      // mixin map insertion
      kernel = kernel.replace("//%MAPINSERT",
                              Parse_Map(mapper.data, mapper.funcname))
                     .replace("//%SCENEINSERT",
                        "void Update_Camera(Camera* camera, float time){}");
      break;
      case KI.FileType.TXT: assert(0);
      case KI.FileType.DTQ:
        scene_file_structure_updated = false;
        string camera     = RScene_Map_Match("CAMERA"),
               update_map = RScene_Map_Match("UPDATEMAP"),
               materials  = RScene_Map_Match("MATERIALS");
        kernel = kernel
          .replace("//%MAPINSERT",
            `Update_Map(a, origin, &res, time, textures, (float3)(0.0f));`)
          .replace("//%SCENEINSERT",
              camera ~ "\n" ~ update_map);
        // materials
        import functional;
        DIMG.Create_Images(materials.split("\n").filter!(n => n != "")
                                    .map!(n => n.split(" ").array[1]).array);
      break;
    }
  }
  // mixin map function definition/declarations
  string[2] parse_info = RParsed_File_Info();
  kernel = kernel
    .replace("//%MAPFUNCDECLARATIONS", parse_info[0])
    .replace("//%MAPFUNCDEFINITIONS" , parse_info[1]);

  //%SCENEINSERT
  alias KV = KI.KernelInfo.Var;
  kernel = kernel
    .replace("//%MARCH_DIST", KI.RVar(KV.March_Dist).to!string)
    .replace("//%MARCH_REPS", KI.RVar(KV.March_Reps).to!string)
    .replace("//%MARCH_ACC",  KI.RVar(KV.March_Acc ).to!string)
    .replace("//%KERNELTYPE", KI.RKernel_Type_String);
  if ( KI.RKernel_Info.flags[KI.KernelInfo.Flag.Show_Normals] ) {
    kernel = kernel.replace("//#define SHOW_NORMALS", "#define SHOW_NORMALS");
  }

  return kernel;
}

string Parse_Texture_Kernel() {
  import globals, opencl_kernel, std.string : replace;
  string kernel = Texture_kernel;
  auto filename = KI.RFilename();
  if ( KI.Should_Reparse_Files(false) ) {
    parsed_files.clear();
  }

  Reparse_Files();

  string[2] parse_info = RParsed_File_Info();
  kernel = kernel
    .replace("//%MAPFUNCDECLARATIONS", parse_info[0])
    .replace("//%MAPFUNCDEFINITIONS" , parse_info[1]);

  if ( std.file.exists(filename) ) {
    Recompile_Files([filename]);
    kernel = kernel
      .replace("//%MAPINSERT",
        "return " ~ filename.To_Funcname ~ "(origin);");
  }
  return kernel;
}

string Parse_Kernel ( ) {
  if ( KI.RFile_Type == KI.FileType.TXT ) {
    return Parse_Texture_Kernel();
  }
  return Parse_DTOADQ_Kernel();
}
