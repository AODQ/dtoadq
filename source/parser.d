module parser;
import globals;
static import KI = kernelinfo;
static import DIMG = dtoadqimage;
static import std.file;
import std.datetime : SysTime;

MapParser[string] parsed_files;
SceneParser       parsed_scene;

private class MapParser {
  string filename, declarations, data;
  SysTime last_modification;

  this ( string filename_ ) {
    assert(std.file.exists(filename_), "FILE '" ~ filename_ ~ "' NOT FOUND");
    filename = filename_;
    Reprocess_File(true);
  }

  bool Reprocess_File ( bool force = false ) {
    if ( !std.file.exists(filename) ) return true;
    auto tlast_modification = RModification_Time(filename);
    if ( !force && tlast_modification == last_modification )
      return false;
    last_modification = tlast_modification;
    data = std.file.read(filename).to!string;
    declarations = RDeclarations();
    return true;
  }

  auto RDeclarations () {
    import std.regex;
    auto rgx = r"(float\d*|int)\s+\w+\s*\([^{]+";
    auto matches = data.matchAll(rgx);
    string declarations = "";
    foreach ( m; matches ) {
      declarations ~= m.hit~";\n";
    }
    return declarations;
  }
}

private class SceneParser : MapParser {
  this ( string filename_ ) {
    super(filename_) ;
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
  if ( !std.file.exists(filename) ) {
    writeln("file ", filename, " not found");
    return;
  }
  if ( filename !in parsed_files )
    parsed_files[filename] = null;
  parsed_files[filename] = new MapParser(filename);
}

private void Recompile_Files ( string[] files ) {
  foreach ( fil; files ) {
    auto location = fil in parsed_files;
    if ( !location ) { Create_Map_Parser(fil); }
    else
      location.Reprocess_File();
  }
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
  bool status = false;
  import functional;
  if ( parsed_files.length == 0 ) {
    "projects/".Reparse_Directory;
    status = true;
  }
  parsed_files.each!(n => status |= n.Reprocess_File);
  /// -- reparse map ---
  if ( KI.RProcedural_Type == KI.ProceduralType.Scene ) {
    // if one does exist we check filename in case scene was swapped
    if ( parsed_scene is null || parsed_scene.filename != KI.RFilename ) {
      if ( !std.file.exists(KI.RFilename) ) {
        assert(false, "File '" ~ KI.RFilename ~ "' does not exist");
      }
      parsed_scene = new SceneParser(KI.RFilename);
    }
    status |= parsed_scene.Reprocess_File;
  }
  return status;
}

private string[2] RParsed_File_Info ( ) {
  string[2] info;
  foreach ( mfile; parsed_files ) {
    info[0] ~= mfile.declarations;
    info[1] ~= mfile.data  ~ "\n";
  }
  return info;
}

private string Parse_DTOADQ_Kernel ( ) {
   import globals, dtoadqkernel, std.string : replace;
   // -- grab kernel (either raytrace or MLT)
   string kernel;
   if ( KI.RKernel_Type == KI.KernelType.Raytrace ) {
     import raytracekernel;
     kernel = Raytrace_kernel;
   } else {
     import dtoadqkernel;
     kernel = DTOADQ_kernel;
   }

   auto filename = KI.RFilename;
   if ( KI.Should_Recompile(false) ) {
     parsed_files.clear();
   }

  Reparse_Files();

  if ( std.file.exists(filename) ) {
    string camera       = RScene_Map_Match("CAMERA"),
           update_map   = RScene_Map_Match("UPDATEMAP"),
           materials    = RScene_Map_Match("MATERIALS"),
           post_process = RScene_Map_Match("POSTPROCESS");
    kernel = kernel
      .replace("//%MAPINSERT",
        `Update_Map(a, origin, &res, time, textures, dval);`)
      .replace("//%SCENEINSERT",
          camera ~ "\n" ~ update_map)
      .replace("//%POSTPROCESS", post_process);
    // materials
    import functional;
    DIMG.Create_Images(materials.split("\n").filter!(n => n != "")
                                .map!(n => n.split(" ").array[1]).array);
  }

  // mixin map function definition/declarations
  string[2] parse_info = RParsed_File_Info;
  kernel = kernel
    .replace("//%MAPFUNCDECLARATIONS", parse_info[0])
    .replace("//%MAPFUNCDEFINITIONS" , parse_info[1]);

  alias KV = KI.KernelVar;
  kernel = kernel
    .replace("//%MARCH_DIST", KI.RVar(KV.March_Dist).to!string)
    .replace("//%MARCH_REPS", KI.RVar(KV.March_Reps).to!string)
    .replace("//%MARCH_ACC",  KI.RVar(KV.March_Acc ).to!string);

  return kernel;
}


private string To_Funcname ( string filename ) {
  string funcname = filename[filename.lastIndexOf('/')+1 .. $];
  funcname        = funcname[0 .. funcname.indexOf('.')];
  return funcname;
}

string Parse_Texture_Kernel() {
  import globals, std.string : replace;
  import texturekernel;
  string kernel = Texture_kernel;
  auto filename = KI.RFilename();
  if ( KI.Should_Recompile(false) ) {
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
  if ( KI.RProcedural_Type == KI.ProceduralType.Texture ) {
    return Parse_Texture_Kernel();
  }
  return Parse_DTOADQ_Kernel();
}
