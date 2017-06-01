module parser;
import globals;
static import KI = kernelinfo;
static import std.file;
import std.datetime : SysTime;

MapParser[string] parsed_files;
bool source_file_has_updated;

private class MapParser {
  string filename, funcname, declaration, definition, data;
  SysTime last_modification;

  this ( string filename_ ) {
    assert(std.file.exists(filename_), "FILE " ~ filename_ ~ " NOT FOUND");
    filename = filename_;
    funcname = To_Funcname(filename);
    Reprocess_File(true);
  }

  void Reprocess_File ( bool force = false ) {
    if ( !force && RModification_Time(filename) == last_modification ) return;
    last_modification = RModification_Time(filename);
    data = std.file.read(filename).to!string;
    declaration = Parse_Declaration(data, funcname);
    definition = Parse_Definition(data);
    Parse_Requirements(data);

    source_file_has_updated = true;
  }
}

private auto RModification_Time ( string filename ) {
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

private string To_Funcname ( string filename ) {
  string funcname = filename[filename.lastIndexOf('/')+1 .. $];
  funcname        = funcname[0 .. funcname.indexOf('.')];
  return funcname;
}
immutable static string
  Extract_require = r"!require\s+([\/\w\.\s]+);",
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

/// Creates the map parsers if they aren't required yet or updates them if
/// they have been modified
private void Parse_Requirements ( inout string data ) {
  import functional, std.regex;
  auto req = data.matchAll(Extract_require);
  if ( !req.empty ) {
    string[] files = req.front[1].split(" ").filter!(n => n != "").array;
    Recompile_Files(files);
  }
}

/// returns the regex match of the declaration "ret fn(..)"
private string Parse_Declaration ( inout string data, inout string funcname ){
  import functional, std.regex;
  string Extract_declaration = r"(\w+)\s*" ~ funcname ~
                               r"\s*\(\s*([a-zA-Z0-9,\s]+)\)";
  auto dec = matchFirst(data, Extract_declaration);
  if ( dec.empty )
    assert(false, "Could not parse declaration from " ~ data);
  return dec.hit ~ ";\n";
}

private string Parse_Definition ( inout string data ) {
  import std.regex, std.string;
  // -- remove !require .. ; from definition --
  auto mall = data.matchAll(Extract_require);
  if ( !mall.empty )
    return data.replace(mall.hit, "");
  return data;
}

private string Parse_Map ( inout string data, inout string funcname ) {
  import functional;
  string output = "res = MapUnionG(a, res, (float2)(%s(origin + %s), 0.0f));";
  string args = KI.model_info.params.map!(n => n.To_String() ~ ",")
                             .joiner.array[0..$-1].to!string;
  return output.format(funcname, args);
}

private void Create_Model_Info ( MapParser model ) {
  import functional, std.regex;
  KI.ModelInfo model_info;
  model_info.name = model.funcname;
  auto params = matchAll(model.declaration, Extract_params);
  while ( !params.empty ) {
    model_info.params ~= KI.ModelInfo.Data(params.front[1], params.front[2]);
    params.popFront();
  }
  KI.model_info = model_info;
}

private string[2] RParsed_File_Info ( ) {
  string[2] info;
  foreach ( mapparser; parsed_files ) {
    info[0] ~= mapparser.declaration ~ "\n";
    info[1] ~= mapparser.definition  ~ "\n";
  }
  return info;
}

string Parse_Kernel ( ) {
   import globals, opencl_kernel, std.string : replace;
   // -- grab default kernel
   string kernel = DTOADQ_kernel;
   auto filename = KI.RFilename;
   if ( KI.Should_Reparse_Files(false) ) {
     parsed_files.clear();
   }

  if ( std.file.exists(filename) ) {
    Recompile_Files([filename]);
    auto mapper = parsed_files[filename];
    if ( source_file_has_updated )
      Create_Model_Info(mapper);
    source_file_has_updated = false;
    // mixin map insertion
    kernel = kernel.replace("//%MAPINSERT",
                            Parse_Map(mapper.data, mapper.funcname));
    // mixin map function definition/declarations
    string[2] parse_info = RParsed_File_Info();
    kernel = kernel
      .replace("//%MAPFUNCDECLARATIONS", parse_info[0])
      .replace("//%MAPFUNCDEFINITIONS" , parse_info[1]);
  } else {
    writeln("Could not find file ", filename);
    kernel = kernel.replace("//%MAPINSERT",
      `res = MapUnionG(a, res, (float2)(length(origin)-1.25f, 0.0f));`);
  }

  alias KV = KI.KernelInfo.Var;
  kernel = kernel
    .replace("//%MARCH_DIST", KI.RVar(KV.March_Dist).to!string)
    .replace("//%MARCH_REPS", KI.RVar(KV.March_Reps).to!string)
    .replace("//%KERNELTYPE", KI.RKernel_Type_String);
  if ( KI.RKernel_Info.flags[KI.KernelInfo.Flag.Show_Normals] ) {
    kernel = kernel.replace("//#define SHOW_NORMALS", "#define SHOW_NORMALS");
  }

  return kernel;
 }
