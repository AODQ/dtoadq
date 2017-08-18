module parser.texture_parser;
static import stl, parser.file;

private string Filename_To_Function ( string filename ) {
  import std.string;
  string funcname = filename[filename.lastIndexOf('/')+1 .. $];
  funcname        = funcname[0 .. funcname.indexOf('.')];
  return funcname;
}

string Parse() {
  import stl : replace;
  string kernel  = parser.file.parsed_kernel.Reprocess_Data.data;

  auto function_data = parser.file.RParsed_Function_Data();

  return kernel
    .replace("//%MAPFUNCDECLARATIONS", function_data.declarations)
    .replace("//%MAPFUNCDEFINITIONS",  function_data.definitions)
    .replace("//%MAPINSERT",
      "return " ~ Filename_To_Function(parser.file.parsed_scene.name)
                ~ "(origin);");
}
