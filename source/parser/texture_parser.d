module parser.texture_parser;
static import stl, file = parser.file;

private string Filename_To_Function ( string filename ) {
  string funcname = filename[filename.lastIndexOf('/')+1 .. $];
  funcname        = funcname[0 .. funcname.indexOf('.')];
  return funcname;
}

string Parse() {
  import stl : replace;
  string data = file.parsed_kernel.Reprocess_Data.data;

  auto function_data = file.RParsed_Function_Data();

  return data
    .replace("//%MAPFUNCDECLARATIONS", function_data.declarations)
    .replace("//%MAPFUNCDEFINITIONS",  function_data.definitions)
    .replace("//%MAPINSERT",
      "return " ~ Filename_To_Function(file.parsed_kernel.name) ~ "(origin)");
}
