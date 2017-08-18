module parser.checker;
static import stl, file = parser.file;

private bool Valid_Function_File (string file) {
  return stl.file.isFile(file) &&
         stl.RMatching_File_Extensions(file, ".dtq") != ".dtq";
}

private string[] RNew_Function_Files ( string directory ) {
  import functional, std.regex;
  return stl.file.dirEntries(directory, stl.file.SpanMode.breadth)
            .map   !(n => n.name)
            .filter!(n => n.Valid_Function_File && n !in file.parsed_functions)
            .array;
}

/// Returns if any files have been modified
bool Recheck_Files ( ) {
  static int recheck_file_counter = 249;
  bool status = false;
  // check new files created
  if ( !(++recheck_file_counter % 250) ) {
    auto files = RNew_Function_Files("projects/globals");
    status = files.length > 0;
    foreach ( f; files ) { file.FunctionFile.Create_File(f); }
  }

  // check files modified
  foreach ( f; file.parsed_functions ) status |= f.Reprocess();
  status |= file.parsed_kernel.Reprocess();
  status |= file.parsed_scene .Reprocess();
  return status;
}
