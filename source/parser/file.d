module parser.file ;
static import std.file;
import std.datetime : SysTime;
import stl;

/*
  needs to provide functionalities . . .

    hot reloading
    insert / replace text
    raytracer , dtoadq , texture files
    raytracer , dtoadq , texture kernels
*/

/// File Type of file (for map, raytracer/dtoadq is agnostic)
enum FileType { Raytracer, DTOADQ, Texture, Function };
/// File structure to maintain file type, definition/declaration and last
/// modification time for updates
struct File {
  string name;
  SysTime last_modification;
  FileType file_type;

  ///
  this ( string filename ) {
    name = stl.strip(filename);
    Enforce_Exist();
    Reprocess(true);
  }

  this ( this ) {
    this.name = name;
    this.last_modification = last_modification; this.file_type = file_type;
  }

  /// Returns the last modification time (according to linux, not file struct)
  auto RModification_Time() {
    Enforce_Exist();
    SysTime access, modif;
    stl.file.getTimes(name, access, modif);
    return modif;
  }

  /// Returns true on file update (including deletion)
  bool Reprocess ( bool force_reprocess = false ) {
    if ( !Exists() ) return true;
    auto curr_modification = RModification_Time();
    if ( !force_reprocess && curr_modification == last_modification )
      return false;
    last_modification = curr_modification;
    Reprocess_Data();
    return true;
  }

  /***
    returns struct { unformatted data, extracted formatted declarations }
    declarations == "" on kernel/map filetype (raytracer, dtoadq, texture)
  ***/
  auto Reprocess_Data ( ) {
    struct Data { string data, declarations; }
    Data data;
    data.data = stl.file.read(name).to!string;
    { // declarations reprocess
      import functional, std.regex;
      auto rgx = ctRegex!r"(float\d*|int)\s+\w+\s*\([^{]+";
      data.declarations = data.data.matchAll(rgx).map!(m => m.hit ~ ";\n")
                              .joiner.to!string;
    }
    return data;
  }

  /// Checks if file still exists
  bool Exists ( ) { return stl.file.exists(name); }
  /// Enforces existence of file
  void Enforce_Exist ( ) {
    stl.enforce(Exists(), "File " ~ name ~ " not found");
  }

  ///
  mixin stl.RValRef;
}

private template File_Template_Functions ( string varname ) {
  File file_;
  alias file_ this;
  this ( File f ) {
    file_ = f.Ref;
    file_type = RFile_Type(name);
  }
  this ( this ) {
    file_ = Ref();
    file_type = file_type;
  }
  static auto New ( string filename ) {
    return typeof(this)(File(filename));
  }
  /// Creates and stores from a given filename (for *File structs)
  static void Create_File ( string filename ) {
    import std.string : format;
    mixin("%s = typeof(this).New(filename);".format(varname));
  }
}

/// Kernel for dtoadq, raytracer and texture
struct KernelFile {
  static FileType RFile_Type(string name) {
    switch ( name ) {
      import parser;
      default: assert(false);
      case DTOADQ_filename:    return FileType.DTOADQ;
      case Raytrace_filename:  return FileType.Raytracer;
      case Raycast_filename:   return FileType.Raytracer;
      case Texture_filename:   return FileType.Texture;
    }
  }
  mixin File_Template_Functions!("parsed_kernel");
}

/// map (DTOADQ) and texture scenes
struct SceneFile {
  static FileType RFile_Type(string name) {
    switch ( stl.RMatching_File_Extensions(name, "dtq", "txt") ) {
      default: assert(false);
      case "dtq": return FileType.DTOADQ;
      case "txt": return FileType.Texture;
    }
  }
  mixin File_Template_Functions!("parsed_scene");
}

/// generic functions
struct FunctionFile {
  static FileType RFile_Type(string name) { return FileType.Function; }
  mixin File_Template_Functions!("parsed_functions[filename]");
}

///
KernelFile           parsed_kernel;
///
SceneFile            parsed_scene;
/// Unlike scene and kernel, multiple function files may exist
FunctionFile[string] parsed_functions;

/// Recompiles files from given array of filenames, if there are any
/// modifications to the file and the file has been stored
void Recompile_Files ( string[] filenames ) {
  foreach ( file; filenames ) {
    if ( file == parsed_kernel.name ) { parsed_kernel.Reprocess(); continue; }
    if ( file == parsed_scene .name ) { parsed_scene.Reprocess(); continue; }
    auto scene_location = file in parsed_functions;
    if ( scene_location ) { scene_location.Reprocess(); continue; }
  }
}

/// returns the data of functions as struct{declarations, definitions}
auto RParsed_Function_Data ( ) {
  import functional;
  struct SceneData { string declarations, definitions; }
  if ( parsed_functions.length == 0 ) return SceneData("", "");
  // [{dec, data}, ..] [[dec, data], ..] -> [[dec, ..], [data, ..]] ->
  // [dec ~ .., data ~ ..]
  auto value = parsed_functions.values.map!((f) {
    auto data = f.Reprocess_Data();
    return [data.declarations, data.data];
  }).array.transposed.map!(both => both.joiner.array.to!string).array;
  return SceneData(value[0], value[1]);
}
