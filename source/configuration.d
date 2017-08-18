module configuration;
static import stl;
static import std.file;

static immutable string Config_filename = "dtoadq.cfg";

/// Checks if config file exists
bool Exists ( ) { return stl.file.exists(Config_filename); }

/// Returns the parsed json value of config file
auto JSON_Value ( ) {
  return stl.json.parseJSON(stl.to!string(stl.file.read(Config_filename)));
}

/// Converts a string of format "1920 1080" to int[]
int[] To_Int_Array ( string str ) { // not to!(int[]) b wary of json syntax
  import functional;
  return str.split.map!(to!int).array;
}

/** Configures all module configs */
auto Configure ( ) {
  static import glfw.config, core.config, emitter.config, gui.config;
  glfw   .config.Configure();
  core   .config.Configure();
  emitter.config.Configure();
  gui    .config.Configure();
}
