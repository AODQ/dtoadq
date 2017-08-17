module glfw.config;
static import conf = configuration, stl;

auto Configure ( ) {
}

int[] RWindow_Size ( ) {
  static immutable auto Default_window_size = [640, 480];
  if ( !conf.Exists ) return Default_window_size.dup;
  stl.json.JSONValue load = conf.JSON_Value;

  try { return conf.To_Int_Array(load["window-size"].str); }
  catch ( Exception e ) { return Default_window_size.dup; }
}

int[] RWindow_Pos ( ) {
  static immutable auto Default_window_pos = [-1, -1];
  if ( !conf.Exists ) return Default_window_pos.dup;
  stl.json.JSONValue load = conf.JSON_Value;

  try { return load.To_Int_Array(load["window-pos"].str); }
  catch ( Exception e ) { return Default_window_pos.dup; }
}
