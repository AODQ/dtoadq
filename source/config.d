module config;
import std.json, globals;
static import std.file;

private string filename = "dtoadq-config.json";

bool Check_File_Exists ( ) {
  return std.file.exists(filename);
}

auto RJSONValue ( ) { return parseJSON(std.file.read(filename).to!string); }
auto To_Int_Array ( string str ) { // to!(int[]) might work but redundant syntax
  import functional;
  return str.split.map!(to!int).array;
}

int[] RWindow_Size ( ) {
  if ( !Check_File_Exists ) return [640, 480];
  JSONValue load = RJSONValue;

  try {
    return load["window-size"].str.To_Int_Array;
  } catch ( Exception e ) {
    writeln("ERROR: Could not find/parse window-size, expecting format: 'x y'");
    return [640, 480];
  }
}

int[] RWindow_Pos ( ) {
  if ( !Check_File_Exists ) return [-1, -1];
  auto load = RJSONValue;
  string str;

  try {
    str = load["window-pos"].str;
  } catch ( Exception e ) {
    return [-1, -1]; // be silent, window pos isn't necessary
  }

  try { return str.To_Int_Array; } catch ( Exception e ) {
    writeln("ERROR: Could not parse window-pos, expecting format: 'x y'");
    return [-1, -1];
  }
}
