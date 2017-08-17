module parser.scene_parser;
static import stl, file = parser.file;

string RData_Subsection (string subsection)(ref string data) {
  auto rgx = stl.regex.ctRegex!(subsection~r"START([\s\S]*)"~subsection~r"END");
  try {
    auto r = stl.regex.matchFirst(data, rgx);
    r.popFront();
    return r.front();
  } catch ( Exception e ) {
    stl.writeln("Label ", subsection, " was not found in scene");
    return "";
  }
}

string Parse ( ) {
  static import kernel.info;
  import stl : replace;
  string data = file.parsed_kernel.Reprocess_Data.data;

  string camera     = RData_Subsection!("CAMERA"    )(data),
         update_map = RData_Subsection!("UPDATEMAP" )(data),
         materials  = RData_Subsection!("MATERIALS" )(data),
         textures   = RData_Subsection!("TEXTURES"  )(data),
         emitter    = RData_Subsection!("EMITTER"   )(data);

  { // textures (extract `0 filename.txt`)
    import functional;
    static import glfw;
    glfw.image.Create_Images(
      textures.split("\n").filter!(n => n != "")
              .map!(n => n.split(" ").array[1]).array);
  }

  { // materials
    static import dtoadq;
    dtoadq.Set_Material(materials);
  }

  // map and scene insert
  data = data
    .replace("//%MAPINSERT",
         `Update_Map(a, origin, &res, si->time, Tx, si->debug_values);`)
    .replace("//%SCENEINSERT",
         emitter ~ "\n" ~ camera ~ "\n" ~ update_map);

  // functions insert
  auto function_data = file.RParsed_Function_Data();
  data = data
    .replace("//%MAPFUNCDECLARATIONS", function_data.declarations)
    .replace("//%MAPFUNCDEFINITIONS",  function_data.definitions);

  // constants insert
  alias KV = kernel.info.KernelVar;
  return data
    .replace("//%MARCH_DIST", KI.RVar(KV.March_Dist).to!string)
    .replace("//%MARCH_REPS", KI.RVar(KV.March_Reps).to!string)
    .replace("//%MARCH_ACC",  KI.RVar(KV.March_Acc ).to!string);
}
