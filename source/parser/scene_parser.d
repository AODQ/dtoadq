module parser.scene_parser;
static import stl, core, parser.file;

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
  import stl : replace;
  string kernel = parser.file.parsed_kernel.Reprocess_Data.data,
         scene  = parser.file.parsed_scene .Reprocess_Data.data;

  string camera     = RData_Subsection!("CAMERA"    )(scene),
         update_map = RData_Subsection!("UPDATEMAP" )(scene),
         materials  = RData_Subsection!("MATERIALS" )(scene),
         textures   = RData_Subsection!("TEXTURES"  )(scene),
         emitter    = RData_Subsection!("EMITTER"   )(scene);

  { // textures (extract `0 filename.txt`)
    import functional;
    static import glfw;
    core.Create_Images(
      textures.split("\n").filter!(n => n != "")
              .map!(n => n.split(" ").array[1]).array);
  }

  { // materials
    static import shared_info = core.shared_info;
    core.Set_Material(stl.json.parseJSON(materials));
    kernel = kernel
      .replace("//%MAT_LENGTH",
          stl.to!string(shared_info.material.length));
  }

  // map and scene insert
  kernel = kernel
    .replace("//%MAPINSERT",
         `Update_Map(a, origin, &res, si->time, Tx, si->debug_values);`)
    .replace("//%SCENEINSERT",
         emitter ~ "\n" ~ camera ~ "\n" ~ update_map);

  // functions insert
  auto function_data = parser.file.RParsed_Function_Data();
  kernel = kernel
    .replace("//%MAPFUNCDECLARATIONS", function_data.declarations)
    .replace("//%MAPFUNCDEFINITIONS",  function_data.definitions);

  // constants insert
  alias KV = core.info.KernelVar;
  return kernel
    .replace("//%MARCH_DIST", stl.to!string(core.RKernel_Var(KV.March_Dist)))
    .replace("//%MARCH_REPS", stl.to!string(core.RKernel_Var(KV.March_Reps)))
    .replace("//%MARCH_ACC",  stl.to!string(core.RKernel_Var(KV.March_Acc )));
}
