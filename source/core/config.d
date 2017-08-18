module core.config;
static import conf = configuration, stl;
static import core.info, core.kernel;

auto Configure ( ) {
  // initialize variables before configuring
  core.kernel.Initialize_Variables();
  core.info.Set_Kernel_Type(core.info.KernelType.DTQ);
  // now config
  core.info.Set_Map_Function(RDefault_Project());
  alias KV = core.info.KernelVar;
  foreach ( tup; RZipped_KVars(KV.March_Dist, KV.March_Reps, KV.March_Acc) ) {
    core.info.Set_Kernel_Var(tup[0], RKernel_Var(tup[1]));
  }
}

private auto RZipped_KVars(KV...)(KV vars) {
  import functional : zip;
  core.info.KernelVar[] kv;
  string[] str;
  foreach ( v; vars ) { kv ~= v; str ~= stl.to!string(v); }
  return zip(kv, str);
}

string RDefault_Project ( ) {
  static immutable auto Default_project_name = "globals/defaultscene.dtq";
  if ( !conf.Exists ) return Default_project_name;
  stl.json.JSONValue load = conf.JSON_Value;

  try { return load["default-project"].str; }
  catch ( Exception e ) { return Default_project_name; }
}

private int RDefault_Kernel_Var ( string var ) {
  switch ( var ) {
    default: assert(false);
    case "March_Dist": return 64;
    case "March_Reps": return 16;
    case "March_Acc":  return 10;
  }
}

int RKernel_Var ( string var ) {
  immutable auto Default_value = RDefault_Kernel_Var(var);
  if ( !conf.Exists ) return Default_value;
  stl.json.JSONValue load = conf.JSON_Value;

  try { return stl.to!int(load["default-kernel-var-"~var].str); }
  catch ( Exception e ) { return Default_value; }
}
