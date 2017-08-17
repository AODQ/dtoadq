module kernel.config;
static import conf = configuration, stl;

auto Configure ( ) {
  import kernel.info;
  kernel.info.Set_Map_Function(RDefault_Project());
}

string RDefault_Project ( ) {
  static immutable auto Default_project_name = "globals/defaultscene.dtq";
  if ( !conf.Exists ) return Default_project_name;
  stl.json.JSONValue load = conf.JSON_Value;

  try { return load["default-project"].str; }
  catch ( Exception e ) { return Default_project_name; }
}
