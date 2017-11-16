
module gui.config;
static import conf = configuration, stl, core, gui.files;

auto Configure ( ) {
  auto proj_dir = core.config.RDefault_Project();
  import stl : lastIndexOf;
  gui.files.Set_Project_Directory(proj_dir[0 .. proj_dir.lastIndexOf("/")]);
  import derelict.imgui.imgui;
  auto g = igGetStyle();
}
