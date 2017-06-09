module gui.files;
import derelict.imgui.imgui, globals;

private string project_dir = "";

void Initialize ( ) {
  import config;
  project_dir = "projects/" ~ RDefault_Project();
}

string Ext (T...) ( string str, T extensions ) {
  foreach ( ext; extensions ) {
    if ( str.length <= ext.length ) continue;
    if (  str[$-ext.length .. $] == ext ) return ext;
  }
  return "";
}

auto RFiles ( string directory ) {
  directory = directory;
  static import std.file;
  import functional, std.regex : matchFirst;
  return std.file.dirEntries(directory, std.file.SpanMode.breadth)
    .filter!(n => n.isFile &&
        !matchFirst(cast(string)(n)[0..2], r".*\.#") && // emacs
        Ext(n.name, ".cl", ".json")) // hand-written or node graph
    .array;
}

auto RProject_Dir ( ) { return project_dir; }
auto RGlobals_Dir ( ) { return "projects/globals"; }

void Update_Directory ( string dir_prestr ) {
  class Dir {
    string dir_name;
    Dir[string] dirs;
    string[] files;

    this ( string dir_name_ ) { dir_name = dir_name_; }

    void Add ( bool is_file, string fname, string name ) {
      auto ind = name.indexOf("/");
      if ( ind == -1 ) {
        if ( is_file ) files ~= fname;
        else dirs[name] = new Dir(name);
      } else {
        string dir = name[0..ind];
        dirs[dir].Add(is_file, fname, name[ind+1..$]);
      }
    }

    void Render ( ) {
      foreach ( dir; dirs ) {
        // put igTreeNode here so the parent directory isn't displayed
        if ( igTreeNode(dir.dir_name.toStringz) ) {
          dir.Render();
          igTreePop();
        }
      }
      foreach ( fil; files ) {
        static import KI = kernelinfo;
        bool open = KI.RFilename == fil;
        if ( igCheckbox(fil[fil.lastIndexOf('/')+1..$].toStringz, &open) ) {
          KI.Set_Map_Function(fil);
        }
      }
    }
  }

  static import Fil = std.file;
  import functional;
  string dir = dir_prestr;
  auto base_dir = new Dir(dir);
  Fil.dirEntries(dir, Fil.SpanMode.breadth)
     .each!(n => base_dir.Add(n.isFile, n.name, n.name.replace(dir~"/", "")));
  base_dir.Render();
}

void Update ( ref bool open ) {
  igBegin("File Browser", &open);
  if ( open ) {
    if ( igCollapsingHeader("globals") )
      Update_Directory(RGlobals_Dir);
    if ( igCollapsingHeader(project_dir.toStringz) )
      Update_Directory(project_dir);
  }
  igEnd();
}
