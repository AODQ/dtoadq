module gui.files;
import derelict.imgui.imgui;
static import stl, core.info;

private string project_dir = "";

void Set_Project_Directory ( string dir ) {
  project_dir = dir;
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
  struct Dir {
    string dir_name;
    Dir[string] dirs;
    string[] files;

    this ( string dir_name_ ) { dir_name = dir_name_; }

    void Add ( bool is_file, string fname, string name ) {
      auto ind = stl.indexOf(name, "/");
      if ( ind == -1 ) {
        if ( is_file ) files ~= fname;
        else dirs[name] = Dir(name);
      } else {
        string dir = name[0..ind];
        dirs[dir].Add(is_file, fname, name[ind+1..$]);
      }
    }

    void Render ( ) {
      foreach ( dir; dirs ) {
        // put igTreeNode here so the parent directory isn't displayed
        if ( igTreeNode(stl.toStringz(dir.dir_name)) ) {
          dir.Render();
          igTreePop();
        }
      }
      foreach ( fil; files ) {
        import std.string;
        bool open = core.info.RFilename == fil;
        if ( igCheckbox(fil[fil.lastIndexOf('/')+1..$].toStringz, &open) ) {
          core.info.Set_Any_Function(fil);
        }
      }
    }
  }

  static import Fil = std.file;
  import functional;
  string dir = dir_prestr;
  auto base_dir = Dir(dir);
  Fil.dirEntries(dir, Fil.SpanMode.breadth)
     .filter!(n => stl.Not_Swap_File(n.name))
     .each!(n => base_dir.Add(n.isFile, n.name, n.name.replace(dir~"/", "")));
  base_dir.Render();
}

void Update ( ref bool open ) {
  igBegin("File Browser", &open);
  if ( open ) {
    import std.string;
    if ( igCollapsingHeader("globals") )
      Update_Directory(RGlobals_Dir);
    if ( igCollapsingHeader(project_dir.toStringz) )
      Update_Directory(project_dir);
  }
  igEnd();
}
