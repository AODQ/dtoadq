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
  dir_prestr ~= "/";
  // I may have overengineered this function . . .
  import functional, std.string;
  static import KI = kernelinfo;
  alias Dir = KI.ProceduralType;
  // -- dir info --
  static string[] directories =
    iota(0, Dir.max+1).map!(n =>
      (cast(Dir)n).to!(string).toLower ~ "s"
  ).array;
  auto dirs = iota(0, Dir.max+1).map!(n => dir_prestr ~ directories[n]);

  // display all files in a directory
  auto Display_Dir(int dir_type) {
    import gui.gui : gdInputText, gdButton;
    static string def_filename = "";
    // gdInputText("New File", def_filename);
    // if ( def_filename != "" && gdButton("Create ", cast(Dir)dir_type) ) {
    //   static import NP = gui.node_parser;
    //   string filename = dir_prestr ~ directories[dir_type] ~ "/" ~
    //                     def_filename ~ ".json";
    //   writeln("Creating filename: ", filename);
    //   import gui.node_parser;
    //   Create_Default_Graph(filename);
    //   Load_Graph(filename);
    // }
    static int open_file = 0, open_dir = 0;
    auto files = dirs[dir_type].RFiles;
    files.sort!("a.name < b.name");
    foreach ( it, fil; files ) {
      // so open_fil index of each entry isn't highlighted
      int t_open_file = open_file;
      if ( dir_type != open_dir ) t_open_file = -1;
      auto name = Util.Truncate_Directory(fil.name);
      if ( igRadioButton(name.toStringz, &t_open_file, cast(int)it) ) {
        open_dir = dir_type;
        auto filetype = KI.String_To_FileType(Ext(name, ".cl", ".dtq", ".txt")
                           [1..$]);
        KI.Set_Map_Function(cast(Dir)dir_type, filetype, fil.name);
      }
      open_file = t_open_file;
    }
  }

  import std.typecons : Tuple;

  alias DirTuple = Tuple!(string, Dir);

  // -- render --
  zip(directories, iota(0, Dir.max+1).map!(n => cast(Dir)n))
  .each!(delegate void(DirTuple dir) {
    if ( igTreeNode((dir_prestr ~ dir[0]).toStringz) ) {
      Display_Dir(dir[1]);
      igTreePop();
    }
  });
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
