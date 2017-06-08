module gui.fuzzyinput;
import globals, gui.gui;
import derelict.imgui.imgui;
import std.typecons : Tuple, tuple;

private string curr_filename;
private size_t cursor;
private Tuple!(size_t, string)[] filenames;

void Initialize_ImGui ( ref bool open ) {
  igSetNextWindowPos(ImVec2(config.RWindow_Size[0]/2 - 200, 200));
  igSetNextWindowSize(ImVec2(400, 400));
  igSetNextWindowFocus();
  if ( !igBegin("Fuzzy File Search", &open,
                ImGuiWindowFlags_NoMove|
                ImGuiWindowFlags_NoScrollbar|
                ImGuiWindowFlags_NoCollapse) ) {
    writeln("WTF?");
    curr_filename = "";
    cursor = 0;
    igEnd();
  }
}

void Deinitialize_ImGui ( ) {
  igEnd();
}

auto Set_Fuzzy_Files ( ) {
  import gui.files;
  auto files = RFiles(RProject_Dir) ~ RFiles(RGlobals_Dir);
  // iterate through all files and record score
  foreach ( fil; files ) {
    string filename = Util.Truncate_DirExt(fil.name);
    import std.algorithm.comparison : levenshteinDistance;
    auto dist = levenshteinDistance!((a, b) => toLower(a) == toLower(b))
                                    (filename, curr_filename);
    if ( dist < 16 ) {
      filenames ~= tuple(dist, fil.name);
    }
  }
  // sort based on score
  import functional;
  filenames.sort!("a[0] < b[0]")();
}

void Reset ( ) {
  filenames = [];
  cursor = 0;
}

string Update ( ref bool open ) {
  if ( !open ) return "";
  Initialize_ImGui(open);
  // -- get input --
  if ( gdInputText("File Name", curr_filename) ) {
    Reset();
    Set_Fuzzy_Files();
  }

  // -- update output --
  foreach ( it, immutable fil; filenames ) {
    if ( it == cursor )
      gdTextColour(fil[1], 0.8f, 0.3f, 0.3f);
    else
      gdText(fil[1]);
  }

  // -- check cursor change --
  bool up   = input.RKey_Input(74, false),
       down = input.RKey_Input(75, false);
  if ( input.RKey_Input(341) || input.RKey_Input(345) )
    cursor += up - down;

  if ( cursor >= filenames.length ) {
    cursor = 0;
    if ( down ) cursor = filenames.length-1;
  }

  Deinitialize_ImGui();

  if ( input.RKey_Input(257) && filenames.length > 0 ) {
    open = false;
    curr_filename = "";
    auto fname = filenames[cursor][1];
    Reset();
    return fname;
  }
  return "";
}
