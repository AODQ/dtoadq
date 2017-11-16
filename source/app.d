static import core, stl, glfw;

bool running = true;

bool Process_Preinit ( string[] arguments ) {
  if ( arguments.length > 1 ) {
    switch ( arguments[1] ) {
      default: break;
      case "--help":
        stl.writeln(
`DTOADQ, a light transport SDF video renderer and prototyping tool

Usage: dtoadq             Run real-time editor
       dtoadq [arguments] Render scene

Arguments:
  -scene       The file of scene to be rendered
  -out         The output file (mp4 format)
  -spp         Samples-per-pixel
  -fps         Frames-per-second
  -time        Seconds of scene render time
  --help       Print help (this message) and exit
  --version    Print version information and exit
        `);
      return false;
      case "--version":
        stl.writeln("DTOADQ - Not A Damn Demo Tool");
        stl.writeln("Version - 0.0");
      return false;
    }
  }
  return true;
}

bool Process_Postinit ( string[] arguments ) {
  import stl : to;
  if ( arguments.length > 1 ) {
    string scene = "", outfile = "out.mp4";
    int spp = 1, fps = 25;
    float time = 2.0f;
    for ( int i = 1; i != arguments.length; ++ i ) {
      switch ( arguments[i] ) {
        default:
          stl.writeln("Unknown command: ", arguments[i]);
        return false;
        case "-scene": scene   = arguments[++i];          break;
        case "-out":   outfile = arguments[++i];          break;
        case "-spp":   spp     = arguments[++i].to!int;   break;
        case "-fps":   fps     = arguments[++i].to!int;   break;
        case "-time":  time    = arguments[++i].to!float; break;
      }
    }
    outfile = stl.replace(outfile, ".mp4", ".y4m");
    static import emitter.video;
    emitter.video.Render(scene, outfile, spp.to!int, time.to!float, fps.to!int);
  }
  return true;
}

void main(string[] arguments) {
  if ( !arguments.Process_Preinit() ) return;

  core.Initialize();
  scope ( exit ) {
    stl.writeln("Exitting");
    core.Clean_Up();
  }

  if ( !arguments.Process_Postinit() ) return;

  while ( !glfw.Should_Close_Window() && running ) {
    try {
      core.Update();
      core.Render();
      // -- close ? --
      import glfw.input;
      running = !RKey_Input(96) & core.RRunning;
    } catch ( Exception e ) {
      stl.writeln("--------------------------------------------");
      stl.writeln("Caught exception: ", e.msg);
      stl.writeln("At: ", e.file, ":", e.line);
      stl.writeln("++++++ callstack +++++++++++++++++++++++++++");
      stl.writeln(e.info);
      stl.writeln("--------------------------------------------");
    }
  }
}
