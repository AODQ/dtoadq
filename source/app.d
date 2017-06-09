import std.stdio;
import globals;
import derelict.imgui.imgui,
       derelict.opengl3.gl3,
       derelict.glfw3,
       imgui_glfw;
static import DTOADQ = dtoadq;
static import GLFW   = glfw;
static import Files  = gui.files;

void Init () {
  GLFW.Initialize();
  DTOADQ.Initialize();
  Files.Initialize();
}

void Update () {
  // -- init frame --
  ImGuiIO* io = igGetIO();
  glfwPollEvents();
  igImplGlfwGL3_NewFrame();
  static float curr_time = 0.0f, prev_time = 0.0f;
  curr_time = glfwGetTime();
  DTOADQ.Add_Time(curr_time - prev_time);
  prev_time = curr_time;
  DTOADQ.Update();
}

void Render ( ) {
  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
  DTOADQ.Render();
  igRender();
  GLFW.Swap_Buffer();
}

bool running = true;

void main(string[] arguments) {

  // check for --version and --help
  if ( arguments.length > 1 ) {
    switch ( arguments[1] ) {
      default: break;
      case "--help":
        writeln("DTOADQ, a MLT Video Renderer. Written by AODQ.");
        writeln();
        writeln("Usage: dtoadq             Run real-time editor");
        writeln("       dtoadq [arguments] Render scene to output");
        writeln();
        writeln("Arguments:");
        writeln("  -scene       The file of the scene to be rendered");
        writeln("  -out         The output file (mp4 format)");
        writeln("  -spp         Samples-Per-Pixel of each frame");
        writeln("  -fps         Frames-Per-Second");
        writeln("  -time        The time, in seconds, to render to");
        writeln("  --help       Print Help (this message) and exit");
        writeln("  --version    Print version information and exit");
      return;
      case "--version":
        writeln("DTOADQ - Not A Damn Demo Tool");
        writeln("Version - 0.0");
      return;
    }
  }

  scope ( exit ) {
    glfwTerminate();
    igImplGlfwGL3_Shutdown();
    DTOADQ.Clean_Up();
  }

  Init();

  if ( arguments.length > 1 ) {
    string scene = "", outfile = "out.mp4";
    int spp = 1, fps = 25;
    float time = 2.0f;
    for ( int i = 1; i != arguments.length; ++ i ) {
      switch ( arguments[i] ) {
        default:
          writeln("Unknown command: ", arguments[i]);
        return;
        case "-scene": scene   = arguments[++i];          break;
        case "-out":   outfile = arguments[++i];          break;
        case "-spp":   spp     = arguments[++i].to!int;   break;
        case "-fps":   fps     = arguments[++i].to!int;   break;
        case "-time":  time    = arguments[++i].to!float; break;
      }
    }
    import VI = videorender;
    VI.Render(scene, outfile, spp.to!int, time.to!float, fps.to!int);
  }

  while ( !GLFW.Should_Close_Window() && running ) {
    Update();
    Render();
    // -- close ? --
    import input;
    running = DTOADQ.RRunning;
    if ( RKey_Input(96) ) running = false;
  }
}
