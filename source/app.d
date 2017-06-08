import std.stdio;
import globals;
import derelict.imgui.imgui,
       derelict.opengl3.gl3,
       derelict.glfw3,
       imgui_glfw;
static import DTOADQ = dtoadq;
static import GLFW   = glfw;
static import Files  = gui.files;

void Init ( bool window = true ) {
  GLFW.Initialize(window);
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
  float[] dumbydata;
  DTOADQ.Update(false, dumbydata);
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
  scope ( exit ) {
    writeln("Terminating glfw...");
    glfwTerminate();
    igImplGlfwGL3_Shutdown();
    writeln("Terminating OCL/DTOADQ");
    DTOADQ.Clean_Up();
    writeln("ended");
  }

  if ( arguments.length > 1 ) {
    Init(false);
    string scene = arguments[1],
           file  = arguments[2],
           spp   = arguments[3],
           time  = arguments[4],
           fps   = arguments[5];
    writeln("SCENE: ", scene);
    writeln("FILE: ", file);
    writeln("SPP: ", spp);
    writeln("TIME: ", time);
    writeln("FPS: ", fps);
    import VI = videorender;
    VI.Render(scene, file, spp.to!int, time.to!float, fps.to!int);
  } else
    Init();

  while ( !GLFW.Should_Close_Window() && running ) {
    Update();
    Render();
    // -- close ? --
    import input;
    if ( RKey_Input(96) ) running = false;
  }
}
