import std.stdio;
import globals;
import derelict.imgui.imgui,
       derelict.opengl3.gl3,
       derelict.glfw3,
       imgui_glfw;
static import DTOADQ = dtoadq;
static import GLFW   = glfw;

void Init ( ) {
  GLFW.Initialize();
  DTOADQ.Initialize();
}

void Update () {
  // -- init frame --
  ImGuiIO* io = igGetIO();
  glfwPollEvents();
  igImplGlfwGL3_NewFrame();
  DTOADQ.Update(glfwGetTime());
}

void Render ( ) {
  glClearColor(0.0f, 1.0f, 0.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
  DTOADQ.Render();
  igRender();
  GLFW.Swap_Buffer();
}

bool running = true;

void main() {
  scope ( exit ) {
    writeln("Terminating glfw...");
    glfwTerminate();
    igImplGlfwGL3_Shutdown();
    writeln("Terminating OCL/DTOADQ");
    DTOADQ.Clean_Up();
    writeln("ended");
  }

  Init();

  while ( !GLFW.Should_Close_Window() && running ) {
    Update();
    Render();
    // -- close ? --
    import input;
    if ( RKey_Input(96) ) running = false;
  }
}
