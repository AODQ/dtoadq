import std.stdio;
import globals;
import derelict.imgui.imgui,
       derelict.opengl3.gl3,
       derelict.glfw3,
       imgui_glfw;
static import raytracer;

GLFWwindow* Init ( ) {
  import derelict.glfw3.glfw3;
  writeln("GL3");
  DerelictGL3.load();
  writeln("GLfw3");
  DerelictGLFW3.load();
  writeln("imgui");
  DerelictImgui.load();

  if (!glfwInit())
		return null;
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
  glfwWindowHint(GLFW_FLOATING, GL_TRUE);
  glfwWindowHint( GLFW_REFRESH_RATE, 0 );
  glfwSwapInterval(30);
	glfwInit();

  auto window = glfwCreateWindow(Win_dim, Win_dim, "Raytracer", null, null);
  glfwSetWindowPos(window, -5, 2);

	glfwMakeContextCurrent(window);
  DerelictGL3.reload();
  igImplGlfwGL3_Init(window, true);
  raytracer.Initialize();

  import input;
  glfwSetCursorPosCallback(window, &Cursor_Position_Callback);
  glfwSetMouseButtonCallback(window, &Cursor_Button_Callback);
  glfwSetKeyCallback(window, &Key_Input_Callback);

  writeln("done i guess");
  return window;
}

void Update ( ref float[3] clear_colour ) {
  // -- init frame --
  ImGuiIO* io = igGetIO();
  glfwPollEvents();
  igImplGlfwGL3_NewFrame();

  // -- base imgui window --
  import gui.gui;
  igColorEdit3("background colour", clear_colour);
  igText("FPS: %.3f ms/frame (%.1f FPS)", 1000.0f / igGetIO().Framerate,
                                                    igGetIO().Framerate);
  igSliderInt("DIM: ", &Img_dim,        8, 1080);

  // -- raytrace update --
  raytracer.Update(glfwGetTime());
}

bool running = true;

import opencl : CLImage;

void main() {
  scope ( exit ) {
    glfwTerminate();
    igImplGlfwGL3_Shutdown();
    raytracer.Remove();
    writeln("Successfully ended");
  }
  auto window = Init();

  float[3] clear_colour = [0.1f, 0.1f, 0.1f];

  while ( !glfwWindowShouldClose(window) && running ) {
    glClearColor(clear_colour[0], clear_colour[1], clear_colour[2], 0);
    glClear(GL_COLOR_BUFFER_BIT);
    Update(clear_colour);
    // -- draw debug --
    igRender();
    glfwSwapBuffers(window);
    // -- close ? --
    import input;
    if ( RKey_Input(96) ) running = false;
  }
}
