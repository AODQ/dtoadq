module glfw;
/***
  A minimalistic implementation of my Art of Dwarfiqorn game framework, but
  trying GLFW instead of SDL and focused on rendering a single cl/gl interop
  texture.
***/
public import glfw.gl_renderer, glfw.input, glfw.image;
import derelict.imgui.imgui, derelict.opengl3.gl3, derelict.glfw3.glfw3;
import glfw.imgui;
static import glfw.config;

GLFWwindow* window;

void Initialize ( ) {
  DerelictGL3.load();
  DerelictGLFW3.load();
  DerelictImgui.load();
  if (!glfwInit()) {
    assert(false, "glfwinit failed");
  }

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE );
  glfwWindowHint(GLFW_RESIZABLE,      GL_FALSE                 );
  glfwWindowHint(GLFW_FLOATING,       GL_TRUE                  );
  glfwWindowHint( GLFW_REFRESH_RATE,  0                        );
  glfwSwapInterval(30);
	glfwInit();

  auto win_size = glfw.config.RWindow_Size();
  window = glfwCreateWindow(win_size[0], win_size[1], "DTOADQ", null, null);
  auto win_pos  = glfw.config.RWindow_Pos();
  if ( win_pos[0] != -1 && win_pos[1] != -1 )
    glfwSetWindowPos(window, win_pos[0], win_pos[1]);

	glfwMakeContextCurrent(window);
  DerelictGL3.reload();
  glClampColor(GL_CLAMP_READ_COLOR, GL_FALSE);
  igImplGlfwGL3_Init(window, true);

  import input;
  glfwSetCursorPosCallback   (window, &Cursor_Position_Callback );
  glfwSetMouseButtonCallback (window, &Cursor_Button_Callback   );
  glfwSetKeyCallback         (window, &Key_Input_Callback       );
}

auto Should_Close_Window ( ) { return glfwWindowShouldClose(window); }
void Swap_Buffer( ) { glfwSwapBuffers(window); }
