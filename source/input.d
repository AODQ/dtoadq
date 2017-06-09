module input;

import derelict.glfw3;

private bool mouse_left, mouse_right, mouse_middle, mouse_x1;
private float mouse_x, mouse_y, mouse_x_stick, mouse_y_stick;

extern(C) static void Cursor_Button_Callback ( GLFWwindow* window,
                                int button, int action, int mods ) nothrow {
  bool mouse_pressed = action == GLFW_PRESS;
  if ( button == GLFW_MOUSE_BUTTON_LEFT   ) mouse_left   = mouse_pressed;
  if ( button == GLFW_MOUSE_BUTTON_RIGHT  ) mouse_right  = mouse_pressed;
  if ( button == GLFW_MOUSE_BUTTON_MIDDLE ) mouse_middle = mouse_pressed;
  if ( button == GLFW_MOUSE_BUTTON_4      ) mouse_x1     = mouse_pressed;
  if ( mouse_pressed ) {
    mouse_x_stick = mouse_x;
    mouse_y_stick = mouse_y;
  }
}

extern(C) static void Cursor_Position_Callback ( GLFWwindow* window,
                                double xpos, double ypos ) nothrow {
  mouse_x = xpos;
  mouse_y = ypos;
}

private bool[400] keystate;

extern(C) static void Key_Input_Callback ( GLFWwindow* window,
                        int key, int scancode, int action, int mods ) nothrow {
  if ( key < 0 || key >= keystate.length ) return;
  keystate[key] = (action == GLFW_PRESS || action == GLFW_REPEAT);
}

bool RMouse_Left   ( ) { return mouse_left;   }
bool RMouse_Right  ( ) { return mouse_right;  }
bool RMouse_Middle ( ) { return mouse_middle; }
bool RMouse_X1     ( ) { return mouse_x1;     }

import std.stdio;
float RMouse_X       ( ) { return mouse_x;       }
float RMouse_Y       ( ) { return mouse_y;       }
float RMouse_X_Stick ( ) { return mouse_x_stick; }
float RMouse_Y_Stick ( ) { return mouse_y_stick; }
void Unstick ( ) {
  mouse_x_stick = mouse_x;
  mouse_y_stick = mouse_y;
}

/** Inputs .. http://www.glfw.org/docs/latest/group__keys.html */
bool RKey_Input ( size_t key, bool silent = true ) in {
  assert(key >= 0 && key <= keystate.length, "Invalid key input");
} body {
  auto ks = keystate[key];
  keystate[key] = silent&&ks;
  return ks;
}
