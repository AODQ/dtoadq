module glfw.input;

import derelict.glfw3;
import stl;

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

private struct CameraInputInfo {
  struct InputAxis {
    float cos_theta, sin_theta;
  }

  InputAxis position, direction;
  bool update_camera;

  static auto New ( ) {
    CameraInputInfo info;
    bool update_camera;
    auto Input_Arr()(int u, int l, int r, int d) {
      float cos_theta = cast(float)(RKey_Input(u) - RKey_Input(d)),
            sin_theta = cast(float)(RKey_Input(l) - RKey_Input(r));
      update_camera |= stl.fabs(cos_theta)+stl.fabs(sin_theta) > float.epsilon;
      return InputAxis(cos_theta, sin_theta);
    }
    return CameraInputInfo(Input_Arr(87, 65, 68, 83), Input_Arr(75, 72, 76, 74),
                           update_camera);
  }
}

float RKey_Velocity() {
  if ( RKey_Input(32) ) return 0.01f;
  if ( RKey_Input(340) ) return 0.001f;
  return 0.005f;
}

import ocl : Camera;
bool Update_Camera ( ref Camera camera ) {
  bool update_camera = false;
  import stl : cos, sin;

  auto input_info = CameraInputInfo.New();
  { // position
    float A = stl.PI - camera.lookat[0]*2.0f*PI;
    float X = input_info.position.sin_theta,
          Y = -input_info.position.cos_theta;
    float vel = RKey_Velocity();
    camera.position[0] +=  (A.cos*X + -A.sin*Y)*vel;
    camera.position[2] += -(A.sin*X + A.cos*Y)*vel;
    camera.position[1] += cast(float)(RKey_Input(81) - RKey_Input(69))*vel;
    input_info.update_camera |= RKey_Input(81) | RKey_Input(69);
  }
  { // lookat
    camera.lookat[0] += input_info.direction.sin_theta*-0.01f;
    camera.lookat[1] += input_info.direction.cos_theta*0.01f;
  }
  return input_info.update_camera;
}
