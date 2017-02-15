/**
  <br><br>
  The main interface to Art of Dwarficorn. Importing this file alone will
  give you access to the majority of the library. The functions defined here
  are only part of the Realm interface.
  <br>
Example:
<br>
---
// This is the "standard" way to initialize the engine. The console is set up
// first so that errors can be received as the AOD Engine is initialized.
// Afterwards the camera is adjusted to the center of the screen, the font is
// loaded , and the console key is assigned. Then we load the key config
void Init () {
  import AOD;
  Console.console_open = false;
  Console.Set_Console_Output_Type(AOD.Console.Type.Debug_In);
  initialize(16, "ART OF DWARFICORN", 640, 480);
  Camera.Set_Size(AOD.Vector(AOD.R_Window_Width(), AOD.R_Window_Height()));
  Camera.Set_Position(AOD.Vector(AOD.R_Window_Width() /2,
                                     AOD.R_Window_Height()/2));
  Text.Set_Default_Font("assets/DejaVuSansMono.ttf", 13);
  Console.initialize(Console.Type.Debug_In);
  Set_BG_Colour(.08, .08, .095);
  Load_Config();
}
---
*/
module AODCore.aod;

import derelict.opengl3.gl3;
import derelict.sdl2.sdl;
import std.stdio;
import std.string;

public import AODCore.entity, AODCore.matrix, AODCore.realm, AODCore.image,
              AODCore.animation, AODCore.render_base,
              AODCore.sound, AODCore.text, AODCore.vector, AODCore.input;
public import AODCore.shader : Shader;
public import UI     = AODCore.imgui;
public import Util   = AODCore.util;
public import CV     = AODCore.clientvars;
public import Camera = AODCore.camera;
public import Serializer = AODCore.serializer;

void GL_Error() {
  GLenum error = glGetError();
  switch ( error ) {
    default: writeln("UNKNOWN ERROR"); assert(false);
    case GL_NO_ERROR: break;
    case GL_INVALID_ENUM:
      writeln("GL_INVALID_ENUM");
    assert(false);
    case GL_INVALID_VALUE:
      writeln("GL_INVALID_VALUE");
    assert(false);
    case GL_INVALID_OPERATION:
      writeln("GL_INVALID_OPERATION");
    assert(false);
    case GL_INVALID_FRAMEBUFFER_OPERATION:
      writeln("GL_INVALID_FRAMEBUFFER_OPERATION");
    assert(false);
    case GL_OUT_OF_MEMORY:
      writeln("GL_OUT_OF_MEMORY");
    assert(false);
    // case GL_STACK_UNDERFLOW:
    //   writeln("GL_STACK_UNDERFLOW");
    // assert(false);
    // case GL_STACK_OVERFLOW:
    //   writeln("GL_STACK_OVERFLOW");
    // assert(false);
  }
}

Realm realm = null;

/** initializes the engine
  Params:
    msdt   = Amount of milliseconds between each update frame call
    name   = Name of the application (for the window title)
    width  = Window X dimension
    height = Window Y dimension
    ico    = File location of the icon to use for the application
*/
void initialize(uint msdt, string name, int width, int height, string ico = "")
in {
  assert(realm is null);
} body {
  if ( name == "" )
    name = "Art of Dwarficorn";
  realm = new Realm(width, height, msdt, name.ptr, ico.ptr);
}
/** Changes the amount of milliseconds between each update frame call */
void Change_MSDT(Uint32 ms_dt) in { assert(realm !is null); } body {
  realm.Change_MSDT(ms_dt);
}
/** Resets the engine */
@disable void Reset() in { assert(realm  is null); } body { /* ... todo ... */ }
/** Ends the engine and deallocates all resources */
void End()   in { assert(realm !is null); } body {
  realm.End();
  realm = null;
}

/** Adds rbase to the engine to be updated and rendered */
Render_Base Add (T...)(T n)
in {
  assert(realm !is null);
  assert(n.length > 0 );
} body {
  foreach ( t; n ) realm.Add(n);
  return n[0];
}
/** Removes rbase from the engine and deallocates it */
void Remove (T...)(T n)
in {
  assert(realm !is null);
  assert(n.length > 0 );
} body {
  foreach ( t; n ) realm.Remove(n);
}
/** Removes entities and playing sounds (not loaded sounds or images)
Params:
  rendereable = Adds a rendereable after cleanup is done (cleanup doesn't take
                 place until after the frame is finished)
*/
void Clean_Up(Render_Base rendereable = null) in {
  assert(realm !is null);
} body { realm.Clean_Up(rendereable); }
/** Sets the background colour when rendering
  Params:
    r = Red
    g = Green
    b = Blue
*/
void Set_BG_Colour(GLfloat r, GLfloat g, GLfloat b) in {
  assert(realm !is null);
} body {
  realm.Set_BG_Colours(r, g, b);
}

/** Runs the engine (won't return until SDL_Quit is called) */
void Run() in { assert(realm !is null); } body {
  // do splash screen
  /* import Splashscreen; */
  /* AOD.Add(new Splash); */
  // now run
  realm.Run();
}

/** Returns the current MS per frame */
float R_MS() in {assert(realm !is null);} body {
  return realm.R_MS();
}
/** Returns current fps */
float R_FPS() in {assert(realm !is null);} body {
  return realm.R_FPS();
}
/** Calculates the minimal amount of frames required for the duration to occur
  Params:
    x = duration
*/
float To_MS(float x) in{assert(realm !is null);} body {
  return realm.To_MS(x);
}

/**
  Return:
    Returns the window width in pixels
*/
int R_Window_Width() in{assert(realm !is null);} body {
  return realm.R_Width();
}
/**
  Return:
    Returns the window height in pixels
*/
int R_Window_Height() { return realm.R_Height(); }
