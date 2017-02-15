/** Check AOD.d instead, this module is reserved for engine use only */
module AODCore.realm;

import AODCore.aod;
import derelict.devil.il;
import derelict.devil.ilu;
import derelict.devil.ilut;
import derelict.freetype.ft;
import derelict.openal.al;
import derelict.opengl3.gl3;
import derelict.sdl2.sdl;
import derelict.vorbis.file;
import derelict.vorbis.vorbis;
import std.stdio : writeln;
/** */
private SDL_Window* screen = null;
/** Returns the SDL Window */
SDL_Window* R_SDL_Window() { return screen; }

/** */
class Realm {
/** objects in realm, index [layer][it]*/
  Render_Base[][] objects;
/** objects to remove at end of each frame */
  Render_Base[] objs_to_rem;
  bool cleanup_this_frame;
  Render_Base add_after_cleanup;

/** colour to clear buffer with */
  GLfloat bg_red, bg_blue, bg_green;

/** if the realm run loop has started yet */
  bool started;
/** */
  bool ended;
/** width/height of window */
  int width, height;
/** delta time (in milliseconds) for a frame */
  uint ms_dt;
/** calculates frames per second */
  float[20] fps = [ 0 ];
public:
/** */
  void Change_MSDT(uint ms_dt_) in {
    assert(ms_dt_ > 0);
  } body {
    ms_dt = ms_dt_;
  }

/** */
  auto R_Object_List() { return objects; }

/** */
  int R_Width ()       { return width;                        }
/** */
  int R_Height()       { return height;                       }
/** */
  float R_MS  ()       { return cast(float)ms_dt;             }
/** */
  float To_MS(float x) { return x/ms_dt; }

/** */
  this(int window_width, int window_height, uint ms_dt_,
       immutable(char)* window_name, immutable(char)* icon = "") {
    Util.Seed_Random();
    width  = window_width;
    height = window_height;
    ended = 0;
    ms_dt = ms_dt_;
    import std.conv : to;
    import derelict.util.exception;
    import std.stdio;

    template Load_Library(string lib, string params) {
      const char[] Load_Library =
        "try { " ~ lib ~ ".load(" ~ params ~ ");" ~
        "} catch ( DerelictException de ) {" ~
            "writeln(\"--------------------------------------------------\");"~
            "writeln(\"Failed to load: " ~ lib ~ ", \" ~ to!string(de));"     ~
        "}";
    }

    mixin(Load_Library!("DerelictGL3"       ,""));
    // mixin(Load_Library!("DerelictGL"        ,""));
    version(linux) {
      mixin(Load_Library!("DerelictSDL2", "SharedLibVersion(2, 0, 2)"));
    } else {
      mixin(Load_Library!("DerelictSDL2",
                          "\"SDL2.dll\",SharedLibVersion(2 ,0 ,2)"));
    }
    mixin(Load_Library!("DerelictIL"        ,""));
    mixin(Load_Library!("DerelictILU"       ,""));
    mixin(Load_Library!("DerelictILUT"      ,""));
    mixin(Load_Library!("DerelictAL"        ,""));
    version (linux) {
      mixin(Load_Library!("DerelictVorbis"    , ""));
      mixin(Load_Library!("DerelictVorbisFile", ""));
    } else {
      mixin(Load_Library!("DerelictVorbis"    ,"\"libvorbis-0.dll\""));
      mixin(Load_Library!("DerelictVorbisFile","\"libvorbisfile-3.dll\""));
    }
    {
      import derelict.imgui.imgui;
      mixin(Load_Library!("DerelictImgui", ""));
    }

    SDL_Init ( SDL_INIT_EVERYTHING );


    screen = SDL_CreateWindow(window_name,
                SDL_WINDOWPOS_CENTERED_DISPLAY(1),
                  SDL_WINDOWPOS_CENTERED_DISPLAY(1),
                window_width, window_height,
                SDL_WINDOW_OPENGL | SDL_WINDOW_FOREIGN );

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,  24);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE,   8);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS,
                        SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
    import std.conv : to;
    if ( screen is null ) {
      throw new Exception("Error SDL_CreateWindow: "
                          ~ to!string(SDL_GetError()));
    }
    int[][string] asdf;

    if ( SDL_GL_CreateContext(screen) is null ) {
      throw new Exception("Error SDL_GL_CreateContext: "
                          ~ to!string(SDL_GetError()));
    }

    try {
      DerelictGL3.reload();
      // DerelictGL.reload();
    } catch ( DerelictException de ) {
      writeln("\n----------------------------------------------------------\n");
      writeln("Failed to reload DerelictGL3: " ~ to!string(de));
      writeln("\n----------------------------------------------------------\n");
    }

    // glEnable        (GL_TEXTURE_2D);
    // glEnable        (GL_BLEND);
    // glBlendEquation (GL_FUNC_ADD);
    // glBlendFunc     (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // glDisable       (GL_CULL_FACE);
    // glDisable       (GL_DEPTH_TEST);
    // glEnable        (GL_SCISSOR_TEST);
    // glEnable(GL_BLEND);
    if ( icon != "" ) {
      SDL_Surface* ico = SDL_LoadBMP(icon);
      if ( ico == null ) {
        writeln("Error loading window icon");
      }
      SDL_SetWindowIcon(screen, ico);
    }

    glEnable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    ilInit();
    iluInit();
    ilutInit();
    if ( !ilutRenderer(ILUT_OPENGL) )
      writeln("Error setting ilut Renderer to ILUT_OPENGL");
    glViewport(0, 0, window_width, window_height);

    { // others
      SoundEng.Set_Up();
      UI.Init();
      static import AODCore.shader;
      AODCore.shader.Create_Default();
      /* objs_to_rem = []; */
      /* bg_red   = 0; */
      /* bg_blue  = 0; */
      /* bg_green = 0; */
      InputEngine.Refresh_Input();
    }
    Camera.Set_Position(Vector(0, 0));
    Camera.Set_Size(Vector(cast(float)window_width,
                           cast(float)window_height));
  }

/** */
  Render_Base Add(Render_Base o) in {
    assert(o !is null);
  } body {
    int l = o.R_Layer();
    if ( objects.length <= l ) objects.length = l+1;
    objects[l] ~= o;
    o.Added_To_Realm();
    return o;
  }
/** */
  void End() {
    Clean_Up(null);
    End_Sound();
    ended = true;
    import std.c.stdlib;
    exit(0);
  }
/** */
  void Remove(Render_Base o) in {
    assert(o !is null);
  } body {
    objs_to_rem ~= o;
  }
/**  */
  void Clean_Up(Render_Base rendereable) {
    cleanup_this_frame = true;
    add_after_cleanup  = rendereable;
  }
/** */
  void Set_BG_Colours(GLfloat r, GLfloat g, GLfloat b) {
    bg_red = r;
    bg_green = g;
    bg_blue = b;
  }

/** */
  void Run() {
    float prev_dt        = 0, // DT from previous frame
          curr_dt        = 0, // DT for beginning of current frame
          elapsed_dt     = 0, // DT elapsed between previous and this frame
          accumulated_dt = 0; // DT needing to be processed
    started = 1;
    SDL_Event _event;
    _event.user.code = 2;
    _event.user.data1 = null;
    _event.user.data2 = null;
    SDL_PushEvent(&_event);

    // so I can set up keys and not have to rely that update is ran first
    SDL_PumpEvents();
    InputEngine.Refresh_Input();
    SDL_PumpEvents();
    InputEngine.Refresh_Input();

    while ( SDL_PollEvent(&_event) ) {
      switch ( _event.type ) {
        case SDL_QUIT:
          if ( !ended )
            End();
        break;
        default: break;
      }
    }
      prev_dt = cast(float)SDL_GetTicks();
    while ( true ) {
      // refresh time handlers
      curr_dt = cast(float)SDL_GetTicks();
      elapsed_dt = curr_dt - prev_dt;
      accumulated_dt += elapsed_dt;
      //-----
      UI.New_Frame();
      //-----
      // refresh calculations
      while ( accumulated_dt >= ms_dt ) {
        // sdl
        SDL_PumpEvents();
        InputEngine.Refresh_Input();
        if ( keystate[ SDL_SCANCODE_ESCAPE ] )
          End();

        // actual update
        accumulated_dt -= ms_dt;
        Update();
        if ( ended ) break;

        string tex;
        string to_handle;
        bool alnum;
        char* chptr = null;

        /* auto input = Console::input->R_Str(), */
        /*      input_after = Console::input_after->R_Str(); */

        while ( SDL_PollEvent(&_event) ) {
          UI.Process_Event(&_event);
          switch ( _event.type ) {
            default: break;
            case SDL_QUIT:
              if ( !ended )
                End();
            return;
          }
        }
      }

      if ( ended ) {
        destroy(this);
        return;
      }

      { // refresh screen
        fps = elapsed_dt ~ fps[0 .. $-1];
        Render();
      }

      { // sleep until temp dt reaches ms_dt
        float temp_dt = accumulated_dt;
        temp_dt = cast(float)(SDL_GetTicks()) - curr_dt;
        while ( temp_dt < ms_dt ) {
          SDL_PumpEvents();
          SDL_Delay(1);
          temp_dt = cast(float)(SDL_GetTicks()) - curr_dt;
        }
      }
      // set current frame timemark
      prev_dt = curr_dt;
    }
  }
/** */
  float R_FPS() {
    float fps_count = 0;
    foreach ( f; fps ) {
      fps_count += f/fps.length;
    }
    return 1000/fps_count;
  }
/** */
  void Update() {
    // update objects
    foreach ( ent_l; objects ) {
      foreach( ent; ent_l ) {
        ent.Update();
      }
    }
    // find duplicates of removals
    if ( objs_to_rem.length > 1 ) {
      for ( int it = 0; it < objs_to_rem.length-1; ++ it ) {
        for ( int ot = it+1; ot < objs_to_rem.length; ++ ot ) {
          if ( objs_to_rem[it] is objs_to_rem[ot] ) {
            objs_to_rem = Util.Remove(objs_to_rem, ot);
          }
        }
      }
    }
    // remove objects
    foreach ( rem_it; 0 .. objs_to_rem.length ) {
      if ( objs_to_rem[rem_it] is null ) continue;
      int layer_it = objs_to_rem[rem_it].R_Layer();
      foreach ( obj_it; 0 .. objects[layer_it].length ) {
        if ( objects[layer_it][obj_it] is objs_to_rem[rem_it] ) {
          writeln("destroying: ");
          writeln(objects[layer_it]);
          writeln(objects[layer_it][obj_it]);
          destroy(objects[layer_it][obj_it]);
          objects[layer_it][obj_it] = null;
          objects[layer_it] = objects[layer_it][0 .. obj_it] ~
                              objects[layer_it][obj_it+1 .. $];
          break;
        }
      }
    }
    objs_to_rem = [];

    // destroy everything this frame?
    if ( cleanup_this_frame ) {
      cleanup_this_frame = false;
      for ( int i = 0; i != objects.length; ++ i ) {
        for ( int j = 0; j != objects[i].length; ++ j ) {
          destroy(objects[i][j]);
          objects[i][j] = null;
        }
      }
      objects = [];
      if ( add_after_cleanup !is null ) {
        Add(add_after_cleanup);
      }
      Clean_Up_Sound();
      objs_to_rem = []; // in case destructors added new objs
    }
  }

  ~this() {
    // todo...
    SDL_DestroyWindow(screen);
    SDL_Quit();
  }

/** */
  void Render() {
    glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);
    glClearColor(bg_red,bg_green,bg_blue,0);
    // TODO: seperate this with func pointers or something? idk
    // static import map;
    // map.Render();
    // --- rendereables ---
    for ( size_t layer = objects.length-1; layer != -1; -- layer ) {
      foreach ( obj ; objects[layer] ) {
        obj.Render();
      }
    }
    import derelict.imgui.imgui;
    igRender(); // render imgui

    SDL_GL_SwapWindow(screen);
  }
}
