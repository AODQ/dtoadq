module softwarerenderer.raytracer;
import softwarerenderer.kernel, softwarerenderer.image,
       softwarerenderer.scene;
import camera, globals;
import std.concurrency : Tid;

import aodheme : ParseTree;

class Raytracer : AOD.Entity {
  immutable(int) Img_dim = 256;
  Camera camera;
  Image image;
  Scene scene;
  Tid[16] threads;
  ParseTree aodheme_map;
public:
  this ( string aodheme_map_ ) {
    super();
    import aodheme : Evaluate;
    aodheme_map = Evaluate(aodheme_map_);
    Set_Position(AOD.R_Window_Width/2.0f, AOD.R_Window_Height/2.0f);
    Set_Size(AOD.Vector(Img_dim, Img_dim), true);
    camera = new Camera(gln.vec3( 1.0f, 1.0f, -1.0f),
                        gln.vec3(-1.0f, 0.0f, 0.0f), Img_dim, Img_dim);
    image = new Image(Img_dim, Img_dim);


    immutable(Scene) Scene_frame = cast(immutable)scene;
    immutable(Camera)     cam  = cast(immutable)camera;
    shared(Image)         img  = cast(shared)image;
    foreach ( i; 0 .. 16 ) {
      size_t dim = Img_dim/4;
      size_t ix = i%4, iy = i/4;
      size_t lx = dim*ix, hx = dim*ix + dim,
             ly = dim*iy, hy = dim*iy + dim;
      import std.concurrency;
      threads[i] = spawn(&Raytrace_Thread_Manager, thisTid, Scene_frame,
                                                   cam, lx, ly, hx, hy);
    }
  }

  bool Update_Control ( ) {
    import derelict.sdl2.sdl;
    auto left     = AOD.RKeystate( SDL_SCANCODE_A ),
         right    = AOD.RKeystate( SDL_SCANCODE_D ),
         forward  = AOD.RKeystate( SDL_SCANCODE_W ),
         backward = AOD.RKeystate( SDL_SCANCODE_S ),
         up       = AOD.RKeystate( SDL_SCANCODE_E ),
         down     = AOD.RKeystate( SDL_SCANCODE_Q ),
         rotleft  = AOD.RKeystate( SDL_SCANCODE_Z ),
         rotright = AOD.RKeystate( SDL_SCANCODE_C );
    bool rc = left || right || forward || backward || up || down
                   || rotleft || rotright;
    float cx = cast(int)(right)   - cast(int)(left),
          cy = cast(int)(up)      - cast(int)(down),
          cz = cast(int)(forward) - cast(int)(backward),
          rx = cast(int)(rotleft) - cast(int)(rotright)*0.00001;
    if ( AOD.RKeystate( SDL_SCANCODE_LSHIFT ) ) {
      cx *= 5.0f; cy *= 5.0f; cz *= 5.0f; rx *= 5.0f;
    }
    camera.position.vector[0] += cx*0.1f;
    camera.position.vector[1] += cy*0.1f;
    camera.position.vector[2] += cz*0.1f;
    static float lmx, lmy;
    import std.math : abs;
    if ( AOD.R_Mouse_Left || abs(rx) >= 0.01 ) {
      rc = true;
      float cmx = AOD.R_Mouse_X(0) - lmx,
            cmy = AOD.R_Mouse_Y(0) - lmy;
      camera.lookat.vector[0] += cmx * 0.5f;
      camera.lookat.vector[1] += rx  * 0.1f;
      camera.lookat.vector[2] -= cmy * 0.5f;
      // normalize
      import std.math : sqrt;
      float mag = sqrt((camera.lookat.vector[0]*camera.lookat.vector[0]) +
                       (camera.lookat.vector[1]*camera.lookat.vector[1]) +
                       (camera.lookat.vector[2]*camera.lookat.vector[2]));
      if ( mag ) {
        camera.lookat.vector[0] /= mag;
        camera.lookat.vector[1] /= mag;
        camera.lookat.vector[2] /= mag;
      }
    }
    lmx = AOD.R_Mouse_X(0);
    lmy = AOD.R_Mouse_Y(0);
    return rc;
  }

  void Update_Thread_Manager ( bool update_camera ) {
    import std.concurrency;
    import core.time;
    import std.variant;
    static size_t img_id;
    size_t counter = 100;
    if ( update_camera ) {
      // writeln("Clearing image");
      ++ img_id;
      foreach ( t; threads ) {
        send(t, cast(immutable)camera);
      }
      image.Clear();
    }
    while (
      receiveTimeout(1.usecs,
        (size_t x, size_t y, PixelInfo info, size_t id) {
          if ( info.hit && img_id == id )
            image.Apply(x, y, info.pixel);
        },
        (Variant variant) {
          assert(false, "Parameter mismatch");
        }
      )
    ) {
      if ( -- counter < 0 )
        break;
    }
  }

  override void Update ( ) @trusted {
    static bool reset_camera;
    bool rupdate = Update_Control;
    writeln("RES: ", reset_camera, " RUPDATE: ", rupdate);
    Update_Thread_Manager(reset_camera && !rupdate);
    reset_camera = rupdate;
  }

  override void Render ( ) @trusted {

    static int counter = 0;
    if ( ++ counter > 60 ) {
      writeln("FPS: ", AOD.R_FPS());
      counter = 0;
    }
    Set_Sprite(image.To_OGL_Sprite());
    Set_Size(AOD.Vector(AOD.R_Window_Height, AOD.R_Window_Height), true);
    super.Render();
  }
}
