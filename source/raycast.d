module raycast;
import opencl;
import opencl_program : Test_raycast_string;
import globals;
import core.time : MonoTime;
import scene;

struct RNG {
  cl_ulong[16] seed;
  cl_ulong p;
};

RNG Generate_New_RNG() {
  RNG rng;
  import std.random;
  foreach ( ref i; rng.seed ) {
    i = uniform(0, cl_ulong.max).to!(cl_ulong);
  }
  return rng;
}

class Raycaster : AOD.Entity {
  immutable(int) Img_dim = 256;
  OpenCLProgram program;
  OpenCLImage img_buffer_write, img_buffer_read, img_buffer_env;
  OpenCLBuffer!Triangle vertice_buffer;
  OpenCLBuffer!Material material_buffer;
  OpenCLSingleton!float timer;
  OpenCLSingleton!RNG rng_buffer;
  OpenCLSingleton!Camera camera_buffer;
  Camera camera;
  MonoTime timer_start;
  bool reset_camera;
public:
  this ( ) {
    super();
    Set_Position(AOD.R_Window_Width/2, AOD.R_Window_Height/2);
    auto scene = Create_Scene("test");
    program = Compile(Test_raycast_string);
    program.Set_Kernel("Kernel_Raycast");
    auto RO = BufferType.read_only,
         WO = BufferType.write_only;
    img_buffer_write = program.Set_Image_Buffer(WO, Img_dim);
    img_buffer_read  = program.Set_Image_Buffer(RO, Img_dim);
    timer  = program.Set_Singleton!float(RO, 0.0f);
    vertice_buffer  = program.Set_Buffer!Triangle(RO, scene.vertices);
    material_buffer = program.Set_Buffer!Material(RO, scene.materials);
    writeln("MATERIALS: ", Material.sizeof);
    auto rng = Generate_New_RNG();
    writeln(rng);
    rng_buffer      = program.Set_Singleton!RNG(RO, rng);
    img_buffer_env  = program.Set_Image_Buffer(RO, Img_dim);
    import image;
    img_buffer_env.data = Read_Image("testenv.tga");
    camera = Construct_Camera([0.0f, 0.0f, 0.0f], [0.0f, 1.0f, 0.0f],
                                                 [Img_dim, Img_dim]);
    camera_buffer = program.Set_Singleton!Camera(RO, camera);
    program.Write(img_buffer_env);
    timer_start = MonoTime.currTime;
  }
  ~this ()  {
    Clean_Up();
  }

  void Clean_Up ( ) {
    writeln("Cleaning up");
    program.Clean_Up();
  }

  override void Update ( ) {
    import derelict.sdl2.sdl;
    auto left     = AOD.keystate[SDL_SCANCODE_A],
         right    = AOD.keystate[SDL_SCANCODE_D],
         forward  = AOD.keystate[SDL_SCANCODE_W],
         backward = AOD.keystate[SDL_SCANCODE_S],
         up       = AOD.keystate[SDL_SCANCODE_E],
         down     = AOD.keystate[SDL_SCANCODE_Q],
         rotleft  = AOD.keystate[SDL_SCANCODE_Z],
         rotright = AOD.keystate[SDL_SCANCODE_C];
    reset_camera = ( left || right || up || down || forward || backward);
    float cx = cast(int)(right)   - cast(int)(left),
          cy = cast(int)(up)      - cast(int)(down),
          cz = cast(int)(forward) - cast(int)(backward),
          rx = cast(int)(rotleft) - cast(int)(rotright)*0.00001;
    if ( AOD.keystate[SDL_SCANCODE_LSHIFT] ) {
      cx *= 5.0f; cy *= 5.0f; cz *= 5.0f; rx *= 5.0f;
    }
    camera.position[0] += cx;
    camera.position[1] += cy*0.01f;
    camera.position[2] += cz;
    static float lmx, lmy;
    import std.math : abs;
    if ( AOD.R_Mouse_Left || abs(rx) >= 0.01 ) {
      reset_camera = true;
      float cmx = AOD.R_Mouse_X(0) - lmx,
            cmy = AOD.R_Mouse_Y(0) - lmy;
      camera.direction[0] += cmx * 0.005f;
      camera.direction[1] += rx;
      camera.direction[2] -= cmy * 0.005f;
      // normalize
      import std.math : sqrt;
      float mag = sqrt((camera.direction[0]*camera.direction[0]) +
                       (camera.direction[1]*camera.direction[1]) +
                       (camera.direction[2]*camera.direction[2]));
      if ( mag ) {
        camera.direction[0] /= mag;
        camera.direction[1] /= mag;
        camera.direction[2] /= mag;
      }
    }
    lmx = AOD.R_Mouse_X(0);
    lmy = AOD.R_Mouse_Y(0);
  }

  CLImage Run_CL() {
    auto duration = (MonoTime.currTime - timer_start).total!"msecs";
    timer.data[0] = duration/1000.0f;
    if ( reset_camera ) {
      import functional;
      float[] buffer;
      buffer.length = 4*Img_dim*Img_dim;
      buffer = buffer.map!(n => n = 0).array;
      img_buffer_read.data = buffer;
      img_buffer_write.data = buffer;
      program.Write(img_buffer_write);
    } else {
      img_buffer_read.data = img_buffer_write.data;
    }
    camera_buffer.data[0] = camera;

    program.Write(timer);
    program.Write(camera_buffer);
    program.Write(img_buffer_read);
    program.Run([Img_dim, Img_dim, 1], [1, 1, 1]);
    program.Read_Image(img_buffer_write);
    reset_camera = false;
    return CLImage(img_buffer_write.data, Img_dim, Img_dim);
  }

  AOD.SheetContainer CLImage_To_Image(CLImage image) {
    import derelict.opengl3.gl3;
    GLuint texture;
    glGetError();
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width, image.height, 0,
        GL_RGBA, GL_FLOAT, cast(void*)image.buffer.ptr
    );
    glBindTexture(GL_TEXTURE_2D, 0);
    return AOD.SheetContainer(texture, image.width, image.height);
  }

  override void Render ( ) {
    writeln("FPS: ", AOD.R_FPS());
    Set_Sprite(CLImage_To_Image(Run_CL()));
    Set_Size(AOD.Vector(512.0f, 512.0f), true);
    super.Render();
  }
}
