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
  immutable(int) Img_dim = 512;
  OpenCLProgram program;
  OpenCLImage img_buffer_write, img_buffer_read, img_buffer_env;
  OpenCLBuffer!Triangle vertice_buffer;
  OpenCLBuffer!Material material_buffer;
  OpenCLSingleton!float timer;
  OpenCLSingleton!RNG rng_buffer;
  MonoTime timer_start;
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

  CLImage Run_CL() {
    auto duration = (MonoTime.currTime - timer_start).total!"msecs";
    timer.data[0] = duration/1000.0f;
    program.Write(timer);
    program.Write(img_buffer_read);
    program.Run([Img_dim, Img_dim, 1], [1, 1, 1]);
    program.Read_Image(img_buffer_write);
    img_buffer_read.data = img_buffer_write.data;
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
    super.Render();
  }
}
