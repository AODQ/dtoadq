module raycast;
import opencl;
import opencl_program : Test_raycast_string;
import globals;
import core.time : MonoTime;

struct float3 {
  float[3] f;
}

struct Triangle {
public:
  float3[3] vertices;
}

Triangle[] Create_Triangles(float[3] a, float[3] b, float[3] c) {
  return [ Triangle ([
      float3 ( a ), float3 ( b ), float3 ( c )
    ]) ];
}


  // auto result = program.Set_Image_Buffer(BufferType.write_only, Img_dim, 0);
class Raycaster : AOD.Entity {
  immutable(int) Img_dim = 512;
  OpenCLProgram program;
  OpenCLImage result;
  OpenCLBuffer!Triangle buffer;
  OpenCLSingleton!float timer;
  MonoTime timer_start;
  Triangle[] vertices;
public:
  this ( ) {
    super();
    float[3] pA = [20.0f,  0.3f,  20.0f],
             pB = [100.0f, 0.3f,  20.0f],
             pC = [ 20.0f, 0.3f, 100.0f];
    vertices = Create_Triangles(pA, pB, pC);
    Set_Position(AOD.R_Window_Width/2, AOD.R_Window_Height/2);
    program = Compile(Test_raycast_string);
    program.Set_Kernel("Kernel_Raycast");
    result = program.Set_Image_Buffer(BufferType.write_only, Img_dim, 0);
    timer  = program.Set_Singleton!float(BufferType.read_only, 0.0f, 1);
    buffer = program.Set_Buffer!Triangle(BufferType.read_only, vertices, 2);
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
    program.Run([Img_dim, Img_dim, 1], [1, 1, 1]);
    program.Read_Image(result);
    return CLImage(result.data, Img_dim, Img_dim);
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
        GL_RGBA, GL_UNSIGNED_BYTE, cast(void*)image.buffer.ptr
    );
    glBindTexture(GL_TEXTURE_2D, 0);
    return AOD.SheetContainer(texture, image.width, image.height);
  }

  override void Render ( ) {
    Set_Sprite(CLImage_To_Image(Run_CL()));
    super.Render();
  }
}
