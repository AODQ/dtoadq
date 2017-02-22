module raycast;
import opencl;
import opencl_program : Test_raycast_string;
import globals;
import core.time : MonoTime;

auto TF3(float[3] a) {
  cl_float3 vec;
  vec.x = a[0];
  vec.y = a[1];
  vec.z = a[2];
  return vec;
}

struct Triangle {
  cl_float3 A, B, C;
  uint material_index;
}

struct Material {
  cl_float3 ambient, diffuse, specular, emission, shininess;
}

Triangle[] Create_Triangles(float[3][][] triangles) {
  Triangle[] results;
  foreach ( sect; triangles ) {
    cl_float3 A = TF3(sect[0]),
              B = TF3(sect[1]),
              C = TF3(sect[2]);
    uint index = sect[3][0].to!uint;
    results ~= Triangle(A, B, C, index);
  }
  return results;
}

Material[] Create_Materials(float[3][][] material) {
  Material[] results;
  foreach ( sect; material ) {
    cl_float3 A = TF3(sect[0]),
              D = TF3(sect[1]),
              S = TF3(sect[2]),
              E = TF3(sect[3]),
              H = TF3(sect[4]);
    results ~= Material(A, D, S, E, H);
  }
  return results;
}

class Raycaster : AOD.Entity {
  immutable(int) Img_dim = 512;
  OpenCLProgram program;
  OpenCLImage result;
  OpenCLBuffer!Triangle vertice_buffer;
  Triangle[] vertices;
  OpenCLBuffer!Material material_buffer;
  Material[] material;
  OpenCLSingleton!float timer;
  MonoTime timer_start;
public:
  this ( ) {
    super();
    float[3] pA = [ 120.0f,  0.3f, 120.0f ],
             pB = [ 200.0f,  0.3f, 120.0f ],
             pC = [ 120.0f,  0.3f, 200.0f ],
             pD = [ 200.0f,  0.3f, 200.0f ],
             pE = [ 110.0f,  0.5f, 110.0f ],
             pF = [ 190.0f,  0.5f, 110.0f ],
             pG = [ 110.0f,  0.5f, 190.0f ],
             pH = [ 190.0f,  0.5f, 190.0f ];
    float[3] b = [0.0f, 0.0f, 0.0f],
             l = [1.0f, 1.0f, 1.0f];
    vertices = Create_Triangles([
      [pC, pA, pB, b], [pB, pD, pC, b],
      [pA, pC, pG, b], [pA, pG, pE, b],
      [pA, pE, pF, b], [pA, pF, pB, b],
    ]);
    writeln(vertices);
    material = Create_Materials([
      [b, b, b, l, b]
    ]);
    Set_Position(AOD.R_Window_Width/2, AOD.R_Window_Height/2);
    program = Compile(Test_raycast_string);
    program.Set_Kernel("Kernel_Raycast");
    result = program.Set_Image_Buffer(BufferType.write_only, Img_dim, 0);
    timer  = program.Set_Singleton!float(BufferType.read_only, 0.0f, 1);
    vertice_buffer  = program.Set_Buffer!Triangle(BufferType.read_only, vertices, 2);
    material_buffer = program.Set_Buffer!Material(BufferType.read_only, material, 4);
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
