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
  cl_uint material_index;
  cl_uint bufsize1, bufsize2, bufsize3;
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
    Triangle t = Triangle(A, B, C);
    t.material_index = index;
    results ~= t;
  }
  return results;
}

Material Create_Material(float[3][] material) {
  Material results;
  cl_float3 A = TF3(material[0]),
            D = TF3(material[1]),
            S = TF3(material[2]),
            E = TF3(material[3]),
            H = TF3(material[4]);
  return Material(A, D, S, E, H);
}

auto Create_Materials ( ) {
  auto wall = Create_Material([
    [0.5f, 0.5f, 0.5f],
    [0.4f, 0.8f, 0.5f],
    [0.0f, 0.0f, 0.0f],
    [0.0f, 0.0f, 0.0f],
    [0.0f, 0.0f, 0.0f],
  ]);
  auto light = Create_Material([
    [0.0f, 0.0f, 0.0f],
    [0.0f, 0.0f, 0.0f],
    [0.0f, 0.0f, 0.0f],
    [1.0f, 1.0f, 1.0f],
    [0.0f, 0.0f, 0.0f],
  ]);
  auto box = Create_Material([
    [0.5f, 0.5f, 0.5f],
    [1.0f, 1.0f, 0.2f],
    [0.0f, 0.0f, 0.0f],
    [0.0f, 0.0f, 0.0f],
    [0.0f, 0.0f, 0.0f],
  ]);
  return [ wall ] ~ [ light ] ~ [ box ];
}

Triangle[] Create_Box(uint material_index) {
  float[3] pA = [ 120.0f,  0.3f, 120.0f ],
            pB = [ 200.0f,  0.3f, 120.0f ],
            pC = [ 120.0f,  0.3f, 200.0f ],
            pD = [ 200.0f,  0.3f, 200.0f ],
            pE = [ 110.0f,  0.5f, 110.0f ],
            pF = [ 190.0f,  0.5f, 110.0f ],
            pG = [ 110.0f,  0.5f, 190.0f ],
            pH = [ 190.0f,  0.5f, 190.0f ];
  float[3] m = [ material_index, 0.0f, 0.0f ];
  return Create_Triangles([
    [pC, pA, pB, m], [pB, pD, pC, m],
    [pA, pC, pG, m], [pA, pG, pE, m],
    [pA, pE, pF, m], [pA, pF, pB, m],
  ]);
}

Triangle[] Create_Light(uint material_index) {
  float[3] pA = [ 100.0f,  0.2f, 250.0f ],
           pB = [ 200.0f,  0.1f, 250.0f ],
           pC = [ 120.0f,  0.2f, 280.0f ];
  float[3] m = [ material_index, 0.0f, 0.0f ];
  return Create_Triangles([
    [pC, pA, pB, m]
  ]);
}

Triangle[] Create_Walls(uint material_index) {
  float[3] pA = [  20.0f,  -5.3f,  20.0f ],
           pB = [ 400.0f,  -5.3f,  20.0f ],
           pC = [  20.0f,  -5.3f, 200.0f ],
           pD = [ 400.0f,  -5.3f, 400.0f ],
           pE = [  10.0f,   5.5f,  10.0f ],
           pF = [ 390.0f,   5.5f,  10.0f ],
           pG = [  10.0f,   5.5f, 390.0f ],
           pH = [ 390.0f,   5.5f, 390.0f ];
  float[3] m = [ material_index, 0.0f, 0.0f ];
  return Create_Triangles([
    [pC, pA, pB, m], [pB, pD, pC, m],
    [pA, pC, pG, m], [pA, pG, pE, m],
    [pA, pE, pF, m], [pA, pF, pB, m],
  ]);
}

struct Scene {
  Material[] materials;
  Triangle[] vertices;
}

Triangle[] Create_Models ( ) {
  auto material = Create_Materials();
  // --- box ---
  return Create_Box(2) ~ Create_Light(1) ~ Create_Walls(0);
}

Scene Create_Scene ( ) {
  return Scene(Create_Materials(), Create_Models());
}

class Raycaster : AOD.Entity {
  immutable(int) Img_dim = 512;
  OpenCLProgram program;
  OpenCLImage result;
  OpenCLBuffer!Triangle vertice_buffer;
  OpenCLBuffer!Material material_buffer;
  OpenCLSingleton!float timer;
  MonoTime timer_start;
public:
  this ( ) {
    super();
    Set_Position(AOD.R_Window_Width/2, AOD.R_Window_Height/2);
    auto scene = Create_Scene();
    program = Compile(Test_raycast_string);
    program.Set_Kernel("Kernel_Raycast");
    result = program.Set_Image_Buffer(BufferType.write_only, Img_dim, 0);
    timer  = program.Set_Singleton!float(BufferType.read_only, 0.0f, 1);
    auto RO = BufferType.read_only;
    vertice_buffer  = program.Set_Buffer!Triangle(RO, scene.vertices, 2);
    material_buffer = program.Set_Buffer!Material(RO, scene.materials, 4);
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
