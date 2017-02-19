module AODCore.shader;
import std.string;
import derelict.opengl3.gl3;

/**
  Allows you to implement GLSL shaders into the rendering of AoD
*/
struct Shader {
  GLuint id;
public:
  /** postblist */
  this ( this ) {
    id = id;
  }
  /** Generates a GL Shader.*/
  this ( string vertex_file,
         string fragment_file  = "",
         string tess_ctrl_file = "",
         string tess_eval_file = "",
         string compute_file   = "" ) {
    bool is_frag = fragment_file   != "",
         is_tesc = tess_ctrl_file != "",
         is_tese = tess_eval_file != "",
         is_comp  = compute_file   != "";
    // --- create program/shaders
    id = glCreateProgram();
    GLuint vertex_ID    = glCreateShader ( GL_VERTEX_SHADER          ) ,
           fragment_ID  = glCreateShader ( GL_FRAGMENT_SHADER        ) ,
           tess_ctrl_ID = glCreateShader ( GL_TESS_CONTROL_SHADER    ) ,
           tess_eval_ID = glCreateShader ( GL_TESS_EVALUATION_SHADER ) ,
           compute_ID   = glCreateShader ( GL_COMPUTE_SHADER         ) ;
    // --- compile shader source & attach to program
    auto Create_Shader = (GLuint ID, GLuint progID, string fname, bool empty) {
      if ( !empty ) return;
      // compile shader
      static import std.file;
      char* fil = cast(char*)(std.file.read(fname).ptr);
      glShaderSource ( ID, 1, &fil, null );
      glCompileShader(ID);
      // -- check for error --
      GLint compile_status;
      glGetShaderiv(ID, GL_COMPILE_STATUS, &compile_status);
      if ( compile_status == GL_FALSE ) {
        import std.stdio : writeln;
        writeln(fname ~ " shader compilation failed");
        writeln("--------------------------------------");
        GLchar[256] error_message;
        glGetShaderInfoLog(ID, 256, null, error_message.ptr);
        writeln(error_message);
        writeln("--------------------------------------");
        return;
      }
      // attach
      glAttachShader(progID, ID);
    };
    Create_Shader( vertex_ID    , id , vertex_file    , true );
    Create_Shader( fragment_ID  , id , fragment_file  , is_frag );
    Create_Shader( tess_ctrl_ID , id , tess_ctrl_file , is_tesc );
    Create_Shader( tess_eval_ID , id , tess_eval_file , is_tese );
    Create_Shader( compute_ID   , id , compute_file   , is_comp );

    // --- link program
    glLinkProgram ( id );

    // --- cleanup
    glDeleteShader(vertex_ID    ) ;
    glDeleteShader(fragment_ID  ) ;
    glDeleteShader(tess_ctrl_ID ) ;
    glDeleteShader(tess_eval_ID ) ;
    glDeleteShader(compute_ID   ) ;
  }

  void Bind() {
    glUseProgram ( id );
  }
  static void Unbind() {
    glUseProgram ( 0 );
  }
  /** Returns the program/shader ID */
  GLuint R_Shader_ID() {
    return id;
  }
}

const GLchar* vertex_shader = `
  #version 330
  layout (location = 0) in vec2 vertices;
  layout (location = 1) in vec2 texcoord;
  out vec2 Texcoord;
  in vec3 colour;
  out vec3 Colour;
  uniform mat4 view;
  void main ( ) {
    Colour      = colour;
    Texcoord    = texcoord;
    gl_Position = view*vec4(vertices, 0.0, 1.0);
  }
`;

const GLchar* fragment_shader =  `
  #version 330
  in vec3 Colour;
  in vec2 Texcoord;
  uniform sampler2D Tex;
  void main() {
      gl_FragColor = texture(Tex, Texcoord);
  }
`;

int def_vert_handle, def_frag_handle;
int def_shader;
int def_attr_tex, def_attr_pos, def_attr_colour, def_attr_view;
GLuint VBO, VAO, EBO, VBO_UV;

GLfloat[] vertices = [
    -1.0f,  1.0f,
     1.0f,  1.0f,
     1.0f, -1.0f,
    -1.0f, -1.0f
];

GLuint[] elements = [
  0, 1, 2,
  2, 3, 0
];

void Check_Shader_Error ( int handle, string type ) {
  GLint compile_status;
  glGetShaderiv(handle, GL_COMPILE_STATUS, &compile_status);
  if ( compile_status == GL_FALSE ) {
    import std.stdio : writeln;
    writeln(type ~ " shader compilation failed");
    writeln("--------------------------------------");
    GLchar[256] error_message;
    glGetShaderInfoLog(handle, 256, null, error_message.ptr);
    writeln(error_message);
    writeln("--------------------------------------");
    assert(0);
  }
}

private void Create_Default_Shaders() {
  def_shader      = glCreateProgram();
  def_vert_handle = glCreateShader(GL_VERTEX_SHADER);
  def_frag_handle = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource     (def_vert_handle, 1, &vertex_shader, null);
  glShaderSource     (def_frag_handle, 1, &fragment_shader, null);
  glCompileShader    (def_frag_handle);
  Check_Shader_Error (def_frag_handle, "fragment");
  glCompileShader    (def_vert_handle);
  Check_Shader_Error (def_vert_handle, "vertex");
  glAttachShader     (def_shader, def_vert_handle);
  glAttachShader     (def_shader, def_frag_handle);
  glLinkProgram      (def_shader);
  glUseProgram       (def_shader);

  GLint compile_status;
  glGetProgramiv(def_shader, GL_LINK_STATUS, &compile_status);
  if ( compile_status == GL_FALSE ) {
    import std.stdio;
    writeln("link shader compilation failed");
    writeln("--------------------------------------");
    GLchar[256] error_message;
    glGetProgramInfoLog(def_shader, 256, null, error_message.ptr);
    writeln(error_message);
    writeln("--------------------------------------");
    assert(0);
  }
  // def_attr_tex    = glGetUniformLocation (def_shader, "texture"  );
  def_attr_tex  = glGetAttribLocation (def_shader, "Tex"   );
  def_attr_view = glGetUniformLocation(def_shader, "view"  );

}

void Create_Default ( ) {
  Create_Default_Shaders();

  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);
  glGenBuffers(1, &VBO_UV);
  glGenBuffers(1, &EBO);

  glBindVertexArray(VAO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  // Yes, you have to use vertices.length*float.sizeof. vertices.sizeof did
  // not work for some weird reason
  glBufferData(GL_ARRAY_BUFFER, vertices.length*float.sizeof,
                                vertices.ptr, GL_STATIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, null);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, elements.length*elements.sizeof,
                                        elements.ptr, GL_STATIC_DRAW);
}
