/** Most the code here from examples in: https://github.com/ocornut/imgui
                          and  https://github.com/Extrawurst/imguid_test */
/** I thought of integrating this more into the AOD engine, as a lot of this
      could be rewritten with what already exists, but might not happen until
      there is scripting support. */
module AODCore.imgui;
import derelict.opengl3.gl;
import derelict.opengl3.gl3;
import derelict.sdl2.sdl;
import derelict.imgui.imgui;
import AODCore.realm : R_SDL_Window;

private GLuint font_texture = 0;
private bool[3] mouse = [ false, false, false ];
private double time = 0.0f;
private float mouse_wheel = 0.0f;
private int shader_handle, vert_handle, frag_handle,
            attrib_location_tex, attrib_location_proj_mtx,
            attrib_location_position, attrib_location_uv,
            attrib_location_colour;
private uint vbo_handle, vao_handle, elements_handle;

void Init() {
  ImGuiIO* io = igGetIO();
  io.KeyMap[ImGuiKey_Tab        ] = SDL_SCANCODE_TAB       - 1;
  io.KeyMap[ImGuiKey_LeftArrow  ] = SDL_SCANCODE_LEFT      - 1;
  io.KeyMap[ImGuiKey_RightArrow ] = SDL_SCANCODE_RIGHT     - 1;
  io.KeyMap[ImGuiKey_UpArrow    ] = SDL_SCANCODE_UP        - 1;
  io.KeyMap[ImGuiKey_DownArrow  ] = SDL_SCANCODE_DOWN      - 1;
  io.KeyMap[ImGuiKey_PageUp     ] = SDL_SCANCODE_PAGEUP    - 1;
  io.KeyMap[ImGuiKey_PageDown   ] = SDL_SCANCODE_PAGEDOWN  - 1;
  io.KeyMap[ImGuiKey_Home       ] = SDL_SCANCODE_HOME      - 1;
  io.KeyMap[ImGuiKey_End        ] = SDL_SCANCODE_END       - 1;
  io.KeyMap[ImGuiKey_Delete     ] = SDL_SCANCODE_DELETE    - 1;
  io.KeyMap[ImGuiKey_Backspace  ] = SDL_SCANCODE_HOME      - 2;
  io.KeyMap[ImGuiKey_Enter      ] = SDL_SCANCODE_RETURN    - 1;
  io.KeyMap[ImGuiKey_Escape     ] = SDL_SCANCODE_ESCAPE    - 1;
  io.KeyMap[ImGuiKey_A          ] = SDL_SCANCODE_A         - 1;
  io.KeyMap[ImGuiKey_C          ] = SDL_SCANCODE_C         - 1;
  io.KeyMap[ImGuiKey_V          ] = SDL_SCANCODE_V         - 1;
  io.KeyMap[ImGuiKey_X          ] = SDL_SCANCODE_X         - 1;
  io.KeyMap[ImGuiKey_Y          ] = SDL_SCANCODE_Y         - 1;
  io.KeyMap[ImGuiKey_Z          ] = SDL_SCANCODE_Z         - 1;

  io.RenderDrawListsFn  = &Render_Draw_Lists;
  io.SetClipboardTextFn = &SClipboard_Text;
  io.GetClipboardTextFn = &RClipboard_Text;
}

extern(C) void SClipboard_Text(const(char)* text) nothrow {
  SDL_SetClipboardText(text);
}

extern(C) const(char)* RClipboard_Text() nothrow {
  return SDL_GetClipboardText();
}

void End() {
  igShutdown();
}

void Process_Event(SDL_Event* event) {
  ImGuiIO* io = igGetIO();
  switch ( event.type ) {
    default: break;
    case SDL_MOUSEWHEEL:
      if ( event.wheel.y > 0 ) mouse_wheel =  1;
      if ( event.wheel.y < 0 ) mouse_wheel = -1;
    break;
    case SDL_MOUSEBUTTONDOWN:
      if ( event.button.button == SDL_BUTTON_LEFT   ) mouse[0] = true;
      if ( event.button.button == SDL_BUTTON_RIGHT  ) mouse[1] = true;
      if ( event.button.button == SDL_BUTTON_MIDDLE ) mouse[2] = true;
    break;
    case SDL_TEXTINPUT:
      ImGuiIO_AddInputCharactersUTF8(event.text.text.ptr);
    break;
    case SDL_KEYDOWN: case SDL_KEYUP:
      int key = event.key.keysym.sym & ~SDLK_SCANCODE_MASK;
      io.KeysDown[key] = (event.type == SDL_KEYDOWN);
      io.KeysDown[key] = (event.type == SDL_KEYDOWN);
      io.KeyShift = ( (SDL_GetModState ( ) & KMOD_SHIFT ) != 0 );
      io.KeyCtrl  = ( (SDL_GetModState ( ) & KMOD_CTRL  ) != 0 );
      io.KeyAlt   = ( (SDL_GetModState ( ) & KMOD_ALT   ) != 0 );
    break;
  }
}

void New_Frame() {
  if ( !font_texture )
    Create_Device_Objects();
  ImGuiIO* io = igGetIO();

  // setup display size
  int w, h;
  int display_w, display_h;
  SDL_GetWindowSize(R_SDL_Window(), &w, &h);
  SDL_GL_GetDrawableSize(R_SDL_Window(), &display_w, &display_h);
  io.DisplaySize = ImVec2(cast(float)w, cast(float)h);
  io.DisplayFramebufferScale = ImVec2(w > 0 ? (cast(float)display_w/w) : 0,
                                      h > 0 ? (cast(float)display_h/h) : 0);

  // Setup time step
  double current_time = SDL_GetTicks() / 1000.0;
  io.DeltaTime = time > 0.0 ? cast(float)(current_time - time) :
                                cast(float)(1.0f / 60.0f);
  time = current_time;

  // Setup inputs
  int mx, my;
  Uint32 mask = SDL_GetMouseState(&mx, &my);
  if (SDL_GetWindowFlags(R_SDL_Window()) & SDL_WINDOW_MOUSE_FOCUS)
      io.MousePos = ImVec2(cast(float)mx, cast(float)my);
  else
      io.MousePos = ImVec2(-1, -1);
  // If a mouse press event came, always pass it as "mouse held this frame",
  // we don't miss click-release events that are shorter than 1 frame.
  io.MouseDown[0] = mouse[0] || (mask & SDL_BUTTON(SDL_BUTTON_LEFT  )) != 0;
  io.MouseDown[1] = mouse[1] || (mask & SDL_BUTTON(SDL_BUTTON_RIGHT )) != 0;
  io.MouseDown[2] = mouse[2] || (mask & SDL_BUTTON(SDL_BUTTON_MIDDLE)) != 0;
  mouse[0] = mouse[1] = mouse[2] = false;

  io.MouseWheel = mouse_wheel;
  mouse_wheel = 0.0f;

  // Hide OS mouse cursor if ImGui is drawing it
  SDL_ShowCursor(io.MouseDrawCursor ? 0 : 1);

  // Start the frame
  igNewFrame();
}



private void Create_Device_Objects() {
  // Backup GL state
  GLint last_texture, last_array_buffer, last_vertex_array;
  glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
  glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &last_array_buffer);
  glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &last_vertex_array);

  const GLchar* vertex_shader =
      "#version 330\n"
      "uniform mat4 ProjMtx;\n"
      "in vec2 Position;\n"
      "in vec2 UV;\n"
      "in vec4 Color;\n"
      "out vec2 Frag_UV;\n"
      "out vec4 Frag_Color;\n"
      "void main() {\n"
      "	Frag_UV = UV;\n"
      "	Frag_Color = Color;\n"
      "	gl_Position = ProjMtx * vec4(Position.xy,0,1);\n"
      "}\n";

  const GLchar* fragment_shader =
      "#version 330\n"
      "uniform sampler2D Texture;\n"
      "in vec2 Frag_UV;\n"
      "in vec4 Frag_Color;\n"
      "out vec4 Out_Color;\n"
      "void main() {\n"
      "	Out_Color = Frag_Color * texture( Texture, Frag_UV.st);\n"
      "}\n";

  shader_handle = glCreateProgram();
  vert_handle   = glCreateShader(GL_VERTEX_SHADER);
  frag_handle   = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource (vert_handle, 1, &vertex_shader,   null);
  glShaderSource (frag_handle, 1, &fragment_shader, null);
  glCompileShader (vert_handle);
  glCompileShader (frag_handle);
  glAttachShader  (shader_handle, vert_handle);
  glAttachShader  (shader_handle, frag_handle);
  glLinkProgram   (shader_handle);

  attrib_location_tex      = glGetUniformLocation (shader_handle, "Texture"  );
  attrib_location_proj_mtx = glGetUniformLocation (shader_handle, "ProjMtx"  );
  attrib_location_position = glGetAttribLocation  (shader_handle, "Position" );
  attrib_location_uv       = glGetAttribLocation  (shader_handle, "UV"       );
  attrib_location_colour   = glGetAttribLocation  (shader_handle, "Color"    );

  glGenBuffers(1, &vbo_handle);
  glGenBuffers(1, &elements_handle);

  glGenVertexArrays         (1, &vao_handle);
  glBindVertexArray         (vao_handle);
  glBindBuffer              (GL_ARRAY_BUFFER, vbo_handle);
  glEnableVertexAttribArray (attrib_location_position );
  glEnableVertexAttribArray (attrib_location_uv       );
  glEnableVertexAttribArray (attrib_location_colour   );

  // from:  #define OFFSETOF(TYPE, ELEMENT) ((size_t)&(((TYPE *)0).ELEMENT))
  template offset_of ( string element ) {
    const char[] offset_of =
      `cast(GLvoid*)(&((cast(ImDrawVert*)0).`~element~`))`;
  }
  import std.stdio;
  glVertexAttribPointer(attrib_location_position, 2, GL_FLOAT, GL_FALSE,
    ImDrawVert.sizeof, mixin(offset_of!("pos")));
  glVertexAttribPointer(attrib_location_uv, 2, GL_FLOAT, GL_FALSE,
    ImDrawVert.sizeof, mixin(offset_of!("uv")));
  glVertexAttribPointer(attrib_location_colour, 4, GL_UNSIGNED_BYTE, GL_TRUE,
    ImDrawVert.sizeof, mixin(offset_of!("col")));

  Create_Fonts_Texture();

  // Restore modified GL state
  glBindTexture     (GL_TEXTURE_2D, last_texture);
  glBindBuffer      (GL_ARRAY_BUFFER, last_array_buffer);
  glBindVertexArray (last_vertex_array);
}

void Create_Fonts_Texture() {
  ImGuiIO* io = igGetIO();
  ubyte* pixels;
  int width, height;
  ImFontAtlas_GetTexDataAsRGBA32( io.Fonts, &pixels, &width, &height, null);

  GLint last_texture;
  // upload texture to graphics system
  glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
  glGenTextures(1, &font_texture);
  glBindTexture(GL_TEXTURE_2D, font_texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA,
                                          GL_UNSIGNED_BYTE, pixels);
  // store identifier
  ImFontAtlas_SetTexID(io.Fonts, cast(void*)font_texture);

  // return state
  glBindTexture(GL_TEXTURE_2D, last_texture);
}

private void Invalidate_Device_Objects() {
  if ( vao_handle      ) glDeleteVertexArrays(1, &vao_handle      );
  if ( vbo_handle      ) glDeleteBuffers     (1, &vbo_handle      );
  if ( elements_handle ) glDeleteBuffers     (1, &elements_handle );
  vao_handle = vbo_handle = elements_handle = 0;

  if ( shader_handle && vert_handle )
    glDetachShader(shader_handle, vert_handle);
  if ( vert_handle ) glDeleteShader(vert_handle);
  vert_handle = 0;

  if ( shader_handle && frag_handle )
    glDetachShader(shader_handle, frag_handle);
  if ( frag_handle ) glDeleteShader(frag_handle);
  frag_handle = 0;

  if ( shader_handle ) glDeleteProgram(shader_handle);
  shader_handle = 0;

  if ( font_texture ) {
    glDeleteTextures(1, &font_texture);
    ImFontAtlas_SetTexID(igGetIO().Fonts, null);
    font_texture = 0;
  }
}

void dwriteln(T...)(T a) nothrow {
  import std.stdio;
  try { writeln(a); } catch ( Exception ) {}
}


/** main rendering function */
extern(C) void Render_Draw_Lists(ImDrawData* data) nothrow {
  // Avoid rendering when minimized, scale coordinates for retina displays
  // (screen coordinates != framebuffer coordinates)
  ImGuiIO* io = igGetIO();
  int fb_width  = cast(int)(io.DisplaySize.x * io.DisplayFramebufferScale.x);
  int fb_height = cast(int)(io.DisplaySize.y * io.DisplayFramebufferScale.y);
  if (fb_width == 0 || fb_height == 0)
      return;
  // draw_data.ScaleClipRects(io.DisplayFramebufferScale);

  // Backup GL state
  GLint last_program, last_texture, last_active_texture, last_array_buffer,
        last_element_array_buffer, last_vertex_array, last_blend_src,
        last_blend_dst, last_blend_equation_rgb, last_blend_equation_alpha;
  GLint[4] last_viewport, last_scissor_box;
  glGetIntegerv(GL_CURRENT_PROGRAM,              &last_program             );
  glGetIntegerv(GL_TEXTURE_BINDING_2D,           &last_texture             );
  glGetIntegerv(GL_ACTIVE_TEXTURE,               &last_active_texture      );
  glGetIntegerv(GL_ARRAY_BUFFER_BINDING,         &last_array_buffer        );
  glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, &last_element_array_buffer);
  glGetIntegerv(GL_VERTEX_ARRAY_BINDING,         &last_vertex_array        );
  glGetIntegerv(GL_BLEND_SRC,                    &last_blend_src           );
  glGetIntegerv(GL_BLEND_DST,                    &last_blend_dst           );
  glGetIntegerv(GL_BLEND_EQUATION_RGB,           &last_blend_equation_rgb  );
  glGetIntegerv(GL_BLEND_EQUATION_ALPHA,         &last_blend_equation_alpha);
  glGetIntegerv(GL_VIEWPORT,                     last_viewport.ptr         );
  glGetIntegerv(GL_SCISSOR_BOX,                  last_scissor_box.ptr      );

  // Setup render state:alpha-blending enabled, no face culling,
  // no depth testing, scissor enabled
  glEnable        (GL_BLEND);
  glBlendEquation (GL_FUNC_ADD);
  glBlendFunc     (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDisable       (GL_CULL_FACE);
  glDisable       (GL_DEPTH_TEST);
  glEnable        (GL_SCISSOR_TEST);
  glActiveTexture (GL_TEXTURE0);

  const float width  = io.DisplaySize.x,
              height = io.DisplaySize.y;
  const float[4][4] ortho_projection = [
      [ 2.0f/io.DisplaySize.x, 0.0f,                   0.0f,  0.0f ],
      [ 0.0f,                  2.0f/-io.DisplaySize.y, 0.0f,  0.0f ],
      [ 0.0f,                  0.0f,                  -1.0f,  0.0f ],
      [-1.0f,                  1.0f,                   0.0f,  1.0f ],
  ];
  glUseProgram(shader_handle);
  glUniform1i(attrib_location_tex, 0);
  glUniformMatrix4fv(attrib_location_proj_mtx, 1, GL_FALSE,
                     &ortho_projection[0][0]);
  glBindVertexArray(vao_handle);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_handle);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elements_handle);

  foreach ( n; 0 .. data.CmdListsCount ) {
    ImDrawList* cmd_list = data.CmdLists[n];
      import std.stdio;
    ImDrawIdx* idx_buffer_offset;

    auto countVertices = ImDrawList_GetVertexBufferSize(cmd_list);
    auto countIndices = ImDrawList_GetIndexBufferSize(cmd_list);

    glBufferData(GL_ARRAY_BUFFER, countVertices * ImDrawVert.sizeof,
           cast(GLvoid*)ImDrawList_GetVertexPtr(cmd_list,0), GL_STREAM_DRAW);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, countIndices * ImDrawIdx.sizeof,
           cast(GLvoid*)ImDrawList_GetIndexPtr(cmd_list,0), GL_STREAM_DRAW);

    auto cmdCnt = ImDrawList_GetCmdSize(cmd_list);

    foreach(i; 0..cmdCnt) {
      auto pcmd = ImDrawList_GetCmdPtr(cmd_list, i);
      glBindTexture(GL_TEXTURE_2D, cast(GLuint)pcmd.TextureId);
      glScissor(cast(int)pcmd.ClipRect.x, cast(int)(height-pcmd.ClipRect.w),
                cast(int)(pcmd.ClipRect.z - pcmd.ClipRect.x),
                cast(int)(pcmd.ClipRect.w - pcmd.ClipRect.y));
      glDrawElements(GL_TRIANGLES, pcmd.ElemCount, GL_UNSIGNED_SHORT,
                      idx_buffer_offset);

      idx_buffer_offset += pcmd.ElemCount;
    }
  }

  // Restore modified GL state
  glUseProgram(last_program);
  glActiveTexture(last_active_texture);
  glBindTexture(GL_TEXTURE_2D, last_texture);
  glBindVertexArray(last_vertex_array);
  glBindBuffer(GL_ARRAY_BUFFER, last_array_buffer);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, last_element_array_buffer);
  glBlendEquationSeparate(last_blend_equation_rgb,last_blend_equation_alpha);
  glBlendFunc(last_blend_src, last_blend_dst);
  glViewport(last_viewport[0],last_viewport[1], cast(GLsizei)last_viewport[2],
             cast(GLsizei)last_viewport[3]);
  glScissor(last_scissor_box[0], last_scissor_box[1],
             cast(GLsizei)last_scissor_box[2],
             cast(GLsizei)last_scissor_box[3]);
}
