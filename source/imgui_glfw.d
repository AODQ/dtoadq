﻿module imgui_glfw;

import derelict.imgui.imgui;
import derelict.opengl3.gl3;
import derelict.glfw3.glfw3;

GLFWwindow* g_window; double g_Time = 0.0f; bool[3] g_MousePressed; float
g_MouseWheel = 0.0f; GLuint g_FontTexture = 0; int g_ShaderHandle = 0,
g_VertHandle = 0, g_FragHandle = 0; int g_AttribLocationTex = 0,
g_AttribLocationProjMtx = 0; int g_AttribLocationPosition = 0,
g_AttribLocationUV = 0, g_AttribLocationColor = 0; uint g_VboHandle,
g_VaoHandle, g_ElementsHandle;

extern(C) nothrow void igImplGlfwGL3_RenderDrawLists(ImDrawData* data)
{
    GLint last_program, last_texture;
    glGetIntegerv(GL_CURRENT_PROGRAM, &last_program);
    glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_SCISSOR_TEST);
    glActiveTexture(GL_TEXTURE0);

	auto io = igGetIO();
	const float width = io.DisplaySize.x;
	const float height = io.DisplaySize.y;
	const float[4][4] ortho_projection =
	[
		[ 2.0f/width,	0.0f,			0.0f,		0.0f ],
		[ 0.0f,			2.0f/-height,	0.0f,		0.0f ],
		[ 0.0f,			0.0f,			-1.0f,		0.0f ],
		[ -1.0f,		1.0f,			0.0f,		1.0f ],
	];
	glUseProgram(g_ShaderHandle);
	glUniform1i(g_AttribLocationTex, 0);
	glUniformMatrix4fv(g_AttribLocationProjMtx, 1, GL_FALSE, &ortho_projection[0][0]);

    glBindVertexArray(g_VaoHandle);
    glBindBuffer(GL_ARRAY_BUFFER, g_VboHandle);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_ElementsHandle);

    foreach (n; 0..data.CmdListsCount)
    {
        ImDrawList* cmd_list = data.CmdLists[n];
        ImDrawIdx* idx_buffer_offset;

        auto countVertices = ImDrawList_GetVertexBufferSize(cmd_list);
        auto countIndices = ImDrawList_GetIndexBufferSize(cmd_list);

        glBufferData(GL_ARRAY_BUFFER, countVertices * ImDrawVert.sizeof,
        cast(GLvoid*)ImDrawList_GetVertexPtr(cmd_list,0), GL_STREAM_DRAW);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, countIndices * ImDrawIdx.sizeof,
        cast(GLvoid*)ImDrawList_GetIndexPtr(cmd_list,0), GL_STREAM_DRAW);

        auto cmdCnt = ImDrawList_GetCmdSize(cmd_list);
        foreach(i; 0..cmdCnt)
        {
            auto pcmd = ImDrawList_GetCmdPtr(cmd_list, i);

            if (pcmd.UserCallback)
            {
                pcmd.UserCallback(cmd_list, pcmd);
            }
            else
            {
                glBindTexture(GL_TEXTURE_2D, cast(GLuint)pcmd.TextureId);
                glScissor(cast(int)pcmd.ClipRect.x, cast(int)(height -
                pcmd.ClipRect.w), cast(int)(pcmd.ClipRect.z - pcmd.ClipRect.x),
                cast(int)(pcmd.ClipRect.w - pcmd.ClipRect.y));
                glDrawElements(GL_TRIANGLES, pcmd.ElemCount, GL_UNSIGNED_SHORT,
                idx_buffer_offset);
            }

            idx_buffer_offset += pcmd.ElemCount;
        }
    }

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glUseProgram(last_program);
    glDisable(GL_SCISSOR_TEST);
    glBindTexture(GL_TEXTURE_2D, last_texture);
}

void igImplGlfwGL3_Init(GLFWwindow* window, bool install_callbacks)
{
	g_window = window;

	ImGuiIO* io = igGetIO();
    io.KeyMap[ImGuiKey_Tab]        = GLFW_KEY_TAB;
    io.KeyMap[ImGuiKey_LeftArrow]  = GLFW_KEY_LEFT;
    io.KeyMap[ImGuiKey_RightArrow] = GLFW_KEY_RIGHT;
    io.KeyMap[ImGuiKey_UpArrow]    = GLFW_KEY_UP;
    io.KeyMap[ImGuiKey_DownArrow]  = GLFW_KEY_DOWN;
    io.KeyMap[ImGuiKey_Home]       = GLFW_KEY_HOME;
    io.KeyMap[ImGuiKey_End]        = GLFW_KEY_END;
    io.KeyMap[ImGuiKey_Delete]     = GLFW_KEY_DELETE;
    io.KeyMap[ImGuiKey_Backspace]  = 259;
    io.KeyMap[ImGuiKey_Enter]      = GLFW_KEY_ENTER;
    io.KeyMap[ImGuiKey_Escape]     = GLFW_KEY_ESCAPE;
    io.KeyMap[ImGuiKey_A]          = GLFW_KEY_A;
    io.KeyMap[ImGuiKey_C]          = GLFW_KEY_C;
    io.KeyMap[ImGuiKey_V]          = GLFW_KEY_V;
    io.KeyMap[ImGuiKey_X]          = GLFW_KEY_X;
    io.KeyMap[ImGuiKey_Y]          = GLFW_KEY_Y;
    io.KeyMap[ImGuiKey_Z]          = GLFW_KEY_Z;
    io.RenderDrawListsFn           = &igImplGlfwGL3_RenderDrawLists;
    io.SetClipboardTextFn          = &igImplGlfwGL3_SetClipboardText;
    io.GetClipboardTextFn          = &igImplGlfwGL3_GetClipboardText;
/+#ifdef _MSC_VER
    io.ImeWindowHandle = glfwGetWin32Window(g_Window);
#endif+/

    if (install_callbacks)
    {
        glfwSetMouseButtonCallback(window, &igImplGlfwGL3_MouseButtonCallback);
        glfwSetScrollCallback(window, &igImplGlfwGL3_ScrollCallback);
        glfwSetKeyCallback(window, &igImplGlfwGL3_KeyCallback);
        glfwSetCharCallback(window, &igImplGlfwGL3_CharCallback);
    }
}

void igImplGlfwGL3_CreateDeviceObjects()
{
	const GLchar *vertex_shader = q{
    #version 330
    uniform mat4 ProjMtx;
    in vec2 Position;
    in vec2 UV;
    in vec4 Color;
    out vec2 Frag_UV;
    out vec4 Frag_Color;
    void main()
    {
    	Frag_UV = UV;
    	Frag_Color = Color;
    	gl_Position = ProjMtx * vec4(Position.xy,0,1);
    }};
	const GLchar* fragment_shader = q{
		#version 330
    uniform sampler2D Texture;
    in vec2 Frag_UV;
    in vec4 Frag_Color;
    out vec4 Out_Color;
    void main()
    {
    	Out_Color = Frag_Color * texture( Texture, Frag_UV.st);
    }};

	g_ShaderHandle = glCreateProgram();
	g_VertHandle = glCreateShader(GL_VERTEX_SHADER);
	g_FragHandle = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(g_VertHandle, 1, &vertex_shader, null);
	glShaderSource(g_FragHandle, 1, &fragment_shader, null);
	glCompileShader(g_VertHandle);
	glCompileShader(g_FragHandle);
	glAttachShader(g_ShaderHandle, g_VertHandle);
	glAttachShader(g_ShaderHandle, g_FragHandle);
	glLinkProgram(g_ShaderHandle);

	g_AttribLocationTex = glGetUniformLocation(g_ShaderHandle, "Texture");
	g_AttribLocationProjMtx = glGetUniformLocation(g_ShaderHandle, "ProjMtx");
	g_AttribLocationPosition = glGetAttribLocation(g_ShaderHandle, "Position");
	g_AttribLocationUV = glGetAttribLocation(g_ShaderHandle, "UV");
	g_AttribLocationColor = glGetAttribLocation(g_ShaderHandle, "Color");

    glGenBuffers(1, &g_VboHandle);
    glGenBuffers(1, &g_ElementsHandle);

    glGenVertexArrays(1, &g_VaoHandle);
    glBindVertexArray(g_VaoHandle);
    glBindBuffer(GL_ARRAY_BUFFER, g_VboHandle);
    glEnableVertexAttribArray(g_AttribLocationPosition);
    glEnableVertexAttribArray(g_AttribLocationUV);
    glEnableVertexAttribArray(g_AttribLocationColor);

	glVertexAttribPointer(g_AttribLocationPosition, 2, GL_FLOAT, GL_FALSE,
    ImDrawVert.sizeof, cast(void*)0); glVertexAttribPointer(g_AttribLocationUV,
    2, GL_FLOAT, GL_FALSE, ImDrawVert.sizeof,
    cast(void*)ImDrawVert.uv.offsetof);
    glVertexAttribPointer(g_AttribLocationColor, 4, GL_UNSIGNED_BYTE, GL_TRUE,
    ImDrawVert.sizeof, cast(void*)ImDrawVert.col.offsetof);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	igImplGlfwGL3_CreateFontsTexture();
}

extern(C) nothrow const(char)* igImplGlfwGL3_GetClipboardText()
{
    return glfwGetClipboardString(g_window);
}

extern(C) nothrow void igImplGlfwGL3_SetClipboardText(const(char)* text)
{
    glfwSetClipboardString(g_window, text);
}

extern(C) nothrow void igImplGlfwGL3_MouseButtonCallback(GLFWwindow*, int
button, int action, int /*mods*/) {
    if (action == GLFW_PRESS && button >= 0 && button < 3)
        g_MousePressed[button] = true;
}

extern(C) nothrow void igImplGlfwGL3_ScrollCallback(GLFWwindow*, double
/*xoffset*/, double yoffset) {
    g_MouseWheel += cast(float)yoffset;
}

extern(C) nothrow void igImplGlfwGL3_KeyCallback(GLFWwindow*, int key, int, int
action, int mods)
{
    auto io = igGetIO();
    if (action == GLFW_PRESS)
        io.KeysDown[key] = true;
    if (action == GLFW_RELEASE)
        io.KeysDown[key] = false;
    io.KeyCtrl = (mods & GLFW_MOD_CONTROL) != 0;
    io.KeyShift = (mods & GLFW_MOD_SHIFT) != 0;
    io.KeyAlt = (mods & GLFW_MOD_ALT) != 0;
}

extern(C) nothrow void igImplGlfwGL3_CharCallback(GLFWwindow*, uint c)
{
    if (c > 0 && c < 0x10000)
    {
        ImGuiIO_AddInputCharacter(cast(ushort)c);
    }
}

void igImplGlfwGL3_CreateFontsTexture()
{
	ImGuiIO* io = igGetIO();
	ubyte* pixels;
	int width, height;
	ImFontAtlas_GetTexDataAsRGBA32(io.Fonts,&pixels,&width,&height,null);
	glGenTextures(1, &g_FontTexture);
	glBindTexture(GL_TEXTURE_2D, g_FontTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA,
	GL_UNSIGNED_BYTE, pixels);

	ImFontAtlas_SetTexID(io.Fonts, cast(void*)g_FontTexture);
}

void igImplGlfwGL3_Shutdown()
{
    if (g_VaoHandle) glDeleteVertexArrays(1, &g_VaoHandle);
    if (g_VboHandle) glDeleteBuffers(1, &g_VboHandle);
    if (g_ElementsHandle) glDeleteBuffers(1, &g_ElementsHandle);
    g_VaoHandle = 0;
    g_VboHandle = 0;
    g_ElementsHandle = 0;
    glDetachShader(g_ShaderHandle, g_VertHandle);
    glDeleteShader(g_VertHandle);
    g_VertHandle = 0;
    glDetachShader(g_ShaderHandle, g_FragHandle);
    glDeleteShader(g_FragHandle);
    g_FragHandle = 0;
    glDeleteProgram(g_ShaderHandle);
    g_ShaderHandle = 0;

	if (g_FontTexture)
	{
		glDeleteTextures(1, &g_FontTexture);
        ImFontAtlas_SetTexID(igGetIO().Fonts, cast(void*)0);
		g_FontTexture = 0;
	}

	igShutdown();
}

void igImplGlfwGL3_NewFrame()
{
	if (!g_FontTexture)
		igImplGlfwGL3_CreateDeviceObjects();

	auto io = igGetIO();

	int w, h;
	int display_w, display_h;
	glfwGetWindowSize(g_window, &w, &h);
	glfwGetFramebufferSize(g_window, &display_w, &display_h);
	io.DisplaySize = ImVec2(cast(float)display_w, cast(float)display_h);

    double current_time =  glfwGetTime();
    io.DeltaTime = g_Time > 0.0 ? cast(float)(current_time - g_Time) :
    cast(float)(1.0f/60.0f);
    g_Time = current_time;

    if (glfwGetWindowAttrib(g_window, GLFW_FOCUSED))
    {
        double mouse_x, mouse_y;
        glfwGetCursorPos(g_window, &mouse_x, &mouse_y);
        mouse_x *= cast(float)display_w / w;
        mouse_y *= cast(float)display_h / h;
        io.MousePos = ImVec2(mouse_x, mouse_y);
    }
    else
    {
        io.MousePos = ImVec2(-1,-1);
    }

    for (int i = 0; i < 3; i++)
    {
        io.MouseDown[i] = g_MousePressed[i] || glfwGetMouseButton(g_window, i)
        != 0;
        g_MousePressed[i] = false;
    }
    io.MouseWheel = g_MouseWheel;
    g_MouseWheel = 0.0f;
    glfwSetInputMode(g_window, GLFW_CURSOR, io.MouseDrawCursor ?
    GLFW_CURSOR_HIDDEN : GLFW_CURSOR_NORMAL);
	igNewFrame();
}
