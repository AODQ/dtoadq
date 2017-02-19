import std.stdio;
import globals;

void Init ( ) {
  AOD.CV.Load_Config;
  AOD.initialize(16, "raycast renderer", 1024, 768);
  AOD.Camera.Set_Size(Vector(AOD.R_Window_Width, AOD.R_Window_Height));
  AOD.Set_BG_Colour(0.0, 0.0, 0.0);
}

void Game_Init ( ) {
}

void OpenCL_Test ( ) {
  import opencl;

  Initialize;
  AOD.Add(new RenderMe(Test_OpenCL));
}

AOD.SheetContainer CLImage_To_Image(CLImage image) {
  import derelict.opengl3.gl3;
  GLuint texture;
  glGetError();
  glGenTextures(1, &texture);
  writeln("error: ", glGetError());
  glBindTexture(GL_TEXTURE_2D, texture);
  writeln("error: ", glGetError());
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  writeln("error: %x".format(glGetError()));
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  writeln("error: ", glGetError());
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width, image.height, 0,
      GL_RGBA, GL_UNSIGNED_BYTE, cast(void*)image.buffer.ptr
  );
  writeln("error: ", glGetError());
  glBindTexture(GL_TEXTURE_2D, 0);
  writeln("error: ", glGetError());
  writeln("TEXTURE: ", texture);
  writeln(image.buffer);
  return AOD.SheetContainer(texture, image.width, image.height);
}

import opencl : CLImage;;
class RenderMe : AOD.Entity {
  void Set_Image(CLImage image) {
  }
  this ( CLImage image ) {
    super();
    Set_Position(AOD.R_Window_Width/2, AOD.R_Window_Height/2);
    Set_Sprite(CLImage_To_Image(image));
  }
}

void main() {
  scope ( exit ) {
    writeln("Successfully ended");
  }
  Init();
  Game_Init();

  writeln("--------------------");
  writeln("opencl wrap test");
  OpenCL_Test;
  writeln("--------------------");

  AOD.Run();
  return;
}
