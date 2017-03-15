import std.stdio;
import globals;

void Init ( ) {
  AOD.CV.Load_Config;
  AOD.initialize(16, "raycast renderer", 512, 512);
  AOD.Set_Camera_Size(Vector(AOD.R_Window_Width, AOD.R_Window_Height));
  AOD.Set_BG_Colour(0.0, 0.0, 0.0);
}

void Game_Init ( ) {
}

import opencl : CLImage;

void main() {
  scope ( exit ) {
    writeln("Successfully ended");
  }
  // Init();
  // Game_Init();
  // import opencl : Initialize;
  // import raycast;
  // Initialize();

  // auto raycaster = new Raycaster();
  // AOD.Add(raycaster);

  // AOD.Run();
  // raycaster.Clean_Up();
  // import std.c.stdlib;
  // exit(0);

  return;
}
