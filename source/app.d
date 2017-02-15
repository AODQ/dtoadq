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

void main() {
    scope ( exit ) {
      writeln("Successfully ended");
    }
    Init();
    Game_Init();
    AOD.Run();
    return;
}
