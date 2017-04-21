import std.stdio;
import globals;

void Init ( ) {
  AOD.CV.Load_Config;
  AOD.initialize(32, "raycast renderer", 512, 512);
  AOD.Set_Camera_Size(Vector(AOD.R_Window_Width, AOD.R_Window_Height));
  AOD.Set_BG_Colour(0.0, 0.0, 0.0);
}

import opencl : CLImage;

void main() {
  scope ( exit ) {
    writeln("Successfully ended");
  }
  import cloctree;
  Init();
  // import opencl : Initialize;
  // import raycast;
  // Initialize();

  // auto raycaster = new Raycaster();
  // AOD.Add(raycaster);

  // import softwarerenderer.raytracer;
  // AOD.Add(new Raytracer());
  {
    static import aodheme;
    aodheme.Initialize();
  }

  import structure.grid;

  Grid grid = new Grid ( 200, 16 );

  grid.Insert(new Voxel(gln.vec3(0.0f), gln.vec3(0.0f)));
  grid.Insert(new Voxel(gln.vec3(1.0f, 0.0f, 0.0f), gln.vec3(0.0f)));
  grid.Insert(new Voxel(gln.vec3(1.0f, 1.0f, 1.0f), gln.vec3(0.0f)));

  grid.writeln;

  AOD.Run();
  // raycaster.Clean_Up();
  import std.c.stdlib;
  exit(0);
}
