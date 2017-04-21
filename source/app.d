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

  import softwarerenderer.raytracer;
  AOD.Add(new Raytracer(q{
    (set ~length (lambda (array) {
      (set ~lhelp (lambda (array t) {
        (writeln (empty array))
        (if (empty array)
          (t)
          (lhelp (cdr array) (+ t (car array))))
      }))
      (lhelp array 0.0f)
    }))

    (set ~Map (lambda (point) {
      (- (length point) 4.0f)
    }))
  }));
  {
    static import aodheme;
    aodheme.Initialize();
  }

  AOD.Run();
  // raycaster.Clean_Up();
  import std.c.stdlib;
  exit(0);
}
