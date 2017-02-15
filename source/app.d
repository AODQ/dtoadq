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
  Print_Device_Information;




  // import opencl;
  // import functional;
  // auto someKernelDef = CLKernelDef!("someKernel",
  //         CLBuffer!(float), "input",
  //         float, "b")
  // (q{
  //     input[gLinId] = exp(cos(input[gLinId] * b));
  // });

  // auto platform = getChosenPlatform();

  // auto devices = platform.getDevices(cl.DEVICE_TYPE_GPU);
  // auto context = createContext(devices);
  // auto queue = context.createCommandQueue(devices[0]);
  // auto kernel = context.createProgram(someKernelDef)
  //     .buildProgram
  //     .createKernel!(someKernelDef.name); // 

  // auto input = iota(90).map!(to!float).array;

  // auto buff = context
  //     .newBuffer(cl.MEM_COPY_HOST_PTR, input);

  // kernel.setArgs(buff, 1.4f);
  // queue.enqueueCLKernel(kernel, [10UL, 9UL]);

  // auto output = new float[](90);
  // queue.read(buff, output);

  // writeln("OUTPUT: ", output);
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
