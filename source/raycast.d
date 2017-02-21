module raycast;
import opencl;
import opencl_program : Test_raycast_string;
import globals;

CLImage Test_Raycast ( ) {
  writeln("RAYCAST TEST");
  immutable(int) Img_dim = 1024;
  auto program = Compile(Test_raycast_string);
  program.Set_Kernel("Kernel_Raycast");
  auto result = program.Set_Image_Buffer(BufferType.write_only, Img_dim, 1);
  auto pA = [20.0f,  0.3f,  20.0f],
       pB = [100.0f, 0.3f,  20.0f],
       pC = [ 20.0f, 0.3f, 100.0f];
  auto buffer = program.Set_Buffer!float(BufferType.read_only, pA~pB~pC, 0);
  program.Run([Img_dim, Img_dim, 1], [1, 1, 1]);
  program.Read_Image(result);
  program.Cleanup();

  return CLImage(result.data, Img_dim, Img_dim);
}
