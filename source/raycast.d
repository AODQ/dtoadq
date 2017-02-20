module raycast;
import opencl;
import opencl_program : Test_raycast_string;

CLImage Test_Raycast ( ) {
  immutable(int) Img_dim = 1024;
  auto program = Compile(Test_raycast_string);
  program.Set_Kernel("Kernel_Raycast");
  program.Set_Image_Buffer(BufferType.write_only, Img_dim, "result_image", 1);
  auto pA = [20.0f,  0.3f,  20.0f],
       pB = [100.0f, 0.3f,  20.0f],
       pC = [ 20.0f, 0.3f, 100.0f];
  program.Set_Buffer(BufferType.read_only, "vertices", pA~pB~pC, 0);
  program.Run([Img_dim, Img_dim, 1], [1, 1, 1]);
  auto image_result = program.Read_Buffer("result_image");
  program.Cleanup();

  return CLImage(image_result, 0, 0);
}
