module opencl_program; immutable(string) Test_raycast_string = q{


  __kernel void Kernel_Raycast(
          __global float* vertex_data,
          __write_only image2d_t output_image) {
  }


};

immutable(string) Test_program_string = q{
      __constant sampler_t sampler =
          CLK_NORMALIZED_COORDS_FALSE |
          CLK_ADDRESS_CLAMP_TO_EDGE   |
          CLK_FILTER_NEAREST;

      __kernel void hello(__read_only image2d_t input_image,
                          __write_only image2d_t output_image) {

      int2 out   = (int2)(get_global_id(0), get_global_id(1));
      float4 val = read_imagef(input_image, out);
      write_imagef(output_image, out, val);
    }
  };
