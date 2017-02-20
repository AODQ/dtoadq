module opencl_program; auto program_string = q{
      __constant sampler_t sampler =
          CLK_NORMALIZED_COORDS_FALSE |
          CLK_ADDRESS_CLAMP_TO_EDGE   |
          CLK_FILTER_NEAREST;

      __kernel void hello(__read_only image2d_t input_image,
                          __write_only image2d_t output_image) {

      int2 out   = (int2)(get_global_id(0), 10);
      int t = get_global_id(0);

      float4 val = read_imagef(input_image, out);
      printf("Writing [%f %f %f %f] at: <%d, %d>\n",
        val.x, val.y, val.z, val.w, out.x, out.y);
      write_imagef(output_image, out, (float4)(1.0, 0.5, 0.5, 1.0));
      // if ( out.x < 16 && out.y < 16 ) {
      //   int weight = 0;
      //   float4 col = (float4)(0.0f, 0.0f, 0.0f, 1.0f);
      //   for ( int y = start.y; y <= end.y; ++ y )
      //     for ( int x = start.x; x <= end.x; ++ x ) {
  //col += read_imagef(input_image, (int2)(x, y))*(kernel_weights[weight]*8);
          //   weight += 1;
          // }
      // }
    }
  };
