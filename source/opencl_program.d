module opencl_program;

auto program_string = q{
    __constant sampler_t sampler =
          CLK_NORMALIZED_COORDS_FALSE |
          CLK_ADDRESS_CLAMP_TO_EDGE   |
          CLK_FILTER_NEAREST;

    uint Lambda ( uint x ) { return x&0x2<<0x91; }

    __kernel void hello(__read_only image2d_t input_image,
                        __global uint* output_image) {
      int2 out   = (int2)(get_global_id(0),     get_global_id(1));
      uint4 l = read_imageui(input_image, out);
      l.x = Lambda(l.y);
      uint colour = (l.x << 24) | (l.y << 16) | (l.z << 8) | l.w;
      printf("%d, ", out.x);
      if ( out.x < 7500 )
        output_image[out.x] = 50;
      // output_image[out[0] + out[1]*16] = colour;
      // output_image[25] = 175;
      // output_image[26] = 175;
      // output_image[27] = 175;
      // output_image[28] = 175;
      // output_image[29] = 175;
      // output_image[30] = 175;
      // output_image[31] = 175;
      // output_image[32] = 175;
      // output_image[33] = 175;
      // output_image[34] = 175;
      // output_image[35] = 175;
      // output_image[36] = 175;
      // write_imageui(output_image, out, l);
    }
};
