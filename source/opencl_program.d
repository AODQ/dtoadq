module opencl_program; immutable(string) Test_raycast_string = q{

  __constant int Max_iters = 15;

  float2 rot ( float2 uv, float a ) {
    return (float2)(uv.x*cos(a) - uv.y*sin(a),
                    uv.y*cos(a) + uv.x*sin(a));
  }

  float2 Abs(float2 a) {
    return (float2)(a.x < 0.0f ? -a.x : a.x,
                    a.y < 0.0f ? -a.y : a.y);
  }

  __kernel void Kernel_Raycast(
          __write_only image2d_t output_image,
          __global float* timer_arr,
          __global float* vertex_data,
          __global uint*  vertex_length
          ) {
    float circle_size = 1.0f/(3.0f*pow(1.9f, ( float )(Max_iters)));
    float timer = timer_arr[0];
    int2 out = (int2)(get_global_id(0), get_global_id(1));
    float4 colour;

    float2 uv = (float2)(0.7f - out.x/512.0f, 0.5f - out.y/512.0f);

    uv = rot(uv, timer);
    uv *= sin(timer*0.25)*-5.5f + 1.5f;

    float s = 0.3f;
    for ( int i = 0; i != Max_iters; ++ i ) {
      uv = Abs(uv) - s;
      uv = rot(uv, timer);
      s = s/2.1f;
    }

    float c = length(uv) > circle_size ? 0.0f : 1.0f;

    colour = (float4)(c, c, c, 1.0f);
    write_imagef(output_image, out, colour);
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
