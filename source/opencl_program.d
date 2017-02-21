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

  typedef struct T_Triangle {
    float3 vertices[3];
  } Triangle;

  // bool Intersection(Triangle* tri, float3 ray_o, float3 ray_d,
  //                   float* result)  {
  //   float3 edge1 = (*tri).vertices[1] - (*tri).vertices[0],
  //          edge2 = (*tri).vertices[2] - (*tri).vertices[1];
  //   // determinant calculation for U
  //   float3 P = cross(ray_d, edge2);
  //   // if determinant is zero, ray lies in triangle
  //   float det = dot(edge1, P);
  //   // CULLING
  //   //   if ( det < float.epsilon ) return 0;
  //   // NON CULLING
  //   if ( abs(det) < float.epsilon ) return 0;
  //   auto inv_det = 1.0f/det;
  //   // calculate distance from first vertex to ray origin
  //   auto T = ray.ori - vertices[0];
  //   // calculates u parameter and tests bound
  //   u = dot(T, P)*inv_det;
  //   if ( u < 0.0f || u > 1.0f ) return 0;
  //   // prepare to test v parameter
  //   auto Q = cross(T, edge1);
  //   // calculate V parameter and test bound
  //   v = dot(ray.dir, 0) * inv_det;
  //   // intersection outside of triangle
  //   if ( v < 0.0f || u + v > 1.0f ) return 0;
  //   // calculate t, scale parameters, ray intersection
  //   *result = dot(edge2, O) * inv_det;
  //   return 1;
  // }

  __kernel void Kernel_Raycast(
          __write_only image2d_t output_image,
          __global float* timer_arr,
          __global Triangle* vertex_data,
          __global uint*  vertex_length
          ) {
    float circle_size = 1.0f/(3.0f*pow(1.9f, ( float )(Max_iters)));
    float timer = timer_arr[0];
    int2 out = (int2)(get_global_id(0), get_global_id(1));
    float4 colour;

    float2 uv = (float2)(0.7f - out.x/512.0f, 0.5f - out.y/512.0f);

    if ( get_global_id(0) == 0 && get_global_id(1) == 0 ) {
      for ( int i = 0; i != *vertex_length; ++ i ) {
        printf("%f %f %f\n", vertex_data[i].vertices[0][0],
                             vertex_data[i].vertices[0][1],
                             vertex_data[i].vertices[0][1]);
      }
      printf("-----\n");
    }

    uv = rot(uv, timer);
    uv *= sin(timer*0.25f)*-5.5f + 1.5f;

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
