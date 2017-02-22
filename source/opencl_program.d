module opencl_program; immutable(string) Test_raycast_string = q{

  __constant int Max_iters = 15;
  __constant float Epsilon = 0.00001f;

  float2 rot ( float2 uv, float a ) {
    return (float2)(uv.x*cos(a) - uv.y*sin(a),
                    uv.y*cos(a) + uv.x*sin(a));
  }

  // float2 Abs(float2 a) {
  //   return (float2)(a.x < 0.0f ? -a.x : a.x,
  //                   a.y < 0.0f ? -a.y : a.y);
  // }

  typedef struct T_Triangle {
    float3 A, B, C;
  } Triangle;

  float Dot(float3 a, float3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
  }

  // mollor trombore for now, will use a quad moller-trumbore
  // in the future
  float Intersection(Triangle tri, float3 ray_o, float3 ray_d)  {
    float3 edge1 = tri.B - tri.A,
           edge2 = tri.C - tri.A;
    // determinant calculation for U
    float3 P = cross(ray_d, edge2);
    // if determinant is zero, ray lies in triangle
    float det = dot(edge1, P);
    // CULLING
    if ( det < Epsilon ) return 0;
    // NON CULLING
    // if ( fabs(det) < Epsilon ) return 0.0f;
    float inv_det = 1.0f/det;
    // calculate distance from first vertex to ray origin
    float3 T = ray_o - tri.A;
    // calculates u parameter and tests bound
    float u = dot(T, P)*inv_det;
    if ( u < 0.0f || u > 1.0f ) return 0.0f;
    // prepare to test v parameter
    float3 Q = cross(T, edge1);
    // calculate V parameter and test bound
    float v = dot(ray_d, Q) * inv_det;
    // intersection outside of triangle
    if ( v < 0.0f || u + v > 1.0f ) return 0.0f;
    // calculate t, scale parameters, ray intersection
    return dot(edge2, Q)*inv_det;
  }

  __kernel void Kernel_Raycast(
          __write_only image2d_t output_image,
          __global float* timer_arr,
          __global Triangle* vertex_data,
          __global uint*  vertex_length
          ) {
    float circle_size = 1.0f/(3.0f*pow(1.9f, ( float )(Max_iters)));
    float timer = timer_arr[0];
    int2 out = (int2)(get_global_id(0), get_global_id(1));
    bool isprint = out.x == 0 && out.y == 0;

    float3 ray_o = (float3)(out.x+timer*50.0f,   0.0f, out.y),
           ray_d = normalize((float3)(0.0f,    1.0f,  0.0f));

    float4 colour = (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    for ( uint i = 0; i != *vertex_length; ++ i ) {
      float result = Intersection(vertex_data[i], ray_o, ray_d);
      if ( result >= Epsilon ) {
        float3 v = ray_o + result*ray_d;
        float3 n = normalize(v);
        colour = (float4)(n.x, n.y, n.z, 1.0f);
        out = (int2)((int)round(v.x), (int)round( v.z ));
        break;
      }
    }

    write_imagef(output_image, out, colour);
  }
};
