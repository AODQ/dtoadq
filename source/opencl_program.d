module opencl_program; immutable(string) Test_raycast_string = q{

  __constant int Max_iters = 15;

  float2 rot ( float2 uv, float a ) {
    return (float2)(uv.x*cos(a) - uv.y*sin(a),
                    uv.y*cos(a) + uv.x*sin(a));
  }

  typedef struct T_Material {
    float3 ambient, diffuse, specular, emission, shininess;
  } Material;

  typedef struct T_Triangle {
    float3 A, B, C;
    uint material_index;
  } Triangle;

  typedef struct T_Intersection_Info {
    bool intersection;
    float3 position;
    float distance;
    Material material;
  } Intersection_Info;

  typedef struct T_Ray {
    float3 o, d;
  } Ray;

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
    // if ( det < FLT_EPSILON ) return 0;
    // NON CULLING
    if ( fabs(det) < FLT_EPSILON ) return 0.0f;
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

  Intersection_Info Raycast_Scene(
              __global Triangle* vertex_data,   const uint vertex_length,
              __global Material* material_data, const uint material_length,
              Ray ray) {
    Intersection_Info info;
    info.intersection = false;

    float closest_dist = FLT_MAX;
    for ( uint i = 0; i != vertex_length; ++ i ) {
      float result = Intersection(vertex_data[i], ray.o, ray.d);
      if ( result > FLT_EPSILON && closest_dist > result ) {
        info.distance = result;
        info.intersection = true;
        info.position = ray.o + ray.d*result;
        info.material = material_data[vertex_data[i].material_index];
      }
    }

    return info;
  }

  __kernel void Kernel_Raycast(
          __write_only image2d_t output_image,
          __global float* timer_arr,
          __global Triangle* vertex_data,   global uint* vertex_lengthp,
          __global Material* material_data, global uint* material_lengthp
          ) {
    float circle_size = 1.0f/(3.0f*pow(1.9f, ( float )(Max_iters)));
    float timer = timer_arr[0];
    uint vertex_length = *vertex_lengthp,
         material_length = *material_lengthp;
    int2 out = (int2)(get_global_id(0), get_global_id(1));
    bool isprint = out.x == 0 && out.y == 0;

    Ray ray;
    ray.o = (float3)(out.x, 0.0f, out.y);
    ray.d = normalize((float3)(0.0f, 1.0f, 0.0f));

    float4 colour = (float4)(0.0f, 0.0f, 0.0f, 1.0f);

    Intersection_Info info = Raycast_Scene(vertex_data, vertex_length,
                                           material_data, material_length, ray);

    if ( isprint ) {
      for ( uint i = 0; i != vertex_length; ++ i ) {
        printf("%f %f %f\n", vertex_data[i].A[0],
                             vertex_data[i].A[0],
                             vertex_data[i].A[0]);
        printf("%f %f %f\n", vertex_data[i].B[0],
                             vertex_data[i].B[0],
                             vertex_data[i].B[0]);
        printf("%f %f %f\n", vertex_data[i].C[0],
                             vertex_data[i].C[0],
                             vertex_data[i].C[0]);
        printf("-----\n");
      }
    }
    if ( info.intersection ) {
      colour.x = info.distance;
      if ( isprint )
        printf("%f \n", info.material.emission.x);
      // colour.xyz = info.material.emission;
    }

    write_imagef(output_image, out, colour);
  }
};
