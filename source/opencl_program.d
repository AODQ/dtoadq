module opencl_program; immutable(string) Test_raycast_string = q{

  __constant int Max_iters = 1;

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
    float3 colour;
  } Intersection_Info;

  typedef struct T_Ray {
    float3 o, d;
  } Ray;

  float Dot(float3 a, float3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
  }

  float Hash ( float2 vec ) {
    return fract( sin(vec.x+cos(vec.y*412.0f))*75138.1942 );
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
    info.colour = (float3)(1.0f, 1.0f, 1.0f);

    for ( uint global_it = 0; global_it != Max_iters; ++ global_it ) {
      // --- find closest triangle ---
      float closest_dist = FLT_MAX;
      for ( uint i = 0; i != vertex_length; ++ i ) {
        float result = Intersection(vertex_data[i], ray.o, ray.d);
        if ( result > FLT_EPSILON && closest_dist > result ) {
          info.distance = result;
          info.intersection = true;
          info.position = ray.o + ray.d*result;
          info.material = material_data[vertex_data[i].material_index];
          if ( info.material.emission.x >= 0.0f ) return info;
        }
      }

      // --- reflect off the triangle at a random point ---

    }

    return info;
  }

  __kernel void Kernel_Raycast(
          __write_only image2d_t output_image,
          __global float* timer_arr,
          __global Triangle* vertex_data,   global uint* vertex_lengthp,
          __global Material* material_data, global uint* material_lengthp
          ) {
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

    if ( info.intersection ) {
      // colour.xyz = info.material.emission;
      colour.x = info.distance;
    }

    write_imagef(output_image, out, colour);
  }
};
