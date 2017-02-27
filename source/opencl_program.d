module opencl_program; immutable(string) Test_raycast_string = q{

  float2 rot ( float2 uv, float a ) {
    return (float2)(uv.x*cos(a) - uv.y*sin(a),
                    uv.y*cos(a) + uv.x*sin(a));
  }


  // --- random generation via xorshift1024star
  // it is under the public domain

  typedef struct RNG {
    ulong seed[16];
    ulong p;
  } RNG;

  ulong RNG_Next(__global RNG* rng) {
    const ulong s0 = (*rng).seed[(*rng).p];
    ulong       s1 = (*rng).seed[(*rng).p = ((*rng).p + 1)&15];
    s1 ^= s1 << 31; // a
    (*rng).seed[(*rng).p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b, c
    return (*rng).seed[(*rng).p] * 1181783497276652981L *
                (get_global_id(0) + 250) * (get_global_id(1) + 250);
  }

  typedef struct T_Material {
    float3 base_colour;
    float metallic, subsurface, specular, roughness, specular_tint,
          anisotropic, sheen, sheen_tint, clearcoat, clearcoat_gloss,
          emission;
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
    float3 colour, normal, angle;
  } Intersection_Info;

  typedef struct T_Ray {
    float3 o, d;
  } Ray;

  float HashFloat(__global RNG* rng, const float min, const float max) {
    return min + ((float)RNG_Next(rng) / (float)(ULONG_MAX/(max-min)));
  }

  float3 HashVec3(__global RNG* rng) {
    return (float3)(
      HashFloat(rng, 0.0f, 1.0f),
      HashFloat(rng, 0.0f, 1.0f),
      HashFloat(rng, 0.0f, 1.0f)
    );
  }

  float3 Random_Hemisphere_Direction(float3 normal, __global RNG* rng) {
    // https://pathtracing.wordpress.com/2011/03/03/cosine-weighted-hemisphere/
    int abort = 0;
    while ( true ) {
      float3 v = normalize(HashVec3(rng) - (float3)(0.5f, 0.5f, 0.5f));

      if ( dot(normal, v) > 0.0f ) return v;
      ++ abort;
      if ( abort > 500 )
        return v; // couldn't find a normal :/
    }
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
    if ( det < FLT_EPSILON ) return 0;
    // NON CULLING
    // if ( fabs(det) < FLT_EPSILON ) return 0.0f;
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

    // --- find closest triangle ---
    float closest_dist = FLT_MAX;
    for ( uint i = 0; i != vertex_length; ++ i ) {
      float result = Intersection(vertex_data[i], ray.o, ray.d);
      if ( result > FLT_EPSILON && closest_dist > result ) {
        info.distance = result;
        info.intersection = true;
        info.position = ray.o + ray.d*result;
        info.material = material_data[vertex_data[i].material_index];
        info.normal = cross(vertex_data[i].B - vertex_data[i].A,
                            vertex_data[i].C - vertex_data[i].A);
        info.angle = dot(ray.d, info.normal);
      }
    }

    return info;
  }

  // --- camera shit ---

  Ray Camera_Ray(__global RNG* rng, uint x, uint y) {
    float3 camera_pos = (float3)(0.0f, -2.0f, 0.0f);
    float3 pixel_pos  = (float3)(x, 0.0f, y);
    // lol antialiasing
    float range = 0.01f;
    float2 antialias  = (float2)(
      HashFloat(rng, -range, range),
      HashFloat(rng, -range, range));
    camera_pos += (float3)(antialias, 0.0f);
    Ray result;
    result.o = camera_pos;
    result.d = normalize(pixel_pos - camera_pos);
    return result;
  }

  // --- BRDF ---

  float3 DisneyBRDF(float3 L, float3 V, float3 N, float3 X, float3 Y) {
    float3 R = dot(L, N)*N - L;
    float D = pow(max(0.0f, dot(R, N)), 1.02f);
    D *= (2.0f+1.02f) / (2.0f*3.14159f);
    return (float3)(D);
    // tangent = normalize(cross(float3(0, 1, 0), normal))
    // bitan   = normalize(cross(normal, tangent))
  }

  // --- kernel ---

  #define DEPTH 16

  __kernel void Kernel_Raycast(
          __read_only image2d_t environment_map,
          __write_only image2d_t output_image,
          __read_only  image2d_t input_image,
          __global float* timer_arr,
          __global Triangle* vertex_data,   global uint* vertex_lengthp,
          __global Material* material_data, global uint* material_lengthp,
          __global RNG* rng
          ) {
    float timer = timer_arr[0];
    uint vertex_length = *vertex_lengthp,
         material_length = *material_lengthp;
    int2 out = (int2)(get_global_id(0), get_global_id(1));
    bool isprint = out.x%10 == 0 && out.y%2 == 0;

    Ray ray = Camera_Ray(rng, out.x, out.y);

    float3 colour = (float3)(0.0f, 0.0f, 0.0f);
    float3 weight = (float3)(1.0f, 1.0f, 1.0f);
    bool hit = false;

    for ( int depth = 0; depth != DEPTH; ++ depth ) {
      float depth_ratio = 1.0f - (1.0f / (DEPTH - depth));
      Intersection_Info info = Raycast_Scene(vertex_data, vertex_length,
                                    material_data, material_length, ray);
      if ( !info.intersection ) {
        // colour = weight * info.material.base_colour;
        break;
      }

      if ( info.material.emission > FLT_EPSILON ) {
        float emittance = (info.material.emission/depth_ratio);
        // if ( isprint )
        //   printf("%f \n", emittance);
        colour = weight * info.material.base_colour * emittance;
        hit = true;
        break;
      }

      float3 brdf = DisneyBRDF(info.position, info.angle, info.normal,
                             (float3)(0.0f), (float3)(0.0f));
      if ( isprint ) printf("%f %f %f\n", brdf.x, brdf.y, brdf.z);
      weight *= info.material.base_colour*brdf;
      info.normal = info.normal * -sign(info.angle);
      float3 new_dir = Random_Hemisphere_Direction(info.normal, rng);
      ray.o = info.position;
      ray.d = new_dir;
    }

    if ( hit ) {
      float3 old_colour = read_imagef(input_image, out).xyz;
      float3 rcolour = mix(colour, old_colour, 0.2f);
      rcolour = fmax(0.0f, fmin(1.0f, rcolour));
      rcolour = pow(rcolour, (float3)(1.0f / 0.5f));
      write_imagef(output_image, out, (float4)(rcolour, 1.0));
    }
  }
};
