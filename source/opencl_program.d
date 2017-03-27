module opencl_program; immutable(string) Test_raycast_string = q{
    // --------------- MATERIAL/VOXEL/RAY --------------------------------------
    typedef struct T_Material {
      float3 colour;
      float metallic, subsurface, specular, roughness, specular_tint,
            anisotropic, sheen, sheen_tint, clearcoat, clearcoat_gloss,
            emission;
    } Material;

    typedef struct T_Ray {
      float3 origin, dir, invdir;
      int3 sign;
    } Ray;

    typedef struct T_Voxel {
      float3 position;
      float3 colour;
    } Voxel;

    Ray New_Ray ( float3 o, float3 d ) {
      Ray ray;
      float3 id = (float3)(1.0f/d.x, 1.0f/d.y, 1.0f/d.z);
      ray.origin = o;
      ray.dir    = d;
      ray.invdir = id;
      ray.sign   = (int3)(id.x < 0, id.y < 0, id.z < 0);
      return ray;
    }

    // --------------- OCTREE --------------------------------------------------
    typedef struct T_OctreeNode {
      int8 child_id;
      float3 origin;
      float3 half_size;
      int voxel_id;
    } OctreeNode;

    typedef struct T_OctreeData {
      __global OctreeNode* node_pool;
      __global Voxel*      voxel_pool;
      int node_pool_size, voxel_pool_size;
    } OctreeData;

    unsigned char ROctant_Mask ( __global OctreeNode* node, float3 point ) {
      unsigned char oct = 0;
      if ( node->origin.x < point.x ) oct |= 4;
      if ( node->origin.y < point.y ) oct |= 2;
      if ( node->origin.z < point.z ) oct |= 1;
      return oct;
    }

    bool Is_Leaf ( __global OctreeNode* node ){ return node->child_id[0] == -1;}
    bool Is_Empty ( __global OctreeNode* node ) {
      return Is_Leaf(node) && node->voxel_id == -1;
    }

    float3 RLowerBound ( __global OctreeNode* node ) {
      return (float3)(node->origin[0] - node->half_size[0],
                      node->origin[1] - node->half_size[1],
                      node->origin[2] - node->half_size[2]);
    }

    float3 RHigherBound ( __global OctreeNode* node ) {
      return (float3)(node->origin[0] + node->half_size[0],
                      node->origin[1] + node->half_size[1],
                      node->origin[2] + node->half_size[2]);
    }

    // --------------- RANDOM --------------------------------------------------
    // --- random generation via xorshift1024star
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

    float Uniform(__global RNG* rng, const float min, const float max) {
      return min + ((float)RNG_Next(rng) / (float)(ULONG_MAX/(max-min)));
    }

    float3 Uniform_Float3(__global RNG* rng, const float min, const float max) {
      return (float3)(
        Uniform(rng, min, max),
        Uniform(rng, min, max),
        Uniform(rng, min, max)
      );
    }

    // https://pathtracing.wordpress.com/2011/03/03/cosine-weighted-hemisphere/
    float3 Random_Hemisphere_Direction(float3 normal, __global RNG* rng) {
      int abort = 0;
      float3 origin = (float3)(0.5f, 0.5f, 0.5f);
      while ( true ) {
        float3 v = normalize(Uniform(rng, 0.0f, 1.0f)) - origin;

        if ( dot(normal, v) > 0.0f ) return v;
        ++ abort;
        if ( abort > 500 )
          return v; // couldn't find a normal :/
      }
    }

    // --------------- INTERSECTION --------------------------------------------
    float Ray_Intersection(float3 bmin, float3 bmax, Ray ray) {
      float t1 = (bmin.x - ray.origin.x)*ray.invdir.x,
            t2 = (bmax.x - ray.origin.x)*ray.invdir.x;
      float tmin = min(t1, t2),
            tmax = max(t1, t2);
      for ( int i = 1; i < 3; ++ i ) {
        t1 = (bmin[i] - ray.origin[i])*ray.invdir[i];
        t2 = (bmax[i] - ray.origin[i])*ray.invdir[i];
        tmin = max(tmin, min(t1, t2));
        tmax = min(tmax, max(t1, t2));
      }

      if ( tmax < max(tmin, 0.0f) ) return -1.0f;
      return tmin;
    }

    int RVoxel_Intersection ( OctreeData* data, Ray ray, float* dist ){
      __global OctreeNode* curr_node = &data->node_pool[0];
      float3 min, max;
      for ( int depth = 0; depth != 5; ++ depth ) {
        if ( Is_Leaf(curr_node) ) {
          if ( curr_node->voxel_id != -1 ) {
            min = RLowerBound(curr_node), max = RHigherBound(curr_node);
            *dist = Ray_Intersection(min, max, ray);
            if ( (*dist) > 0.0f ) {
              return curr_node->voxel_id;
            }
          }
        } else {
          float node_dist = FLT_MAX;
          unsigned char  node_mask = 50;
          for ( unsigned char i = 0; i != 8; ++ i ) {
            // get child info
            int child_id = curr_node->child_id[i];
            if ( child_id == -1 )
              continue;
            __global OctreeNode* child_node = &data->node_pool[child_id];
            if ( Is_Empty(child_node) ) continue;

            // get child bounding box and calculate its ray intersection
            min = RLowerBound(child_node), max = RHigherBound(child_node);
            float results = Ray_Intersection(min, max, ray);
            if ( results > 0.0f && results < node_dist ) {
              node_dist = results;
              node_mask = i;
            }
          }
          if ( node_mask == 50 ) {
            return -1;
          }
          curr_node = &data->node_pool[curr_node->child_id[node_mask]];
        }
      }
      return -1;
    }

    typedef struct T_IntersectionInfo {
      bool intersection;
      float3 position;
      float distance;
      Voxel voxel;
      float3 normal, angle;
      float3 colour;
    } IntersectionInfo;

    IntersectionInfo Raycast_Scene(
                OctreeData* data,
                Ray ray) {
      IntersectionInfo info;
      info.intersection = false;

      // --- find closest voxel ---
      float dist;
      int voxel_id = RVoxel_Intersection(data, ray, &dist);

      if ( voxel_id >= 0 ) {
        info.intersection = true;
        info.distance = dist;
        info.position = ray.origin + ray.dir*dist;
        info.colour = data->voxel_pool[voxel_id].colour;
        // info.material = material_data[vertex_data[i].material_index];
        // info.normal = cross(vertex_data[i].B - vertex_data[i].A,
        //                     vertex_data[i].C - vertex_data[i].A);
        // info.angle = dot(ray.d, info.normal);
      }

      return info;
    }

    // --------------- CAMERA --------------------------------------------------
    typedef struct T_Camera {
      float3 position, lookat, up;
      int2 dimensions;
      float fov;
    } Camera;

    Ray Camera_Ray(__global RNG* rng, __global Camera* camera) {
      float2 dim   = (float2)(camera->dimensions.x/2.0f,
                              camera->dimensions.y/2.0f);

      float3 z_axis = normalize(camera->lookat - camera->position);
      float3 x_axis = cross(z_axis, camera->up);
      float3 y_axis = cross(x_axis, z_axis);

      float fov_rad = camera->fov * 3.14159f/180.0f;
      float focus_dist = 1.0f/tan(fov_rad/2.0f);

      float aspect_ratio = dim.x/dim.y;

      float2 pixel = (float2)((float)get_global_id(0), (float)get_global_id(1));

      float x_proj = 2.0f * aspect_ratio * (pixel.x/dim.x - 0.5f),
            y_proj = 2.0f * (0.5f - pixel.y/dim.y);
      float3 direction = x_axis*x_proj + y_axis*y_proj + z_axis*focus_dist;
      return New_Ray(camera->position, direction);
   }







     // viewDir*= float3(0.9,1.0,1.0);
     // Lpos = float3(0.85*cos(time*0.55),1.5, 0.85*sin(time*1.0));
     // // Compute camera properties
     // float  camDist = 10.0f;
     // float3 camUp = float3(0,1.0f,0);
     // float3 camTarget = float3(0,3.0,0);
     // // And from them evaluted ray direction in world space
     // float3 forward = normalize(camTarget - camPos);
     // float3 left = normalize(cross(forward, camUp));
     // float3 up = cross(left, forward);
     // float3 worldDir = viewDir.x*left + viewDir.y*up + viewDir.z*forward;
    // }
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
            __write_only image2d_t output_image,
            __read_only  image2d_t input_image,
            __global OctreeNode* node_pool,  __global uint* node_length,
            __global Voxel*      voxel_pool, __global uint* voxel_length,
            __global RNG* rng,
            // __read_only  image2d_t environment_map,
            __global Camera* camera
            ) {
      int2 out = (int2)(get_global_id(0), get_global_id(1));
      bool isprint  = out.x%10 == 0 && out.y%2 == 0,
           isprintg = out.x == 0 && out.y == 0;

      Ray ray = Camera_Ray(rng, camera);
      if ( isprintg ) {
        printf("camera: position <%f, %f, %f> lookat <%f, %f, %f>, "
               "up <%f, %f, %f>, fov %f\n",
               camera->position.x, camera->position.y, camera->position.z,
               camera->lookat.x, camera->lookat.y, camera->lookat.z,
               camera->up.x, camera->up.y, camera->up.z,
               camera->fov
               );
        printf("Ray: origin <%f, %f, %f> ; dir <%f, %f, %f> ; idir <%f, %f, %f>"
                 ", sign <%d, %d, %d>\n",
                 ray.origin.x, ray.origin.y, ray.origin.z,
                 ray.dir.x, ray.dir.y, ray.dir.z,
                 ray.invdir.x, ray.invdir.y, ray.invdir.z,
                 ray.sign.x, ray.sign.y, ray.sign.z);
      }

      float3 colour = (float3)(0.0f, 0.0f, 0.0f);
      float3 weight = (float3)(1.0f, 1.0f, 1.0f);
      bool hit = false;

      OctreeData data;
      data.node_pool  = node_pool;
      data.voxel_pool = voxel_pool;
      data.node_pool_size = *node_length;
      data.voxel_pool_size = *voxel_length;

      int depth = 0;
      while ( true ) {
        float depth_ratio = 1.0f - (1.0f / (DEPTH - depth)); // 
        if ( Uniform(rng, 0.0f, 1.0f) >= depth_ratio ) break; // 
        IntersectionInfo info = Raycast_Scene(&data, ray);
        if ( !info.intersection ) {
          colour = (float3)(0.0f, 0.0f, 0.0f);
          hit = true;
          break;
        }

        colour = info.colour;
        hit = true;
        // hit = true;
        // if ( info.material.emission > FLT_EPSILON ) {
        //   float emittance = (info.material.emission/depth_ratio);
        //   colour = weight * info.material.base_colour * emittance;
        //   hit = true;
        //   break;
        // }

        // float3 brdf = DisneyBRDF(info.position, info.angle, info.normal,
        //                        (float3)(0.0f), (float3)(0.0f));
        // weight *= info.material.base_colour;
        // info.normal = info.normal * -sign(info.angle);
        // float3 new_dir = Random_Hemisphere_Direction(info.normal, rng);
        // ray.o = info.position;
        // ray.d = new_dir;
        // colour = weight * info.material.base_colour + info.normal;
        // hit = true;
        break;
      }

      if ( hit ) {
        float3 old_colour = read_imagef(input_image, out).xyz;
        float3 rcolour = mix(colour, old_colour, 0.2f);
        rcolour = fmax(0.0f, fmin(1.0f, rcolour));
        // rcolour = pow(rcolour, (float3)(1.0f / 0.5f));
        write_imagef(output_image, out, (float4)(rcolour, 1.0));
      }
    }
};
