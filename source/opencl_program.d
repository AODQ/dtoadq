module opencl_program; immutable(string) Test_raycast_string = q{
// --------------- DEBUG -------------------------------------------------------
bool Is_Debug_Print ( ) {
return get_global_id(0) == 1 && get_global_id(1) == 1;
}
// --------------- MATRICES ----------------------------------------------------
typedef struct T_mat3 {
  float3 x, y, z;
} mat3;


float3 mat3_mul ( mat3 mat, float3 vec ) {
  return (float3)(
    mat.x.x*vec.x + mat.x.y*vec.y + mat.x.z*vec.z,
    mat.y.x*vec.x + mat.y.y*vec.y + mat.y.z*vec.z,
    mat.z.x*vec.x + mat.z.y*vec.y + mat.z.z*vec.z
  );
}


mat3 nmat3 ( float3 x_, float3 y_, float3 z_ ) {
  mat3 mat;
  mat.x = x_;
  mat.y = y_;
  mat.z = z_;
  return mat;
}

mat3 rotate_y ( float fi ) {
  return nmat3(
    (float3)(cos(fi), 0.0f, sin(fi)),
    (float3)(0.0f, 1.0f, 0.0f),
    (float3)(-sin(fi), 0.0f, cos(fi))
  );
}

mat3 rotate_x ( float fi ) {
  return nmat3(
    (float3)(1.0f, 0.0f, 0.0f),
    (float3)(0.0f, cos(fi), -sin(fi)),
    (float3)(0.0f, sin(fi),  cos(fi))
  );
}
// --------------- BASIC GRAPHIC FUNCTIONS -------------------------------------
float3 reflect ( float3 V, float3 N ) {
  return V - 2.0f*dot(V, N)*N;
}

float3 refract(float3 V, float3 N, float refraction) {
  float cosI = -dot(N, V);
  float cosT = 1.0f - refraction*refraction*(1.0f - cosI*cosI);
  return (refraction*V) + (refraction*cosI - sqrt(cosT))*N;
}

float length2(float2 p, int n) {
  return pow(pown(p.x, n) + pown(p.y, n), 1.0f/n);
}
// --------------- MATERIAL/VOXEL/RAY ------------------------------------------
// typedef struct T_Material {
//   float3 colour;
//   float metallic, subsurface, specular, roughness, specular_tint,
//         anisotropic, sheen, sheen_tint, clearcoat, clearcoat_gloss,
//         emission;
// } Material;

typedef struct T_Ray {
  float3 origin, dir, invdir;
  int3 sign;
} Ray;

Ray New_Ray ( float3 o, float3 d ) {
  Ray ray;
  ray.origin = o;
  ray.dir    = d;
  return ray;
}


// --------------- RANDOM ------------------------------------------------------
// // --- random generation via xorshift1024star
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


float2 fract2f ( float2 vec ) {
  float2 itptr;
  return fract(vec, &itptr);
}

float fract1f ( float vec ) {
  float itptr;
  return fract(vec, &itptr);
}


float rand ( float2 n ) {
  return fract1f(sin(dot(n, (float2)(19.9898, 4.1414))) * 43758.5453);
}

float noise ( float2 n ) {
  const float2 d = (float2)(0.0f, 1.0f);
  float2 b = floor(n),
         f = smoothstep((float2)(0.0f), (float2)(1.0f), fract2f(n));
  return mix(mix(rand(b), rand(b + d.yx), f.x),
             mix(rand(b + d.xy), rand(b + d.yy), f.x),
             f.y);
}


// // https://pathtracing.wordpress.com/2011/03/03/cosine-weighted-hemisphere/
// float3 Random_Hemisphere_Direction(float3 normal, __global RNG* rng) {
//   int abort = 0;
//   float3 origin = (float3)(0.5f, 0.5f, 0.5f);
//   while ( true ) {
//     float3 v = normalize(Uniform(rng, 0.0f, 1.0f)) - origin;

//     if ( dot(normal, v) > 0.0f ) return v;
//     ++ abort;
//     if ( abort > 500 )
//       return v; // couldn't find a normal :/
//   }
// }

float PrimitiveA ( float3 point ) {
  return length(max(fabs(point) - (float3)(0.9f), 0.0f));
}


float PrimitiveB ( float3 point ) {
  float2 q = (float2)(length(point.xy) - 1.5f, point.z);
  return length(q) - 0.9f;
}

float Terrain ( float3 pos ) {
  return length(pos) - 2.0f;
}

float Boundary2 ( float3 p, float3 wal, float dist ) {
  return dot(p, wal) + dist;
}

float Boundary ( float3 p ) {
  return fmin(fmin(fmin(fmin(fmin(
    Boundary2(p, (float3)(1.0f, 0.0f, 0.0f), 5.0f),
    Boundary2(p, (float3)(0.0f, 1.0f, 0.0f), 5.0f)),
    Boundary2(p, (float3)(0.0f, 0.0f, 1.0f), 5.0f)),
    Boundary2(p, (float3)(-1.0f, 0.0f, 0.0f), 5.0f)),
    Boundary2(p, (float3)(0.0f, -1.0f, 0.0f), 5.0f)),
    Boundary2(p, (float3)(0.0f, 0.0f, -1.0f), 5.0f));
}

float Sun ( float3 p ) {
  return length(p - (float3)(1.0f, 5.0f, 0.0f)) - 1.0f;
}

float smin ( float a, float b ) {
  float k = -0.1f;
  float h = clamp(0.5f + 0.5f*(b - a)/k, 0.0f, 1.0f);
  return mix(b, a, h) - k*h*(1.0f - h);
}

// --------------- INTERSECTION ------------------------------------------------
float Map ( float3 p ) {
  float r = PrimitiveB(p);
  r = smin(r, PrimitiveA(p));
  r = smin(r, Terrain(p + (float3)(0.0f, 1.2f, 0.0f)));
  r = fmin(r, Boundary(p));
  r = fmin(r, Sun(p));
  return r;
}

float3 Normal ( float3 p ) {
  const float Delta = 0.01f;

  float3 X = (float3)(Delta, 0.0f, 0.0f),
         Y = (float3)(0.0f, Delta, 0.0f),
         Z = (float3)(0.0f, 0.0f, Delta);

  return normalize((float3)(
    Map(p + X) - Map(p - X),
    Map(p + Y) - Map(p - Y),
    Map(p + Z) - Map(p - Z)
  ));
}


float3 Phong_Contribution ( float3 diffuse, float3 specular, float shininess,
                            float3 pos, float3 cam, float3 lpos, float3 lint ) {
  float3 N = Normal(pos);
  float3 L = normalize(lpos - pos);
  float3 V = normalize(cam - pos);
  float3 R = normalize(reflect(-L, N));

  float dotLN = dot(L, N);
  float dotRV = dot(R, V);

  // check if light is visible
  if ( dotLN < 0.0f ) return (float3)(0.0f);

  if ( dotRV < 0.0f ) { // reflection is opposite direction, apply only diffuse
    return lint * (diffuse * dotLN);
  }

  return lint * (diffuse*dotLN + specular * pow(dotRV, shininess));
}

float3 Phong ( float3 ambient, float3 diffuse, float3 specular,
               float shininess, float3 pos, float3 cam, float timer ) {
  float3 ambient_light = 0.5f * (float3)(1.0f, 1.0f, 1.0f);
  float3 colour = 0.5f * (float3)(1.0f, 1.0f, 1.0f) * ambient;

  float3 light_pos = (float3)(4.0f, 2.0f, 4.0f + sin(timer)*-8.0f);
  float3 light_intensity = (float3)(0.4f, 0.2f, 0.4f);

  colour += Phong_Contribution(diffuse, specular, shininess, pos, cam,
                               light_pos, light_intensity);

  light_pos = (float3)(4.0f, 2.0f - sin(timer*0.25f)*0.5f, 4.0f + cos(timer)*-8.0f);
  light_intensity = (float3)(0.0f, 0.1f, 0.0f);

  colour += Phong_Contribution(diffuse, specular, shininess, pos, cam,
                               light_pos, light_intensity);

  return colour;
}

float March ( Ray ray ) {
  const float max_dist = 8.0f;
  float distance = 0.0f;
  float interval;
  for ( int i = 0; i < 32; ++ i ) {
    interval = Map(ray.origin + ray.dir*distance);
    if ( interval < 0.00001f || interval > max_dist ) break;
    distance += interval;
  }
  if ( interval > max_dist ) return 0.0f;
  return distance;
}

typedef struct T_IntersectionInfo {
  bool intersection;
  float3 position;
  float distance;
  float3 normal, angle;
  float3 colour;
} IntersectionInfo;

// --------------- CAMERA ------------------------------------------------------
typedef struct T_Camera {
  float3 position, lookat, up;
  int2 dimensions;
  float fov;
} Camera;



float3 rayDirection(float fieldOfView, float2 size, float2 fragCoord) {
  float2 xy = fragCoord - size / 2.0f;
  float z = size.y / tan(radians(fieldOfView) / 2.0f);
  return normalize((float3)(xy, -z));
}

Ray Camera_Ray(__global RNG* rng, __global Camera* cam, float time) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  float2 outf = (float2)((float)out.x, (float)out.y);
  float2 dim = (float2)((float)cam->dimensions.x, (float)cam->dimensions.y);
  float2 uv = ((2.0f * outf) - dim)/fmin(dim.x, dim.y);

  mat3 rotatey = rotate_y(sin(time)),
       rotatex = rotate_x(sin(time*0.9f)*0.2);

  float3 eye_pos = mat3_mul(rotatex,
                   mat3_mul(rotatey, (float3)(0.0f, 0.0f, -4.0f)));
  float3 up = (float3)(0.0f, 1.0f, 0.0f);
  float3 forward = mat3_mul(rotatey, (float3)(0.0f, 0.0f, 1.0f));
  float3 right = cross(up, forward);

  return New_Ray(eye_pos,
                 normalize(up * uv.y + right * uv.x + forward)
  );
}

// ------- kernel ---

__kernel void Kernel_Raycast(
        __write_only image2d_t output_image,
        __read_only  image2d_t input_image,
        __global RNG* rng,
        __global Camera* camera,
        __global float* time
      ) {
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  Ray ray = Camera_Ray(rng, camera, *time);


  // ray = Camera_Ray(rng, cam, ray.dir);

  if ( Is_Debug_Print() ) {
    printf("ORIGIN: <%f, %f, %f> DIR: <%f, %f, %f>\n",
              ray.origin.x, ray.origin.y, ray.origin.z,
              ray.dir.x, ray.dir.y, ray.dir.z);
  }

  float dist = March(ray);

  if ( dist > 0.0f ) {
    float3 point = ray.origin + ray.dir*dist;

    float3 colour = Phong(
      (float3)(0.7f, 0.2f, 0.0f),
      (float3)(0.5f, 0.1f, 0.0f),
      (float3)(0.5f, 0.5f, 0.4f),
       1.0f,
      point,
      ray.origin,
      *time
    );

    write_imagef(output_image, out, (float4)(colour, 1.0f));
  } else {
    write_imagef(output_image, out, (float4)(0.0f, 0.0f, 0.0f, 1.0f));
  }

  // float3 colour = (float3)(0.0f, 0.0f, 0.0f);
  // float3 weight = (float3)(1.0f, 1.0f, 1.0f);
  // bool hit = false;

  // int depth = 0;
  // while ( true ) {
  //   float depth_ratio = 1.0f - (1.0f / (4 - depth));
  //   if ( Uniform(rng, 0.0f, 1.0f) >= depth_ratio ) break;
  //   IntersectionInfo info = Raycast_Scene(ray);
  //   if ( !info.intersection ) {
  //     hit = false;
  //     break;
  //   }

  //   hit = true;
  //   float dist = info.distance;
  //   weight *= info.colour;
  //   colour = weight;
  //   break;
  // }

  // if ( hit ) {
  //   float3 old_colour = read_imagef(input_image, out).xyz;
  //   float3 rcolour = mix(colour, old_colour, 0.6f);
  //   // rcolour = fmax(0.0f, fmin(1.0f, rcolour));
  //   // rcolour = pow(rcolour, (float3)(1.0f / 0.5f));
  //   write_imagef(output_image, out, (float4)(rcolour, 1.0));
  // }
}};

