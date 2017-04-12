module camera;
import globals, ray;

class Camera {
  gln.vec3 position, lookat, up;
  int[2] dimensions;
  float fov;

  this ( gln.vec3 position_, gln.vec3 lookat_, int x, int y ) {
    dimensions[0] = x;
    dimensions[1] = y;
    position      = position_;
    lookat        = lookat_;
    up            = gln.vec3(0.0f, 1.0f, 0.0f);
    fov           = 60.0f;
  }

  Ray Generate_Ray ( int x, int y ) {
    float[2] dim = [ dimensions[0]/2.0f, dimensions[1]/2.0f ];
    import std.random : uniform;
    import std.math : tan;
    gln.vec3 z_axis = TNormalize(lookat - position),
             x_axis = gln.cross(z_axis, up),
             y_axis = gln.cross(x_axis, z_axis);

    float fov_rad = fov * 3.14159f/180.0f;
    float focus_dist = 1.0f/tan(fov_rad/2.0f);
    float aspect_ratio = dim[0]/dim[1];

    float[2] pixel = [cast(float)x + uniform(-0.5f, 0.5f),
                      cast(float)y + uniform(-0.5f, 0.5f)];
    float x_proj = 2.0f * aspect_ratio * (pixel[0]/dim[0] - 0.5f),
          y_proj = 2.0f * (0.5f - pixel[1]/dim[1]);
    return Ray(position, x_axis*x_proj + y_axis*y_proj + z_axis*focus_dist);
  }
}


auto TNormalize ( gln.vec3 vec ) {
  import std.math : abs;
  float mag = abs(vec.x) + abs(vec.y) + abs(vec.z);
  return vec/mag;
}
