module camera;
import globals;

class Camera {
public:
  gln.vec3 position, lookat, up;
  size_t[2] dimensions;
  float fov;

  this ( gln.vec3 position_, gln.vec3 lookat_, size_t[2] dimensions_ ) {
    dimensions = dimensions_.dup;
    position      = position_;
    lookat        = lookat_;
    up            = gln.vec3(0.0f, 1.0f, 0.0f);
    fov           = 100.0f;
  }

  this ( gln.vec3 position_, gln.vec3 lookat_, size_t x, size_t y ) {
    this ( position_, lookat_, [x, y] );
  }

  this ( inout(Camera) cam ) {
    this ( cam.position, cam.lookat, cam.dimensions );
  }

}


auto TNormalize ( gln.vec3 vec ) {
  import std.math : abs;
  float mag = abs(vec.x) + abs(vec.y) + abs(vec.z);
  return vec/mag;
}
