module ray;
import globals;

struct Ray {
  gln.vec3 origin, dir, invdir;
  int[3] sign;

  this ( gln.vec3 origin_, gln.vec3 dir_ ) {
    origin = origin_;
    dir = dir_;
    invdir = 1.0f/dir;
    sign = [invdir.x < 0.0f, invdir.y < 0.0f, invdir.z < 0.0f];
  }
}
