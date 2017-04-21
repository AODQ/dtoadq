module softwarerenderer.scene;
import globals;

class Scene {
  Primitive[] primitives;
}

class Primitive {
  string distfunc;
public:
  this ( string distfunc_ ) {
    distfunc = distfunc_;
  }

  float RDist ( gln.vec3 point ) {
    import aodheme;
    auto tstr = distfunc.format(point.x, point.y, point.z);
    return tstr.Evaluate().RFloat;
  }
}
