module structure.octree;
import globals;
import ray;

class OctreeNode {
  OctreeNode[8] children;
  gln.vec3 origin, dim;
  Primitive primitive;
public:
  this ( gln.vec3 origin_, gln.vec3 dim_ ) {
    origin = origin_; dim = dim_;
  }

  void Insert ( Primitive prim ) {
    if ( Is_Leaf ) {
      if ( !primitive ) {
        primitive = prim;
        return;
      }
      foreach ( i; 0 .. 8 ) {
        gln.vec3 new_origin = origin, new_dim = dim/2.0f;
        foreach ( p; 0 .. 3 )
          new_origin.vector[p] += dim.vector[p]*(i&(4/(1+p)) ? 0.5f : -0.5f);
        children[i] = new OctreeNode(new_origin, new_dim);
      }

      auto old = primitive;
      primitive = null;
      Insert(old);
    } else {
      children[ROctant_Mask(prim.origin)].Insert(prim);
    }
  }

  void RBounds ( out gln.vec3 lower, out gln.vec3 upper ) inout {
    lower = origin - dim;
    upper = origin + dim;
  }

  immutable(Primitive) Ray_Intersection ( inout(Ray) ray ) inout {
    gln.vec3 min, max;
    if ( Is_Leaf ) {
      return cast(immutable)primitive;
    }
    Primitive prim = null;
    float dist = float.max;
    foreach ( i; 0 .. 8 ) {
      auto chnode = children[i];
      chnode.RBounds(min, max);
      auto res = Ray_Intersection_AABB(min, max, ray);
      if ( res.hit ) {
        auto tprim = chnode.Ray_Intersection(ray);
        if ( tprim is null ) continue;
        return tprim;
        // tprim.RBounds(min, max);
        // auto tdist = Ray_Intersection_AABB(min, max, ray);
        // if ( tdist < dist ) {
        //   tdist = dist;
        //   prim = cast(Primitive)tprim;
        // }
      }
    }
    return cast(immutable)prim;
  }

  bool Is_Leaf ( ) inout { return children[0] is null; }
  size_t ROctant_Mask ( gln.vec3 point ) {
    size_t oct = 0;
    if ( origin.x < point.x ) oct |= 4;
    if ( origin.y < point.y ) oct |= 2;
    if ( origin.z < point.z ) oct |= 1;
    return oct;
  }
}

auto Construct_Octree ( gln.vec3 origin, gln.vec3 dim, Primitive[] primitives ) {
  auto node = new OctreeNode(origin, dim);
  foreach ( prim; primitives ) {
    node.Insert(prim);
  }
  return node;
}

class Primitive {
public:
  gln.vec3 origin, dim, colour;

  this ( gln.vec3 origin_, gln.vec3 dim_, gln.vec3 colour_ ) {
    origin = origin_; dim = dim_; colour = colour_;
  }
  void RBounds ( out gln.vec3 lower, out gln.vec3 upper ) inout {
    lower = origin - dim;
    upper = origin + dim;
  }
}

struct IntInfo {
  bool hit;
  float dist;
}

private auto Ray_Intersection_AABB ( gln.vec3 min, gln.vec3 max, Ray ray ) {
  import std.math : fmin, fmax;
  float t1 = (min.x - ray.origin.x)*ray.invdir.x,
        t2 = (max.x - ray.origin.x)*ray.invdir.x;
  float tmin = fmin(t1, t2),
        tmax = fmax(t1, t2);
  foreach ( i; 1 .. 3 ) {
    t1 = (min.vector[i] - ray.origin.vector[i])*ray.invdir.vector[i];
    t1 = (max.vector[i] - ray.origin.vector[i])*ray.invdir.vector[i];
    tmin = fmax(tmin, fmin(t1, t2));
    tmax = fmin(tmax, fmax(t1, t2));
  }

  if ( tmin < fmax(tmax, 0.0f) ) return IntInfo(false, -1.0f);
  return IntInfo(true, tmin);
}