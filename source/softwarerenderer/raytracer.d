module softwarerenderer.raytracer;
import softwarerenderer.kernel, softwarerenderer.image;
import structure.octree;
import camera, globals;


OctreeNode RNew_Octree ( ) {
  gln.vec3 Rand ( ) {
    float RFl ( ) {
      import std.random : uniform;
      return uniform(-1.0f, 1.0f);
    }
    return gln.vec3(RFl, RFl, RFl);
  }

  immutable size_t Amt_pts = 100;
  Primitive[] primitives;
  foreach ( i; 0 .. Amt_pts ) {
    primitives ~= new Primitive(Rand, Rand/10.0f, Rand);
  }

  gln.vec3 origin = gln.vec3(0.0f, 0.0f, 0.0f),
           dim    = gln.vec3(10.0f, 10.0f, 10.0f);
  return Construct_Octree(origin, dim, primitives);
}


class Raytracer : AOD.Entity {
  immutable(int) Img_dim = 256;
  Camera camera;
  Image image;
  OctreeNode octree;
public:
  this ( ) {
    super();
    Set_Position(AOD.R_Window_Width/2.0f, AOD.R_Window_Height/2.0f);
    Set_Size(AOD.Vector(Img_dim, Img_dim), true);
    octree = RNew_Octree();
    camera = new Camera(gln.vec3(1.0f, 1.0f, -1.0f),
                        gln.vec3(-1.0f, 0.0f, 0.0f), Img_dim, Img_dim);
    image = new Image(Img_dim, Img_dim);
  }

  void Run ( ) {
    // foreach ( i; 0 .. Img_dim )
    //   foreach ( j; 0 .. Img_dim )
    //     Kernel_Raycast(image, octree, camera, i, j);
  }

  override void Render ( ) @trusted {
    static int counter = 0;
    if ( ++ counter > 60 ) {
      writeln("FPS: ", AOD.R_FPS());
      counter = 0;
    }
    Run();
    Set_Sprite(image.To_OGL_Sprite());
  }
}
