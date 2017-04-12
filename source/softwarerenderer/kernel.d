module softwarerenderer.kernel;
import softwarerenderer.image;
import structure.octree;
import camera;
import globals;
import ray;


private struct IntersectionInfo {
  bool intersection;
  gln.vec3 position;
  float distance;
  Primitive primitive;
  gln.vec3 normal, angle, colour;
}

private auto Raycast_Scene ( OctreeNode node, Ray ray )  {
  IntersectionInfo info;
  info.intersection = false;

  float dist;
  Primitive primitive = node.Ray_Intersection(ray);

  if ( primitive ) {
    info.intersection = true;
    info.colour = primitive.colour;
  }

  return info;
}

void Kernel_Raycast ( ref Image image, OctreeNode node,
                      Camera camera, int x, int y ) {
  Ray ray = camera.Generate_Ray(x, y);

  gln.vec3 colour = gln.vec3(0.0f, 0.0f, 0.0f),
           weight = gln.vec3(1.0f, 1.0f, 1.0f);
  bool hit = false;

  while ( true ) {
    auto info = Raycast_Scene(node, ray);
    if ( !info.intersection ) {
      hit = true; break;
    }

    hit = true;

    colour = info.colour;
    break;
  }

  if ( hit ) {
    image.Apply(x, y, colour);
  }
}
