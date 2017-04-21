module softwarerenderer.kernel;
import softwarerenderer.image, softwarerenderer.scene;
import camera;
import globals;
import ray;


private struct IntersectionInfo {
  bool intersection;
  gln.vec3 position;
  float distance;
  Primitive primitive;
}

private auto Raytrace_Scene ( inout Scene scene, inout Ray ray )  {
  IntersectionInfo info;
  info.intersection = false;

  // float dist;
  // immutable(Primitive) primitive = node.Ray_Intersection(ray);

  // if ( primitive ) {
  //   info.intersection = true;
  //   info.colour = primitive.colour;
  // }

  return info;
}

struct PixelInfo {
  bool hit;
  gln.vec3 pixel;
}

auto Kernel_Raytrace ( inout Scene scene,
                       inout Camera camera, size_t x, size_t y ) {
  Ray ray = camera.Generate_Ray(x, y);

  gln.vec3 colour = gln.vec3(0.0f, 0.0f, 0.0f),
           weight = gln.vec3(1.0f, 1.0f, 1.0f);
  bool hit = false;

  while ( true ) {
    auto info = Raytrace_Scene(scene, ray);
    if ( !info.intersection ) {
      hit = true; break;
    }

    hit = true;

    // colour = info.colour;
    break;
  }

  return PixelInfo(hit, colour);
}


import std.concurrency;
void Raytrace_Thread_Manager ( Tid parent_id,
                              inout Scene scene,
                              inout Camera base_camera,
                              size_t lx, size_t ly, size_t hx, size_t hy ) {
  Camera camera = new Camera(base_camera);
  import core.thread, core.time, std.variant : Variant;
  size_t img_id;
  while ( true ) {
    foreach ( x; lx .. hx ) {
      foreach ( y; ly .. hy ) {
        auto result = Kernel_Raytrace(scene, camera, x, y);
        send(parent_id, x, y, result, img_id);
        // -- check camera --
        bool exit;
        receiveTimeout(dur!("nsecs")(1),
          (immutable(Camera) camera_) {
            ++ img_id;
            camera = new Camera(camera_);
            exit = true;
          },
          (Variant variant) {
            "Parameter mismatch!".writeln;
            assert(false, "Parameter mismatch");
          }
        );

        if ( exit ) goto END_PIXEL_CALC;
      }
    }

    END_PIXEL_CALC:
    Thread.sleep(dur!("msecs")(25));
  }
}
