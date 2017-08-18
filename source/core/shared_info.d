module core.shared_info;
import derelict.opencl.cl : cl_event;
static import stl, ocl, glfw, core;

// --- kernel arguments ----
core.Image image;
core.GLImage gl_image;
ocl.ImageMetaData image_metadata;
ocl.Camera        camera;
ocl.Material[]    material;
ubyte[] rw_image;
int[] rng;
float timer = 0.0f;

//// Sets image buffer (camera, gl images, etc) based off resolution
void Set_Image_Buffer ( core.Resolution resolution, bool force ) {
  if ( force || image.resolution != resolution ) {
    image_metadata.clear_img = true;
    image = core.RImage(resolution);
    rw_image.length = image.x*image.y*4;
    camera.dimensions.x = image.x;
    camera.dimensions.y = image.y;
  }
}

/// Updates all kernel arguments, and Returns if should reset progressive pixel
/// sampling
bool Update_Buffer ( ) {
  static import core.info;
  static float previous_timer = 0.0f;
  bool reset =  glfw.Update_Camera(camera)     |
                (stl.abs(previous_timer - timer) > 0.0f);
  previous_timer = timer;

  {
    import functional, std.random;
    rng = iota(0, 32).map!(n => uniform(-int.max, int.max)).array;
  }

  return reset;
}
