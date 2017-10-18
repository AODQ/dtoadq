module core.shared_info;
import derelict.opencl.cl : cl_event;
static import stl, ocl, glfw, core;
import ocl.opencl;

// --- kernel arguments ----
core.Image image;
core.GLImage gl_image;
ocl.ImageMetaData image_metadata;
ocl.Camera        camera;
ocl.Material[]    material;
ubyte[] rw_image;
cl_uint2[] rng_states;
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
  // update rng
  rng_states.length = core.RResolution_Length(resolution);
  import std.random;
  auto gen = () => uniform(0, uint.max);
  foreach ( ref r; rng_states ) {
    r = To_CLUint2([gen(), gen()]);
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

  return reset;
}
