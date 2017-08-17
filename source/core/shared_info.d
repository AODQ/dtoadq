module core.shared_info;
import derelict.opencl.cl : cl_event;
static import core.image, ocl, glfw;

// --- kernel arguments ----
core.image.Image image;
core.image.GLImage gl_image;
ocl.structs.ImageMetaData image_meta_data;
ocl.structs.Camera        camera;
ocl.structs.Material[]    material;
ubyte[] rw_image;
int[] rng;
float timer = 0.0f;

void Set_Image_Buffer ( core.image.Resolution resolution, bool force ) {
    if ( force || image.resolution != resolution ) {
      shared_info.clear_img = true;
      image = core.image.RImage(resolution);
      rw_image.length = image.x*image.y*4;
      camera.dimensions.x = image.x;
      camera.dimensions.y = image.y;
    }
}

/// Returns if should reset progressive pixel sampling
bool Update_Buffer ( ) {
  static import kernel.info;
  static float previous_timer = 0.0f;
  bool reset =  glfw.Update_Camera(camera)     |
                kernel.info.Should_Recompile() |
                (stl.abs(previous_timer - timer) > 0.0f);
  previous_timer = timer;

  {
    import functional;
    rng = iota(0, 32).map!(n => uniform(-int.max, int.max)).array;
  }

  return reset;
}
