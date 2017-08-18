module core.renderer;
static import core.info, ocl, gl_render = glfw.gl_renderer,
              shared_info = core.shared_info;
private ocl.cl_event image_unlock_event;

void Render ( ) {
  import gui : Imgui_Render;
  // Render GUI
  shared_info.image_metadata.clear_img |=
       Imgui_Render(shared_info.material, shared_info.camera);
  // Render texture or raytrace/dtq opengl texture
  if ( core.info.RType == core.info.KernelType.Texture ) {
    if ( shared_info.gl_image !is null )
      gl_render.Render(shared_info.gl_image.gl_texture, image_unlock_event);
  } else {
    gl_render.Render(shared_info.image.RRender, image_unlock_event);
  }
}
