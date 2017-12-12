module core.kernel;
import derelict.opencl.cl : cl_event;
static import stl, ocl, shared_info = core.shared_info, core.image, parser,
              core.info;

float[] debug_values;
bool update_timer  = false,
     update_kernel = true;
bool running = true;

void Initialize_Variables ( ) {
  debug_values = [ 0.0f, 0.0f, 0.0f ];
  shared_info.camera = ocl.Camera([0.0f, 0.0f, 0.0f],
                                  [0.0f, 0.0f, 0.0f],
                                  [shared_info.image.x, shared_info.image.y]);
  shared_info.image_metadata.finished_samples = 0;
  shared_info.material = [];
  shared_info.Set_Image_Buffer(core.image.Resolution.r160_140, true);
  stl.writeln("CAMERA SI: ", shared_info.camera);
}

bool Compile ( string kernel_name ) {
  import parser;
  return ocl.Compile(parser.Parse_Kernel(), kernel_name);
}

void Update_DTOADQ ( ) {
  bool reset_spp = shared_info.Update_Buffer,
       recompile = parser.Recheck_Files() | core.info.Should_Recompile(false);
  ocl.cl_mem imgs = core.image.RImages();
  // -- kernel --
  // stl.writeln("CAMERA: ", shared_info.camera);
  if ( recompile ) {
    import std.datetime;
    stl.writeln("Compile ", (Compile("DTOADQ_Kernel")?"Success":"Failure"),
                ": ", Clock.currTime);
    shared_info.image_metadata.clear_img = true;
    return;
  }
  if ( reset_spp ) {
    shared_info.image_metadata.clear_img = true;
  }
  // update colour
  foreach ( ind; 0 .. shared_info.material.length ) {
    shared_info.ocl_material[ind].albedo =
        ocl.To_CLFloat3(shared_info.material[ind].albedo);
  }
  // stl.writeln(ocl.OCLMaterial.sizeof);
  { // Run kernel
    shared_info.image.Lock(); // Acquire OCL Mem from OGL image object
    ocl.Run(
      ocl.CLStoreMem(shared_info.rw_image),
      shared_info.image.RWrite(),
      ocl.CLStoreMem(shared_info.image_metadata),
      shared_info.camera,
      shared_info.timer,
      ocl.CLPredefinedMem(imgs),
      ocl.CLStoreMem(shared_info.ocl_material),
      debug_values,
      ocl.CLStoreMem(shared_info.rng_states),
      // kernel size
      shared_info.image.x, shared_info.image.y
    );
    shared_info.image.Unlock();
    shared_info.image_metadata.clear_img = false;
  }
}

void Update_Video_Render ( ) {
  Update_DTOADQ();
  static import VI = emitter.video;
  static int cnt = 0;
  auto len = core.image.RResolution_Length(VI.RRender_Resolution);
  if ( shared_info.image_metadata.finished_samples >= len*0.85f ||
       ++cnt > cast(int)(shared_info.image_metadata.spp*1.25f) ) {
    cnt = 0;
    running = !VI.Update(shared_info.rw_image);
    shared_info.image_metadata.finished_samples = 0;
    shared_info.image_metadata.clear_img = true;
  }
}

///
float[] Run_Texture ( string filename ) {
  // save current state
  auto prev_filename  = core.info.RFilename,
       prev_recompile = core.info.Should_Recompile,
       prev_type      = core.info.RType;
  // set map function and run
  core.info.Set_Texture_Function(filename);
  {
    import parser;
    ocl.Compile(Parse_Kernel(), "Texture_kernel");
  }
  float[] t_img;
  t_img.length = 4*1024*1024;
  core.image.GLImage gl_image = new core.image.GLImage(1024, 1024);

  ocl.Run(ocl.CLStoreMem(t_img), 1024.0f, 1024, 1024);
  // restore previous state
  core.info.Set_Map_Function(prev_filename);
  core.info.Set_Kernel_Type(prev_type);
  if ( !prev_recompile ) {
    core.info.Clear_Recompile();
  }
  return t_img;
}
void Update_Texture ( ) {
  if ( shared_info.gl_image is null ) {
    shared_info.gl_image = new core.image.GLImage(1024, 1024);
  }
  // we only want to run once per change
  if ( core.info.Should_Recompile(false) ) {
    shared_info.gl_image.Update(Run_Texture(core.info.RFilename));
  }
}
