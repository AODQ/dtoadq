module glfw.image;
import globals, opencl;
static import OCL = opencl;
static import DIMG = dtoadqimage;
static import KI = kernelinfo;

private void Allocate_GL_Texture ( ref uint gl_texture, int width, int height ){
  import derelict.opengl3.gl3;
  glGenTextures(1, &gl_texture);
  glBindTexture(GL_TEXTURE_2D, gl_texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height,
                              0, GL_RGBA, GL_FLOAT, null);
  glBindTexture(GL_TEXTURE_2D, 0);
}

private struct CLGLImage {
  OCL.CLPredefinedMem cl_handle;
  uint gl_texture;
  this ( int width, int height ) {
    foreach ( it; 0 .. 2 ) {
      Allocate_GL_Texture(gl_texture, width, height);
      cl_handle = OCL.Create_CLGL_Texture(gl_texture);
    }
  }

  ~this ( ) {
  }
}

enum Resolution {
  r160_140,
  r640_360,
  r960_540,
  r1366_768,
  r1920_1080,
  r4096_2304,
  r8192_4608,
}

private struct ImageInfo {
  Resolution resolution;
  string value;
  uint x, y, total;
}
struct Image {
  ImageInfo _image_info;
  alias _image_info this;

  private CLGLImage image;
  private int rw_access = 0;

  this ( ImageInfo _image_info_, CLGLImage image_ ) {
    _image_info = _image_info_;
    image = image_;
  }

  void Lock ( ) {
    import derelict.opengl3.gl3;
    static import OCL = opencl;
    OCL.Lock_CLGL_Image(image.cl_handle);
  }

  void Unlock ( ) {
    static import OCL = opencl;
    OCL.Unlock_CLGL_Image(image.cl_handle);
  }

  auto RWrite  ( ) { return image.cl_handle ; }
  auto RRender ( ) { return image.gl_texture; }

  auto RResolution ( ) { return _image_info.resolution; }
}


public class GLImage {
  int width, height;
  uint gl_texture;
  this ( int width_, int height_ ) {
    width = width_; height = height_;
    Allocate_GL_Texture(gl_texture, width, height);
  }

  void Update(inout float[] img) {
    import derelict.opengl3.gl3;
    glBindTexture(GL_TEXTURE_2D, gl_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0,
                GL_RGBA, GL_FLOAT, cast(void*)img.ptr);
    glBindTexture(GL_TEXTURE_2D, 0);
  }
}

public struct CLImage {
  cl_mem memory;
  this ( int width, int height, int amt, void* data ) {
    memory = OCL.Create_CL_Image(CL_MEM_READ_ONLY,
      cl_image_format(CL_RGBA, CL_FLOAT),
      cl_image_desc  (CL_MEM_OBJECT_IMAGE2D, width, height, 1,
                       0,  0, 0, 0, 0, null), data);
  }
}

private cl_mem cl_image_array;
private string[] cl_image_names;

/// Sets cl_image_array
private void Create_Images(int width, int height, int len, inout float[] data){
  if ( cl_image_array !is null ) {
    clReleaseMemObject(cl_image_array);
  }
  cl_image_array = OCL.Create_CL_Image(CL_MEM_READ_ONLY,
    cl_image_format(CL_RGBA, CL_FLOAT),
    cl_image_desc  (CL_MEM_OBJECT_IMAGE2D_ARRAY, width, height, 1, len,
                    0, 0, 0, 0, null), cast(void*)data.ptr);
}

void Create_Images ( string[] filenames ) {
  if ( filenames.length == 0 ) return;
  // check if change in files (there's order, whatever)
  if ( filenames.length == cl_image_names.length ) {
    bool same_files = true;
    for ( size_t i = 0; i != filenames.length; ++ i )
      same_files &= cl_image_names[i] == filenames[i];
    if ( same_files ) return;
  }
  cl_image_names = filenames;
  // recreate images... might take some time :-[
  float[] data;
  foreach ( fil; filenames ) {
    if ( !std.file.exists(fil) ) {
      writeln("ERROR: TEXTURE FILE '", fil, "' DOES NOT EXIST");
      return;
    }
    data ~= Kernel_Run_Texture(1024, fil);
  }

  Create_Images(1024, 1024, cast(int)filenames.length, data);
}

/// Compiles, runs kernel & returns the texture as a DIMG.GLImage
float[] Kernel_Run_Texture ( int dim, string filename ) {
  static import DIMG = dtoadqimage;
  // save current state
  auto prev_filename  = KI.RFilename,
       prev_recompile = KI.Should_Recompile;
  KI.Set_Map_Function(filename);
  {
    import parser;
    OCL.Compile(Parse_Kernel(), "Texture_kernel");
  }
  float[] t_img;
  t_img.length = 4*dim*dim;
  DIMG.GLImage gl_image = new DIMG.GLImage(dim, dim);

  float fl_dim = dim;
  OCL.Run(OCL.CLStoreMem(t_img), fl_dim, dim, dim);
  // restore previous state
  KI.Set_Map_Function(prev_filename);
  if ( !prev_recompile ) {
    KI.Clear_Recompile();
  }
  return t_img;
}

cl_mem RImages ( ) {
  // sort of like a singleton, we always want cl_image_array to be a valid
  // cl_mem object
  if ( cl_image_array is null ) {
    Create_Images(["projects/globals/textures/txCheckerboard.txt"]);
  }
  return cl_image_array;
}

private Resolution [string] resolution_string;

private ImageInfo[Resolution] image_info;
private Image    [Resolution] images;

auto RImage_Info ( Resolution resolution ) {
  return image_info[resolution];
}

auto RImage ( Resolution resolution ) {
  return images[resolution];
}

ulong RResolution_Length(Resolution resolution) {
  return image_info[resolution].x * image_info[resolution].y;
}

string[] RResolution_Strings ( ) {
  return [ "160x140",
           "640x360",
           "960x540",
           "1366x768",
           "1920x1080",
           "4096x2304",
           "8192x4608",
  ];
}

static this ( ) {
  with ( Resolution ) {
    alias II = ImageInfo;
    image_info = [
      r160_140   : II(r160_140,   "160x140",   160,  140,  22_400    ),
      r640_360   : II(r640_360,   "640x360",   640,  360,  230_400    ),
      r960_540   : II(r960_540,   "960x540",   960,  540,  518_400    ),
      r1366_768  : II(r1366_768,  "1366x768",  1366, 768,  1_049_088  ),
      r1920_1080 : II(r1920_1080, "1920x1080", 1920, 1080, 2_073_600  ),
      r4096_2304 : II(r4096_2304, "4096x2304", 4096, 2304, 9_437_184  ),
      r8192_4608 : II(r8192_4608, "8192x4608", 8192, 4608, 37_748_736 ),
    ];
    resolution_string = [
      "160x140"   : r160_140,
      "640x360"   : r640_360,
      "960x540"   : r960_540,
      "1366x768"  : r1366_768,
      "1920x1080" : r1920_1080,
      "4096x2304" : r4096_2304,
      "8192x4608" : r8192_4608,
    ];
  }
}


void Initialize() {
  foreach ( res; image_info )
    images[res.resolution] = Image(res, CLGLImage(res.x, res.y));
}

void Clean_Up ( ) {
  foreach ( img; images ) {
    import derelict.opencl.cl, derelict.opengl3.gl3;
    clReleaseMemObject(img.image.cl_handle);
    glDeleteTextures(1, &img.image.gl_texture);
  }
}
