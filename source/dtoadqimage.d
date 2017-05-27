module dtoadqimage;
import globals, opencl, scene;
static import OCL = opencl;

private struct CLImage {
  OCL.CLPredefinedMem[2] cl_handle;
  uint[2] gl_texture;
  this ( int width, int height ) {
    foreach ( it; 0 .. 2 ) {
      import derelict.opengl3.gl3;
      glGenTextures(1, &gl_texture[it]);
      glBindTexture(GL_TEXTURE_2D, gl_texture[it]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height,
                                  0, GL_RGBA, GL_FLOAT, null);
      glBindTexture(GL_TEXTURE_2D, 0);

      cl_handle[it] = OCL.Create_CLGL_Texture(gl_texture[it]);
      writeln("GL_HANDLE: ", gl_texture[it], " CL_HANDLE: ", cl_handle[it]._mem,
              "  ", width, "x", height);
    }
  }

  ~this ( ) {
  }
}

enum Resolution {
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

  private CLImage image;
  private int rw_access = 0;

  this ( ImageInfo _image_info_, CLImage image_ ) {
    _image_info = _image_info_;
    image = image_;
  }

  void Lock ( ) {
    import derelict.opengl3.gl3;
    glFlush();
    static import OCL = opencl;
    OCL.Lock_CLGL_Images(image.cl_handle);
  }

  auto Unlock ( ) {
    static import OCL = opencl;
    auto event = OCL.Unlock_CLGL_Images(image.cl_handle);
    OCL.Flush();
    rw_access ^= 1;
    return event;
  }

  auto RRead   ( ) { return image.cl_handle [rw_access  ]; }
  auto RWrite  ( ) { return image.cl_handle [rw_access^1]; }
  auto RRender ( ) { return image.gl_texture[rw_access  ]; }
}

private Resolution [string] resolution_string;

private ImageInfo[Resolution] image_info;
private Image    [Resolution] images;

auto RImage ( Resolution resolution ) {
  return images[resolution];
}

string[] RResolution_Strings ( ) {
  return [ "640x360",
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
      r640_360   : II(r640_360,   "640x360",   640,  360,  230_400    ),
      r960_540   : II(r960_540,   "960x540",   960,  540,  518_400    ),
      r1366_768  : II(r1366_768,  "1366x768",  1366, 768,  1_049_088  ),
      r1920_1080 : II(r1920_1080, "1920x1080", 1920, 1080, 2_073_600  ),
      r4096_2304 : II(r4096_2304, "4096x2304", 4096, 2304, 9_437_184  ),
      r8192_4608 : II(r8192_4608, "8192x4608", 8192, 4608, 37_748_736 ),
    ];
    resolution_string = [
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
    images[res.resolution] = Image(res, CLImage(res.x, res.y));
}

void Clean_Up ( ) {
  foreach ( img; images ) {
    import derelict.opencl.cl, derelict.opengl3.gl3;
    clReleaseMemObject(img.image.cl_handle[0].cl_mem);
    clReleaseMemObject(img.image.cl_handle[1].cl_mem);
    glDeleteTextures(2, &img.image.gl_texture[0]);
  }
}
