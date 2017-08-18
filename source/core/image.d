module core.image;
static import stl, ocl;
import core.kernel : Run_Texture;
/// TODO: split the opencl related components into ocl module

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
  ocl.CLPredefinedMem cl_handle;
  uint gl_texture;
  this ( int width, int height ) {
    foreach ( it; 0 .. 2 ) {
      Allocate_GL_Texture(gl_texture, width, height);
      cl_handle = ocl.Create_CLGL_Texture(gl_texture);
    }
  }

  ~this ( ) {
  }
}

/// Only have to introduce a new enum element here to add a resolution
/// (You can't have run-time resolutions as all OGL-ocl image buffers must be
///  constructed before they share contexts)
enum Resolution {
  r160_140,   r640_360,   r960_540,
  r1366_768,  r1920_1080, r4096_2304,
  r8192_4608,
}

/// Returns struct{enum_ "r160_140", str "160x140", x 160, y 140}
/// Should only be used during compile-time, use map_info for run-time
private auto Resolution_Info ( Resolution r ) {
  struct Info { string enum_, str; int x, y; }
  string enum_ = stl.to!string(r),
         str   = stl.replace(enum_[1 .. $], "_", "x");
  auto dims = stl.array(stl.map!(n => stl.to!int(n))(stl.splitter(str, "x")));
  return Info(enum_, str, dims[0], dims[1]);
}

/// Returns strings of all resolutions, indexeable by enum
string[] Resolution_Strings ( ) {
  import std.traits;
  string[] results;
  foreach ( res; [EnumMembers!Resolution] )
    results ~= stl.to!string(res);
  return results;
}

string RResolution_String ( Resolution r ) {
  string RResolution_String_Switch ( string res_var ) {
    import std.traits, stl : format;
    string results = `final switch ( %s ) {`.format(res_var);
    foreach ( immutable res; [EnumMembers!Resolution] ) {
      auto info = Resolution_Info(res);
      results ~= `case Resolution.%s: return "%s";`
                  .format(info.enum_, info.str);
    }
    return results ~ `}`;
  }

  mixin(RResolution_String_Switch("r"));
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
    ocl.Lock_CLGL_Image(image.cl_handle);
  }

  void Unlock ( ) {
    ocl.Unlock_CLGL_Image(image.cl_handle);
  }

  ocl.CLPredefinedMem RWrite  ( ) { return image.cl_handle ; }
  uint                RRender ( ) { return image.gl_texture; }

  Resolution RResolution ( ) { return _image_info.resolution; }
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
  ocl.cl_mem memory;
  this ( int width, int height, int amt, void* data ) {
    memory = ocl.Create_CL_Image(ocl.CL_MEM_READ_ONLY,
      ocl.cl_image_format(ocl.CL_RGBA, ocl.CL_FLOAT),
      ocl.cl_image_desc  (ocl.CL_MEM_OBJECT_IMAGE2D, width, height, 1,
                       0,  0, 0, 0, 0, null), data);
  }
}

private ocl.cl_mem cl_image_array;
private string[] cl_image_names;

/// Sets cl_image_array
private void Create_Images(int width, int height, int len, inout float[] data){
  import ocl;
  if ( cl_image_array !is null ) {
    clReleaseMemObject(cl_image_array);
  }
  cl_image_array = ocl.Create_CL_Image(CL_MEM_READ_ONLY,
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
    if ( !stl.file.exists(fil) ) {
      stl.writeln("ERROR: TEXTURE FILE '", fil, "' DOES NOT EXIST");
      return;
    }
    data ~= Run_Texture(fil);
  }

  Create_Images(1024, 1024, cast(int)filenames.length, data);
}

ocl.cl_mem RImages ( ) {
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

ImageInfo RImage_Info ( Resolution resolution ) {
  return image_info[resolution];
}

Image RImage ( Resolution resolution ) {
  return images[resolution];
}

ulong RResolution_Length(Resolution resolution) {
  return image_info[resolution].x * image_info[resolution].y;
}

static this ( ) {
  string RImage_Info_Constructor ( string res_var ) {
    import std.traits, stl : format;
    string results = `%s = [`.format(res_var);
    foreach ( immutable res; [EnumMembers!Resolution] ) {
      auto info = Resolution_Info(res);
      // r160_140 : ImageInfo(r160_140, "160x140", 160, 140, 22_400)
      results ~= `%s : ImageInfo(%s, "%s", %d, %d, %d),`
              .format(info.enum_, info.enum_, info.str, info.x, info.y,
                      info.x*info.y);
    }
    return results ~ `];`;
  }

  with ( Resolution ) {
    mixin(RImage_Info_Constructor("image_info"));
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
