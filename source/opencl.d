module opencl;
public import derelict.opencl.cl;
import std.stdio;
import std.conv : to;
import openclmisc;

alias cl_context              CLContext;
alias cl_context_properties[3] CLContextProperties;
alias cl_command_queue        CLCommandQueue;
alias cl_program              CLProgram;
alias cl_device_id            CLDeviceID;
alias cl_platform_id          CLPlatformID;
alias cl_kernel               CLKernel;

class Device {
  CLPlatformID  platform_id;
  CLDeviceID    device_id;
  this ( CLDeviceID _platform_id, CLPlatformID _device_id ) {
    device_id   = _device_id;
    platform_id = _platform_id;
  }
}

private Device device;

struct OpenCLProgram {
  CLContext            context;
  CLContextProperties  properties;
  CLCommandQueue       command_queue;
  CLProgram            program;
  CLKernel             kernel;

  this ( string source ) in {
    assert(device !is null);
  } body {
    properties = [CL_CONTEXT_PLATFORM, cast(int)device.platform_id, 0];
    int err;
    context = clCreateContext(properties.ptr, 1, &device.device_id,
                                null, null, &err);
    CLAssert(err, "Create context");
    command_queue = clCreateCommandQueue(context, device.device_id, 0, &err);
    CLAssert(err, "Create command queue");

    import std.string;
    auto e = source.toStringz;
    program = clCreateProgramWithSource(context, 1, &e, null, &err);
    CLAssert(err, "Create program with source");
    if ( clBuildProgram(program, 0, null, null, null, null) != CL_SUCCESS ) {
      size_t len;
      clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, 0,
                            null, &len);
      char[] log; log.length = len;
      clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG,
                            len, log.ptr, null);
      assert(false, "Error building program:\n" ~ log);
    }
  }

}

void CLAssert(int cond, string err) {
  if ( cond == CL_SUCCESS ) {
    writeln(err, ": SUCCESS");
  } else {
    assert(cond == CL_SUCCESS, err ~ ": " ~ CL_Error_String(cond));
  }
}

auto RPlatforms() {
  cl_platform_id[] platforms;
  cl_uint platform_amt = 10;
  platforms.length = platform_amt;
  CLAssert(clGetPlatformIDs(platform_amt, platforms.ptr, &platform_amt),
           "Get platform id");
  platforms.length = platform_amt;
  return platforms;
}

auto Set_Current_Platform() {
  // --- grab platform from user
  auto platforms = RPlatforms;
  assert(platforms.length >= 1, "No OpenCL platform found");
  int index;
  { // grab index
    import functional;
    foreach ( it; 0 .. platforms.length ) {
      writeln("Platform looper");
      writeln("----- Index: ", it, "\n", RPlatform_Info(platforms[it]));
    }
    write("Please enter a platform to choose: ");
    index = 0;
    // index = readln.chomp.to!int;
    // assert(index >= 0 && index < platforms.length, " invalid platform index");
  }
  return cast(void*)platforms[index];
}

auto Set_Current_Device(CLPlatformID platform_id) {
  // --- get device information
  CLDeviceID device_id;
  CLAssert(clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, null),
           "Get device id");
  return device_id;
}

void Initialize ( ) {
  DerelictCL.load();
  auto platform_id = Set_Current_Platform();
  auto device_id = Set_Current_Device(platform_id);
  DerelictCL.reload(RCL_Version(platform_id));
  // DerelictCL.loadEXT(platform_id);
  device = new Device(platform_id, device_id);
  // writeln(RDevice_Info(device_id));
}
auto Compile(string source) {
  return OpenCLProgram(source);
}

enum PlatformInfo {
  Profile    = 0x0901, Version    = 0x0902,
  Name       = 0x0903, Vendor     = 0x0904,
  Extensions = 0x0905
}

auto RPlatform_Info(PlatformInfo info) in {
  assert(device !is null);
} body {
  return RPlatform_Info(info, device.platform_id);
}
auto RPlatform_Info(PlatformInfo info, CLPlatformID platform) {
  size_t pv_size = 50;
  void[] pv;
  pv.length = 50;
  clGetPlatformInfo(platform, info, pv_size, pv.ptr, &pv_size);
  pv.length = pv_size;
  return pv.to!string;
}

auto RPlatform_Info() in {
  assert(device !is null);
} body {
  return RPlatform_Info(device.platform_id);
}

auto RSupported_Image_Formats(OpenCLProgram platform) {
  import std.string;
  uint if_size;
  clGetSupportedImageFormats(platform.context, CL_MEM_WRITE_ONLY,
      CL_MEM_OBJECT_IMAGE2D, 0, null, &if_size);
  cl_image_format[] imgforms;
  imgforms.length = if_size;
  clGetSupportedImageFormats(platform.context, CL_MEM_WRITE_ONLY,
      CL_MEM_OBJECT_IMAGE2D, if_size, imgforms.ptr, null);
  string output;
  string[int] channel_order_map = [
    CL_R         : "CL_R",         CL_A         : "CL_A",
    CL_INTENSITY : "CL_INTENSITY", CL_LUMINANCE : "CL_LUMINANCE",
    CL_RG        : "CL_RG",        CL_RA        : "CL_RA",
    CL_RGB       : "CL_RGB",       CL_RGBA      : "CL_RGBA",
    CL_ARGB      : "CL_ARGB",      CL_BGRA      : "CL_BGRA"
  ];
  string[int] channel_data_map = [
    CL_FLOAT :            "CL_FLOAT",
    CL_HALF_FLOAT :       "CL_HALF_FLOAT",
    CL_SIGNED_INT32 :     "CL_SIGNED_INT32",
    CL_SIGNED_INT8 :      "CL_SIGNED_INT8",
    CL_SNORM_INT16 :      "CL_SNORM_INT16",
    CL_SNORM_INT8 :       "SNORM_INT8",
    CL_UNORM_INT_101010 : "CL_UNORM_INT_101010",
    CL_UNORM_INT16 :      "CL_UNORM_INT16",
    CL_UNORM_INT8 :       "CL_UNORM_INT8",
    CL_UNORM_SHORT_555 :  "UNORM_SHORT_555",
    CL_UNORM_SHORT_565 :  "CL_UNORM_SHORT_565",
    CL_UNSIGNED_INT16 :   "UNSIGNED_INT16",
    CL_UNSIGNED_INT32 :   "CL_UNSIGNED_INT32",
    CL_UNSIGNED_INT8:     "CL_UNSIGNED_INT8",
  ];
  foreach ( format; imgforms ) {
    auto forder = format.image_channel_order;
    auto fdata  = format.image_channel_data_type;
    auto order = forder in channel_order_map;
    if ( order !is null ) output ~= *order;
    else                  output ~= "Unknown " ~ forder.to!string;
    auto data  = fdata in channel_data_map;
    output ~= " : ";
    if ( data  !is null ) output ~= *data;
    else                  output ~= "Unknown " ~ fdata.to!string;
    output ~= ", ";
  }
  return output;
}

auto RPlatform_Info(CLPlatformID platform) {
  import std.string;
  return
    `Platform ID: %s
     OpenCL Profile: %s
     Platform Version: %s
     Platform Name: %s
     Platform Vendor: %s
     Platform Extensions: %s
    `.format(
       (cast(int)platform).to!string,
       RPlatform_Info(PlatformInfo.Profile,    platform),
       RPlatform_Info(PlatformInfo.Version,    platform),
       RPlatform_Info(PlatformInfo.Name,       platform),
       RPlatform_Info(PlatformInfo.Vendor,     platform),
       RPlatform_Info(PlatformInfo.Extensions, platform)
    );
}

CLVersion RCL_Version(CLPlatformID platform_id) {
  uint major, minor;
  import std.format;
  char[] info = RPlatform_Info(PlatformInfo.Profile, platform_id).to!(char[]);
  formattedRead(info, "OpenCL %d.%d", &major, &minor);
  CLVersion cl_version;
  string version_error = "OpenCL version " ~ major.to!string ~ "." ~
                         minor.to!string ~ " unsupported";
  switch ( major ) {
    default: assert(0, version_error);
    case 1:
      switch ( minor ) {
        default: assert(0, version_error);
        case 0: return CLVersion.CL10;
        case 1: return CLVersion.CL11;
        case 2: return CLVersion.CL12;
      }
    case 2:
      assert(0, version_error ~ " (no D bindings :[)");
  }
}

/*
khronos.org/registry/OpenCL/sdk/2.1/docs/man/xhtml/clCreateFromGLTexture.html
describes corresponding opencl and opengl image formats
*/
struct CLImage {
  ubyte[] buffer;
  int width, height;
}

CLImage Test_OpenCL ( ) {
  int err;
  immutable(int) dim = 64;

  auto program = Compile(q{
    __kernel void hello(__read_only image2d_t input_image,
                        __write_only image2d_t output_image) {
      size_t tid = get_global_id(0);
      const sampler_t sample = CLK_FILTER_NEAREST;
      uint4 col = read_imageui(input_image, sample, (int2)(tid, 20));
      write_imageui(output_image, (int2)(tid, 20), col);
    }
  });

  import std.string;
  auto e = "hello".toStringz;
  auto kernel = clCreateKernel(program.program, "hello", &err);
  CLAssert(err, "Create kernel");
  ubyte[] input_image_buffer;
  input_image_buffer.length = dim*dim;
  import std.random;
  foreach ( ref i; input_image_buffer ) {
    i = uniform(0, 255).to!ubyte;
  }
  writeln("Supported image formats: ", RSupported_Image_Formats(program));
  auto Img_format = cl_image_format(CL_RGBA, CL_UNORM_INT8);
  auto Img_descriptor = cl_image_desc(CL_MEM_OBJECT_IMAGE2D, dim, dim,
                                      1, 0, 0, 0, 0, 0, null);
  auto input  = clCreateImage(program.context, CL_MEM_READ_ONLY,
                    &Img_format, &Img_descriptor, null, &err);
  CLAssert(err, "Create image input");
  auto image  = clCreateImage(program.context, CL_MEM_WRITE_ONLY,
                  &Img_format, &Img_descriptor, null, &err);
  CLAssert(err, "Create image output");

  // load data into input buffer
  size_t[3] origin = [0, 0, 0],
            region = [1, 1, 1];
  ubyte[] image_buffer;
  CLAssert(clEnqueueWriteImage(program.command_queue, input, CL_TRUE,
                   origin.ptr, region.ptr, 0, 0, input_image_buffer.ptr, 0,
                   null, null),
    "enqueue write image");
  // set arg list for kernel command
  CLAssert(clSetKernelArg(kernel, 0, cl_mem.sizeof, &input),
           "Setting kernel of input");
  CLAssert(clSetKernelArg(kernel, 1, cl_mem.sizeof, &image),
           "Setting kernel of image");
  size_t global = input_image_buffer.length;

  import std.datetime;
  writeln("Starting");
  StopWatch watch;
  watch.start;
  // enqueue kernel command for execution
  CLAssert(clEnqueueNDRangeKernel(program.command_queue, kernel, 1, null,
                                  &global, null, 0, null, null),
           "Enqueue kernel command for execution");
  CLAssert(clFinish(program.command_queue),
           "Finish execution");
  watch.stop;
  import core.time;
  writeln("Finished, duration: ", watch.peek().msecs, " milliseconds");
  // copy image from image buffer
  image_buffer.length = dim*dim;
  foreach ( ref i; image_buffer ) {
    i = 255;
  }
  CLAssert(clEnqueueReadImage(program.command_queue, image, CL_TRUE,
            origin.ptr, region.ptr, 0, 0, image_buffer.ptr, 0, null, null),
           "Enqueue image");
  auto img = CLImage(image_buffer.dup, dim, dim);
  // cleanup opencl resources
  clReleaseMemObject(input);
  clReleaseMemObject(image);
  clReleaseProgram(program.program);
  clReleaseKernel(kernel);
  clReleaseCommandQueue(program.command_queue);
  clReleaseContext(program.context);
  writeln("asdf");
  return img;
}
