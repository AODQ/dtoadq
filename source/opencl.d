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
alias cl_mem                  CLMem;

class Device {
  CLPlatformID  platform_id;
  CLDeviceID    device_id;
  this ( CLDeviceID _platform_id, CLPlatformID _device_id ) {
    device_id   = _device_id;
    platform_id = _platform_id;
  }
}

private Device device;

enum BufferType { write_only = CL_MEM_WRITE_ONLY,
                  read_only  = CL_MEM_READ_ONLY };

struct BaseBuffer(BType) {
  BufferType type;
  CLMem cl_handle;
  BType[] data;
}

struct OpenCLImage {
  BaseBuffer!float _mybuffer;
  alias _mybuffer this;
  cl_image_format  format;
  cl_image_desc    description;
  size_t width, height;
  this ( BufferType _type, float[] _data, cl_image_format _format,
         cl_image_desc _description, size_t _width, size_t _height) {
    _mybuffer = BaseBuffer!float(_type, CLMem(), _data);
    format = _format;
    description = _description;
    width = _width; height = _height;
  }
}

struct OpenCLBuffer(BType) {
  BaseBuffer!BType _mybuffer;
  alias _mybuffer this;
  uint length;
  CLMem cl_length_handle;
  this ( BufferType _type, BType[] _data ) {
    _mybuffer = BaseBuffer!BType(_type, CLMem(), _data);
    length = cast(uint)data.length;
  }
}

class OpenCLProgram {
  CLContext            context;
  CLContextProperties  properties;
  CLCommandQueue       command_queue;
  CLProgram            program;
  CLKernel             kernel;
  CLMem[]              mem_objects;
  int err;
public:
  this ( string source ) in {
    assert(device !is null);
  } body {
    err = CL_SUCCESS;
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
  ~this() {
    Cleanup();
  }

  void Set_Kernel(string kernel_name) {
    kernel = clCreateKernel(program, kernel_name.ptr, &err);
    CLAssert(err, "clCreateKernel");
  }

  OpenCLImage Set_Image_Buffer(BufferType type, int dim, int ind) {
    return Set_Image_Buffer(type, dim, dim, ind);
  }
  OpenCLImage Set_Image_Buffer(BufferType type, int width, int height, int ind){
    OpenCLImage image = OpenCLImage(
      type, [],
      cl_image_format(CL_RGBA, CL_FLOAT),
      cl_image_desc (CL_MEM_OBJECT_IMAGE2D, width, height, 1,
                    0, 0, 0, 0, 0, null),
      width, height
    );
    image.cl_handle = clCreateImage(context, type, &image.format,
                                    &image.description, null, &err);
    mem_objects ~= image.cl_handle;
    CLAssert(err, "Creating image ");
    image.data.length = width*height*4;
    CLAssert(clSetKernelArg(kernel, ind, CLMem.sizeof, &image.cl_handle),
             "Set kernel");
    return image;
  }

  OpenCLBuffer!T Set_Buffer(T)(BufferType type, T[] data, int ind) {
    OpenCLBuffer!T buffer = OpenCLBuffer!T(
      type, data.dup
    );
    buffer.cl_handle = clCreateBuffer(context, type | CL_MEM_COPY_HOST_PTR,
                             data.length, buffer.data.ptr, &err);
    buffer.cl_length_handle = clCreateBuffer(context,
                             CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             buffer.length.sizeof, &buffer.length, &err);
    mem_objects ~= buffer.cl_handle;
    mem_objects ~= buffer.cl_length_handle;
    CLAssert(err, "Creating buffer");
    CLAssert(clSetKernelArg(kernel, ind, CLMem.sizeof, &buffer.cl_handle),
             "Set kernel");
    CLAssert(clSetKernelArg(kernel, ind, CLMem.sizeof, &buffer.cl_length_handle),
             "Set kernel length");
    return buffer;
  }

  void Run(size_t[] global, size_t[] local) {
    assert(global.length == local.length);
    CLAssert(clEnqueueNDRangeKernel(command_queue, kernel,
             cast(uint)global.length, null, global.ptr, local.ptr,
             0, null, null), "EnqueueNDRangeKernel");
  }

  void Read_Buffer(BType)(ref OpenCLBuffer buffer) {
    CLAssert(clEnqueueReadBuffer(command_queue, buffer.cl_handle, CL_TRUE, 0,
                  BType.sizeof*buffer.data.length, buffer.data.ptr,
                  0, null, null), "Enqueue read buffer");
  }

  void Read_Image(ref OpenCLImage image) {
    static size_t[3] origin = [0, 0, 0];
    size_t[3] region = [
      image.width, image.height, 1
    ];
    CLAssert(clEnqueueReadImage(command_queue, image.cl_handle, CL_TRUE,
                origin.ptr, region.ptr, 0, 0, image.data.ptr, 0, null, null),
                "Enqueue image");
  }

  void Cleanup() {
    if ( program == null ) return;
    foreach ( mem; mem_objects )
      clReleaseMemObject(mem);
    mem_objects = [];
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(command_queue);
    clReleaseContext(context);
    program = null;
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
  CLAssert(clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, 1, &device_id, null),
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
  return new OpenCLProgram(source);
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
      return CLVersion.CL12;
      // assert(0, version_error ~ " (no D bindings :[)");
  }
}

/*
khronos.org/registry/OpenCL/sdk/2.1/docs/man/xhtml/clCreateFromGLTexture.html
describes corresponding opencl and opengl image formats
*/
struct CLImage {
  ubyte[] buffer;
  int width, height;
  this ( float[] _buffer, int _width, int _height )  {
    width = _width; height = _height;
    foreach ( i; _buffer ) {
      buffer ~= cast(ubyte)(i*255.0f);
    }
  }
}

CLImage Test_OpenCL ( ) {
 return CLImage();
  // CLAssert(err, "Create image input");
  // auto image  = clCreateImage(program.context, CL_MEM_WRITE_ONLY,
  //                 &Img_format, &Img_descriptor, null, &err);
  // CLAssert(err, "Create image output");

  // // load data into input buffer
  // size_t[3] origin = [0, 0, 0],
  //           region = [dim, dim, 1];
  // CLAssert(clEnqueueWriteImage(program.command_queue, input, CL_TRUE,
  //                  origin.ptr, region.ptr, 0, 0, input_image_buffer.ptr, 0,
  //                  null, null),
  //   "enqueue write image");
  // // set arg list for kernel command
  // CLAssert(clSetKernelArg(kernel, 0, cl_mem.sizeof, &input),
  //          "Setting kernel of input");
  // CLAssert(clSetKernelArg(kernel, 1, cl_mem.sizeof, &image),
  //          "Setting kernel of image");
  // size_t[3] global = [dim, dim, 1];
  // size_t[3] local  = [1,  1, 1 ];

  // import std.datetime;
  // writeln("Starting");
  // StopWatch watch;
  // watch.start;
  // // enqueue kernel command for execution
  // CLAssert(clEnqueueNDRangeKernel(program.command_queue, kernel, 3, null,
  //                                 global.ptr, local.ptr, 0, null, null),
  //          "Enqueue kernel command for execution");
  // CLAssert(clFinish(program.command_queue),
  //          "Finish execution");
  // watch.stop;
  // import core.time;
  // writeln("Finished, duration: ", watch.peek().msecs, " milliseconds");
  // float[] image_buffer;
  // image_buffer.length = dim*dim*4;
  // CLAssert(clEnqueueReadImage(program.command_queue, image, CL_TRUE,
  //           origin.ptr, region.ptr, 0, 0, image_buffer.ptr, 0, null, null),
  //          "Enqueue image");
  // auto img = CLImage(image_buffer.dup, dim, dim);
  // // writeln("RESULTS: ", image_buffer);
  // // cleanup opencl resources
  // clReleaseMemObject(input);
  // clReleaseMemObject(image);
  // clReleaseProgram(program.program);
  // clReleaseKernel(kernel);
  // clReleaseCommandQueue(program.command_queue);
  // clReleaseContext(program.context);
  // return img;
}
