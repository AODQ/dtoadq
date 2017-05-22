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
  int index;
}

struct OpenCLImage {
  BaseBuffer!float _mybuffer;
  alias _mybuffer this;
  cl_image_format  format;
  cl_image_desc    description;
  size_t width, height;
  this ( BufferType _type, float[] _data, int _ind, cl_image_format _format,
         cl_image_desc _description, size_t _width, size_t _height) {
    _mybuffer = BaseBuffer!float(_type, CLMem(), _data, _ind);
    format = _format;
    description = _description;
    width = _width; height = _height;
  }
}

struct OpenCLSingleton(BType) {
  BaseBuffer!BType _mybuffer;
  alias _mybuffer this;
  this ( BufferType _type, BType _data, int _ind ) {
    _mybuffer = BaseBuffer!BType(_type, CLMem(), [_data], _ind);
  }
}


OpenCLImage Create_CL_Image(BufferType type, int width, int height,
                                                int param_index) {
  return OpenCLImage(
    type, [], param_index,
    cl_image_format(CL_RGBA, CL_FLOAT),
    cl_image_desc (CL_MEM_OBJECT_IMAGE2D, width, height, 1,
                  0, 0, 0, 0, 0, null),
    width, height
  );
}

struct OpenCLBuffer(BType) {
  BaseBuffer!BType _mybuffer;
  alias _mybuffer this;
  uint length;
  CLMem cl_length_handle;
  this ( BufferType _type, BType[] _data, int _ind ) {
    _mybuffer = BaseBuffer!BType(_type, CLMem(), _data, _ind);
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
  int param_count;
public:
  this ( string source ) in {
    assert(device !is null);
  } body {
    err = CL_SUCCESS;
    properties = [CL_CONTEXT_PLATFORM, cast(int)device.platform_id, 0];
    int err;
    writeln("create context");
    writeln(device.platform_id, " 00 ", device.device_id);
    context = clCreateContext(properties.ptr, 1, &device.device_id,
                                null, null, &err);
    CLAssert(err, "Create context");
    command_queue = clCreateCommandQueue(context, device.device_id, 0, &err);
    CLAssert(err, "Create command queue");

    import std.string;
    auto e = source.toStringz;
    program = clCreateProgramWithSource(context, 1, &e, null, &err);
    CLAssert(err, "Create program with source");
    writeln("building");
    if ( clBuildProgram(program, 0, null, null, null, null) != CL_SUCCESS ) {
      writeln("Error building program");
      size_t len;
      clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, 0,
                            null, &len);
      char[] log; log.length = len;
      clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG,
                            len, log.ptr, null);
      writeln("LOG: ", log);
      assert(false);
    }
    writeln("PROGRAM BUILT .. .");
  }
  ~this() {
    Clean_Up();
  }

  void Set_Kernel(string kernel_name) {
    kernel = clCreateKernel(program, kernel_name.ptr, &err);
    CLAssert(err, "clCreateKernel");
  }

  OpenCLImage Set_Image_Buffer(BufferType type, int dim) {
    return Set_Image_Buffer(type, dim, dim);
  }
  OpenCLImage Set_Image_Buffer(BufferType type, int width, int height){
    OpenCLImage image = Create_CL_Image(type, width, height, param_count);
    image.cl_handle = clCreateImage(context, type, &image.format,
                                    &image.description, null, &err);
    mem_objects ~= image.cl_handle;
    CLAssert(err, "Creating image ");
    image.data.length = width*height*4;
    {import functional; image.data.each!((ref n) => n = 0.0f);}
    CLAssert(clSetKernelArg(kernel, param_count, CLMem.sizeof,&image.cl_handle),
             "clSetKernelArg");
    image.index = param_count;
    ++param_count;
    return image;
  }

  void Modify_Image_Size ( ref OpenCLImage image, int width, int height ) {
    clReleaseMemObject(image.cl_handle);
    image.width = width;
    image.height = height;
    image.description = cl_image_desc(CL_MEM_OBJECT_IMAGE2D, width, height, 1,
                                      0, 0, 0, 0, 0, null);
    image.cl_handle = clCreateImage(context, image.type, &image.format,
                                    &image.description, null, &err);
    CLAssert(err, "Modifying image size ");
    image.data.length = width*height*4;

    CLAssert(clSetKernelArg(kernel, image.index, CLMem.sizeof,&image.cl_handle),
            "clSetKernelArg");
  }

  OpenCLSingleton!T Set_Singleton(T)(BufferType type, inout(T) data) {
    auto singleton = OpenCLSingleton!T(
      type, data, param_count
    );
    singleton.cl_handle = clCreateBuffer(context, type | CL_MEM_COPY_HOST_PTR,
                                         T.sizeof, singleton.data.ptr, &err);
    CLAssert(err, "clCreateBuffer");
    CLAssert(clSetKernelArg(kernel, param_count, CLMem.sizeof,
                            &singleton.cl_handle),
             "clSetKernelArg");
    ++param_count;
    return singleton;
  }

  void Write(OpenCLImage image) {
    size_t[3] origin = [0, 0, 0];
    size_t[3] region = [image.width, image.height, 1];
    CLAssert(clEnqueueWriteImage(command_queue, image.cl_handle,
                  CL_TRUE, origin.ptr, region.ptr, 0, 0,
                  cast(void*)image.data.ptr, 0, null, null),
             "clEnqueueWriteImage");
  }
  void Write(T)(OpenCLBuffer!T buffer) {
    int len = cast(int)buffer.data.length;
    CLAssert(clEnqueueWriteBuffer(command_queue, buffer.cl_handle,
              CL_TRUE, 0, T.sizeof*len, &buffer.data[0],
              0, null, null), "clEnqueueWriteBuffer");
    CLAssert(clEnqueueWriteBuffer(command_queue, buffer.cl_length_handle,
              CL_TRUE, 0, int.sizeof, &len,
              0, null, null), "clEnqueueWriteBuffer length");
  }
  void Write(T)(OpenCLSingleton!T singleton) in {
    assert(singleton.data.length == 1,
            "Length mismatch on singleton: " ~ singleton.data.length.to!string);
  } body {
    CLAssert(clEnqueueWriteBuffer(command_queue, singleton.cl_handle,
                  CL_TRUE, 0, T.sizeof, &singleton.data[0],
                  0, null, null), "clEnqueueWriteBuffer singleton");
  }

  OpenCLBuffer!T Set_Buffer(T)(BufferType type, inout(T[]) data) {
    OpenCLBuffer!T buffer = OpenCLBuffer!T(
      type, data.dup, param_count
    );

    auto flags = type | (data is null? 0 : CL_MEM_COPY_HOST_PTR);
    buffer.cl_handle = clCreateBuffer(context, flags, data.length*T.sizeof,
                                      buffer.data.ptr, &err);
    int temp = cast(int)data.length;
    buffer.cl_length_handle = clCreateBuffer(context,
                             CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             int.sizeof, &temp, &err);
    mem_objects ~= buffer.cl_handle;
    mem_objects ~= buffer.cl_length_handle;
    CLAssert(err, "clCreateBufferl");
    CLAssert(clSetKernelArg(kernel, param_count, CLMem.sizeof,
                            &buffer.cl_handle),
             "clSetKernelArg");
    ++ param_count;
    CLAssert(clSetKernelArg(kernel, param_count, CLMem.sizeof,
             &buffer.cl_length_handle), "Set kernel length");
    ++ param_count;
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

  void Clean_Up() {
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
  assert(cond == CL_SUCCESS, err ~ ": " ~ CL_Error_String(cond));
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
    write("CHOOSE A PLATFORM: ");
    index = readln.chomp.to!int;
    writeln();
    assert(index >= 0 && index < platforms.length, " invalid platform index");
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
  writeln("det load");
  DerelictCL.load();
  writeln("det platform");
  auto platform_id = Set_Current_Platform();
  writeln("det device id");
  auto device_id = Set_Current_Device(platform_id);
  writeln("det reload");
  DerelictCL.reload(RCL_Version(platform_id));
  // DerelictCL.loadEXT(platform_id);
  writeln("det device");
  device = new Device(platform_id, device_id);
  // writeln(RDevice_Info(device_id));
  writeln("done init");
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
  float[] buffer;
  int width, height;
  this ( float[] _buffer, int _width, int _height )  {
    width = _width; height = _height;
    import functional;
    buffer = _buffer;
  }
}


/** conversion helpers */

/** */
auto To_CLFloat3(float[3] a) {
  cl_float3 vec;
  vec.x = a[0]; vec.y = a[1]; vec.z = a[2];
  return vec;
}

auto To_CLInt8(uint[8] a) {
  cl_int8 vec;
  foreach ( i; 0 .. 8 )
    vec[i] = a[i];
  return vec;
}
