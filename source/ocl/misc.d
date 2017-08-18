module ocl.misc;
import derelict.opencl.cl;
static import stl;
import std.conv : to;

auto RDevice_Info(cl_device_info info, cl_device_id device_id) {
  size_t dv_size = 5000;
  char[] dv;
  dv.length = 5000;
  clGetPlatformInfo(device_id, info, dv_size, dv.ptr, &dv_size);
  dv.length = dv_size;
  return stl.to!string(dv);
}

enum PlatformInfo {
  Profile    = 0x0901, Version    = 0x0902,
  Name       = 0x0903, Vendor     = 0x0904,
  Extensions = 0x0905
}

auto RPlatform_Info(PlatformInfo info, cl_platform_id platform) {
  size_t pv_size = 50;
  void[] pv;
  pv.length = 50;
  clGetPlatformInfo(platform, info, pv_size, pv.ptr, &pv_size);
  pv.length = pv_size;
  return stl.to!string(pv);
}

auto RSupported_Image_Formats(ref cl_context context) {
  import std.string;
  uint if_size;
  clGetSupportedImageFormats(context, CL_MEM_WRITE_ONLY,
      CL_MEM_OBJECT_IMAGE2D, 0, null, &if_size);
  cl_image_format[] imgforms;
  imgforms.length = if_size;
  clGetSupportedImageFormats(context, CL_MEM_WRITE_ONLY,
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
    else                  output ~= "Unknown " ~ stl.to!string(forder);
    auto data  = fdata in channel_data_map;
    output ~= " : ";
    if ( data  !is null ) output ~= *data;
    else                  output ~= "Unknown " ~ stl.to!string(fdata);
    output ~= ", ";
  }
  return output;
}

auto RPlatform_Info(cl_platform_id platform) {
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

CLVersion RCL_Version(cl_platform_id platform_id) {
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

/**
  from:
  stackoverflow.com/questions/24326432/convenient-way-to-show-opencl-error-codes
*/
string CL_Error_String(int error) {
  string[int] err_map = [
    // run-time and JIT compiler errors
    0: "SUCCESS",                        -1: "DEVICE_NOT_FOUND",
    -2: "DEVICE_NOT_AVAILABLE",          -3: "COMPILER_NOT_AVAILABLE",
    -4: "MEM_OBJECT_ALLOCATION_FAILURE", -5: "OUT_OF_RESOURCES",
    -6: "OUT_OF_HOST_MEMORY",            -7: "PROFILING_INFO_NOT_AVAILABLE",
    -8: "MEM_COPY_OVERLAP",              -9: "IMAGE_FORMAT_MISMATCH",
    -10: "IMAGE_FORMAT_NOT_SUPPORTED",   -11: "BUILD_PROGRAM_FAILURE",
    -12: "MAP_FAILURE",                  -13: "MISALIGNED_SUB_BUFFER_OFFSET",
    -14: "EVENTS_IN_WAIT_LIST",          -15: "COMPILE_PROGRAM_FAILURE",
    -16: "LINKER_NOT_AVAILABLE",         -17: "LINK_PROGRAM_FAILURE",
    -18: "DEVICE_PARTITION_FAILED",      -19: "KERNEL_ARG_INFO_NOT_AVAILABLE",

    // compile-time errors
    -30: "INVALID_VALUE",              -31: "INVALID_DEVICE_TYPE",
    -32: "INVALID_PLATFORM",           -33: "INVALID_DEVICE",
    -34: "INVALID_CONTEXT",            -35: "INVALID_QUEUE_PROPERTIES",
    -36: "INVALID_COMMAND_QUEUE",      -37: "INVALID_HOST_PTR",
    -38: "INVALID_MEM_OBJECT",         -39: "INVALID_IMAGE_FORMAT_DESCRIPTOR",
    -40: "INVALID_IMAGE_SIZE",         -41: "INVALID_SAMPLER",
    -42: "INVALID_BINARY",             -43: "INVALID_BUILD_OPTIONS",
    -44: "INVALID_PROGRAM",            -45: "INVALID_PROGRAM_EXECUTABLE",
    -46: "INVALID_KERNEL_NAME",        -47: "INVALID_KERNEL_DEFINITION",
    -48: "INVALID_KERNEL",             -49: "INVALID_ARG_INDEX",
    -50: "INVALID_ARG_VALUE",          -51: "INVALID_ARG_SIZE",
    -52: "INVALID_KERNEL_ARGS",        -53: "INVALID_WORK_DIMENSION",
    -54: "INVALID_WORK_GROUP_SIZE",    -55: "INVALID_WORK_ITEM_SIZE",
    -56: "INVALID_GLOBAL_OFFSET",      -57: "INVALID_EVENT_WAIT_LIST",
    -58: "INVALID_EVENT",              -59: "INVALID_OPERATION",
    -60: "INVALID_GL_OBJECT",          -61: "INVALID_BUFFER_SIZE",
    -62: "INVALID_MIP_LEVEL",          -63: "INVALID_GLOBAL_WORK_SIZE",
    -64: "INVALID_PROPERTY",           -65: "INVALID_IMAGE_DESCRIPTOR",
    -66: "INVALID_COMPILER_OPTIONS",   -67: "INVALID_LINKER_OPTIONS",
    -68: "INVALID_DEVICE_PARTITION_COUNT",
    // extension errors
    -1000: "INVALID_GL_SHAREGROUP_REFERENCE_KHR",
    -1001: "PLATFORM_NOT_FOUND_KHR",
    -1002: "INVALID_D3D10_DEVICE_KHR",
    -1003: "INVALID_D3D10_RESOURCE_KHR",
    -1004: "D3D10_RESOURCE_ALREADY_ACQUIRED_KHR",
    -1005: "D3D10_RESOURCE_NOT_ACQUIRED_KHR",
  ];
  if ( error !in err_map ) return "Unknown OpenCL error";
  return err_map[error];
}

string RDevice_Info(cl_device_id device_id) {
  size_t size;
  char[] str;
  int err;
  err = clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, 0, null, &size);
  if ( err != 0 ) {
    stl.writeln("clGetDeviceInfo for size ERROR: ", err.CL_Error_String);
    assert(0);
  }
  str.length = size;
  err = clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, size, str.ptr, null);
  if ( err != 0 ) {
    stl.writeln("clGetDeviceInfo ERROR: ", err.CL_Error_String);
    assert(0);
  }
  return stl.to!string(str);
}
