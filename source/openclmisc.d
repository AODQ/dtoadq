module openclmisc;
import derelict.opencl.cl;
import globals;
import std.conv : to;

auto RDevice_Info(cl_device_info info, cl_device_id device_id) {
  size_t dv_size = 5000;
  char[] dv;
  dv.length = 5000;
  clGetPlatformInfo(device_id, info, dv_size, dv.ptr, &dv_size);
  dv.length = dv_size;
  return dv.to!string;
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
  return pv.to!string;
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
    else                  output ~= "Unknown " ~ forder.to!string;
    auto data  = fdata in channel_data_map;
    output ~= " : ";
    if ( data  !is null ) output ~= *data;
    else                  output ~= "Unknown " ~ fdata.to!string;
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
  switch ( error ) {
    // run-time and JIT compiler errors
    case 0:   return "CL_SUCCESS"; case -1:  return "CL_DEVICE_NOT_FOUND"; case -2:  return "CL_DEVICE_NOT_AVAILABLE"; case -3:  return "CL_COMPILER_NOT_AVAILABLE"; case -4:  return "CL_MEM_OBJECT_ALLOCATION_FAILURE"; case -5:  return "CL_OUT_OF_RESOURCES"; case -6:  return "CL_OUT_OF_HOST_MEMORY"; case -7:  return "CL_PROFILING_INFO_NOT_AVAILABLE"; case -8:  return "CL_MEM_COPY_OVERLAP"; case -9:  return "CL_IMAGE_FORMAT_MISMATCH"; case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED"; case -11: return "CL_BUILD_PROGRAM_FAILURE"; case -12: return "CL_MAP_FAILURE"; case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET"; case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST"; case -15: return "CL_COMPILE_PROGRAM_FAILURE"; case -16: return "CL_LINKER_NOT_AVAILABLE"; case -17: return "CL_LINK_PROGRAM_FAILURE"; case -18: return "CL_DEVICE_PARTITION_FAILED"; case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE"; 
    // compile-time errors
    case -30: return "CL_INVALID_VALUE"; case -31: return "CL_INVALID_DEVICE_TYPE"; case -32: return "CL_INVALID_PLATFORM"; case -33: return "CL_INVALID_DEVICE"; case -34: return "CL_INVALID_CONTEXT"; case -35: return "CL_INVALID_QUEUE_PROPERTIES"; case -36: return "CL_INVALID_COMMAND_QUEUE"; case -37: return "CL_INVALID_HOST_PTR"; case -38: return "CL_INVALID_MEM_OBJECT"; case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR"; case -40: return "CL_INVALID_IMAGE_SIZE"; case -41: return "CL_INVALID_SAMPLER"; case -42: return "CL_INVALID_BINARY"; case -43: return "CL_INVALID_BUILD_OPTIONS"; case -44: return "CL_INVALID_PROGRAM"; case -45: return "CL_INVALID_PROGRAM_EXECUTABLE"; case -46: return "CL_INVALID_KERNEL_NAME"; case -47: return "CL_INVALID_KERNEL_DEFINITION"; case -48: return "CL_INVALID_KERNEL"; case -49: return "CL_INVALID_ARG_INDEX"; case -50: return "CL_INVALID_ARG_VALUE"; case -51: return "CL_INVALID_ARG_SIZE"; case -52: return "CL_INVALID_KERNEL_ARGS"; case -53: return "CL_INVALID_WORK_DIMENSION"; case -54: return "CL_INVALID_WORK_GROUP_SIZE"; case -55: return "CL_INVALID_WORK_ITEM_SIZE"; case -56: return "CL_INVALID_GLOBAL_OFFSET"; case -57: return "CL_INVALID_EVENT_WAIT_LIST"; case -58: return "CL_INVALID_EVENT"; case -59: return "CL_INVALID_OPERATION"; case -60: return "CL_INVALID_GL_OBJECT"; case -61: return "CL_INVALID_BUFFER_SIZE"; case -62: return "CL_INVALID_MIP_LEVEL"; case -63: return "CL_INVALID_GLOBAL_WORK_SIZE"; case -64: return "CL_INVALID_PROPERTY"; case -65: return "CL_INVALID_IMAGE_DESCRIPTOR"; case -66: return "CL_INVALID_COMPILER_OPTIONS"; case -67: return "CL_INVALID_LINKER_OPTIONS"; case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";
    // extension errors
    case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR"; case -1001: return "CL_PLATFORM_NOT_FOUND_KHR"; case -1002: return "CL_INVALID_D3D10_DEVICE_KHR"; case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR"; case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR"; case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
    default: return "Unknown OpenCL error";
  }
}

string RDevice_Info(cl_device_id device_id) {
  size_t size;
  char[] str;
  int err;
  err = clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, 0, null, &size);
  if ( err != 0 ) {
    writeln("clGetDeviceInfo for size ERROR: ", err.CL_Error_String);
    assert(0);
  }
  str.length = size;
  err = clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, size, str.ptr, null);
  if ( err != 0 ) {
    writeln("clGetDeviceInfo ERROR: ", err.CL_Error_String);
    assert(0);
  }
  return str.to!string;
}
