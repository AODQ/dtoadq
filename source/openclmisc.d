module openclmisc;
import opencl;
import std.conv : to;

auto RDevice_Info(cl_device_info info, CLDeviceID device_id) {
  size_t dv_size = 5000;
  char[] dv;
  dv.length = 5000;
  clGetPlatformInfo(device_id, info, dv_size, dv.ptr, &dv_size);
  dv.length = dv_size;
  return dv.to!string;
}

/**
  from:
  stackoverflow.com/questions/24326432/convenient-way-to-show-opencl-error-codes
*/
string CL_Error_String(int error) {
  switch ( error ) {
    // run-time and JIT compiler errors
    case 0:   return "CL_SUCCESS";
    case -1:  return "CL_DEVICE_NOT_FOUND";
    case -2:  return "CL_DEVICE_NOT_AVAILABLE";
    case -3:  return "CL_COMPILER_NOT_AVAILABLE";
    case -4:  return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case -5:  return "CL_OUT_OF_RESOURCES";
    case -6:  return "CL_OUT_OF_HOST_MEMORY";
    case -7:  return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case -8:  return "CL_MEM_COPY_OVERLAP";
    case -9:  return "CL_IMAGE_FORMAT_MISMATCH";
    case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case -11: return "CL_BUILD_PROGRAM_FAILURE";
    case -12: return "CL_MAP_FAILURE";
    case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
    case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
    case -15: return "CL_COMPILE_PROGRAM_FAILURE";
    case -16: return "CL_LINKER_NOT_AVAILABLE";
    case -17: return "CL_LINK_PROGRAM_FAILURE";
    case -18: return "CL_DEVICE_PARTITION_FAILED";
    case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

    // compile-time errors
    case -30: return "CL_INVALID_VALUE";
    case -31: return "CL_INVALID_DEVICE_TYPE";
    case -32: return "CL_INVALID_PLATFORM";
    case -33: return "CL_INVALID_DEVICE";
    case -34: return "CL_INVALID_CONTEXT";
    case -35: return "CL_INVALID_QUEUE_PROPERTIES";
    case -36: return "CL_INVALID_COMMAND_QUEUE";
    case -37: return "CL_INVALID_HOST_PTR";
    case -38: return "CL_INVALID_MEM_OBJECT";
    case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case -40: return "CL_INVALID_IMAGE_SIZE";
    case -41: return "CL_INVALID_SAMPLER";
    case -42: return "CL_INVALID_BINARY";
    case -43: return "CL_INVALID_BUILD_OPTIONS";
    case -44: return "CL_INVALID_PROGRAM";
    case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
    case -46: return "CL_INVALID_KERNEL_NAME";
    case -47: return "CL_INVALID_KERNEL_DEFINITION";
    case -48: return "CL_INVALID_KERNEL";
    case -49: return "CL_INVALID_ARG_INDEX";
    case -50: return "CL_INVALID_ARG_VALUE";
    case -51: return "CL_INVALID_ARG_SIZE";
    case -52: return "CL_INVALID_KERNEL_ARGS";
    case -53: return "CL_INVALID_WORK_DIMENSION";
    case -54: return "CL_INVALID_WORK_GROUP_SIZE";
    case -55: return "CL_INVALID_WORK_ITEM_SIZE";
    case -56: return "CL_INVALID_GLOBAL_OFFSET";
    case -57: return "CL_INVALID_EVENT_WAIT_LIST";
    case -58: return "CL_INVALID_EVENT";
    case -59: return "CL_INVALID_OPERATION";
    case -60: return "CL_INVALID_GL_OBJECT";
    case -61: return "CL_INVALID_BUFFER_SIZE";
    case -62: return "CL_INVALID_MIP_LEVEL";
    case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
    case -64: return "CL_INVALID_PROPERTY";
    case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
    case -66: return "CL_INVALID_COMPILER_OPTIONS";
    case -67: return "CL_INVALID_LINKER_OPTIONS";
    case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";

    // extension errors
    case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
    case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
    case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
    case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
    case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
    case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
    default: return "Unknown OpenCL error";
  }
}

auto RDevice_Info(CLDeviceID device_id) {
  import std.string;
  return `
    Device Address Bits: %s
    Device Available: %s
    Device Compiler Available: %s
    Device Double Fp Config: %s
    Device Endian Little: %s
    Device Error Correction Support: %s
    Device Execution Capabilities: %s
    Device Extensions: %s
    Device Global Mem Cache Size: %s
    Device Global Mem Cache Type: %s
    Device Global Mem Cacheline Size: %s
    Device Global Mem Size: %s
    Device Half Fp Config: %s
    Device Image Support: %s
    Device Image2D Max Height: %s
    Device Image2D Max Width: %s
    Device Image3D Max Depth: %s
    Device Image3D Max Height: %s
    Device Image3D Max Width: %s
    Device Local Mem Size: %s
    Device Local Mem Type: %s
    Device Max Clock Frequency: %s
    Device Max Compute Units: %s
    Device Max Constant Args: %s
    Device Max Constant Buffer Size: %s
    Device Max Mem Alloc Size: %s
    Device Max Parameter Size: %s
    Device Max Read Image Args: %s
    Device Max Samplers: %s
    Device Max Work Group Size: %s
    Device Max Work Item Dimensions: %s
    Device Max Work Item Sizes: %s
    Device Max Write Image Args: %s
    Device Mem Base Addr Align: %s
    Device Min Data Type Align Size: %s
    Device Name: %s
    Device Platform: %s
    Device Preferred Vector Width Char: %s
    Device Preferred Vector Width Short: %s
    Device Preferred Vector Width Int: %s
    Device Preferred Vector Width Long: %s
    Device Preferred Vector Width Float: %s
    Device Preferred Vector Width Double: %s
    Device Profile: %s
    Device Profiling Timer Resolution: %s
    Device Queue Properties: %s
    Device Single Fp Config: %s
    Device Type: %s
    Device Vendor: %s
    Device Vendor ID: %s
    Device Version: %s
    Driver Version %s
  `.format(
    RDevice_Info(CL_DEVICE_ADDRESS_BITS,                  device_id),
    RDevice_Info(CL_DEVICE_AVAILABLE,                     device_id),
    RDevice_Info(CL_DEVICE_COMPILER_AVAILABLE,            device_id),
    RDevice_Info(CL_DEVICE_DOUBLE_FP_CONFIG,              device_id),
    RDevice_Info(CL_DEVICE_ENDIAN_LITTLE,                 device_id),
    RDevice_Info(CL_DEVICE_ERROR_CORRECTION_SUPPORT,      device_id),
    RDevice_Info(CL_DEVICE_EXECUTION_CAPABILITIES,        device_id),
    RDevice_Info(CL_DEVICE_EXTENSIONS,                    device_id),
    RDevice_Info(CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,         device_id),
    RDevice_Info(CL_DEVICE_GLOBAL_MEM_CACHE_TYPE,         device_id),
    RDevice_Info(CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE,     device_id),
    RDevice_Info(CL_DEVICE_GLOBAL_MEM_SIZE,               device_id),
    RDevice_Info(CL_DEVICE_HALF_FP_CONFIG,                device_id),
    RDevice_Info(CL_DEVICE_IMAGE_SUPPORT,                 device_id),
    RDevice_Info(CL_DEVICE_IMAGE2D_MAX_HEIGHT,            device_id),
    RDevice_Info(CL_DEVICE_IMAGE2D_MAX_WIDTH,             device_id),
    RDevice_Info(CL_DEVICE_IMAGE3D_MAX_DEPTH,             device_id),
    RDevice_Info(CL_DEVICE_IMAGE3D_MAX_HEIGHT,            device_id),
    RDevice_Info(CL_DEVICE_IMAGE3D_MAX_WIDTH,             device_id),
    RDevice_Info(CL_DEVICE_LOCAL_MEM_SIZE,                device_id),
    RDevice_Info(CL_DEVICE_LOCAL_MEM_TYPE,                device_id),
    RDevice_Info(CL_DEVICE_MAX_CLOCK_FREQUENCY,           device_id),
    RDevice_Info(CL_DEVICE_MAX_COMPUTE_UNITS,             device_id),
    RDevice_Info(CL_DEVICE_MAX_CONSTANT_ARGS,             device_id),
    RDevice_Info(CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,      device_id),
    RDevice_Info(CL_DEVICE_MAX_MEM_ALLOC_SIZE,            device_id),
    RDevice_Info(CL_DEVICE_MAX_PARAMETER_SIZE,            device_id),
    RDevice_Info(CL_DEVICE_MAX_READ_IMAGE_ARGS,           device_id),
    RDevice_Info(CL_DEVICE_MAX_SAMPLERS,                  device_id),
    RDevice_Info(CL_DEVICE_MAX_WORK_GROUP_SIZE,           device_id),
    RDevice_Info(CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,      device_id),
    RDevice_Info(CL_DEVICE_MAX_WORK_ITEM_SIZES,           device_id),
    RDevice_Info(CL_DEVICE_MAX_WRITE_IMAGE_ARGS,          device_id),
    RDevice_Info(CL_DEVICE_MEM_BASE_ADDR_ALIGN,           device_id),
    RDevice_Info(CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE,      device_id),
    RDevice_Info(CL_DEVICE_NAME,                          device_id),
    RDevice_Info(CL_DEVICE_PLATFORM,                      device_id),
    RDevice_Info(CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,   device_id),
    RDevice_Info(CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,  device_id),
    RDevice_Info(CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT,    device_id),
    RDevice_Info(CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,   device_id),
    RDevice_Info(CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,  device_id),
    RDevice_Info(CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, device_id),
    RDevice_Info(CL_DEVICE_PROFILE,                       device_id),
    RDevice_Info(CL_DEVICE_PROFILING_TIMER_RESOLUTION,    device_id),
    RDevice_Info(CL_DEVICE_QUEUE_PROPERTIES,              device_id),
    RDevice_Info(CL_DEVICE_SINGLE_FP_CONFIG,              device_id),
    RDevice_Info(CL_DEVICE_TYPE,                          device_id),
    RDevice_Info(CL_DEVICE_VENDOR,                        device_id),
    RDevice_Info(CL_DEVICE_VENDOR_ID,                     device_id),
    RDevice_Info(CL_DEVICE_VERSION,                       device_id),
    RDevice_Info(CL_DRIVER_VERSION,                       device_id));
}
