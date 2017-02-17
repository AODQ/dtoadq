module opencl;

public import derelict.opencl.cl;
public static import wrap = clwrap.clwrap;
import clwrap.clwrap, std.stdio;
import std.conv : to;

void OpenCLAssert(bool condition, string error) {
  assert(condition, error);
}

auto RPlatforms() {
  platform_id[] platforms;
  cl_uint platform_amt = 10;
  platforms.length = platform_amt;
  getPlatformIDs(platform_amt, platforms.ptr, &platform_amt);
  platforms.length = platform_amt;
  return platforms;
}

struct Device {
  cl_platform_id platform_id;
  cl_device_id   device_id;
}

auto RPlatform() {
  auto platforms = RPlatforms;
  OpenCLAssert(platforms.length >= 1, "No OpenCL platform found");
  cl_platform_id platform_id = cast(void*)platforms[0];
  cl_device_id device_id;
  writeln("GETTING device ID stuff");
  clGetDeviceIDs(&platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, null);
  writeln("PLATFORM ID: ", platform_id);
  writeln("DEVICE ID: ", device_id);
  return Device(platform_id, device_id);
}

struct CLProgram {
  cl_context context;
  cl_context_properties[] properties;
  cl_command_queue command_queue;
  cl_program program;
  cl_device_id device_id;
  cl_platform_id platform_id;
  cl_kernel kernel;
}

auto Compile_Program ( string source, Device device) {
  cl_context_properties[] properties = [
    CL_CONTEXT_PLATFORM,
    cast(cl_context_properties) device.platform_id,
    0
  ];

  auto context =
    clCreateContext(properties.ptr, 1, &device.device_id, null, null, null);

  auto command_queue = clCreateCommandQueue(context, device.device_id, 0, null);

  import std.string;
  auto e = source.toStringz;
  auto program = clCreateProgramWithSource(context, 1, &e, null, null);
  OpenCLAssert(
    clBuildProgram(program, 0, null, null, null, null) != CL_SUCCESS,
    "Error building program");
  return CLProgram(
    context, properties, command_queue, program,
    device.device_id, device.platform_id
  );
}

auto Compile(string source) {
  return Compile_Program(source, RPlatform);
}

alias PlatformID = platform_id;

auto RPlatform_Info(PlatformID platform, platform_info info) {
  size_t pv_size = 50;
  void[] pv;
  pv.length = 50;
  getPlatformInfo(platform, info, pv_size, pv.ptr, &pv_size);
  pv.length = pv_size;
  return pv.to!string;
}

auto RPlatform_Info(PlatformID platform) {
  import std.string;
  return
    `Platform ID: %s
     OpenCL Profile: %s
     Platform Version: %s
     Platforn Name: %s
     Platform Vendor: %s
     Platform Extensions: %s
    `
     .format(platform.to!string,
       RPlatform_Info(platform, platform_info(0x0900)),
       RPlatform_Info(platform, platform_info(0x0901)),
       RPlatform_Info(platform, platform_info(0x0902)),
       RPlatform_Info(platform, platform_info(0x0903)),
       RPlatform_Info(platform, platform_info(0x0904))
       );
}

void Print_Device_Information ( ) {
  auto platforms = RPlatforms;
  writeln("PLATFORMS: ", platforms);
  import std.algorithm;
  platforms.each!(n => n.RPlatform_Info.writeln);
}


void Initialize ( ) {
  DerelictCL.load();
}


void Test_OpenCL ( ) {
  float[10] input_data = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ];

  auto program = Compile(`
    __kernel void hello(__global float* input, __global float* output) {
      size_t id = get_global_id(0);
      output[id] = input[id] * input[id];
    }
  `);

  import std.string;
  auto e = "hello".toStringz;
  auto size = float.sizeof*input_data.length;
  auto kernel = clCreateKernel(program.program, "hello", null);
  auto input  = clCreateBuffer(program.context, CL_MEM_READ_ONLY,
                    size, null, null);
  auto output = clCreateBuffer(program.context, CL_MEM_WRITE_ONLY,
                    size, null, null);

  // load data into input buffer
  clEnqueueWriteBuffer(program.command_queue, input, CL_TRUE, 0,
                    size, input_data.ptr, 0, null, null);
  // set arg list for kernel command
  clSetKernelArg(kernel, 0, cl_mem.sizeof, &input);
  clSetKernelArg(kernel, 1, cl_mem.sizeof, &output);
  size_t global = input_data.length;

  // enqueue kernel command for execution
  clEnqueueNDRangeKernel(program.command_queue, kernel, 1, null, &global, null,
                         0, null, null);
  clFinish(program.command_queue);

  // copy results from out of output buffer
  float[10] results;
  clEnqueueReadBuffer(program.command_queue, output, CL_TRUE, 0, size,
                      results.ptr, 0, null, null);
  writeln("Output: ");
  import std.algorithm, std.range;
  results.each!writeln;

  // cleanup opencl resources
  clReleaseMemObject(input);
  clReleaseMemObject(output);
  clReleaseProgram(program.program);
  clReleaseKernel(kernel);
  clReleaseCommandQueue(program.command_queue);
  clReleaseContext(program.context);
}
