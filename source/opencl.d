module opencl;
public import derelict.opencl.cl;
import std.stdio;
import openclmisc;

alias cl_context              CLContext;
alias cl_context_properties[] CLContextProperties;
alias cl_command_queue        CLCommandQueue;
alias cl_program              CLProgram;
alias cl_device_id            CLDeviceID;
alias cl_platform_id          CLPlatformID;
alias cl_kernel               CLKernel;
alias cl_mem                  CLMem;

struct CLPredefinedMem { // so memory can be allocated to OpenCL before Run
  CLMem _mem;
  alias _mem this;

  this ( CLMem mem_ ) {
    _mem = mem_;
  }
}

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

private CLContext           context;
private CLContextProperties properties;
private CLCommandQueue      command_queue;
private CLProgram           program;
private CLKernel            kernel;
private CLMem[]             mem_objects;
private int                 err;
private string              kernel_name;


private void Initialize_Kernel ( string kernel_name_ ) {
  kernel_name = kernel_name_;
  {
    import derelict.opengl3.glx;
    properties = [CL_CONTEXT_PLATFORM, cast(int)device.platform_id,
                  CL_GL_CONTEXT_KHR,   cast(int)glXGetCurrentContext(),
                  CL_GLX_DISPLAY_KHR,  cast(int)glXGetCurrentDisplay(),
                  0];
  }
  context = clCreateContext(properties.ptr, 1, &device.device_id,
                            null, null, &err);
  CLAssert(err, "clCreateContext");
  command_queue = clCreateCommandQueue(context, device.device_id, 0, &err);
  CLAssert(err, "clCreateCommandQueue");
}

void Compile ( string kernel ) {
  import std.string : toStringz;
  auto cstr_kernel = kernel.toStringz;
  program = clCreateProgramWithSource(context, 1, &cstr_kernel, null, &err);
  CLAssert(err, "clCreateProgramWithSource");
  if ( clBuildProgram(program, 0, null, null, null, null) != CL_SUCCESS ) {
      writeln("clBuildProgram");
      size_t len;
      clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, 0,
                            null, &len);
      char[] log; log.length = len;
      clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG,
                            len, log.ptr, null);
      writeln("LOG: ", log);
      assert(false);
  }
}

auto Create_CLGL_Texture ( uint gl_texture ) {
  import derelict.opengl3.gl3;
  auto cl_handle = clCreateFromGLTexture(context, CL_MEM_READ_WRITE,
                                GL_TEXTURE_2D, 0, gl_texture, &err);
  CLAssert(err, "clCreateFromGLTexture");
  return cl_handle;
}
void Lock_CLGL_Image ( CLPredefinedMem cl_handle ) {
  clEnqueueAcquireGLObjects(command_queue, 1, &image.cl_handle,
                            0, null, null);
}

void Unlock_CLGL_Image ( CLPredefinedMem  cl_handle ) {
    clEnqueueReleaseGLObjects(command_queue, 1, &image.cl_handle,
                              0, null, null);
}

void Flush ( ) { clFlush(null); }

private auto Load_Argument(T)(T data, int it) {
  import std.traits;
  static auto flgs = CL_MEM_READ_WRITE;

  CLMem cl_handle;
  static if ( isArray!(T) )
    cl_handle = clCreateBuffer(context, flgs, data.length*T.sizeof, &data[0],
                               &err);
  else static if ( is(T == CLPredefinedMem) ) {
    cl_handle = data;
    err = 0;
  } else
    cl_handle = clCreateBuffer(context, flgs, T.sizeof, &data, &err);
  CLAssert(err, "clCreateBuffer");

  static if ( isArray!(T) )
    CLAssert(clEnqueueWriteBuffer(command_queue, cl_handle, CL_TRUE, 0,
                        data.length*T.sizeof, &data[0], 0, null, null),
             "clEnqueueWriteBuffer");
  else static if ( !is(T == CLPredefinedMem) )
    CLAssert(clEnqueueWriteBuffer(command_queue, cl_handle, CL_TRUE, 0,
                                   T.sizeof, &data,     0, null, null),
             "clEnqueueWriteBuffer");
  // don't do writing for a predefinedclmem
  CLAssert(clSetKernelArg(kernel, it, CLMem.sizeof, cl_handle),
           "clSetKernelArg");
  return cl_handle;
}

void Run ( T... ) ( T params, size_t X, size_t Y ) {
  CLMem[] cl_handles_dealloc;

  foreach ( it, ref p; params ) {
    CLMem cl_handle;
    cl_handle = Load_Argument(p, it);
    static if ( !is(T == CLPredefinedMem) )
      cl_handles_dealloc ~= cl_handle;
  }

  ulong[] global = [cast(ulong)X, cast(ulong)Y, 1];
  ulong[] local  = [1, 1, 1];
  CLAssert(clEnqueueNDRangeKernel(command_queue, kernel,
            cast(uint)global.length, null, global.ptr, local.ptr,
            0, null, null), "clEnqueueNDRangeKernel");

  foreach ( mem; cl_handles_dealloc ) {
    clReleaseMemObject(mem);
  }
}

void Clean_Up() {
  clReleaseProgram(program);
  clReleaseKernel(kernel);
  clReleaseCommandQueue(command_queue);
  clReleaseContext(context);
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
      writeln("----- Index: ", it, "\n", RPlatform_Info(platforms[it]));
    }
    write("CHOOSE A PLATFORM: ");
    // index = readln.chomp.to!int;
    index = 0;
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

void Initialize ( string kernel_name ) {
  DerelictCL.load();
  auto platform_id = Set_Current_Platform();
  auto device_id = Set_Current_Device(platform_id);
  DerelictCL.reload(RCL_Version(platform_id));
  DerelictCL.loadEXT(platform_id);
  device = new Device(platform_id, device_id);
  auto device_extensions = RDevice_Info(device_id);
  import std.string : count;
  if ( device_extensions.count("cl_khr_gl_sharing") == 0 ) {
    writeln("Requires cl_khr_gl_sharing (for the moment)!");
    assert(0);
  }

  Initialize_Kernel(kernel_name);
}

// --- conversion helpers ---
private string To_CL_Mixin ( string name ) {
  import std.conv, std.string;
  string res;
  foreach ( i; 2 .. 9 ) {
    string si = name ~ i.to!string;
    res ~= q{
      auto To_CL%s(%s[3] a) {
        cl_%s vec;
        foreach ( i; 0 .. %s ) vec[i] = a[i];
        return vec;
      }
    }.format(si[0].to!(string).toUpper ~ si[1..$], name, si, i.to!string);
  }
  return res;
}

// --- To_CLFloat#/To_CLInt#
mixin(To_CL_Mixin("float"));
mixin(To_CL_Mixin("int"  ));
