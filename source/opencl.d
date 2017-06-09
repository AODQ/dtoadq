module opencl;
/**
  A general-purpose opencl module, though most parts are specific to my
    project. For example, to run the kernel you may only use two-dimensional
    global group, most functions only exist on a per-demand-basis, etc
*/
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

static import std.traits;

struct CLStoreMem { // store buffer into this after run
  CLMem _mem;
  alias _mem this;
  void* data;
  size_t size;
  this (T) ( ref T data_ ) {
    static if ( std.traits.isArray!T ) {
      data = cast(void*)data_.ptr;
      size = (data_[0]).sizeof*data_.length;
    } else {
      data = cast(void*)&data_;
      size = T.sizeof;
    }
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

public abstract class CL {
static:
  CLContext           context;
  CLContextProperties properties;
  CLCommandQueue      command_queue;
  CLProgram           program;
  CLKernel            kernel;
}
int err;

private void Initialize_Kernel ( ) {
  {
    import derelict.opengl3.glx;
    alias CLCP = cl_context_properties;
    CL.properties = [CL_CONTEXT_PLATFORM, cast(CLCP)device.platform_id,
                     CL_GL_CONTEXT_KHR,   cast(CLCP)glXGetCurrentContext(),
                     CL_GLX_DISPLAY_KHR,  cast(CLCP)glXGetCurrentDisplay(),
                     0];
  }
  CL.context = clCreateContext(CL.properties.ptr, 1, &device.device_id,
                            null, null, &err);
  CLAssert(err, "clCreateContext");
  CL.command_queue = clCreateCommandQueue(CL.context, device.device_id, 0,&err);
  CLAssert(err, "clCreateCommandQueue");
}

static bool working_program = true;
void Compile ( inout string source, string kernel_name ) {
  import std.string : toStringz;
  auto sourcestr = source.toStringz;
  CL.program = clCreateProgramWithSource(CL.context, 1, &sourcestr, null, &err);
  CLAssert(err, "clCreateProgramWithSource");
  static import std.file;
  std.file.write("BUILD", source);
  working_program =
    clBuildProgram(CL.program, 0, null, null, null, null) == CL_SUCCESS;
  if ( !working_program ) {
    size_t len;
    clGetProgramBuildInfo(CL.program, device.device_id, CL_PROGRAM_BUILD_LOG,
                          0, null, &len);
    char[] log; log.length = len;
    clGetProgramBuildInfo(CL.program, device.device_id, CL_PROGRAM_BUILD_LOG,
                          len, log.ptr, null);
    writeln("LOG: ", log);
    writeln("-- could not compile --");
    return;
  }
  if ( CL.kernel !is null )
    clReleaseKernel(CL.kernel);
  CL.kernel = clCreateKernel(CL.program, kernel_name.toStringz, &err);
  CLAssert(err, "clCreateKernel");
}

auto Create_CLGL_Texture ( uint gl_texture ) {
  import derelict.opengl3.gl3;
  auto cl_handle = clCreateFromGLTexture(CL.context, CL_MEM_READ_WRITE,
                                   GL_TEXTURE_2D, 0, gl_texture, &err);
  CLAssert(err, "clCreateFromGLTexture");
  return cl_handle;
}

auto Create_CL_Image ( cl_mem_flags flags, cl_image_format format,
                       cl_image_desc desc, void* data ) {
  auto mem = clCreateImage(CL.context, flags, &format, &desc, data, &err);
  CLAssert(err, "clCreateImage");
  return mem;
}

private auto RHandles ( CLPredefinedMem[] cl_predefined_mem ) {
  import functional;
  return cl_predefined_mem.map!(n => n._mem).array;
}
void Lock_CLGL_Image ( CLPredefinedMem cl_predefined_mem ) {
  if ( !working_program ) return;
  import derelict.opengl3.gl3;
  auto handle = cl_predefined_mem._mem;
  CLAssert(clEnqueueAcquireGLObjects(CL.command_queue, 1, &handle,
            0, null, null), "clEnqueueAcquireGLObjects");
}

void Unlock_CLGL_Image ( CLPredefinedMem cl_predefined_mem ) {
  if ( !working_program ) return;
  import derelict.opengl3.gl3;
  auto handle = cl_predefined_mem._mem;
  CLAssert(clEnqueueReleaseGLObjects(CL.command_queue, 1, &handle,
            0, null, null), "clEnqueueReleaseGLObjects");
}

void Sync_GL_Event ( cl_event event ) {
  import derelict.opengl3.gl3;
  // event = glCreateSyncFromCLeventARB(cast(_cl_context)CL.context, event, 0);
  // glWaitSynx(event, 0, GL_TIMEOUT_IGNORE);
  glFinish(); // TODO remove this if possible
  clFinish(CL.command_queue); // and this
}

void Flush ( ) {
  CLAssert(clFlush(CL.command_queue), "clFlush");
  CLAssert(clFinish(CL.command_queue), "clFinish");
}

private auto Load_Argument(T)(ref T data, int it) {
  import std.traits;
  static auto flgs = CL_MEM_READ_WRITE;

  CLMem cl_handle;
  static if ( isArray!(T) )
    cl_handle = clCreateBuffer(CL.context, flgs, data.length*T.sizeof, &data[0],
                               &err);
  else static if ( is(T == CLPredefinedMem) ) {
    cl_handle = data._mem;
    err = 0;
  } else static if ( is(T == CLStoreMem) ) {
    data._mem = cl_handle =
      clCreateBuffer(CL.context, flgs, data.size, data.data, &err);
  } else {
    cl_handle = clCreateBuffer(CL.context, flgs, T.sizeof, &data, &err);
  }
  CLAssert(err, "clCreateBuffer");

  static if ( isArray!(T) ) {
    CLAssert(clEnqueueWriteBuffer(CL.command_queue, cl_handle, CL_TRUE, 0,
        data.length*T.sizeof, &data[0], 0, null, null), "clEnqueueWriteBuffer");
  } else static if ( !is(T == CLPredefinedMem) ) {
    void* vp_data = &data;
    auto length = T.sizeof;
    static if ( is(T == CLStoreMem) ) {
      vp_data = data.data;
      length  = data.size;
    }
    CLAssert(clEnqueueWriteBuffer(CL.command_queue, cl_handle, CL_TRUE, 0,
                    length, vp_data, 0, null, null), "clEnqueueWriteBuffer");
  }
  // don't do writing for a predefinedclmem
  CLAssert(clSetKernelArg(CL.kernel, it, CLMem.sizeof, &cl_handle),
                                                "clSetKernelArg");
  return cl_handle;
}

void Run ( T... ) ( T params, size_t X, size_t Y ) {
  if ( !working_program ) return;
  CLMem[] cl_handles_dealloc;
  CLStoreMem[] cl_store_mem;

  foreach ( it, ref p; params ) {
    CLMem cl_handle = Load_Argument(p, it);
    static if ( !is(typeof(p) == CLPredefinedMem) )
      cl_handles_dealloc ~= cl_handle;
    static if (  is(typeof(p) == CLStoreMem) )
      cl_store_mem ~= p;
  }

  ulong[] global = [cast(ulong)X, cast(ulong)Y];
  cl_event kernel_event;
  CLAssert(clEnqueueNDRangeKernel(CL.command_queue, CL.kernel, 2, null,
          global.ptr, null, 0, null, &kernel_event), "clEnqueueNDRangeKernel");

  CLAssert(clWaitForEvents(1, &kernel_event), "clWaitForEvents");

  foreach ( mem; cl_store_mem ) {
    CLAssert(clEnqueueReadBuffer(CL.command_queue, mem._mem, CL_TRUE, 0,
            mem.size, mem.data, 0, null, null), "clEnqueueReadBuffer");
  }

  foreach ( mem; cl_handles_dealloc ) {
    clReleaseMemObject(mem);
  }

}

void Clean_Up() {
  clReleaseProgram(CL.program);
  clReleaseKernel(CL.kernel);
  clReleaseCommandQueue(CL.command_queue);
  clReleaseContext(CL.context);
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

void Initialize ( ) {
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

  Initialize_Kernel();
}

// --- conversion helpers ---
string To_OpenCL_Float ( float float_valf ) {
  import std.conv : to;
  return float_valf.to!string.To_OpenCL_Float;
}

string To_OpenCL_Float ( string float_val ) {
  import functional;
  // -- "3"|"3."|"3.0" -> "3.0f" (there are no f in imgui floats)
  if ( float_val.filter!(n => n == '.').empty ) float_val ~= ".";
  if ( float_val[$-1] == '.' ) float_val ~= "0"; // to include "3." & "3"
  return float_val ~ "f";
}

string To_OpenCL_Float_Array ( float[] float_vals ) {
  import functional;
  auto v = float_vals.map!(n => n.To_OpenCL_Float).array;
  if ( v.length == 2 ) return "(float2)(" ~ v[0] ~ ", " ~ v[1] ~ ")";
  return "(float3)(" ~ v[0] ~ ", " ~ v[1] ~ ", " ~ v[2] ~ ")";
}

private string To_CL_Mixin ( string name ) {
  import std.conv, std.string;
  string res;
  foreach ( i; 2 .. 5 ) {
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
