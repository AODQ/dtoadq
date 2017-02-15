module clwrap.functions;

import clwrap.types, derelict.opencl.cl;
nothrow:

auto getPlatformIDs(uint a, platform_id* b, uint* c)
{
    debug assert(clGetPlatformIDs);
    auto ret = cast(int)clGetPlatformIDs(cast(cl_uint)a, cast(cl_platform_id*)b, cast(cl_uint*)c);
    return ret;
}

auto getPlatformInfo(platform_id a, platform_info b, size_t c, void* d, size_t* e)
{
    debug assert(clGetPlatformInfo);
    auto ret = cast(int)clGetPlatformInfo(cast(cl_platform_id)a, cast(cl_platform_info)b, cast(size_t)c, cast(void*)d, cast(size_t*)e);
    return ret;
}

auto getDeviceIDs(platform_id a, device_type b, uint c, device_id* d, uint* e)
{
    debug assert(clGetDeviceIDs);
    auto ret = cast(int)clGetDeviceIDs(cast(cl_platform_id)a, cast(cl_device_type)b, cast(cl_uint)c, cast(cl_device_id*)d, cast(cl_uint*)e);
    return ret;
}

auto getDeviceInfo(device_id a, device_info b, size_t c, void* d, size_t* e)
{
    debug assert(clGetDeviceInfo);
    auto ret = cast(int)clGetDeviceInfo(cast(cl_device_id)a, cast(cl_device_info)b, cast(size_t)c, cast(void*)d, cast(size_t*)e);
    return ret;
}

extern(System) alias createContext_FuncAlias = void function(const(char*), const(void*), size_t, void*);
auto createContext(const(context_properties*) a, uint b, const(device_id*) c, createContext_FuncAlias d, void* e, int* f)
{
    debug assert(clCreateContext);
    auto ret = cast(context)clCreateContext(cast(const(cl_context_properties*))a, cast(cl_uint)b, cast(const(cl_device_id*))c, cast(createContext_FuncAlias)d, cast(void*)e, cast(cl_int*)f);
    return ret;
}

extern(System) alias createContextFromType_FuncAlias = void function(const(char*), const(void*), size_t, void*);
auto createContextFromType(const(context_properties*) a, device_type b, createContextFromType_FuncAlias c, void* d, int* e)
{
    debug assert(clCreateContextFromType);
    auto ret = cast(context)clCreateContextFromType(cast(const(cl_context_properties*))a, cast(cl_device_type)b, cast(createContextFromType_FuncAlias)c, cast(void*)d, cast(cl_int*)e);
    return ret;
}

auto retainContext(context a)
{
    debug assert(clRetainContext);
    auto ret = cast(int)clRetainContext(cast(cl_context)a);
    return ret;
}

auto releaseContext(context a)
{
    debug assert(clReleaseContext);
    auto ret = cast(int)clReleaseContext(cast(cl_context)a);
    return ret;
}

auto getContextInfo(context a, context_info b, size_t c, void* d, size_t* e)
{
    debug assert(clGetContextInfo);
    auto ret = cast(int)clGetContextInfo(cast(cl_context)a, cast(cl_context_info)b, cast(size_t)c, cast(void*)d, cast(size_t*)e);
    return ret;
}

auto createCommandQueue(context a, device_id b, command_queue_properties c, int* d)
{
    debug assert(clCreateCommandQueue);
    auto ret = cast(command_queue)clCreateCommandQueue(cast(cl_context)a, cast(cl_device_id)b, cast(cl_command_queue_properties)c, cast(cl_int*)d);
    return ret;
}

auto retainCommandQueue(command_queue a)
{
    debug assert(clRetainCommandQueue);
    auto ret = cast(int)clRetainCommandQueue(cast(cl_command_queue)a);
    return ret;
}

auto releaseCommandQueue(command_queue a)
{
    debug assert(clReleaseCommandQueue);
    auto ret = cast(int)clReleaseCommandQueue(cast(cl_command_queue)a);
    return ret;
}

auto getCommandQueueInfo(command_queue a, command_queue_info b, size_t c, void* d, size_t* e)
{
    debug assert(clGetCommandQueueInfo);
    auto ret = cast(int)clGetCommandQueueInfo(cast(cl_command_queue)a, cast(cl_command_queue_info)b, cast(size_t)c, cast(void*)d, cast(size_t*)e);
    return ret;
}

auto createBuffer(context a, mem_flags b, size_t c, void* d, int* e)
{
    debug assert(clCreateBuffer);
    auto ret = cast(mem)clCreateBuffer(cast(cl_context)a, cast(cl_mem_flags)b, cast(size_t)c, cast(void*)d, cast(cl_int*)e);
    return ret;
}

auto retainMemObject(mem a)
{
    debug assert(clRetainMemObject);
    auto ret = cast(int)clRetainMemObject(cast(cl_mem)a);
    return ret;
}

auto releaseMemObject(mem a)
{
    debug assert(clReleaseMemObject);
    auto ret = cast(int)clReleaseMemObject(cast(cl_mem)a);
    return ret;
}

auto getSupportedImageFormats(context a, mem_flags b, mem_object_type c, uint d, image_format* e, uint* f)
{
    debug assert(clGetSupportedImageFormats);
    auto ret = cast(int)clGetSupportedImageFormats(cast(cl_context)a, cast(cl_mem_flags)b, cast(cl_mem_object_type)c, cast(cl_uint)d, cast(cl_image_format*)e, cast(cl_uint*)f);
    return ret;
}

auto getMemObjectInfo(mem a, mem_info b, size_t c, void* d, size_t* e)
{
    debug assert(clGetMemObjectInfo);
    auto ret = cast(int)clGetMemObjectInfo(cast(cl_mem)a, cast(cl_mem_info)b, cast(size_t)c, cast(void*)d, cast(size_t*)e);
    return ret;
}

auto getImageInfo(mem a, image_info b, size_t c, void* d, size_t* e)
{
    debug assert(clGetImageInfo);
    auto ret = cast(int)clGetImageInfo(cast(cl_mem)a, cast(cl_image_info)b, cast(size_t)c, cast(void*)d, cast(size_t*)e);
    return ret;
}

auto createSampler(context a, bool b, addressing_mode c, filter_mode d, int* e)
{
    debug assert(clCreateSampler);
    auto ret = cast(sampler)clCreateSampler(cast(cl_context)a, cast(cl_bool)b, cast(cl_addressing_mode)c, cast(cl_filter_mode)d, cast(cl_int*)e);
    return ret;
}

auto retainSampler(sampler a)
{
    debug assert(clRetainSampler);
    auto ret = cast(int)clRetainSampler(cast(cl_sampler)a);
    return ret;
}

auto releaseSampler(sampler a)
{
    debug assert(clReleaseSampler);
    auto ret = cast(int)clReleaseSampler(cast(cl_sampler)a);
    return ret;
}

auto getSamplerInfo(sampler a, sampler_info b, size_t c, void* d, size_t* e)
{
    debug assert(clGetSamplerInfo);
    auto ret = cast(int)clGetSamplerInfo(cast(cl_sampler)a, cast(cl_sampler_info)b, cast(size_t)c, cast(void*)d, cast(size_t*)e);
    return ret;
}

auto createProgramWithSource(context a, uint b, const(char*)* c, const(size_t*) d, int* e)
{
    debug assert(clCreateProgramWithSource);
    auto ret = cast(program)clCreateProgramWithSource(cast(cl_context)a, cast(cl_uint)b, cast(const(char*)*)c, cast(const(size_t*))d, cast(cl_int*)e);
    return ret;
}

auto createProgramWithBinary(context a, uint b, const(device_id*) c, const(size_t*) d, const(ubyte*)* e, int* f, int* g)
{
    debug assert(clCreateProgramWithBinary);
    auto ret = cast(program)clCreateProgramWithBinary(cast(cl_context)a, cast(cl_uint)b, cast(const(cl_device_id*))c, cast(const(size_t*))d, cast(const(ubyte*)*)e, cast(cl_int*)f, cast(cl_int*)g);
    return ret;
}

auto createProgramWithBuiltInKernels(context a, uint b, const(device_id*) c, const(char*) d, int* e)
{
    debug assert(clCreateProgramWithBuiltInKernels);
    auto ret = cast(program)clCreateProgramWithBuiltInKernels(cast(cl_context)a, cast(cl_uint)b, cast(const(cl_device_id*))c, cast(const(char*))d, cast(cl_int*)e);
    return ret;
}

auto retainProgram(program a)
{
    debug assert(clRetainProgram);
    auto ret = cast(int)clRetainProgram(cast(cl_program)a);
    return ret;
}

auto releaseProgram(program a)
{
    debug assert(clReleaseProgram);
    auto ret = cast(int)clReleaseProgram(cast(cl_program)a);
    return ret;
}

extern(System) alias buildProgram_FuncAlias = void function(cl_program, void*);
auto buildProgram(program a, uint b, const(device_id*) c, const(char*) d, buildProgram_FuncAlias e, void* f)
{
    debug assert(clBuildProgram);
    auto ret = cast(int)clBuildProgram(cast(cl_program)a, cast(cl_uint)b, cast(const(cl_device_id*))c, cast(const(char*))d, cast(buildProgram_FuncAlias)e, cast(void*)f);
    return ret;
}

auto getProgramInfo(program a, program_info b, size_t c, void* d, size_t* e)
{
    debug assert(clGetProgramInfo);
    auto ret = cast(int)clGetProgramInfo(cast(cl_program)a, cast(cl_program_info)b, cast(size_t)c, cast(void*)d, cast(size_t*)e);
    return ret;
}

auto getProgramBuildInfo(program a, device_id b, program_build_info c, size_t d, void* e, size_t* f)
{
    debug assert(clGetProgramBuildInfo);
    auto ret = cast(int)clGetProgramBuildInfo(cast(cl_program)a, cast(cl_device_id)b, cast(cl_program_build_info)c, cast(size_t)d, cast(void*)e, cast(size_t*)f);
    return ret;
}

auto createKernel(program a, const(char*) b, int* c)
{
    debug assert(clCreateKernel);
    auto ret = cast(kernel)clCreateKernel(cast(cl_program)a, cast(const(char*))b, cast(cl_int*)c);
    return ret;
}

auto createKernelsInProgram(program a, uint b, kernel* c, uint* d)
{
    debug assert(clCreateKernelsInProgram);
    auto ret = cast(int)clCreateKernelsInProgram(cast(cl_program)a, cast(cl_uint)b, cast(cl_kernel*)c, cast(cl_uint*)d);
    return ret;
}

auto retainKernel(kernel a)
{
    debug assert(clRetainKernel);
    auto ret = cast(int)clRetainKernel(cast(cl_kernel)a);
    return ret;
}

auto releaseKernel(kernel a)
{
    debug assert(clReleaseKernel);
    auto ret = cast(int)clReleaseKernel(cast(cl_kernel)a);
    return ret;
}

auto setKernelArg(kernel a, uint b, size_t c, const(void*) d)
{
    debug assert(clSetKernelArg);
    auto ret = cast(int)clSetKernelArg(cast(cl_kernel)a, cast(cl_uint)b, cast(size_t)c, cast(const(void*))d);
    return ret;
}

auto getKernelInfo(kernel a, kernel_info b, size_t c, void* d, size_t* e)
{
    debug assert(clGetKernelInfo);
    auto ret = cast(int)clGetKernelInfo(cast(cl_kernel)a, cast(cl_kernel_info)b, cast(size_t)c, cast(void*)d, cast(size_t*)e);
    return ret;
}

auto getKernelArgInfo(kernel a, uint b, kernel_arg_info c, size_t d, void* e, size_t* f)
{
    debug assert(clGetKernelArgInfo);
    auto ret = cast(int)clGetKernelArgInfo(cast(cl_kernel)a, cast(cl_uint)b, cast(cl_kernel_arg_info)c, cast(size_t)d, cast(void*)e, cast(size_t*)f);
    return ret;
}

auto getKernelWorkGroupInfo(kernel a, device_id b, kernel_work_group_info c, size_t d, void* e, size_t* f)
{
    debug assert(clGetKernelWorkGroupInfo);
    auto ret = cast(int)clGetKernelWorkGroupInfo(cast(cl_kernel)a, cast(cl_device_id)b, cast(cl_kernel_work_group_info)c, cast(size_t)d, cast(void*)e, cast(size_t*)f);
    return ret;
}

auto waitForEvents(uint a, const(event*) b)
{
    debug assert(clWaitForEvents);
    auto ret = cast(int)clWaitForEvents(cast(cl_uint)a, cast(const(cl_event*))b);
    return ret;
}

auto getEventInfo(event a, event_info b, size_t c, void* d, size_t* e)
{
    debug assert(clGetEventInfo);
    auto ret = cast(int)clGetEventInfo(cast(cl_event)a, cast(cl_event_info)b, cast(size_t)c, cast(void*)d, cast(size_t*)e);
    return ret;
}

auto retainEvent(event a)
{
    debug assert(clRetainEvent);
    auto ret = cast(int)clRetainEvent(cast(cl_event)a);
    return ret;
}

auto releaseEvent(event a)
{
    debug assert(clReleaseEvent);
    auto ret = cast(int)clReleaseEvent(cast(cl_event)a);
    return ret;
}

auto getEventProfilingInfo(event a, profiling_info b, size_t c, void* d, size_t* e)
{
    debug assert(clGetEventProfilingInfo);
    auto ret = cast(int)clGetEventProfilingInfo(cast(cl_event)a, cast(cl_profiling_info)b, cast(size_t)c, cast(void*)d, cast(size_t*)e);
    return ret;
}

auto flush(command_queue a)
{
    debug assert(clFlush);
    auto ret = cast(int)clFlush(cast(cl_command_queue)a);
    return ret;
}

auto finish(command_queue a)
{
    debug assert(clFinish);
    auto ret = cast(int)clFinish(cast(cl_command_queue)a);
    return ret;
}

auto enqueueReadBuffer(command_queue a, mem b, bool c, size_t d, size_t e, void* f, uint g, const(event*) h, event* i)
{
    debug assert(clEnqueueReadBuffer);
    auto ret = cast(int)clEnqueueReadBuffer(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_bool)c, cast(size_t)d, cast(size_t)e, cast(void*)f, cast(cl_uint)g, cast(const(cl_event*))h, cast(cl_event*)i);
    return ret;
}

auto enqueueWriteBuffer(command_queue a, mem b, bool c, size_t d, size_t e, const(void*) f, uint g, const(event*) h, event* i)
{
    debug assert(clEnqueueWriteBuffer);
    auto ret = cast(int)clEnqueueWriteBuffer(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_bool)c, cast(size_t)d, cast(size_t)e, cast(const(void*))f, cast(cl_uint)g, cast(const(cl_event*))h, cast(cl_event*)i);
    return ret;
}

auto enqueueCopyBuffer(command_queue a, mem b, mem c, size_t d, size_t e, size_t f, uint g, const(event*) h, event* i)
{
    debug assert(clEnqueueCopyBuffer);
    auto ret = cast(int)clEnqueueCopyBuffer(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_mem)c, cast(size_t)d, cast(size_t)e, cast(size_t)f, cast(cl_uint)g, cast(const(cl_event*))h, cast(cl_event*)i);
    return ret;
}

auto enqueueReadImage(command_queue a, mem b, bool c, const(size_t*) d, const(size_t*) e, size_t f, size_t g, void* h, uint i, const(event*) j, event* k)
{
    debug assert(clEnqueueReadImage);
    auto ret = cast(int)clEnqueueReadImage(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_bool)c, cast(const(size_t*))d, cast(const(size_t*))e, cast(size_t)f, cast(size_t)g, cast(void*)h, cast(cl_uint)i, cast(const(cl_event*))j, cast(cl_event*)k);
    return ret;
}

auto enqueueWriteImage(command_queue a, mem b, bool c, const(size_t*) d, const(size_t*) e, size_t f, size_t g, const(void*) h, uint i, const(event*) j, event* k)
{
    debug assert(clEnqueueWriteImage);
    auto ret = cast(int)clEnqueueWriteImage(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_bool)c, cast(const(size_t*))d, cast(const(size_t*))e, cast(size_t)f, cast(size_t)g, cast(const(void*))h, cast(cl_uint)i, cast(const(cl_event*))j, cast(cl_event*)k);
    return ret;
}

auto enqueueCopyImage(command_queue a, mem b, mem c, const(size_t*) d, const(size_t*) e, const(size_t*) f, uint g, const(event*) h, event* i)
{
    debug assert(clEnqueueCopyImage);
    auto ret = cast(int)clEnqueueCopyImage(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_mem)c, cast(const(size_t*))d, cast(const(size_t*))e, cast(const(size_t*))f, cast(cl_uint)g, cast(const(cl_event*))h, cast(cl_event*)i);
    return ret;
}

auto enqueueCopyImageToBuffer(command_queue a, mem b, mem c, const(size_t*) d, const(size_t*) e, size_t f, uint g, const(event*) h, event* i)
{
    debug assert(clEnqueueCopyImageToBuffer);
    auto ret = cast(int)clEnqueueCopyImageToBuffer(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_mem)c, cast(const(size_t*))d, cast(const(size_t*))e, cast(size_t)f, cast(cl_uint)g, cast(const(cl_event*))h, cast(cl_event*)i);
    return ret;
}

auto enqueueCopyBufferToImage(command_queue a, mem b, mem c, size_t d, const(size_t*) e, const(size_t*) f, uint g, const(event*) h, event* i)
{
    debug assert(clEnqueueCopyBufferToImage);
    auto ret = cast(int)clEnqueueCopyBufferToImage(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_mem)c, cast(size_t)d, cast(const(size_t*))e, cast(const(size_t*))f, cast(cl_uint)g, cast(const(cl_event*))h, cast(cl_event*)i);
    return ret;
}

auto enqueueMapBuffer(command_queue a, mem b, bool c, map_flags d, size_t e, size_t f, uint g, const(event*) h, event* i, int* j)
{
    debug assert(clEnqueueMapBuffer);
    auto ret = cast(void*)clEnqueueMapBuffer(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_bool)c, cast(cl_map_flags)d, cast(size_t)e, cast(size_t)f, cast(cl_uint)g, cast(const(cl_event*))h, cast(cl_event*)i, cast(cl_int*)j);
    return ret;
}

auto enqueueMapImage(command_queue a, mem b, bool c, map_flags d, const(size_t*) e, const(size_t*) f, size_t* g, size_t* h, uint i, const(event*) j, event* k, int* l)
{
    debug assert(clEnqueueMapImage);
    auto ret = cast(void*)clEnqueueMapImage(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_bool)c, cast(cl_map_flags)d, cast(const(size_t*))e, cast(const(size_t*))f, cast(size_t*)g, cast(size_t*)h, cast(cl_uint)i, cast(const(cl_event*))j, cast(cl_event*)k, cast(cl_int*)l);
    return ret;
}

auto enqueueUnmapMemObject(command_queue a, mem b, void* c, uint d, const(event*) e, event* f)
{
    debug assert(clEnqueueUnmapMemObject);
    auto ret = cast(int)clEnqueueUnmapMemObject(cast(cl_command_queue)a, cast(cl_mem)b, cast(void*)c, cast(cl_uint)d, cast(const(cl_event*))e, cast(cl_event*)f);
    return ret;
}

auto enqueueNDRangeKernel(command_queue a, kernel b, uint c, const(size_t*) d, const(size_t*) e, const(size_t*) f, uint g, const(event*) h, event* i)
{
    debug assert(clEnqueueNDRangeKernel);
    auto ret = cast(int)clEnqueueNDRangeKernel(cast(cl_command_queue)a, cast(cl_kernel)b, cast(cl_uint)c, cast(const(size_t*))d, cast(const(size_t*))e, cast(const(size_t*))f, cast(cl_uint)g, cast(const(cl_event*))h, cast(cl_event*)i);
    return ret;
}

auto enqueueTask(command_queue a, kernel b, uint c, const(event*) d, event* e)
{
    debug assert(clEnqueueTask);
    auto ret = cast(int)clEnqueueTask(cast(cl_command_queue)a, cast(cl_kernel)b, cast(cl_uint)c, cast(const(cl_event*))d, cast(cl_event*)e);
    return ret;
}

extern(System) alias enqueueNativeKernel_FuncAlias = void function(void*);
auto enqueueNativeKernel(command_queue a, enqueueNativeKernel_FuncAlias b, void* c, size_t d, uint e, const(mem*) f, const(void*)* g, uint h, const(event*) i, event* j)
{
    debug assert(clEnqueueNativeKernel);
    auto ret = cast(int)clEnqueueNativeKernel(cast(cl_command_queue)a, cast(enqueueNativeKernel_FuncAlias)b, cast(void*)c, cast(size_t)d, cast(cl_uint)e, cast(const(cl_mem*))f, cast(const(void*)*)g, cast(cl_uint)h, cast(const(cl_event*))i, cast(cl_event*)j);
    return ret;
}

auto setCommandQueueProperty(command_queue a, command_queue_properties b, bool c, command_queue_properties* d)
{
    debug assert(clSetCommandQueueProperty);
    auto ret = cast(int)clSetCommandQueueProperty(cast(cl_command_queue)a, cast(cl_command_queue_properties)b, cast(cl_bool)c, cast(cl_command_queue_properties*)d);
    return ret;
}

auto createSubBuffer(mem a, mem_flags b, buffer_create_type c, const(void*) d, int* e)
{
    debug assert(clCreateSubBuffer);
    auto ret = cast(mem)clCreateSubBuffer(cast(cl_mem)a, cast(cl_mem_flags)b, cast(cl_buffer_create_type)c, cast(const(void*))d, cast(cl_int*)e);
    return ret;
}

extern(System) alias setMemObjectDestructorCallback_FuncAlias = void function(cl_mem, void*);
auto setMemObjectDestructorCallback(mem a, setMemObjectDestructorCallback_FuncAlias b, void* c)
{
    debug assert(clSetMemObjectDestructorCallback);
    auto ret = cast(int)clSetMemObjectDestructorCallback(cast(cl_mem)a, cast(setMemObjectDestructorCallback_FuncAlias)b, cast(void*)c);
    return ret;
}

auto createUserEvent(context a, int* b)
{
    debug assert(clCreateUserEvent);
    auto ret = cast(event)clCreateUserEvent(cast(cl_context)a, cast(cl_int*)b);
    return ret;
}

auto setUserEventStatus(event a, int b)
{
    debug assert(clSetUserEventStatus);
    auto ret = cast(int)clSetUserEventStatus(cast(cl_event)a, cast(cl_int)b);
    return ret;
}

extern(System) alias setEventCallback_FuncAlias = void function(cl_event, cl_int, void*);
auto setEventCallback(event a, int b, setEventCallback_FuncAlias c, void* d)
{
    debug assert(clSetEventCallback);
    auto ret = cast(int)clSetEventCallback(cast(cl_event)a, cast(cl_int)b, cast(setEventCallback_FuncAlias)c, cast(void*)d);
    return ret;
}

auto enqueueReadBufferRect(command_queue a, mem b, bool c, const(size_t*) d, const(size_t*) e, const(size_t*) f, size_t g, size_t h, size_t i, size_t j, void* k, uint l, const(event*) m, event* n)
{
    debug assert(clEnqueueReadBufferRect);
    auto ret = cast(int)clEnqueueReadBufferRect(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_bool)c, cast(const(size_t*))d, cast(const(size_t*))e, cast(const(size_t*))f, cast(size_t)g, cast(size_t)h, cast(size_t)i, cast(size_t)j, cast(void*)k, cast(cl_uint)l, cast(const(cl_event*))m, cast(cl_event*)n);
    return ret;
}

auto enqueueWriteBufferRect(command_queue a, mem b, bool c, const(size_t*) d, const(size_t*) e, const(size_t*) f, size_t g, size_t h, size_t i, size_t j, const(void*) k, uint l, const(event*) m, event* n)
{
    debug assert(clEnqueueWriteBufferRect);
    auto ret = cast(int)clEnqueueWriteBufferRect(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_bool)c, cast(const(size_t*))d, cast(const(size_t*))e, cast(const(size_t*))f, cast(size_t)g, cast(size_t)h, cast(size_t)i, cast(size_t)j, cast(const(void*))k, cast(cl_uint)l, cast(const(cl_event*))m, cast(cl_event*)n);
    return ret;
}

auto enqueueCopyBufferRect(command_queue a, mem b, mem c, const(size_t*) d, const(size_t*) e, const(size_t*) f, size_t g, size_t h, size_t i, size_t j, uint k, const(event*) l, event* m)
{
    debug assert(clEnqueueCopyBufferRect);
    auto ret = cast(int)clEnqueueCopyBufferRect(cast(cl_command_queue)a, cast(cl_mem)b, cast(cl_mem)c, cast(const(size_t*))d, cast(const(size_t*))e, cast(const(size_t*))f, cast(size_t)g, cast(size_t)h, cast(size_t)i, cast(size_t)j, cast(cl_uint)k, cast(const(cl_event*))l, cast(cl_event*)m);
    return ret;
}

auto createImage2D(context a, mem_flags b, const(image_format*) c, size_t d, size_t e, size_t f, void* g, int* h)
{
    debug assert(clCreateImage2D);
    auto ret = cast(mem)clCreateImage2D(cast(cl_context)a, cast(cl_mem_flags)b, cast(const(cl_image_format*))c, cast(size_t)d, cast(size_t)e, cast(size_t)f, cast(void*)g, cast(cl_int*)h);
    return ret;
}

auto createImage3D(context a, mem_flags b, const(image_format*) c, size_t d, size_t e, size_t f, size_t g, size_t h, void* i, int* j)
{
    debug assert(clCreateImage3D);
    auto ret = cast(mem)clCreateImage3D(cast(cl_context)a, cast(cl_mem_flags)b, cast(const(cl_image_format*))c, cast(size_t)d, cast(size_t)e, cast(size_t)f, cast(size_t)g, cast(size_t)h, cast(void*)i, cast(cl_int*)j);
    return ret;
}

auto enqueueMarker(command_queue a, event* b)
{
    debug assert(clEnqueueMarker);
    auto ret = cast(int)clEnqueueMarker(cast(cl_command_queue)a, cast(cl_event*)b);
    return ret;
}

auto enqueueWaitForEvents(command_queue a, uint b, const(event*) c)
{
    debug assert(clEnqueueWaitForEvents);
    auto ret = cast(int)clEnqueueWaitForEvents(cast(cl_command_queue)a, cast(cl_uint)b, cast(const(cl_event*))c);
    return ret;
}

auto enqueueBarrier(command_queue a)
{
    debug assert(clEnqueueBarrier);
    auto ret = cast(int)clEnqueueBarrier(cast(cl_command_queue)a);
    return ret;
}

auto unloadCompiler()
{
    debug assert(clUnloadCompiler);
    auto ret = cast(int)clUnloadCompiler();
    return ret;
}

auto getExtensionFunctionAddress(const(char*) a)
{
    debug assert(clGetExtensionFunctionAddress);
    auto ret = cast(void*)clGetExtensionFunctionAddress(cast(const(char*))a);
    return ret;
}

auto createSubDevices(device_id a, const(device_partition_property*) b, uint c, device_id* d, uint* e)
{
    debug assert(clCreateSubDevices);
    auto ret = cast(int)clCreateSubDevices(cast(cl_device_id)a, cast(const(cl_device_partition_property*))b, cast(cl_uint)c, cast(cl_device_id*)d, cast(cl_uint*)e);
    return ret;
}

auto retainDevice(device_id a)
{
    debug assert(clRetainDevice);
    auto ret = cast(int)clRetainDevice(cast(cl_device_id)a);
    return ret;
}

auto releaseDevice(device_id a)
{
    debug assert(clReleaseDevice);
    auto ret = cast(int)clReleaseDevice(cast(cl_device_id)a);
    return ret;
}

auto createImage(context a, mem_flags b, const(image_format*) c, const(image_desc*) d, void* e, int* f)
{
    debug assert(clCreateImage);
    auto ret = cast(mem)clCreateImage(cast(cl_context)a, cast(cl_mem_flags)b, cast(const(cl_image_format*))c, cast(const(cl_image_desc*))d, cast(void*)e, cast(cl_int*)f);
    return ret;
}

extern(System) alias compileProgram_FuncAlias = void function(cl_program, void*);
auto compileProgram(program a, uint b, const(device_id*) c, const(char*) d, uint e, const(program*) f, const(char*)* g, compileProgram_FuncAlias h, void* i)
{
    debug assert(clCompileProgram);
    auto ret = cast(int)clCompileProgram(cast(cl_program)a, cast(cl_uint)b, cast(const(cl_device_id*))c, cast(const(char*))d, cast(cl_uint)e, cast(const(cl_program*))f, cast(const(char*)*)g, cast(compileProgram_FuncAlias)h, cast(void*)i);
    return ret;
}

extern(System) alias linkProgram_FuncAlias = void function(cl_program, void*);
auto linkProgram(context a, uint b, const(device_id*) c, const(char*) d, uint e, const(program*) f, linkProgram_FuncAlias g, void* h, int* i)
{
    debug assert(clLinkProgram);
    auto ret = cast(program)clLinkProgram(cast(cl_context)a, cast(cl_uint)b, cast(const(cl_device_id*))c, cast(const(char*))d, cast(cl_uint)e, cast(const(cl_program*))f, cast(linkProgram_FuncAlias)g, cast(void*)h, cast(cl_int*)i);
    return ret;
}

auto unloadPlatformCompiler(platform_id a)
{
    debug assert(clUnloadPlatformCompiler);
    auto ret = cast(int)clUnloadPlatformCompiler(cast(cl_platform_id)a);
    return ret;
}

auto enqueueFillBuffer(command_queue a, mem b, const(void*) c, size_t d, size_t e, size_t f, uint g, const(event*) h, event* i)
{
    debug assert(clEnqueueFillBuffer);
    auto ret = cast(int)clEnqueueFillBuffer(cast(cl_command_queue)a, cast(cl_mem)b, cast(const(void*))c, cast(size_t)d, cast(size_t)e, cast(size_t)f, cast(cl_uint)g, cast(const(cl_event*))h, cast(cl_event*)i);
    return ret;
}

auto enqueueFillImage(command_queue a, mem b, const(void*) c, const(size_t*) d, const(size_t*) e, uint f, const(event*) g, event* h)
{
    debug assert(clEnqueueFillImage);
    auto ret = cast(int)clEnqueueFillImage(cast(cl_command_queue)a, cast(cl_mem)b, cast(const(void*))c, cast(const(size_t*))d, cast(const(size_t*))e, cast(cl_uint)f, cast(const(cl_event*))g, cast(cl_event*)h);
    return ret;
}

auto enqueueMigrateMemObjects(command_queue a, uint b, const(mem*) c, mem_migration_flags d, uint e, const(event*) f, event* g)
{
    debug assert(clEnqueueMigrateMemObjects);
    auto ret = cast(int)clEnqueueMigrateMemObjects(cast(cl_command_queue)a, cast(cl_uint)b, cast(const(cl_mem*))c, cast(cl_mem_migration_flags)d, cast(cl_uint)e, cast(const(cl_event*))f, cast(cl_event*)g);
    return ret;
}

auto enqueueMarkerWithWaitList(command_queue a, uint b, const(event*) c, event* d)
{
    debug assert(clEnqueueMarkerWithWaitList);
    auto ret = cast(int)clEnqueueMarkerWithWaitList(cast(cl_command_queue)a, cast(cl_uint)b, cast(const(cl_event*))c, cast(cl_event*)d);
    return ret;
}

auto enqueueBarrierWithWaitList(command_queue a, uint b, const(event*) c, event* d)
{
    debug assert(clEnqueueBarrierWithWaitList);
    auto ret = cast(int)clEnqueueBarrierWithWaitList(cast(cl_command_queue)a, cast(cl_uint)b, cast(const(cl_event*))c, cast(cl_event*)d);
    return ret;
}

auto getExtensionFunctionAddressForPlatform(platform_id a, const(char*) b)
{
    debug assert(clGetExtensionFunctionAddressForPlatform);
    auto ret = cast(void*)clGetExtensionFunctionAddressForPlatform(cast(cl_platform_id)a, cast(const(char*))b);
    return ret;
}

