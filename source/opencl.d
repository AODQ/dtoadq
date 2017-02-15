module opencl;

public import derelict.opencl.cl;
public static import wrap = clwrap.clwrap;
import clwrap.clwrap, std.stdio;
import std.conv : to;

auto RPlatforms() {
  platform_id[] platforms;
  cl_uint platform_amt = 10;
  platforms.length = platform_amt;
  getPlatformIDs(platform_amt, platforms.ptr, &platform_amt);
  platforms.length = platform_amt;
  return platforms;
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
  auto platforms = RPlatforms
  writeln("PLATFORMS: ", platforms);
  import std.algorithm;
  platforms.each!(n => n.RPlatform_Info.writeln);
}


void Initialize ( ) {
  DerelictCL.load();
}
