module emitter.video;
import progress : Progress;
static import core.info;
static import std.file;
import std.conv : to;
static import stl, y4md, core;

Progress progress;
y4md.Y4MWriter output;
string output_name;

float time, end_time;
int fps, spp;
int frame;
private core.Resolution resolution;

auto RRender_Resolution ( ) { return resolution; }

void Render (string scene, string out_, int spp_, float end_time_, int fps_) {
  time = 0.0f;
  frame = 0;
  end_time = end_time_;
  fps = fps_;
  spp = spp_;
  assert(spp_ >= 0 && spp_ <= 255, "SPP incorrect");
  {
    static import core.shared_info;
    core.shared_info.image_metadata.spp = cast(ubyte)spp;
  }
  resolution = core.Resolution.r960_540;
  static import core.shared_info, core.image;
  core.shared_info.Set_Image_Buffer(core.image.Resolution.r960_540, true);

  // -- set up video emitter --
  output_name = out_;
  try {
    std.file.mkdir("gifout");
  } catch ( Exception e ) {}
  // output = new y4md.Y4MWriter(out_,
  //       core.RImage(resolution).x, core.RImage(resolution).y,
  //       y4md.Rational(cast(int)fps, 1), y4md.Rational(0, 0),
  //       y4md.Interlacing.Progressive, y4md.Subsampling.C420);

  // -- set up progress bar --
  progress = new Progress(cast(int)(end_time*fps));
  progress.title = "Rendering";

  // -- set up scene/renderer --
  core.Set_Image_Buffer(resolution);
  core.info.Set_Map_Function(scene);
  core.info.Set_Kernel_Type(core.info.KernelType.VideoRender);
}

auto Pad(int f) {
  import std.string : rightJustify;
  string t = f.to!string;
  return rightJustify(t, 7, '0');
}

// Returns true when render is finished
bool Update ( ubyte[] data ) {
  time += 1.0f/fps;
  core.Set_Time(time);
  progress.next();
  // output.writeFrame(data);
  ++frame;
  string filnam    = "gifout/tempout-"~Pad(frame)~".pbm",
         filnampng = "gifout/tempout-"~Pad(frame)~".png";
  int width  = core.RImage(resolution).x,
      height = core.RImage(resolution).y;
  std.file.write(filnam, "P3\n"~width.to!string ~ " " ~
                                height.to!string~"\n255\n");
  foreach ( j; 0 .. height ) {
    foreach ( i; 0 .. width ) {
      auto pos = (height-j-1)*width*4 + i*4;
      std.file.append(filnam, data[pos  ].to!string ~ " " ~
                              data[pos+1].to!string ~ " " ~
                              data[pos+2].to!string ~ " ");
    }
    std.file.append(filnam, "\n");
  }
  import std.process : spawnShell, wait;
  {//convert to png
    {
      auto pid = spawnShell(`convert "` ~ filnam ~ `" "` ~ filnampng ~ `"`);
      wait(pid);
    } {
      auto pid = spawnShell(`rm "` ~ filnam ~ `"`);
      wait(pid);
    }
  }
  if ( time >= end_time ) {
    {
      auto pid = spawnShell(`convert "gifout/tempout-*.png" "`~output_name~`"`);
      wait(pid);
    }
    // {
    //   auto pid = spawnShell(`rm -rf gifout/`);
    //   wait(pid);
    // }
    return true;
  }
  return false;
}
