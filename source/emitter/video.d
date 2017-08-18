module emitter.video;
import progress : Progress;
static import stl, y4md, core;

Progress progress;
y4md.Y4MWriter output;
string output_name;

float time, end_time;
int fps, spp;
private core.Resolution resolution;

auto RRender_Resolution ( ) { return resolution; }

void Render (string scene, string out_, int spp_, float end_time_, int fps_) {
  time = 0.0f;
  end_time = end_time_;
  fps = fps_;
  spp = spp_;
  resolution = core.Resolution.r1920_1080;

  // -- set up video emitter --
  output_name = out_;
  output = new y4md.Y4MWriter(out_,
        core.RImage(resolution).x, core.RImage(resolution).y,
        y4md.Rational(cast(int)fps, 1), y4md.Rational(32, 27),
        y4md.Interlacing.Progressive, y4md.Subsampling.C444);

  // -- set up progress bar --
  progress = new Progress(cast(int)(end_time*fps));
  progress.title = "Rendering";

  // -- set up scene/renderer --
  core.Set_Image_Buffer(resolution);
  core.info.Set_Map_Function(scene);
  core.info.Set_Kernel_Type(core.info.KernelType.VideoRender);
}

// Returns true when render is finished
bool Update ( ubyte[] data ) {
  time += 1.0f/fps;
  core.Set_Time(time);
  progress.next();
  output.writeFrame(data);
  if ( time >= end_time ) {
    import std.process;
    auto mp4proc = spawnShell("ffmpeg -i " ~ output_name ~ " " ~
                              stl.replace(output_name, "y4m", "mp4") ~ " -y");
    wait(mp4proc);
    stl.file.remove(output_name);
    return true;
  }
  return false;
}
