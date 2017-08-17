module emitter.video;
import progress : Progress;
import y4md, dtoadq;

Progress progress;
YFM.Y4MWriter output;
string output_name;

float time, end_time;
int fps, spp;
private DIMG.Resolution resolution;

auto RRender_Resolution ( ) { return resolution; }

void Render (string scene, string out_, int spp_, float end_time_, int fps_) {
  time = 0.0f;
  end_time = end_time_;
  fps = fps_;
  spp = spp_;
  resolution = DIMG.Resolution.r1920_1080;

  // -- set up video emitter --
  output_name = out_;
  output = new YFM.Y4MWriter(out_,
        DIMG.RImage(resolution).x, DIMG.RImage(resolution).y,
        YFM.Rational(cast(int)fps, 1), YFM.Rational(32, 27),
        YFM.Interlacing.Progressive, YFM.Subsampling.C444);

  // -- set up progress bar --
  progress = new Progress(cast(int)(end_time*fps));
  progress.title = "Rendering";

  // -- set up scene/renderer --
  DTOADQ.Set_Image_Buffer(resolution);
  KI.Set_Map_Function(scene);
  KI.Set_Kernel_Type(KI.KernelType.VideoRender);
}

// Returns true when render is finished
bool Update ( ubyte[] data ) {
  time += 1.0f/fps;
  DTOADQ.Set_Time(time);
  progress.next();
  output.writeFrame(data);
  if ( time >= end_time ) {
    import std.process;
    auto mp4proc = spawnShell("ffmpeg -i " ~ output_name ~ " " ~
        output_name.replace("y4m", "mp4") ~ " -y");
    wait(mp4proc);
    static import std.file;
    std.file.remove(output_name);
    return true;
  }
  return false;
}
