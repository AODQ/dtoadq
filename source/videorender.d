module videorender;
import globals;
import progress : Progress;
import YFM = y4md;
import DTOADQ = dtoadq;
import KI = kernelinfo;
import DIMG = dtoadqimage;

void Render (string scene, string outname, int spp, float end_time, int fps) {
  float time = 0.0f;

  // -- set up video emitter --
  auto output = new YFM.Y4MWriter(outname, 1920, 1080,
                  YFM.Rational(cast(int)fps, 1));
  float[] floatdata;
  floatdata.length = 1920*1080*4;

  // -- set up progress bar --
  Progress p = new Progress(cast(int)(end_time*fps));
  p.title = "Rendering";

  // -- set up scene/renderer --
  DTOADQ.Initialize(false);
  DTOADQ.Set_Image_Buffer(DIMG.Resolution.r1920_1080);
  KI.Set_Map_Function(KI.ProceduralType.Scene, KI.FileType.DTQ, scene);


  // -- render loop --
  while ( time < end_time ) {
    DTOADQ.Set_Time(time);
    DTOADQ.Update(true, floatdata);
    { // write results
      import functional;
      output.writeFrame(floatdata.map!(n => cast(ubyte)(n*255.0f)).array);
    }
    { // done with frame, set time, update progress, sleep
      time += 1.0f/fps;
      p.next();
      import core.time, core.thread;
      Thread.sleep(dur!("msecs")(50));
    }
  }

  writeln();
  writeln("Converting to mp4");
  import std.process;
  auto mp4proc = spawnShell("ffmpeg -i "~outname~" out.mp4");
  wait(mp4proc);
  writeln("FINISHED! :)");
}
