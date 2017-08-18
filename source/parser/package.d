module parser;
public import parser.checker;
import parser.file : KernelFile;
static import stl;

static immutable string DTOADQ_filename   = "dtoadqkernel.c",
                        Raytrace_filename = "raytracekernel.c",
                        Texture_filename  = "texturekernel.c";

/// Parses kernel and returns resulting kernel source code string
string Parse_Kernel ( ) {
  static import core.info;
  string function() kernel_func;
  if ( core.info.RType() == core.info.KernelType.Texture ) {
    static import parser.texture_parser;
    kernel_func = &parser.texture_parser.Parse;
  } else {
    static import parser.scene_parser;
    kernel_func = &parser.scene_parser.Parse;
  }
  return kernel_func();
}

static import core.info;
/// Changes Kernel
void Change_Kernel ( core.info.KernelType kernel_type ) {
  alias KT = core.info.KernelType;
  final switch ( kernel_type ) {
    case KT.VideoRender: assert(false, "Video Render kernel must be DTQ or RT");
    case KT.DTQ:      KernelFile.Create_File(DTOADQ_filename); break;
    case KT.Raytrace: KernelFile.Create_File(Raytrace_filename); break;
    case KT.Texture:  KernelFile.Create_File(Texture_filename); break;
  }
}

/// Changes Scene
void Change_Scene ( string filename ) {
  import parser.file : SceneFile;
  SceneFile.Create_File(filename);
}
