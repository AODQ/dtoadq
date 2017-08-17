module parser;

static immutable string DTOADQ_filename   = "dtoadqkernel.cl",
                        Raytrace_filename = "raytracekernel.cl",
                        Texture_filename  = "texturekernel.cl";

/// Parses kernel and returns resulting kernel source code string
string Parse_Kernel ( ) {
  static import kernel.info;
  string function() kernel_func;
  if ( kernel.info.RKernelType() == kernel.info.KernelType.Texture ) {
    static import parser.texture_parser;
    kernel_func = &parser.texture_parser.Parse;
  } else {
    static import parser.scene_parser;
    kernel_func = &parser.scene_parser.Parse;
  }
  return kernel_func();
}

static import kernel.info;
/// Changes Kernel
void Change_Kernel ( kernel.info.KernelType kernel_type ) {
  import parser.file : KernelFile;
  alias KT = kernel.info.KernelType;
  final switch ( kernel_type ) {
    case KT.VideoRender: assert(false, "Video Render kernel must be DTQ or RT");
    case KT.DTQ:      KernelFile.Create_File("dtoadqkernel.cl"  ); break;
    case KT.Raytrace: KernelFile.Create_File("raytracekernel.cl"); break;
    case KT.Texture:  KernelFile.Create_file("texturekernel.cl" ); break;
  }
}

/// Changes Scene
void Change_Scene ( string filename ) {
  import parser.file : SceneFile;
  SceneFile.Create_File(filename);
}
