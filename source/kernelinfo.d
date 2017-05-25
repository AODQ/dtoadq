module kernelinfo;

// -- kernel compile shit --
struct Kernel_Info {
  enum Type { Raycast, Raytrace, MLT, }
  enum Flag { Show_Normals,           }
  enum Var  { March_Dist, March_Reps, }

  Type type;
  bool[Flag] flags;
  int [Var ] vars;
  string map_function;
}

private Kernel_Info kernel_info;
private bool        recompile = true;

auto Set_Kernel_Type ( Kernel_Info.Type type ) {
  recompile = true;
  kernel_info.type = type;
}

auto Set_Kernel_Flag ( Kernel_Info.Flag flag, bool status ) {
  recompile = true;
  kernel_info.flags[flag] = status;
}

auto Set_Kernel_Var  ( Kernel_Info.Var var, int value ) {
  recompile = true;
  kernel_info.vars[var] = value;
}

void Set_Map_Function ( string value ) {
  recompile = true;
  kernel_info.map_function = value;
}

auto Should_Recompile ( bool silent = true ) {
  recompile ^= !silent;
  return recompile^(!silent);
}

immutable(Kernel_Info) RKernel_Info ( ) { return cast(immutable)kernel_info; }

auto RKernel_Type_String ( ) {
  import opencl_kernel;
  final switch ( kernel_info.type ) with ( Kernel_Info.Type ) {
    case Raycast:  return Raycast_kernel_function;
    case Raytrace: return Raytrace_kernel_function;
    case MLT:      return MLT_kernel_function;
  }
}

static this ( ) {
  with ( Kernel_Info ) {
    kernel_info = Kernel_Info(
      Type.Raycast,
      [ Flag.Show_Normals : false ],
      [ Var.March_Dist : 64, Var.March_Reps : 128 ]
    );
  }
}
